#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "TRandom3.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"

#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TPad.h"
#include "TLegend.h"
#include "TString.h"
#include "TColor.h"

//////////////////////////////////////////////////////////////////////////// class

const int bins_theta = 60;
const int bins_dm2   = 60;

class TClass_CLs {
public:
  TClass_CLs() {
    cout<<endl<<" ---> The beginning ... "<<endl<<endl;
    
    h2d_space = new TH2D("h2d_space", "h2d_space", bins_theta, -4, 0, bins_dm2, -1, 2);
    
    dm2_low = 0.61;/// hhack, ue = 0.15, ee = 0.65
    dm2_hgh = 94;
    dm2_step = 0.05;
  }

  /////////////////////////////// data memeber

  double dm2_low;
  double dm2_hgh;
  double dm2_step;

  TH2D *h2d_space;

  /// distribution
  vector<double> array_vec_dchi2_4vPseudo[bins_dm2+1][bins_theta+1];
  vector<double> array_vec_dchi2_3vPseudo[bins_dm2+1][bins_theta+1];

  /// sensitivity and error band
  vector<double> array_vec_dchi2_sm_toys[bins_dm2+1][bins_theta+1];

  /// data result
  double array_dchi2_data[bins_dm2+1][bins_theta+1];

  ///
  vector<TGraph*>vec_graph_95_toys;
  vector<TGraph*>vec_graph_95_toys_T;// Transpose

  ///
  TGraph *gh_sensitivity_median;
  TGraph *gh_sensitivity_1s_LL;
  TGraph *gh_sensitivity_1s_RR;
  TGraph *gh_sensitivity_1s_band;
  TGraph *gh_sensitivity_2s_LL;
  TGraph *gh_sensitivity_2s_RR;  
  TGraph *gh_sensitivity_2s_band;

  ///
  TGraph *gh_data_result;
  
  /////////////////////////////// memeber function

  void func_read_distribution(TString roofile);
  void plot_distribution(int idm2, int itheta, double dchi2, int index);  
  
  void func_read_toys(TString roofile);
  void func_CLs_toys();
  void plot_CLs_toys();
  
  double func_cal_CLs(int idm2, int itheta, double dchi2);
  void func_get_contours(TH2D *h2_CL, TGraph *gh_CL[3], int index);

  double func_quantiles(double qq, double dm2_val, bool flag_plot);
  void func_sensitivity();

  void func_data_CLs(TString roofile);
};

/////// ccc

void TClass_CLs::func_data_CLs(TString roofile)
{  
  TString roostr = "";
  
  TFile *inputfile = new TFile(roofile, "read");
  TTree *tree = (TTree*)inputfile->Get("tree");

  // Declaration of leaf types
  vector<int>     *vec_dm2;
  vector<int>     *vec_t14;
  vector<double>  *vec_chi2_4v;
  vector<double>  *vec_chi2_3v;

  // List of branches
  TBranch        *b_vec_dm2;   //!
  TBranch        *b_vec_t14;   //!
  TBranch        *b_vec_chi2_4v;   //!
  TBranch        *b_vec_chi2_3v;   //!
  
  // Set object pointer
  vec_dm2 = 0;
  vec_t14 = 0;
  vec_chi2_4v = 0;
  vec_chi2_3v = 0;
  
  // Set branch addresses and branch pointers
  tree->SetBranchAddress("vec_dm2", &vec_dm2, &b_vec_dm2);
  tree->SetBranchAddress("vec_t14", &vec_t14, &b_vec_t14);
  tree->SetBranchAddress("vec_chi2_4v", &vec_chi2_4v, &b_vec_chi2_4v);
  tree->SetBranchAddress("vec_chi2_3v", &vec_chi2_3v, &b_vec_chi2_3v);

  int entries = tree->GetEntries();
  
  cout<<" ---> func_data_CLs, entries "<<entries<<endl<<endl;

  tree->GetEntry(0);

  int vc_size = vec_dm2->size();

  for(int idx=0; idx<vc_size; idx++) {
    int idm2 = vec_dm2->at(idx);
    int itheta = vec_t14->at(idx);
    double dchi2 = vec_chi2_4v->at(idx) - vec_chi2_3v->at(idx);
    
    double val_CLs = func_cal_CLs(idm2, itheta, dchi2);	
    h2d_space->SetBinContent(itheta, idm2, val_CLs);
  }

  ///////
    
  const int Ncontour = 3;
  TGraph *gh_CLs_data[Ncontour];
  for(int index=0; index<Ncontour; index++) {
    roostr = TString::Format("gh_CLs_data_%d", index);
    gh_CLs_data[index] = new TGraph(); gh_CLs_data[index]->SetName(roostr);
  }
  func_get_contours( h2d_space, gh_CLs_data, 9876543);

  //gh_data_result = (TGraph*)gh_CLs_data[1]->Clone("gh_data_result");// 1 ---> 95%, 2 ---> 99.73%
  
  gh_data_result = new TGraph(); gh_data_result->SetName("gh_data_result");
  for(int idx=0; idx<gh_CLs_data[1]->GetN(); idx++) {
    double xx, yy;
    gh_CLs_data[1]->GetPoint(idx, xx, yy);
    if( yy<dm2_low ) continue;
    gh_data_result->SetPoint(gh_data_result->GetN(), xx, yy);
  }
  
  ///////

  if( 0 ) {
    roostr = "canv_data";
    TCanvas *canv_data = new TCanvas(roostr, roostr, 900, 650);
    canv_data->SetBottomMargin(0.15); canv_data->SetTopMargin(0.1); canv_data->SetLeftMargin(0.15); canv_data->SetRightMargin(0.1);
    canv_data->SetLogx();
    canv_data->SetLogy();

    roostr = "h2d_data";
    TH2D *h2d_data = new TH2D(roostr, "", 10, 1e-2, 1, 10, 0.1, 100);
    h2d_data->Draw();
    h2d_data->SetXTitle("sin^{2}2#theta");
    h2d_data->SetYTitle("#Deltam^{2}_{41} (eV^{2})");
  
    h2d_data->GetXaxis()->CenterTitle(); h2d_data->GetXaxis()->SetTitleSize(0.05); h2d_data->GetXaxis()->SetLabelSize(0.05);
    h2d_data->GetYaxis()->CenterTitle(); h2d_data->GetYaxis()->SetTitleSize(0.05); h2d_data->GetYaxis()->SetLabelSize(0.05);  
    h2d_data->GetXaxis()->SetTitleOffset(1.2);
    h2d_data->GetYaxis()->SetTitleOffset(1.0);
  
    gh_data_result->Draw("same l"); gh_data_result->SetLineColor(kBlack); gh_data_result->SetLineWidth(3);
  }
  
}

/////// ccc

void TClass_CLs::func_sensitivity()
{
  gh_sensitivity_median = new TGraph(); gh_sensitivity_median->SetName("gh_sensitivity_median");
  
  gh_sensitivity_1s_LL = new TGraph(); gh_sensitivity_1s_LL->SetName("gh_sensitivity_1s_LL");
  gh_sensitivity_1s_RR = new TGraph(); gh_sensitivity_1s_RR->SetName("gh_sensitivity_1s_RR");
  
  gh_sensitivity_2s_LL = new TGraph(); gh_sensitivity_2s_LL->SetName("gh_sensitivity_2s_LL");
  gh_sensitivity_2s_RR = new TGraph(); gh_sensitivity_2s_RR->SetName("gh_sensitivity_2s_RR");

  ///////
  
  double val_1sigma = 0.6827;
  double val_2sigma = 0.9545;

  
  for(int idx=1; idx<=10000; idx++) {
    double dm2_val = idx * dm2_step;
    if( dm2_val < dm2_low ) continue;
    if( dm2_val > dm2_hgh ) break;

    double qq = 0;
    double tt_val = 0;

    qq = 0.5;
    tt_val = func_quantiles( qq, dm2_val, 0 );
    gh_sensitivity_median->SetPoint(gh_sensitivity_median->GetN(), tt_val, dm2_val);

    qq = 0.5 - val_1sigma/2;
    tt_val = func_quantiles( qq, dm2_val, 0 );
    gh_sensitivity_1s_LL->SetPoint(gh_sensitivity_1s_LL->GetN(), tt_val, dm2_val);
 
    qq = 0.5 + val_1sigma/2;
    tt_val = func_quantiles( qq, dm2_val, 0 );
    gh_sensitivity_1s_RR->SetPoint(gh_sensitivity_1s_RR->GetN(), tt_val, dm2_val);
 
    qq = 0.5 - val_2sigma/2;
    tt_val = func_quantiles( qq, dm2_val, 0 );
    gh_sensitivity_2s_LL->SetPoint(gh_sensitivity_2s_LL->GetN(), tt_val, dm2_val);
 
    qq = 0.5 + val_2sigma/2;
    tt_val = func_quantiles( qq, dm2_val, 0 );
    gh_sensitivity_2s_RR->SetPoint(gh_sensitivity_2s_RR->GetN(), tt_val, dm2_val);
    
  }// for(int idx=1; idx<=10000; idx++)

  ///////
  gh_sensitivity_1s_band = new TGraph(); gh_sensitivity_1s_band->SetName("gh_sensitivity_1s_band");
  
  for(int idx=0; idx<gh_sensitivity_1s_LL->GetN(); idx++) {
    double xx, yy;
    gh_sensitivity_1s_LL->GetPoint(idx, xx, yy);
    gh_sensitivity_1s_band->SetPoint(gh_sensitivity_1s_band->GetN(), xx, yy);
  }
  
  for(int idx=gh_sensitivity_1s_RR->GetN()-1; idx>=0; idx--) {
    double xx, yy;
    gh_sensitivity_1s_RR->GetPoint(idx, xx, yy);
    gh_sensitivity_1s_band->SetPoint(gh_sensitivity_1s_band->GetN(), xx, yy);
  }
  
  ///////
  gh_sensitivity_2s_band = new TGraph(); gh_sensitivity_2s_band->SetName("gh_sensitivity_2s_band");
  
  for(int idx=0; idx<gh_sensitivity_2s_LL->GetN(); idx++) {
    double xx, yy;
    gh_sensitivity_2s_LL->GetPoint(idx, xx, yy);
    gh_sensitivity_2s_band->SetPoint(gh_sensitivity_2s_band->GetN(), xx, yy);
  }
  
  for(int idx=gh_sensitivity_2s_RR->GetN()-1; idx>=0; idx--) {
    double xx, yy;
    gh_sensitivity_2s_RR->GetPoint(idx, xx, yy);
    gh_sensitivity_2s_band->SetPoint(gh_sensitivity_2s_band->GetN(), xx, yy);
  }
  
  
  //////////////////////////////////////////////// plotting

  if( 1 ) {
    TString roostr = "";

    roostr = "canv_sens";
    TCanvas *canv_sens = new TCanvas(roostr, roostr, 900, 650);
    canv_sens->SetBottomMargin(0.15); canv_sens->SetTopMargin(0.1); canv_sens->SetLeftMargin(0.15); canv_sens->SetRightMargin(0.1);
    canv_sens->SetLogx();
    canv_sens->SetLogy();

    roostr = "h2d_sens";
    TH2D *h2d_sens = new TH2D(roostr, "", 10, 1e-2, 1, 10, 0.1, 100);
    h2d_sens->Draw();
    h2d_sens->SetXTitle("sin^{2}2#theta");
    h2d_sens->SetYTitle("#Deltam^{2}_{41} (eV^{2})");
  
    h2d_sens->GetXaxis()->CenterTitle(); h2d_sens->GetXaxis()->SetTitleSize(0.05); h2d_sens->GetXaxis()->SetLabelSize(0.05);
    h2d_sens->GetYaxis()->CenterTitle(); h2d_sens->GetYaxis()->SetTitleSize(0.05); h2d_sens->GetYaxis()->SetLabelSize(0.05);  
    h2d_sens->GetXaxis()->SetTitleOffset(1.2);
    h2d_sens->GetYaxis()->SetTitleOffset(1.0);
  
    gh_sensitivity_2s_band->Draw("same f"); gh_sensitivity_2s_band->SetFillColor(kYellow); gh_sensitivity_2s_band->SetLineWidth(0); 
    gh_sensitivity_1s_band->Draw("same f"); gh_sensitivity_1s_band->SetFillColor(kGreen); gh_sensitivity_1s_band->SetLineWidth(0); 
    gh_sensitivity_median->Draw("same l"); gh_sensitivity_median->SetLineColor(kGray+2); gh_sensitivity_median->SetLineWidth(3);

    if( gh_data_result!=NULL ) {
      gh_data_result->Draw("same l"); gh_data_result->SetLineColor(kBlack); gh_data_result->SetLineWidth(3);
    }
  }
  
}

/////// ccc

double TClass_CLs::func_quantiles(double qq, double dm2_val, bool flag_plot)
{
  double val_median = 0;

  TH1D *h1d_temp = new TH1D("h1d_temp", "h1d_temp", 1000, 0, 1);

  vector<double>vc_exact_cal;
  
  int vc_size = vec_graph_95_toys_T.size();

  for(int idx=0; idx<vc_size; idx++) {
    int nps = vec_graph_95_toys_T.at(idx)->GetN();

    double dm2_edge1, tt_edge1, dm2_edge2, tt_edge2;
    vec_graph_95_toys_T.at(idx)->GetPoint(0, dm2_edge1, tt_edge1);
    vec_graph_95_toys_T.at(idx)->GetPoint(nps-1, dm2_edge2, tt_edge2);
    bool flag_in = false;
    if( (dm2_val>=dm2_edge1 && dm2_val<=dm2_edge2) || (dm2_val>=dm2_edge2 && dm2_val<=dm2_edge1) ) flag_in = true;
    if( !flag_in ) continue;

    double tt_val = vec_graph_95_toys_T.at(idx)->Eval(dm2_val);
    h1d_temp->Fill(tt_val);

    vc_exact_cal.push_back(tt_val);
  }
  
  double xx;
  h1d_temp->GetQuantiles(1, &xx, &qq);
  val_median = xx;

  ///////

  if( 1 ) {// exact value
    sort( vc_exact_cal.begin(), vc_exact_cal.end() );
    int size_vc = vc_exact_cal.size();
    for(int idx=0; idx<size_vc; idx++) {
      if( idx+1>=size_vc*qq ) {
	val_median = vc_exact_cal.at(idx);
	break;
      }
    }    
  }
  
  ///////
  
  if( flag_plot ) {
    TCanvas *canv_temp = new TCanvas("canv_temp", "canv_temp", 900, 650);
    h1d_temp->Draw("hist");
  }

  if( !flag_plot ) {
    delete h1d_temp;
  }
  
  return val_median;
}

/////// ccc

void TClass_CLs::plot_CLs_toys()
{
  TString roostr = "";

  roostr = "canv_CLs_toys";
  TCanvas *canv_CLs_toys = new TCanvas(roostr, roostr, 900, 650);
  canv_CLs_toys->SetBottomMargin(0.15); canv_CLs_toys->SetTopMargin(0.1); canv_CLs_toys->SetLeftMargin(0.15); canv_CLs_toys->SetRightMargin(0.1);
  canv_CLs_toys->SetLogx();
  canv_CLs_toys->SetLogy();

  roostr = "h2d_CLs_toys";
  TH2D *h2d_CLs_toys = new TH2D(roostr, "", 10, 1e-2, 1, 10, 0.1, 100);
  h2d_CLs_toys->Draw();
  h2d_CLs_toys->SetXTitle("sin^{2}2#theta");
  h2d_CLs_toys->SetYTitle("#Deltam^{2}_{41} (eV^{2})");
  
  h2d_CLs_toys->GetXaxis()->CenterTitle(); h2d_CLs_toys->GetXaxis()->SetTitleSize(0.05); h2d_CLs_toys->GetXaxis()->SetLabelSize(0.05);
  h2d_CLs_toys->GetYaxis()->CenterTitle(); h2d_CLs_toys->GetYaxis()->SetTitleSize(0.05); h2d_CLs_toys->GetYaxis()->SetLabelSize(0.05);  
  h2d_CLs_toys->GetXaxis()->SetTitleOffset(1.2);
  h2d_CLs_toys->GetYaxis()->SetTitleOffset(1.0);
  
  int vc_size = vec_graph_95_toys.size();
  
  for(int idx=0; idx<vc_size; idx++) {    
    vec_graph_95_toys.at(idx)->Draw("same l");
    vec_graph_95_toys.at(idx)->SetLineColor(kGray+1);
  }
}

/////// ccc

void TClass_CLs::func_CLs_toys()
{
  TString roostr = "";
  
  int vc_size = array_vec_dchi2_sm_toys[1][1].size();
  cout<<" ---> func_CLs_toys, size "<<vc_size<<endl<<endl;

  for(int idx=0; idx<vc_size; idx++) {
    if( (idx+1)%50==0 ) cout<<"      ---> processing toy "<<idx+1<<endl;

    for(int idm2=1; idm2<=bins_dm2; idm2++) {
      for(int itheta=1; itheta<=bins_theta; itheta++) {

	double dchi2_toy = array_vec_dchi2_sm_toys[idm2][itheta].at(idx);
	double val_CLs = func_cal_CLs(idm2, itheta, dchi2_toy);	
	h2d_space->SetBinContent(itheta, idm2, val_CLs);
	  
      }// for(int itheta=1; itheta<=bins_theta; itheta++)
    }// for(int idm2=1; idm2<=bins_dm2; idm2++)

    ///////
    
    const int Ncontour = 3;
    TGraph *gh_CLs_toys[Ncontour];
    for(int index=0; index<Ncontour; index++) {
      roostr = TString::Format("gh_CLs_toys_%d", index);
      gh_CLs_toys[index] = new TGraph(); gh_CLs_toys[index]->SetName(roostr);
    }
    func_get_contours( h2d_space, gh_CLs_toys, idx);

    roostr = TString::Format("graph_95_toys_%06d", idx);
    TGraph *graph_95_toys = (TGraph*)gh_CLs_toys[1]->Clone(roostr);// 1 ---> 95%, 2 ---> 99.73%
    vec_graph_95_toys.push_back(graph_95_toys);

    ///////
    
    roostr = TString::Format("graph_95_toys_T_%06d", idx);
    TGraph *graph_95_toys_T = new TGraph(); graph_95_toys_T->SetName(roostr);
    int nps = graph_95_toys->GetN();
    for(int ip=0; ip<nps; ip++) {
      double xx, yy;
      graph_95_toys->GetPoint(ip, xx, yy);
      graph_95_toys_T->SetPoint(ip, yy, xx);
    }
    vec_graph_95_toys_T.push_back(graph_95_toys_T);
    
  }// for(int idx=0; idx<vc_size; idx++)
  
}

/////// ccc

double TClass_CLs::func_cal_CLs(int idm2, int itheta, double dchi2)
{
  double CL_sb = 0;
  double CL_b  = 0;  
  double CL_s = 1;
  
  double dchi2_data = dchi2;
  int vec_size = array_vec_dchi2_4vPseudo[idm2][itheta].size();
  
  int num_sb = 0;
  for(int idx=0; idx<vec_size; idx++) {
    if( dchi2_data < array_vec_dchi2_4vPseudo[idm2][itheta].at(idx) ) num_sb++;
  }
  
  int num_b = 0;
  for(int idx=0; idx<vec_size; idx++) {
    if( dchi2_data < array_vec_dchi2_3vPseudo[idm2][itheta].at(idx) ) num_b++;
  }

  CL_sb = num_sb*1./vec_size;
  CL_b  = num_b*1./vec_size;
  
  if( CL_b!=0 ) CL_s = CL_sb/CL_b;
  return 1 - CL_s;  
}

/////// ccc

void TClass_CLs::func_get_contours(TH2D *h2_CL, TGraph *gh_CL[3], int index)
{
  gROOT->SetBatch( 1 );
  
  TString roostr = "";
  
  const int Ncontour = 3;
  double contours[Ncontour] = {0};
  contours[0] = 0.9;
  contours[1] = 0.95;
  contours[2] = 0.9973;
  
  ///////
  roostr = TString::Format("canv_h2_CL_%d", index);
  TCanvas *canv_h2_CL = new TCanvas(roostr, roostr, 800, 600);
  h2_CL->SetStats(0);
  h2_CL->SetContour(Ncontour, contours);
  h2_CL->Draw("cont z list");
  canv_h2_CL->Update(); // Needed to force the plotting and retrieve the contours in TGraphs
     
  // Get Contours
  TObjArray *conts_mc = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
  TList* contLevel_mc = NULL;
  TGraph* gh_curv_mc[10] = {0};

  Int_t nGraphs_mc    = 0;
  Int_t TotalConts_mc = 0;

  if (conts_mc == NULL){
    printf("*** No Contours Were Extracted!\n");
    TotalConts_mc = 0;
    exit(1);
  } else {
    TotalConts_mc = conts_mc->GetSize();
  }

  //printf("TotalConts_mc = %d\n", TotalConts_mc);

  for(int i = 0; i < TotalConts_mc; i++){
    contLevel_mc = (TList*)conts_mc->At(i);
    //printf("Contour %d has %d Graphs\n", i, contLevel_mc->GetSize());
    nGraphs_mc += contLevel_mc->GetSize();
  }

  nGraphs_mc = 0;
  for(int i = 0; i < TotalConts_mc; i++){
    contLevel_mc = (TList*)conts_mc->At(i);

    // Get first graph from list on curves on this level
    gh_curv_mc[i] = (TGraph*)contLevel_mc->First();     
  }

  ///////  
  for(int ic=0; ic<Ncontour; ic++) {
    int np_curv = gh_curv_mc[ic]->GetN();
    for(int idx=0; idx<np_curv; idx++) {
      double t14 = 0;
      double m41 = 0;
      gh_curv_mc[ic]->GetPoint(idx, t14, m41);
      
      double xx_value = pow(10., t14);
      
      /// user convert from sin2_theta_14 to sin2_2theta_14
      //if( xx_value>0.5 ) continue;      
      //xx_value = 4*xx_value*(1-xx_value);

      /// user convert from sin2_theta_24 to sin2_2theta_ue with fix t14
      //double fix_sin2_theta_14 = 0.99;
      //xx_value = 4*fix_sin2_theta_14*(1-fix_sin2_theta_14)*xx_value;
      
      gh_CL[ic]->SetPoint(gh_CL[ic]->GetN(), xx_value, pow(10., m41) );
      
    }    
  }

  delete canv_h2_CL;

  gROOT->SetBatch( 0 );
}

/////// ccc

void TClass_CLs::func_read_toys(TString roofile)
{
  TFile *inputfile = new TFile(roofile, "read");
  TTree *tree = (TTree*)inputfile->Get("tree");

  // Declaration of leaf types
  vector<int>     *vec_dm2;
  vector<int>     *vec_t14;
  vector<double>  *vec_chi2_4v;
  vector<double>  *vec_chi2_3v;

  // List of branches
  TBranch        *b_vec_dm2;   //!                                                                                                                                       
  TBranch        *b_vec_t14;   //!                                                                                                                                       
  TBranch        *b_vec_chi2_4v;   //!                                                                                                                                   
  TBranch        *b_vec_chi2_3v;   //!
   
  // Set object pointer
  vec_dm2 = 0;
  vec_t14 = 0;
  vec_chi2_4v = 0;
  vec_chi2_3v = 0;
   
  // Set branch addresses and branch pointers
  tree->SetBranchAddress("vec_dm2", &vec_dm2, &b_vec_dm2);
  tree->SetBranchAddress("vec_t14", &vec_t14, &b_vec_t14);
  tree->SetBranchAddress("vec_chi2_4v", &vec_chi2_4v, &b_vec_chi2_4v);
  tree->SetBranchAddress("vec_chi2_3v", &vec_chi2_3v, &b_vec_chi2_3v);

  int entries = tree->GetEntries();

  cout<<" ---> func_read_toys, entries "<<entries<<endl<<endl;

  for(int ientry=0; ientry<entries; ientry++) {
    tree->GetEntry(ientry);

    int vc_size = vec_dm2->size();

    for(int idx=0; idx<vc_size; idx++) {
      array_vec_dchi2_sm_toys[vec_dm2->at(idx)][vec_t14->at(idx)].push_back(vec_chi2_4v->at(idx)-vec_chi2_3v->at(idx));
    }
  }

  delete tree;
  delete inputfile;
  
}

/////// ccc

void TClass_CLs::plot_distribution(int idm2, int itheta, double dchi2, int index)
{
  TString roostr = "";

  double dm2_vv = pow(10, h2d_space->GetYaxis()->GetBinCenter(idm2));
  double theta_vv = pow(10, h2d_space->GetXaxis()->GetBinCenter(itheta));
  
  int vc_size = array_vec_dchi2_3vPseudo[idm2][itheta].size();
  
  double xmin_3v = array_vec_dchi2_3vPseudo[idm2][itheta].at(0);
  double xmax_3v = array_vec_dchi2_3vPseudo[idm2][itheta].at(vc_size-1);  
  double xmin_4v = array_vec_dchi2_4vPseudo[idm2][itheta].at(0);
  double xmax_4v = array_vec_dchi2_4vPseudo[idm2][itheta].at(vc_size-1);
  double xmin_vv = (xmin_3v>xmin_4v)?xmin_4v:xmin_3v;
  double xmax_vv = (xmax_3v>xmax_4v)?xmax_3v:xmax_4v;

  roostr = TString::Format("h1d_dchi2_3v_%02d_%02d", idm2, itheta); TH1D *h1d_dchi2_3v = new TH1D(roostr, "", 50, xmin_vv-1, xmax_vv+1);
  roostr = TString::Format("h1d_dchi2_4v_%02d_%02d", idm2, itheta); TH1D *h1d_dchi2_4v = new TH1D(roostr, "", 50, xmin_vv-1, xmax_vv+1);

  int lg_d4v = 0;
  int lg_d3v = 0;
  
  for(int idx=0; idx<vc_size; idx++) {
    h1d_dchi2_3v->Fill( array_vec_dchi2_3vPseudo[idm2][itheta].at(idx) );
    h1d_dchi2_4v->Fill( array_vec_dchi2_4vPseudo[idm2][itheta].at(idx) );

    if(  array_vec_dchi2_4vPseudo[idm2][itheta].at(idx) > dchi2 ) lg_d4v++;
    if(  array_vec_dchi2_3vPseudo[idm2][itheta].at(idx) > dchi2 ) lg_d3v++;    
  }

  
  if( 1 ) {/// hhack
    if( dchi2>0 ) {
      if( lg_d3v>0 ) {
	cout<<"  ---> CLs "<< 1 - lg_d4v*1./lg_d3v<<endl;
      }
    }
  }
  
  ///////

  double ymax_3v = h1d_dchi2_3v->GetMaximum();
  double ymax_4v = h1d_dchi2_4v->GetMaximum();
  double ymax_vv = (ymax_3v>ymax_4v)?ymax_3v:ymax_4v;
  h1d_dchi2_3v->SetMaximum( ymax_vv*1.2 );
  
  roostr = TString::Format("canv_dchi2_%02d_%02d", idm2, itheta);
  TCanvas *canv_dchi2 = new TCanvas(roostr, roostr, 900, 650);
  canv_dchi2->SetBottomMargin(0.15); canv_dchi2->SetTopMargin(0.1); canv_dchi2->SetLeftMargin(0.15); canv_dchi2->SetRightMargin(0.05);
  
  h1d_dchi2_3v->Draw("hist");
  h1d_dchi2_3v->SetLineColor(kBlue);
  h1d_dchi2_3v->SetXTitle("#Delta#chi^{2} = #chi^{2}_{4#nu} - #chi^{2}_{3#nu}");
  h1d_dchi2_3v->SetYTitle("Events");
  
  h1d_dchi2_3v->GetXaxis()->CenterTitle(); h1d_dchi2_3v->GetXaxis()->SetTitleSize(0.05); h1d_dchi2_3v->GetXaxis()->SetLabelSize(0.05);
  h1d_dchi2_3v->GetYaxis()->CenterTitle(); h1d_dchi2_3v->GetYaxis()->SetTitleSize(0.05); h1d_dchi2_3v->GetYaxis()->SetLabelSize(0.05);

  h1d_dchi2_4v->Draw("same hist");
  h1d_dchi2_4v->SetLineColor(kRed);

  h1d_dchi2_3v->Draw("same axis");

  if( 0 ) {
    TLine *line_dchi2 = new TLine(dchi2, 0, dchi2, ymax_vv*1.2*0.5); line_dchi2->Draw();
    line_dchi2->SetLineColor(kBlack); line_dchi2->SetLineWidth(3);
  }

  TLegend *lg = new TLegend(0.65, 0.62, 0.88, 0.85);
  roostr = TString::Format("lg_%02d_%02d", idm2, itheta); lg->SetName(roostr);
  lg->SetHeader(TString::Format("#splitline{#Deltam^{2}_{41} = %4.2f eV^{2}}{sin^{2}2#theta = %6.4f}", dm2_vv, theta_vv));
  lg->AddEntry( "", "", "");
  lg->AddEntry( h1d_dchi2_3v, TString::Format("#color[%d]{3#nu H is true}", kBlue), "l");
  lg->AddEntry( h1d_dchi2_4v, TString::Format("#color[%d]{4#nu H is true}", kRed), "l");
  lg->Draw();
  lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.05);
}

/////// ccc

void TClass_CLs::func_read_distribution(TString roofile)
{
  TFile *inputfile = new TFile(roofile, "read");
  TTree *tree = (TTree*)inputfile->Get("tree");

  // Declaration of leaf types
  Int_t           idm2;
  Int_t           it14;
  vector<double>  *vec_chi2_4v_on_4vPseudo;
  vector<double>  *vec_chi2_3v_on_4vPseudo;
  vector<double>  *vec_chi2_4v_on_3vPseudo;
  vector<double>  *vec_chi2_3v_on_3vPseudo;

  // List of branches
  TBranch        *b_idm2;   //!
  TBranch        *b_it14;   //!
  TBranch        *b_vec_chi2_4v_on_4vPseudo;   //!
  TBranch        *b_vec_chi2_3v_on_4vPseudo;   //!
  TBranch        *b_vec_chi2_4v_on_3vPseudo;   //!
  TBranch        *b_vec_chi2_3v_on_3vPseudo;   //!
 
  // Set object pointer
  vec_chi2_4v_on_4vPseudo = 0;
  vec_chi2_3v_on_4vPseudo = 0;
  vec_chi2_4v_on_3vPseudo = 0;
  vec_chi2_3v_on_3vPseudo = 0;
  
  // Set branch addresses and branch pointers
  tree->SetBranchAddress("idm2", &idm2, &b_idm2);
  tree->SetBranchAddress("it14", &it14, &b_it14);
  tree->SetBranchAddress("vec_chi2_4v_on_4vPseudo", &vec_chi2_4v_on_4vPseudo, &b_vec_chi2_4v_on_4vPseudo);
  tree->SetBranchAddress("vec_chi2_3v_on_4vPseudo", &vec_chi2_3v_on_4vPseudo, &b_vec_chi2_3v_on_4vPseudo);
  tree->SetBranchAddress("vec_chi2_4v_on_3vPseudo", &vec_chi2_4v_on_3vPseudo, &b_vec_chi2_4v_on_3vPseudo);
  tree->SetBranchAddress("vec_chi2_3v_on_3vPseudo", &vec_chi2_3v_on_3vPseudo, &b_vec_chi2_3v_on_3vPseudo); 

  int entries = tree->GetEntries();

  cout<<" ---> func_read_distribution, entries "<<entries<<endl<<endl;

  for(int ientry=0; ientry<entries; ientry++) {
    tree->GetEntry(ientry);

    int vc_size = vec_chi2_4v_on_4vPseudo->size();
    for(int idx=0; idx<vc_size; idx++) {
      double dchi2_4vPseudo = vec_chi2_4v_on_4vPseudo->at(idx) - vec_chi2_3v_on_4vPseudo->at(idx);
      double dchi2_3vPseudo = vec_chi2_4v_on_3vPseudo->at(idx) - vec_chi2_3v_on_3vPseudo->at(idx);

      array_vec_dchi2_4vPseudo[idm2][it14].push_back( dchi2_4vPseudo );
      array_vec_dchi2_3vPseudo[idm2][it14].push_back( dchi2_3vPseudo );      
    }// for(int idx=0; idx<vc_size; idx++)

    sort(array_vec_dchi2_4vPseudo[idm2][it14].begin(), array_vec_dchi2_4vPseudo[idm2][it14].end());
    sort(array_vec_dchi2_3vPseudo[idm2][it14].begin(), array_vec_dchi2_3vPseudo[idm2][it14].end());
    
  }// for(int ientry=0; ientry<entries; ientry++)

  delete tree;
  delete inputfile;
}

////////////////////////////////////////////////////////////////////////////
//////////////////////////// MAIN //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void ea_frequentis_CLs()
{
  //////////////////////////////////////////////////////////////////////////////////////// Draw style
  
  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kBird);

  double snWidth = 2;

  // use medium bold lines and thick markers
  gStyle->SetLineWidth(snWidth);
  gStyle->SetFrameLineWidth(snWidth);
  gStyle->SetHistLineWidth(snWidth);
  gStyle->SetFuncWidth(snWidth);
  gStyle->SetGridWidth(snWidth);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetEndErrorSize(4);
  gStyle->SetEndErrorSize(0);
 
  TString roostr = "";
  
  //////////////////////////////////////////////////////////////////////////////////////// hhack

  TString roofile_distribution = "./ub_new_numuDISonly_CLs/total_distribution.root";
  TString roofile_toys = "./ub_new_numuDISonly_CLs/total_toys.root";
  TString roofile_data = "./ub_new_numuDISonly_CLs/sub_sm_CLs_000001.root";
  //roofile_data = "./ub_new_numuDISonly_CLs/sub_sm_CLs_000001_3vAsimov.root";

  // roofile_distribution = "./ub_new_nueDISonly_CLs/total_distribution.root";
  // roofile_toys = "./ub_new_nueDISonly_CLs/total_toys.root";
  // roofile_data = "./ub_new_nueDISonly_CLs/sub_sm_CLs_000001.root";

  // roofile_distribution = "./ub_new_nueAPPonly_CLs/total_distribution.root";
  // roofile_toys = "./ub_new_nueAPPonly_CLs/total_toys.root";
  // roofile_data = "./ub_new_nueAPPonly_CLs/sub_sm_CLs_000001.root";
  
  ///////
  
  TClass_CLs *test_CLs = new TClass_CLs();

  /////// distribution
  test_CLs->func_read_distribution(roofile_distribution);
  //test_CLs->plot_distribution(46, 50, 64, 1);
  //test_CLs->plot_distribution(20, 10, 0, 2);

  
  /////// data
  test_CLs->func_data_CLs(roofile_data);

  /////// toys: sensitivity
  test_CLs->func_read_toys(roofile_toys);
  test_CLs->func_CLs_toys();
  ////test_CLs->plot_CLs_toys();
  test_CLs->func_sensitivity();

  //////////////////////////////////////////////

  TFile *outfile = new TFile("outfile_CLs.root", "recreate");
  
  test_CLs->gh_sensitivity_median->Write();
  test_CLs->gh_sensitivity_1s_band->Write();
  test_CLs->gh_sensitivity_2s_band->Write();
  test_CLs->gh_data_result->Write();

  outfile->Close();
 
}
