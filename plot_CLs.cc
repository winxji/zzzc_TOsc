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

#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

////////////////////////////////////////////////////////////////////// ccc

void func_get_contours(TH2D *h2_CL, TGraph *gh_CL[3], int index)
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

  printf("TotalConts_mc = %d\n", TotalConts_mc);

  for(int i = 0; i < TotalConts_mc; i++){
    contLevel_mc = (TList*)conts_mc->At(i);
    printf("Contour %d has %d Graphs\n", i, contLevel_mc->GetSize());
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
      if( xx_value>0.5 ) continue;      
      xx_value = 4*xx_value*(1-xx_value);

      /// user convert from sin2_theta_24 to sin2_2theta_ue with fix t14
      //double fix_sin2_theta_14 = 0.99;
      //xx_value = 4*fix_sin2_theta_14*(1-fix_sin2_theta_14)*xx_value;
      
      gh_CL[ic]->SetPoint(gh_CL[ic]->GetN(), xx_value, pow(10., m41) );
      
    }    
  }

  delete canv_h2_CL;
}

//////////////////////////////////////////////////////////////////
///////////////////////// MAIN ///////////////////////////////////
//////////////////////////////////////////////////////////////////

void plot_CLs()
{
 
  //////////////////////////////////////////////////////////////////////////////////////// Draw style
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

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

  ////////////////////////////////////////////////////////////////////////////////////////

  TString roostr = "";
       
  ///////
  const int Ncontour = 3;
  //int colors[Ncontour] = {kBlue, kGreen+3, kRed};
  int colors[Ncontour] = {kBlue, kBlack, kRed};
  
  ////////////////////////////////////////////////////////////////////////////////////////

  
  ////////////////////////////////////////////////////////////////////////////////////////
  
  TString file_roostr = "./zb_nuedisapp_NuMIBNB/sum_NuMIBNB_both_scalePOT_NuMIrun123_8d4_BNBrun123_1d9.dat";
  file_roostr = "./zd_numudisapp_NuMIBNB_numuCC_related7chs/sum_NuMIBNB_BNBonly.dat";
  file_roostr = "./zd_numudisapp_NuMIBNB_numuCC_related7chs/sum_NuMIBNB_NuMIonly.dat";
  file_roostr = "./zd_numudisapp_NuMIBNB_numuCC_related7chs/sum_NuMIBNB_both.dat";

  file_roostr = "./zm_nueapp/sum_both_t24_fix0d0.dat";
  file_roostr = "./zm_nueapp/sum_both_t24_fix0d1.dat";


  
  file_roostr = "./tresult_note/sum_BNBonly_nueDisonly.dat";
  file_roostr = "./tresult_note/sum_both_nueDisonly.dat";

  file_roostr = "./tresult_note/sum_BNBonly_nueAPPonly.dat";
  file_roostr = "./tresult_note/sum_both_nueAPPonly.dat";

  file_roostr = "./tresult_note/sum_BNBonly_nueDisonly_270m.dat";
  file_roostr = "./tresult_note/sum_BNBonly_nueDisonly_670m.dat";
  file_roostr = "./tresult_note/sum_BNBonly_nueDisonly_470m.dat";



  file_roostr = "./tresult_note/sum_BNBonly_nueDisonly.dat";
  file_roostr = "./tresult_note/sum_BNBonly_nueAPPonly.dat";


  
  file_roostr = "./tresult_note/sum_BNBonly_APPonly_670m.dat";


  
  file_roostr = "./tresult_note/sum_BNBonly_nueAPPonly.dat";
  //file_roostr = "./tresult_note/sum_both_nueAPPonly.dat";


  file_roostr = "./tresult_note/sum_BNBonly_numuDISonly_fix470m.dat";
  file_roostr = "./tresult_note/sum_BNBonly_numuDISonly.dat";
  
  //file_roostr = "./tresult_note/sum_BNBonly_nueDisonly.dat";

  file_roostr = "./xc_GaussCLs_numuDIS_rate/sum.txt";


  file_roostr = "./wa_note_nueDIS/sum.dat";
  file_roostr = "./wb_new_nueDIS_mcstat/sum.dat";


  
  file_roostr = "./ua_new_numuDISonly_GaussCLs/sum.dat";  
  file_roostr = "./ua_new_nueDISonly_GaussCLs/sum.dat";  
  file_roostr = "./ua_new_nueAPPonly_GaussCLs/sum.dat";  


  
  file_roostr = "sum_fullosc_fixed_t24_0d0126.dat";
  file_roostr = "sum_fullosc_fixed_t24_0d0126_noNCpi0Channel.dat";
  file_roostr = "sum_fullosc_fixed_t24_0d0126_NCpi0Channel_1bin.dat";
  file_roostr = "sum_fullosc_fixed_t24_0d0126_ct34_0d5.dat";

  
  
  int bins_theta = 60;
  int bins_dm2   = 60;
  
  ////// X: sin22t14, 1e-2 -> 1   ---> "log10()" ---> -2 -> 0
  ////// Y: m41^2,    1e-1 -> 20  ---> "log10()" ---> -1 -> 1.30103

  double xlow = -4;
  double xhgh = 0;
  
  roostr = "h2_space_data";
  TH2D *h2_space_data = new TH2D(roostr, roostr, bins_theta, xlow, xhgh, bins_dm2, -1, 2);

  roostr = "h2_space_pred";
  TH2D *h2_space_pred = new TH2D(roostr, roostr, bins_theta, xlow, xhgh, bins_dm2, -1, 2);

  roostr = "h2_space_pred_1sigma_plus";
  TH2D *h2_space_pred_1sigma_plus = new TH2D(roostr, roostr, bins_theta, xlow, xhgh, bins_dm2, -1, 2);

  roostr = "h2_space_pred_1sigma_minus";
  TH2D *h2_space_pred_1sigma_minus = new TH2D(roostr, roostr, bins_theta, xlow, xhgh, bins_dm2, -1, 2);

  roostr = "h2_space_pred_2sigma_plus";
  TH2D *h2_space_pred_2sigma_plus = new TH2D(roostr, roostr, bins_theta, xlow, xhgh, bins_dm2, -1, 2);

  roostr = "h2_space_pred_2sigma_minus";
  TH2D *h2_space_pred_2sigma_minus = new TH2D(roostr, roostr, bins_theta, xlow, xhgh, bins_dm2, -1, 2);


  //////////////////////////////////////////////

  roostr = "h2_wilk_chi2_data";
  TH2D *h2_wilk_chi2_data = new TH2D(roostr, roostr, bins_theta, xlow, xhgh, bins_dm2, -1, 2);
  roostr = "h2_wilk_chi2_pred";
  TH2D *h2_wilk_chi2_pred = new TH2D(roostr, roostr, bins_theta, xlow, xhgh, bins_dm2, -1, 2);
  
  roostr = "h2_wilk_CL_data";
  TH2D *h2_wilk_CL_data = new TH2D(roostr, roostr, bins_theta, xlow, xhgh, bins_dm2, -1, 2);
  roostr = "h2_wilk_CL_pred";
  TH2D *h2_wilk_CL_pred = new TH2D(roostr, roostr, bins_theta, xlow, xhgh, bins_dm2, -1, 2);
  
  //////////////////////////////////////////////
  
  ifstream InputFile(file_roostr, ios::in);
  if(!InputFile) { cerr<<" No input-list"<<endl; exit(1); }

  for(int idx=1; idx<=bins_theta*bins_dm2; idx++) {
    int theta(0), dm2(0), t24(0);    
    double chi2_4v_4vAsimov(0), chi2_3v_4vAsimov(0),
      chi2_4v_3vAsimov(0), chi2_3v_3vAsimov(0),
      chi2_4v_data(0), chi2_3v_data(0),
      CL_data(0), CL_pred(0),
      CL_pred_1sigma_plus(0), CL_pred_1sigma_minus(0),
      CL_pred_2sigma_plus(0), CL_pred_2sigma_minus(0);

    InputFile >> theta >> dm2 >> t24
	      >> chi2_4v_4vAsimov >> chi2_3v_4vAsimov
	      >> chi2_4v_3vAsimov >> chi2_3v_3vAsimov
	      >> chi2_4v_data >> chi2_3v_data
	      >> CL_data >> CL_pred
	      >> CL_pred_1sigma_plus >> CL_pred_1sigma_minus
	      >> CL_pred_2sigma_plus >> CL_pred_2sigma_minus;

    CL_data *= 0.01;
    CL_pred *= 0.01;
    CL_pred_1sigma_plus  *= 0.01;
    CL_pred_1sigma_minus *= 0.01;
    CL_pred_2sigma_plus  *= 0.01;
    CL_pred_2sigma_minus *= 0.01;
    
    h2_space_data->SetBinContent( theta, dm2, CL_data );
    h2_space_pred->SetBinContent( theta, dm2, CL_pred );
    h2_space_pred_1sigma_plus->SetBinContent( theta, dm2, CL_pred_1sigma_plus );
    h2_space_pred_1sigma_minus->SetBinContent( theta, dm2, CL_pred_1sigma_minus );
    h2_space_pred_2sigma_plus->SetBinContent( theta, dm2, CL_pred_2sigma_plus );
    h2_space_pred_2sigma_minus->SetBinContent( theta, dm2, CL_pred_2sigma_minus );
  
    
    h2_wilk_chi2_pred->SetBinContent(theta, dm2, chi2_3v_4vAsimov);
    //double pvalue = TMath::Prob(chi2_3v_4vAsimov, 2); // wrong
    double pvalue = TMath::Prob(chi2_4v_3vAsimov, 2);// right
    double CL = 1 - pvalue;
    h2_wilk_CL_pred->SetBinContent(theta, dm2, CL);
  }

  ////////////////////////////////////////////////////////////

  int index = 0;
  
  // h2_space_pred_1sigma_minus->SaveAs("h2_space_pred_1sigma_minus.root");
  // h2_space_pred_1sigma_plus->SaveAs("h2_space_pred_1sigma_plus.root");
    
  ////////////////////////////////////////////////////////////
  
  index++;
  TGraph *gh_CLs_data[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_data_%d", idx);
    gh_CLs_data[idx] = new TGraph();
    gh_CLs_data[idx]->SetName(roostr);
    gh_CLs_data[idx]->SetLineColor( colors[idx] );
  }
  func_get_contours( h2_space_data, gh_CLs_data, index);

  ////////////////////////////////////////////////////////////
  
  index++;
  TGraph *gh_CLs_pred[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_pred_%d", idx);
    gh_CLs_pred[idx] = new TGraph();
    gh_CLs_pred[idx]->SetName(roostr);
    gh_CLs_pred[idx]->SetLineColor( colors[idx] );
  }
  func_get_contours( h2_space_pred, gh_CLs_pred, index);

  ////////////////////////////////////////////////////////////
  
  index++;
  TGraph *gh_CLs_pred_1sigma_plus[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_pred_1sigma_plus_%d", idx);
    gh_CLs_pred_1sigma_plus[idx] = new TGraph();
    gh_CLs_pred_1sigma_plus[idx]->SetName(roostr);
    gh_CLs_pred_1sigma_plus[idx]->SetLineColor( colors[idx] );
  }
  func_get_contours( h2_space_pred_1sigma_plus, gh_CLs_pred_1sigma_plus, index);

  index++;
  TGraph *gh_CLs_pred_1sigma_minus[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_pred_1sigma_minus_%d", idx);
    gh_CLs_pred_1sigma_minus[idx] = new TGraph();
    gh_CLs_pred_1sigma_minus[idx]->SetName(roostr);
    gh_CLs_pred_1sigma_minus[idx]->SetLineColor( colors[idx] );
  }
  func_get_contours( h2_space_pred_1sigma_minus, gh_CLs_pred_1sigma_minus, index);

  ////////////////////////////////////////////////////////////
  
  index++;
  TGraph *gh_CLs_pred_2sigma_plus[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_pred_2sigma_plus_%d", idx);
    gh_CLs_pred_2sigma_plus[idx] = new TGraph();
    gh_CLs_pred_2sigma_plus[idx]->SetName(roostr);
    gh_CLs_pred_2sigma_plus[idx]->SetLineColor( colors[idx] );
  }
  func_get_contours( h2_space_pred_2sigma_plus, gh_CLs_pred_2sigma_plus, index);

  index++;
  TGraph *gh_CLs_pred_2sigma_minus[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_CLs_pred_2sigma_minus_%d", idx);
    gh_CLs_pred_2sigma_minus[idx] = new TGraph();
    gh_CLs_pred_2sigma_minus[idx]->SetName(roostr);
    gh_CLs_pred_2sigma_minus[idx]->SetLineColor( colors[idx] );
  }
  func_get_contours( h2_space_pred_2sigma_minus, gh_CLs_pred_2sigma_minus, index);
  
  ////////////////////////////////////////////////////////////

  TGraph *gh_contour_1sigma[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_contour_1sigma_%d", idx);
    gh_contour_1sigma[idx] = new TGraph();
    gh_contour_1sigma[idx]->SetName(roostr);
    gh_contour_1sigma[idx]->SetFillColor(kGreen);
    
    int num_plus = gh_CLs_pred_1sigma_plus[idx]->GetN();
    for(int ip=0; ip<num_plus; ip++) {
      double xx(0), yy(0);
      gh_CLs_pred_1sigma_plus[idx]->GetPoint(ip, xx, yy);
      gh_contour_1sigma[idx]->SetPoint( gh_contour_1sigma[idx]->GetN(), xx, yy );
    }// ip
    
    int num_minus = gh_CLs_pred_1sigma_minus[idx]->GetN();
    for(int ip=num_minus-1; ip>=0; ip--) {
      double xx(0), yy(0);
      gh_CLs_pred_1sigma_minus[idx]->GetPoint(ip, xx, yy);
      gh_contour_1sigma[idx]->SetPoint( gh_contour_1sigma[idx]->GetN(), xx, yy );
    }// ip          
  }// idx
  
  TGraph *gh_contour_2sigma[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_contour_2sigma_%d", idx);
    gh_contour_2sigma[idx] = new TGraph();
    gh_contour_2sigma[idx]->SetName(roostr);
    gh_contour_2sigma[idx]->SetFillColor(kYellow);
    
    int num_plus = gh_CLs_pred_2sigma_plus[idx]->GetN();
    for(int ip=0; ip<num_plus; ip++) {
      double xx(0), yy(0);
      gh_CLs_pred_2sigma_plus[idx]->GetPoint(ip, xx, yy);
      gh_contour_2sigma[idx]->SetPoint( gh_contour_2sigma[idx]->GetN(), xx, yy );
    }// ip
    
    int num_minus = gh_CLs_pred_2sigma_minus[idx]->GetN();
    for(int ip=num_minus-1; ip>=0; ip--) {
      double xx(0), yy(0);
      gh_CLs_pred_2sigma_minus[idx]->GetPoint(ip, xx, yy);
      gh_contour_2sigma[idx]->SetPoint( gh_contour_2sigma[idx]->GetN(), xx, yy );
    }// ip          
  }// idx

  //////////////////////////////////////////////////////////// Wilk

  cout<<endl<<" ---> Wilk "<<endl<<endl;
  
  index++;
  TGraph *gh_wilk_CL_pred[Ncontour];
  for(int idx=0; idx<Ncontour; idx++) {
    roostr = TString::Format("gh_wilk_CL_pred_%d", idx);
    gh_wilk_CL_pred[idx] = new TGraph();
    gh_wilk_CL_pred[idx]->SetName(roostr);
    gh_wilk_CL_pred[idx]->SetLineColor( colors[idx] );
  }
  func_get_contours( h2_wilk_CL_pred, gh_wilk_CL_pred, index);

  
  //////////////////////////////////////////////////////////// plotting  
  //////////////////////////////////////////////////////////// plotting
  
  gROOT->SetBatch( 0 );  
  
  int index_90 = 0;
  int index_95 = 1;
  int index_99 = 2;
  
  double xxlow = h2_space_data->GetXaxis()->GetBinLowEdge(1);
  double xxhgh = h2_space_data->GetXaxis()->GetBinUpEdge(bins_theta);
    
  double yylow = h2_space_data->GetYaxis()->GetBinLowEdge(1);
  double yyhgh = h2_space_data->GetYaxis()->GetBinUpEdge(bins_dm2);

  
  TGraph *gh_n4 = new TGraph();
  gh_n4->SetPoint(0, 0.26, 7.2);
  gh_n4->SetMarkerStyle(22);
  gh_n4->SetMarkerSize(1.8);
  gh_n4->SetMarkerColor(kRed);
 
  TGraph *gh_n4new = new TGraph();
  gh_n4new->SetPoint(0, 0.36, 7.3);
  gh_n4new->SetMarkerStyle(22);
  gh_n4new->SetMarkerSize(1.8);
  gh_n4new->SetMarkerColor(kRed);

  
  roostr = "h2_basic_CLs_data";
  // TH2D *h2_basic_CLs_data = new TH2D(roostr, "",
  // 				     bins_theta, pow(10, xxlow), pow(10, xxhgh),
  // 				     bins_dm2, pow(10, yylow), pow(10, yyhgh));

  /// hhack
  TH2D *h2_basic_CLs_data = new TH2D(roostr, "",
				     10, 1e-2, 1,
				     10, 1e-1, 100);

  /////////////////////////////////////////////////////// 95
  
  roostr = "canv_h2_basic_CLs_data_95";
  TCanvas *canv_h2_basic_CLs_data_95 = new TCanvas(roostr, roostr, 800, 600);
  canv_h2_basic_CLs_data_95->SetLeftMargin(0.15);
  canv_h2_basic_CLs_data_95->SetRightMargin(0.1);
  canv_h2_basic_CLs_data_95->SetBottomMargin(0.18);
  canv_h2_basic_CLs_data_95->SetLogx();
  canv_h2_basic_CLs_data_95->SetLogy();
  
  h2_basic_CLs_data->Draw();
  h2_basic_CLs_data->SetXTitle("sin^{2}2#theta_{ee}");
  h2_basic_CLs_data->SetYTitle("#Deltam^{2}_{41} (eV^{2})");
  h2_basic_CLs_data->GetXaxis()->CenterTitle(1);
  h2_basic_CLs_data->GetYaxis()->CenterTitle(1);
  h2_basic_CLs_data->GetYaxis()->SetLabelSize(0.045);
  h2_basic_CLs_data->GetYaxis()->SetTitleSize(0.045);    
  h2_basic_CLs_data->GetXaxis()->SetLabelSize(0.045);
  h2_basic_CLs_data->GetXaxis()->SetTitleSize(0.045);
  h2_basic_CLs_data->GetYaxis()->SetTitleOffset(1.4);
  h2_basic_CLs_data->GetXaxis()->SetTitleOffset(1.4);

  //////////////
  
  //gh_contour_2sigma[index_95]->Draw("same f");
  //gh_contour_1sigma[index_95]->Draw("same f");
  
  //gh_CLs_data[index_95]->Draw("same l");
  gh_CLs_data[index_95]->SetLineColor(kBlack);
  gh_CLs_data[index_95]->SetLineWidth(3);
    
  gh_CLs_pred[index_95]->Draw("same l");
  gh_CLs_pred[index_95]->SetLineStyle(1);
  gh_CLs_pred[index_95]->SetLineColor(kRed);
  gh_CLs_pred[index_95]->SetLineWidth(3);

  gh_wilk_CL_pred[index_95]->Draw("same");
  gh_wilk_CL_pred[index_95]->SetLineStyle(1);
  gh_wilk_CL_pred[index_95]->SetLineColor(kBlue);  
  gh_wilk_CL_pred[index_95]->SetLineWidth(3);
 

  //gh_n4new->Draw("same p");
 
  /////////////////////////////////////////////////////// 99
    
  roostr = "canv_h2_basic_CLs_data_99";
  TCanvas *canv_h2_basic_CLs_data_99 = new TCanvas(roostr, roostr, 800, 600);
  canv_h2_basic_CLs_data_99->SetLeftMargin(0.15);
  canv_h2_basic_CLs_data_99->SetRightMargin(0.1);
  canv_h2_basic_CLs_data_99->SetBottomMargin(0.18);
  canv_h2_basic_CLs_data_99->SetLogx();
  canv_h2_basic_CLs_data_99->SetLogy();
  
  h2_basic_CLs_data->Draw();
  h2_basic_CLs_data->SetXTitle("sin^{2}2#theta_{ee}");
  h2_basic_CLs_data->SetYTitle("#Deltam^{2}_{41} (eV^{2})");
  h2_basic_CLs_data->GetXaxis()->CenterTitle(1);
  h2_basic_CLs_data->GetYaxis()->CenterTitle(1);
  h2_basic_CLs_data->GetYaxis()->SetLabelSize(0.045);
  h2_basic_CLs_data->GetYaxis()->SetTitleSize(0.045);    
  h2_basic_CLs_data->GetXaxis()->SetLabelSize(0.045);
  h2_basic_CLs_data->GetXaxis()->SetTitleSize(0.045);
  h2_basic_CLs_data->GetYaxis()->SetTitleOffset(1.4);
  h2_basic_CLs_data->GetXaxis()->SetTitleOffset(1.4);

  //////////////
  
  //gh_contour_2sigma[index_99]->Draw("same f");
  //gh_contour_1sigma[index_99]->Draw("same f");
  
  //gh_CLs_data[index_99]->Draw("same l");
  gh_CLs_data[index_99]->SetLineColor(kBlack);
  gh_CLs_data[index_99]->SetLineWidth(3);
    
  gh_CLs_pred[index_99]->Draw("same l");
  gh_CLs_pred[index_99]->SetLineStyle(1);
  gh_CLs_pred[index_99]->SetLineColor(kRed);
  gh_CLs_pred[index_99]->SetLineWidth(3);

  gh_wilk_CL_pred[index_99]->Draw("same");
  gh_wilk_CL_pred[index_99]->SetLineStyle(1);
  gh_wilk_CL_pred[index_99]->SetLineColor(kBlue);  
  gh_wilk_CL_pred[index_99]->SetLineWidth(3);
 

  //gh_n4new->Draw("same p");
 
  ////////////////////////////////////////////////////////////

  // roostr = "canv_h2_space_pred_1sigma_minus";
  // TCanvas *canv_h2_space_pred_1sigma_minus = new TCanvas(roostr, roostr, 900, 650);
  // h2_space_pred_1sigma_minus->Draw("colz");
  
  // roostr = "canv_h2_space_pred_2sigma_minus";
  // TCanvas *canv_h2_space_pred_2sigma_minus = new TCanvas(roostr, roostr, 900, 650);
  // h2_space_pred_2sigma_minus->Draw("colz");

  // cout<<h2_space_pred_1sigma_minus->GetBinContent(1,1)<<endl;
  // cout<<h2_space_pred_1sigma_minus->GetBinContent(60,60)<<endl;
  
  //////////////////////////////////////////////////////

  TFile *roofile = new TFile("za_roofile_GaussCLs.root", "recreate");

  h2_basic_CLs_data->Write();
  
  gh_CLs_data[index_95]->Write();
  gh_CLs_pred[index_95]->Write();
  gh_contour_1sigma[index_95]->Write();
  gh_contour_2sigma[index_95]->Write();  
  gh_wilk_CL_pred[index_95]->Write();  
    
  gh_CLs_data[index_99]->Write();
  gh_CLs_pred[index_99]->Write();
  gh_contour_1sigma[index_99]->Write();
  gh_contour_2sigma[index_99]->Write();  
  gh_wilk_CL_pred[index_99]->Write();  

  roofile->Close();
}
