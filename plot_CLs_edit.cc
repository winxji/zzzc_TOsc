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

#include "TPolyLine3D.h"

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
      //if( xx_value>0.5 ) continue;      
      //xx_value = 4*xx_value*(1-xx_value);

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

void plot_CLs_edit()
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

  ////////////////////////////////////////////////////////////////////////////////////////

  TString roostr = "";
       
  ///////
  const int Ncontour = 3;
  //int colors[Ncontour] = {kBlue, kGreen+3, kRed};
  int colors[Ncontour] = {kBlack, kBlue, kRed};
  
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  
  int bins_theta = 60;
  int bins_dm2   = 60;
  
  /////// X: sin2_theta, 1e-3 -> 1   ---> "log10()" ---> -3 -> 0
  /////// Y: m41^2,      1e-1 -> 20  ---> "log10()" ---> -1 -> 1.30103      
  TH2D *h2_space_dm2_theta = new TH2D("h2_space_dm2_theta", "h2_space_dm2_theta", bins_theta, -4, 0, bins_dm2, -1, 2);
  TH2D *h2_space_t24_t14 = new TH2D("h2_space_t24_t14", "h2_space_t24_t14", bins_theta, -4, 0, bins_theta, -4, 0);


  //roostr = "./tresult_note/sum_BNBonly_APPDIS.dat";
  roostr = "./tresult_note/sum_both_APPDIS.dat";
  
  
  map<int, map<int, map<int, double> > >map_t14_dm2_t24_predCL;
  
  ifstream InputFile(roostr, ios::in);
  if(!InputFile) { cerr<<" No input-list"<<endl; exit(1); }

  
  for(int idx=1; idx<=bins_theta; idx++)
    for(int jdx=1; jdx<=bins_dm2; jdx++)
      for(int kdx=1; kdx<=bins_theta; kdx++)
	{
	  map_t14_dm2_t24_predCL[idx][jdx][kdx] = -1;
	}

  
  for(int idx=1; idx<=bins_theta*bins_dm2*bins_theta; idx++) {
    //for(int idx=1; idx<=215969; idx++) {
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

    map_t14_dm2_t24_predCL[theta][dm2][t24] = CL_pred;
    map_t14_dm2_t24_predCL[theta][dm2][t24] = CL_data;
  }

  for(int idx=1; idx<=60; idx++)
    for(int jdx=1; jdx<=60; jdx++)
      for(int kdx=1; kdx<=60; kdx++)
	{
	  if( map_t14_dm2_t24_predCL[idx][jdx][kdx]==-1 ) map_t14_dm2_t24_predCL[idx][jdx][kdx] = 100;
	}
  
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  TFile *outfile = new TFile("outfile_CL.root", "recreate");

  
  if( 1 ) {// dm2 vs. t14  <--- fixing t24
    cout<<endl<<" ---> processing dm2 vs. t14"<<endl<<endl;
    
    for(int it24=1; it24<=bins_theta; it24++) {
      cout<<TString::Format(" ---> %2d", it24)<<endl;
      
      TGraph *gh_CLs_pred[Ncontour];
      for(int idx=0; idx<Ncontour; idx++) {
	roostr = TString::Format("gh_CLs_pred_dm2VSt14_%d_%d", it24, idx+1);
	gh_CLs_pred[idx] = new TGraph(); gh_CLs_pred[idx]->SetName(roostr); gh_CLs_pred[idx]->SetLineColor( colors[idx] );
      }
      
      for(int it14=1; it14<=bins_theta; it14++) {
	for(int idm2=1; idm2<=bins_theta; idm2++) {	  
	  h2_space_dm2_theta->SetBinContent(it14, idm2, map_t14_dm2_t24_predCL[it14][idm2][it24]);	  	  
	}// for(int idm2=1; idm2<=bins_theta; idm2++)
      }// for(int it14=1; it14<=bins_theta; it14++)

      func_get_contours( h2_space_dm2_theta, gh_CLs_pred, it24);
      for(int idx=0; idx<Ncontour; idx++) gh_CLs_pred[idx]->Write();      
    }// for(int it24=1; it24<=bins_theta; it24++)    
  }

  
  if( 1 ) {// dm2 vs. t24  <--- fixing t14
    cout<<endl<<" ---> processing dm2 vs. t24"<<endl<<endl;
    
    for(int it14=1; it14<=bins_theta; it14++) {
      cout<<TString::Format(" ---> %2d", it14)<<endl;
      
      TGraph *gh_CLs_pred[Ncontour];
      for(int idx=0; idx<Ncontour; idx++) {
	roostr = TString::Format("gh_CLs_pred_dm2VSt24_%d_%d", it14, idx+1);
	gh_CLs_pred[idx] = new TGraph(); gh_CLs_pred[idx]->SetName(roostr); gh_CLs_pred[idx]->SetLineColor( colors[idx] );
      }
      
      for(int it24=1; it24<=bins_theta; it24++) {
	for(int idm2=1; idm2<=bins_theta; idm2++) {	  
	  h2_space_dm2_theta->SetBinContent(it24, idm2, map_t14_dm2_t24_predCL[it14][idm2][it24]);
	}// for(int idm2=1; idm2<=bins_theta; idm2++)
      }// for(int it24=1; it24<=bins_theta; it24++)

      func_get_contours( h2_space_dm2_theta, gh_CLs_pred, it14*100);
      for(int idx=0; idx<Ncontour; idx++) gh_CLs_pred[idx]->Write();      
    }// for(int it14=1; it14<=bins_theta; it14++)    
  }

  
  if( 1 ) {// t24 vs. t14  <--- fixing dm2
    cout<<endl<<" ---> processing t24 vs. t14"<<endl<<endl;
    
    for(int idm2=1; idm2<=bins_theta; idm2++) {
      cout<<TString::Format(" ---> %2d", idm2)<<endl;
      
      TGraph *gh_CLs_pred[Ncontour];
      for(int idx=0; idx<Ncontour; idx++) {
	roostr = TString::Format("gh_CLs_pred_t24VSt14_%d_%d", idm2, idx+1);
	gh_CLs_pred[idx] = new TGraph(); gh_CLs_pred[idx]->SetName(roostr); gh_CLs_pred[idx]->SetLineColor( colors[idx] );
      }
      
      for(int it14=1; it14<=bins_theta; it14++) {
	for(int it24=1; it24<=bins_theta; it24++) {	  
	  h2_space_dm2_theta->SetBinContent(it14, it24, map_t14_dm2_t24_predCL[it14][idm2][it24]);
	}// for(int it24=1; it24<=bins_theta; it24++)
      }// for(int it14=1; it14<=bins_theta; it14++)

      func_get_contours( h2_space_dm2_theta, gh_CLs_pred, idm2*10000);
      for(int idx=0; idx<Ncontour; idx++) gh_CLs_pred[idx]->Write();      
    }// for(int idm2=1; idm2<=bins_theta; idm2++)    
  }
  
 

  
  outfile->Close();
  cout<<" ---> Finish writing files"<<endl;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  
  
}
