#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<vector>
#include<map>
#include<set>

#include "WCPLEEANA/TOsc.h"

#include "WCPLEEANA/Configure_Osc.h"

#include "TApplication.h"

//#include <chrono> // timer
//auto time_start = chrono::high_resolution_clock::now();
//auto time_stop = chrono::high_resolution_clock::now();
//auto time_duration = chrono::duration_cast<chrono::seconds>(time_stop - time_start);
//cout<<endl<<" ---> check time duration "<<time_duration.count()<<endl<<endl;
//milliseconds, minutes

//
// When do minimization, the initial values should not be at the exact edge.
// Otherwise, the inital values may be automatically changed.
//

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  TString roostr = "";

  cout<<endl<<" ---> A story ..."<<endl<<endl;

  int ifile = 1;
  double scaleF_POT_BNB  = 1;
  double scaleF_POT_NuMI = 1;
  int display = 0;

  int it14 = 0;
  int idm2 = 0;
  int it24 = 0;
  
  for(int i=1; i<argc; i++) {
    if( strcmp(argv[i],"-f")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>ifile ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-pbnb")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>scaleF_POT_BNB ) ) { cerr<<" ---> Error scaleF_POT_BNB !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-pnumi")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>scaleF_POT_NuMI ) ) { cerr<<" ---> Error scaleF_POT_NuMI !"<<endl; exit(1); }
    }    
    if( strcmp(argv[i],"-d")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>display ) ) { cerr<<" ---> Error display !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-it14")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>it14 ) ) { cerr<<" ---> Error it14 !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-idm2")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>idm2 ) ) { cerr<<" ---> Error idm2 !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-it24")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>it24 ) ) { cerr<<" ---> Error it24 !"<<endl; exit(1); }
    }    
  }

  ///////////////////////////////////////////////////////////
  
  if( !display ) {
    gROOT->SetBatch( 1 );
  }
  
  TApplication theApp("theApp",&argc,argv);
  
  /////////////////////////////////////////////////////////// Draw style

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

  ///////////////////////////////////////////////////////////

  TOsc *osc_test = new TOsc();

  ///////

  osc_test->scaleF_POT_BNB  = scaleF_POT_BNB;
  osc_test->scaleF_POT_NuMI = scaleF_POT_NuMI;
     
  ///////

  osc_test->flag_syst_dirt   = Configure_Osc::flag_syst_dirt;
  osc_test->flag_syst_mcstat = Configure_Osc::flag_syst_mcstat;
  osc_test->flag_syst_flux   = Configure_Osc::flag_syst_flux;
  osc_test->flag_syst_geant  = Configure_Osc::flag_syst_geant;
  osc_test->flag_syst_Xs     = Configure_Osc::flag_syst_Xs;
  osc_test->flag_syst_det    = Configure_Osc::flag_syst_det;
  
  ///////

  osc_test->flag_NuMI_nueCC_from_intnue       = Configure_Osc::flag_NuMI_nueCC_from_intnue;
  osc_test->flag_NuMI_nueCC_from_overlaynumu  = Configure_Osc::flag_NuMI_nueCC_from_overlaynumu;
  osc_test->flag_NuMI_nueCC_from_appnue       = Configure_Osc::flag_NuMI_nueCC_from_appnue;
  osc_test->flag_NuMI_nueCC_from_appnumu      = Configure_Osc::flag_NuMI_nueCC_from_appnumu;
  osc_test->flag_NuMI_nueCC_from_overlaynueNC = Configure_Osc::flag_NuMI_nueCC_from_overlaynueNC;
  osc_test->flag_NuMI_nueCC_from_overlaynumuNC= Configure_Osc::flag_NuMI_nueCC_from_overlaynumuNC;
  
  osc_test->flag_NuMI_numuCC_from_overlaynumu   = Configure_Osc::flag_NuMI_numuCC_from_overlaynumu;
  osc_test->flag_NuMI_numuCC_from_overlaynue    = Configure_Osc::flag_NuMI_numuCC_from_overlaynue;
  osc_test->flag_NuMI_numuCC_from_appnue        = Configure_Osc::flag_NuMI_numuCC_from_appnue;
  osc_test->flag_NuMI_numuCC_from_appnumu       = Configure_Osc::flag_NuMI_numuCC_from_appnumu;
  osc_test->flag_NuMI_numuCC_from_overlaynumuNC = Configure_Osc::flag_NuMI_numuCC_from_overlaynumuNC;
  osc_test->flag_NuMI_numuCC_from_overlaynueNC  = Configure_Osc::flag_NuMI_numuCC_from_overlaynueNC;  
  
  osc_test->flag_NuMI_CCpi0_from_overlaynumu  = Configure_Osc::flag_NuMI_CCpi0_from_overlaynumu;
  osc_test->flag_NuMI_CCpi0_from_appnue       = Configure_Osc::flag_NuMI_CCpi0_from_appnue;
  osc_test->flag_NuMI_CCpi0_from_overlaynumuNC= Configure_Osc::flag_NuMI_CCpi0_from_overlaynumuNC;
  osc_test->flag_NuMI_CCpi0_from_overlaynueNC = Configure_Osc::flag_NuMI_CCpi0_from_overlaynueNC;
  
  osc_test->flag_NuMI_NCpi0_from_overlaynumu  = Configure_Osc::flag_NuMI_NCpi0_from_overlaynumu;
  osc_test->flag_NuMI_NCpi0_from_appnue       = Configure_Osc::flag_NuMI_NCpi0_from_appnue;
  osc_test->flag_NuMI_NCpi0_from_overlaynumuNC= Configure_Osc::flag_NuMI_NCpi0_from_overlaynumuNC;
  osc_test->flag_NuMI_NCpi0_from_overlaynueNC = Configure_Osc::flag_NuMI_NCpi0_from_overlaynueNC;


  ///////
  
  osc_test->flag_BNB_nueCC_from_intnue       = Configure_Osc::flag_BNB_nueCC_from_intnue;
  osc_test->flag_BNB_nueCC_from_overlaynumu  = Configure_Osc::flag_BNB_nueCC_from_overlaynumu;
  osc_test->flag_BNB_nueCC_from_appnue       = Configure_Osc::flag_BNB_nueCC_from_appnue;
  osc_test->flag_BNB_nueCC_from_appnumu      = Configure_Osc::flag_BNB_nueCC_from_appnumu;
  osc_test->flag_BNB_nueCC_from_overlaynueNC = Configure_Osc::flag_BNB_nueCC_from_overlaynueNC;
  osc_test->flag_BNB_nueCC_from_overlaynumuNC= Configure_Osc::flag_BNB_nueCC_from_overlaynumuNC;
  
  osc_test->flag_BNB_numuCC_from_overlaynumu   = Configure_Osc::flag_BNB_numuCC_from_overlaynumu;
  osc_test->flag_BNB_numuCC_from_overlaynue    = Configure_Osc::flag_BNB_numuCC_from_overlaynue;
  osc_test->flag_BNB_numuCC_from_appnue        = Configure_Osc::flag_BNB_numuCC_from_appnue;
  osc_test->flag_BNB_numuCC_from_appnumu       = Configure_Osc::flag_BNB_numuCC_from_appnumu;
  osc_test->flag_BNB_numuCC_from_overlaynumuNC = Configure_Osc::flag_BNB_numuCC_from_overlaynumuNC;
  osc_test->flag_BNB_numuCC_from_overlaynueNC  = Configure_Osc::flag_BNB_numuCC_from_overlaynueNC;  
  
  osc_test->flag_BNB_CCpi0_from_overlaynumu  = Configure_Osc::flag_BNB_CCpi0_from_overlaynumu;
  osc_test->flag_BNB_CCpi0_from_appnue       = Configure_Osc::flag_BNB_CCpi0_from_appnue;
  osc_test->flag_BNB_CCpi0_from_overlaynumuNC= Configure_Osc::flag_BNB_CCpi0_from_overlaynumuNC;
  osc_test->flag_BNB_CCpi0_from_overlaynueNC = Configure_Osc::flag_BNB_CCpi0_from_overlaynueNC;
  
  osc_test->flag_BNB_NCpi0_from_overlaynumu  = Configure_Osc::flag_BNB_NCpi0_from_overlaynumu;
  osc_test->flag_BNB_NCpi0_from_appnue       = Configure_Osc::flag_BNB_NCpi0_from_appnue;
  osc_test->flag_BNB_NCpi0_from_overlaynumuNC= Configure_Osc::flag_BNB_NCpi0_from_overlaynumuNC;
  osc_test->flag_BNB_NCpi0_from_overlaynueNC = Configure_Osc::flag_BNB_NCpi0_from_overlaynueNC;
  
  /////// set only one time
  
  osc_test->Set_default_cv_cov(Configure_Osc::default_cv_file,
                               Configure_Osc::default_dirtadd_file,
                               Configure_Osc::default_mcstat_file,
                               Configure_Osc::default_fluxXs_dir,
                               Configure_Osc::default_detector_dir);
  
  osc_test->Set_oscillation_base();
  
  /////// Set_oscillation_pars(double val_dm2_41, double val_sin2_2theta_14, double val_sin2_theta_24, double val_sin2_theta_34)
  
  double val_dm2_41         = 7.3;
  double val_sin2_2theta_14 = 0.1;
  double val_sin2_theta_24  = 0.01;
  double val_sin2_theta_34  = 0;

  /// standard order
  val_dm2_41         = 7.3;
  val_sin2_2theta_14 = 0.1;
  val_sin2_theta_24  = 0.01;  
  osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();// meas, CV, COV: all ready
  
  osc_test->Set_meas2fitdata();
  //osc_test->Set_asimov2fitdata();

  //osc_test->Set_toy_variations( 1 );
  //osc_test->Set_toy2fitdata( 1 );  
    
  ///////
  //osc_test->Plot_user();
  //osc_test->Minimization_OscPars_FullCov(7.2, 0.12, 0.008, 0, "str_flag_fixpar");
  //osc_test->Minimization_OscPars_FullCov(0.11, 0.1, 0.1, 0, "str_flag_fixpar");
  //osc_test->Minimization_OscPars_FullCov(100000, 0.2, 0, 0, "dm2_t24");

  ///////

  if( 0 ) {
    //osc_test->Minimization_OscPars_FullCov(1.29532, 0.935754, 0, 0, "str_flag_fixpar");
    osc_test->Minimization_OscPars_FullCov(2, 0.06, 0, 0, "abc");
    
    double pars_4v[4] = {osc_test->minimization_dm2_41_val, osc_test->minimization_sin2_2theta_14_val, osc_test->minimization_sin2_theta_24_val, 0};
    double chi2_min = osc_test->FCN( pars_4v );
    cout<<" ---> chi2_min: "<<chi2_min<<endl;

    double pars_3v[4] = {0};
    double chi2_null = osc_test->FCN( pars_3v );
    cout<<" ---> chi2_null: "<<chi2_null<<endl;
    cout<<endl;
    cout<<" ---> diff: "<<chi2_null - chi2_min<<endl;
  }

  if( 0 ) {    
    double t24_low  = 0;
    double t24_hgh  = 0.99;
    double t24_step = 0.01;

    for(int idx=1; idx<=10000; idx++) {
      double val_t24 = idx*t24_step;
      if( val_t24<t24_low ) continue;
      if( val_t24>t24_hgh ) break;
      osc_test->Minimization_OscPars_FullCov(0.10592537, 1.1e-4, val_t24, 0, "t24");

      int min_status               = osc_test->minimization_status;
      double min_chi2              = osc_test->minimization_chi2;
      double min_dm2_41_val        = osc_test->minimization_dm2_41_val;
      double min_sin2_theta_14_val = osc_test->minimization_sin2_2theta_14_val;
      double min_sin2_theta_24_val = osc_test->minimization_sin2_theta_24_val;
      
      cout<<" abcdefg "<<min_status<<"\t"<<min_chi2<<"\t"<<min_dm2_41_val<<"\t"<<min_sin2_theta_14_val<<"\t"<<min_sin2_theta_24_val<<endl;
    }
    
  }
  
  ///////

  if( 0 ) {
    TFile *roofile = new TFile("cov_cov.root", "recreate");
    osc_test->matrix_eff_newworld_abs_syst_total.Write("mat_cov");
    roofile->Close();
  }

  if( 0 ) {

    osc_test->Set_oscillation_pars(0, 0, 0, 0);
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    TMatrixD matrix_noosc = osc_test->matrix_eff_newworld_pred;
    TMatrixD matrix_noosc_cov = osc_test->matrix_eff_newworld_abs_syst_total;
    
    val_dm2_41         = 1.29532;
    val_sin2_2theta_14 = 0.935754;
    val_sin2_theta_24  = 0;    
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, 0);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready

    TFile *out_spectra = new TFile("out_spectra.root", "recreate");
    osc_test->matrix_eff_newworld_meas.Write("matrix_meas");    
    matrix_noosc.Write("matrix_pred_noosc");
    matrix_noosc_cov.Write("matrix_systcov_noosc");
    osc_test->matrix_eff_newworld_pred.Write("matrix_pred_bestfit");
    osc_test->matrix_eff_newworld_abs_syst_total.Write("matrix_systcov_bestfit");
    out_spectra->Close();
    
  }
  
  // BNB fixed baseline = 470m
  // chi2 86.6479
  // dm2 1.1885
  // t14 0.926119
  // t24 0.000107978

  // FVAL  = 86.6089426492080605
  // Edm   = 1.94853705477949532e-07
  // Nfcn  = 94
  // dm2_41	  = 1.23581	 +/-  0.643853	(limited)
  // sin2_theta_14	  = 0.935806	 +/-  0.0549119	(limited)
  // sin2_theta_24	  = 5.53627e-09	 +/-  0.0295546	(limited)

  // chi2_sm = 89.1536
  if( 0 ) {
    double pars_3v[4] = {0};
    cout<<osc_test->FCN( pars_3v )<<endl;
  }

  if( 0 ) {

    if( 1 ) {
      int bins_theta = 60;
      int bins_dm2   = 60;      
      TH2D *h2_space = new TH2D("h2_space_whole", "h2_space_whole", bins_theta, -4, 0, bins_dm2, -1, 2);
	
      cout<<endl;
      for(int idx=1; idx<=60; idx++) {
	double center_dm2 = h2_space->GetYaxis()->GetBinCenter(idx); center_dm2 = pow(10, center_dm2);
	double center_tt = h2_space->GetXaxis()->GetBinCenter(idx); center_tt = pow(10, center_tt);
	cout<<" ---> idx  "<<idx<<"\t"<<center_dm2<<"\t"<<center_tt<<endl;
      }
      cout<<endl;
    }
    
    /////// sin22tt = 0.107978, index 46
    /////// dm2 = 18.8365 eV2, index 46
    
  }
  
  ///////////////////////////////////////////////////////////
  
  if( 0 ) {
   
    val_dm2_41         = 18.8365;// index 46
    val_sin2_2theta_14 = 0.316228;// index 53
    val_sin2_theta_24  = 0;

    double pars_4v[4] = {val_dm2_41, val_sin2_2theta_14, 0, 0};
    double pars_3v[4] = {0};

    //////
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, 0, 0);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    osc_test->Set_asimov2fitdata();

    double chi2_4v_4v = osc_test->FCN( pars_4v );
    double chi2_3v_4v = osc_test->FCN( pars_3v );

    /////
    osc_test->Set_oscillation_pars(0, 0, 0, 0);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    osc_test->Set_asimov2fitdata();
    
    double chi2_4v_3v = osc_test->FCN( pars_4v );
    double chi2_3v_3v = osc_test->FCN( pars_3v );

    //////
    osc_test->Set_meas2fitdata();

    // for(int idx=1; idx<=26; idx++) {
    //   if( idx<=10 ) continue;      
    //   double fc = osc_test->matrix_fitdata_newworld(0, 26*2 + idx-1);
    //   osc_test->matrix_fitdata_newworld(0, 26*2 + idx-1) = fc / (1-0.05*(idx-11));      
    //   double pc = osc_test->matrix_fitdata_newworld(0, 26*3 + idx-1);
    //   osc_test->matrix_fitdata_newworld(0, 26*3 + idx-1) = pc / (1-0.05*(idx-11));
    //   //cout<<idx<<"\t"<<osc_test->matrix_fitdata_newworld(0, 26*2 + idx-1)<<"\t"<<osc_test->matrix_eff_newworld_pred(0, 26*2 + idx-1)<<endl;
    // }    
    
    // for(int idx=1; idx<=26; idx++) {
    //   if( idx>10 ) continue;      
    //   double fc = osc_test->matrix_fitdata_newworld(0, 26*2 + idx-1);
    //   osc_test->matrix_fitdata_newworld(0, 26*2 + idx-1) = fc / (1-0.04*(idx-11));      
    //   double pc = osc_test->matrix_fitdata_newworld(0, 26*3 + idx-1);
    //   osc_test->matrix_fitdata_newworld(0, 26*3 + idx-1) = pc / (1-0.04*(idx-11));
    //   //cout<<idx<<"\t"<<osc_test->matrix_fitdata_newworld(0, 26*2 + idx-1)<<"\t"<<osc_test->matrix_eff_newworld_pred(0, 26*2 + idx-1)<<endl;
    // }    
    
    double chi2_4v_dd = osc_test->FCN( pars_4v );
    double chi2_3v_dd = osc_test->FCN( pars_3v );
    
    //////
    cout<<endl;
    cout<<" chi2_4v_4v "<<chi2_4v_4v<<endl;
    cout<<" chi2_3v_4v "<<chi2_3v_4v<<endl;
    cout<<" chi2_4v_3v "<<chi2_4v_3v<<endl;
    cout<<" chi2_3v_3v "<<chi2_3v_3v<<endl;
    cout<<" chi2_4v_dd "<<chi2_4v_dd<<endl;
    cout<<" chi2_3v_dd "<<chi2_3v_dd<<endl;   
    cout<<endl;

    double delta_4v = chi2_4v_4v - chi2_3v_4v;
    double delta_3v = chi2_4v_3v - chi2_3v_3v;
    double delta_dd = chi2_4v_dd - chi2_3v_dd;
    
    double GaussCLs_ss = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v) * 100.;
    double GaussCLs_dd = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_dd) * 100.;

    cout<<endl<<" ---> "<<delta_4v<<"\t"<<delta_3v<<"\t"<<delta_dd<<"\t"<<GaussCLs_ss<<"\t"<<GaussCLs_dd<<endl<<endl;
  }
 
  /////////////////////////////////////////////////////////// 

  if( 0 ) {
    cout<<endl<<" ---------------> Gaussian CLs" <<endl<<endl;

    val_dm2_41         = 59.5662;
    val_sin2_2theta_14 = 0.217068;
    val_sin2_theta_24  = 0;
    
    double pars_4v[4] = {val_dm2_41, val_sin2_2theta_14, 0, 0};
    double pars_3v[4] = {0};     
    
    double chi2_4v_on_4vPseudo = 0;
    double chi2_3v_on_4vPseudo = 0;
    double chi2_4v_on_3vPseudo = 0;
    double chi2_3v_on_3vPseudo = 0;

    roostr = TString::Format("fc_%06d.root", ifile);
    TFile *subroofile = new TFile(roostr, "recreate");
    TTree *tree = new TTree("tree", "tree");
    
    tree->Branch( "chi2_4v_on_4vPseudo", &chi2_4v_on_4vPseudo, "chi2_4v_on_4vPseudo/D" );
    tree->Branch( "chi2_3v_on_4vPseudo", &chi2_3v_on_4vPseudo, "chi2_3v_on_4vPseudo/D" );
    tree->Branch( "chi2_4v_on_3vPseudo", &chi2_4v_on_3vPseudo, "chi2_4v_on_3vPseudo/D" );
    tree->Branch( "chi2_3v_on_3vPseudo", &chi2_3v_on_3vPseudo, "chi2_3v_on_3vPseudo/D" );

    for(int idx=1; idx<=10; idx++) {
      cout<<" ---> processing "<<idx<<endl;
      
      //////
      osc_test->Set_oscillation_pars(0, 0, 0, 0);
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT();// meas, CV, COV: all ready
      osc_test->Set_toy_variations( 1 );
      osc_test->Set_toy2fitdata( 1 );

      chi2_4v_on_3vPseudo = osc_test->FCN( pars_4v );
      chi2_3v_on_3vPseudo = osc_test->FCN( pars_3v );

      ///////
      osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, 0, 0);
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT();// meas, CV, COV: all ready
      osc_test->Set_toy_variations( 1 );
      osc_test->Set_toy2fitdata( 1 );
      
      chi2_4v_on_4vPseudo = osc_test->FCN( pars_4v );      
      chi2_3v_on_4vPseudo = osc_test->FCN( pars_3v );

      //////      
      cout<<" chi2_4v_on_4vPseudo "<<chi2_4v_on_4vPseudo<<endl;
      cout<<" chi2_3v_on_4vPseudo "<<chi2_3v_on_4vPseudo<<endl;
      cout<<" chi2_4v_on_3vPseudo "<<chi2_4v_on_3vPseudo<<endl;
      cout<<" chi2_3v_on_3vPseudo "<<chi2_3v_on_3vPseudo<<endl;

      tree->Fill();
    }
    
    tree->Write();
    subroofile->Close();
  }
  
  /////////////////////////////////////////////////////////// 

  if( 0 ) {
    val_dm2_41         = 7.3;
    val_sin2_2theta_14 = 0.5;
    val_sin2_theta_24  = 0.01;
    
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, 0);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready

    TFile *file = new TFile("cov_asimov.root", "recreate");
    osc_test->matrix_eff_newworld_abs_syst_total.Write("covmat");
    file->Close();
  }
   
  /////////////////////////////////////////////////////////// 
  
  if( 0 ) {
    cout<<endl<<" ---> Asimov scan grids"<<endl<<endl;
    
    val_dm2_41         = 7.3;
    val_sin2_2theta_14 = 0.1;
    val_sin2_theta_24  = 0.01;

    int    min_status             = 10;
    int    flag_negative          = 0;
    double min_chi2               = 0;
    double min_dm2_41_val         = 0;
    double min_sin2_2theta_14_val = 0;
    double min_sin2_theta_24_val  = 0;
    double min_sin2_theta_34_val  = 0;

    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, 0);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready   

    roostr = TString::Format("sub_fit_%06d.root", ifile);
    TFile *subroofile = new TFile(roostr, "recreate");
    TTree *tree = new TTree("tree", "tree");
    
    tree->Branch( "min_status",             &min_status,             "min_status/I" );
    tree->Branch( "flag_negative",          &flag_negative,          "flag_negative/I" );
    tree->Branch( "min_chi2",               &min_chi2,               "min_chi2/D" );
    tree->Branch( "min_dm2_41_val",         &min_dm2_41_val,         "min_dm2_41_val/D" );
    tree->Branch( "min_sin2_2theta_14_val", &min_sin2_2theta_14_val, "min_sin2_2theta_14_val/D" );
    tree->Branch( "min_sin2_theta_24_val",  &min_sin2_theta_24_val,  "min_sin2_theta_24_val/D" );
    tree->Branch( "min_sin2_theta_34_val",  &min_sin2_theta_34_val,  "min_sin2_theta_34_val/D" );

    double dm2_low  = 0;
    double dm2_hgh  = 10;
    double dm2_step = 0.1;
    int num_dm2 = (dm2_hgh-dm2_low)/dm2_step;
    num_dm2 = num_dm2;
    
    double t14_low  = 0;
    double t14_hgh  = 1;
    double t14_step = 0.01;
    int num_t14 = (t14_hgh-t14_low)/t14_step;
    num_t14 = num_t14;
    
    double t24_low  = 0;
    double t24_hgh  = 0.5;
    double t24_step = 0.005;
    int num_t24 = (t24_hgh-t24_low)/t24_step;
    num_t24 = num_t24;
              
    /////// X: sin22t14, 1e-4 -> 1   ---> "log10()" ---> -4 -> 0
    /////// Y: m41^2,    1e-1 -> 1e2  ---> "log10()" ---> -1 -> 2
    int bins_theta = 60;
    int bins_dm2   = 60;      
    TH2D *h2_space = new TH2D("h2_space_whole", "h2_space_whole", bins_theta, -4, 0, bins_dm2, -1, 2);
    h2_space->SetBinContent(1,1,1);

    
    for(int itoy=1; itoy<=1; itoy++) {

      cout<<endl<<" ---> proceesing toy "<<itoy<<endl;
      
      osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, 0);  
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT();// meas, CV, COV: all ready
      osc_test->Set_meas2fitdata();
      //osc_test->Set_asimov2fitdata();
      //osc_test->Set_toy_variations( 1 );
      //osc_test->Set_toy2fitdata( 1 );
      
      double pars_4v[4] = {0};

      //for(int jdm2=1; jdm2<=num_dm2; jdm2++) {
      for(int jdm2=ifile; jdm2<=ifile; jdm2++) {
	cout<<" ---> processing "<<jdm2<<endl;
	
	//for(int jt14=0; jt14<num_t14; jt14++) {
	for(int jt14=1; jt14<=bins_theta; jt14++) {
	  cout<<TString::Format("      jdm2 %3d, jt14 %3d", jdm2, jt14)<<endl;
	  
	  //for(int jt24=0; jt24<num_t24; jt24++) {
	  for(int jt24=1; jt24<=bins_theta; jt24++) {
	    double dm2_val = dm2_low + jdm2*dm2_step;
	    double t14_val = t14_low + jt14*t14_step;
	    double t24_val = t24_low + jt24*t24_step;
	    
	    dm2_val = pow(10, h2_space->GetYaxis()->GetBinCenter(jdm2) );
	    t14_val = pow(10, h2_space->GetXaxis()->GetBinCenter(jt14) );
	    t24_val = pow(10, h2_space->GetXaxis()->GetBinCenter(jt24) );
	      
	    pars_4v[0] = dm2_val;
	    pars_4v[1] = t14_val;
	    pars_4v[2] = t24_val;
	      
	    double chi2_val = osc_test->FCN( pars_4v );// 1000 times ~ 3 min
	      
	    min_status            = 0;
	    min_chi2              = chi2_val;
	    min_dm2_41_val        = dm2_val;
	    min_sin2_2theta_14_val= t14_val;
	    min_sin2_theta_24_val = t24_val;
	    min_sin2_theta_34_val = 0;
	      
	    tree->Fill();    
	      
	  }
	}
      }
    }

    tree->Write();
    subroofile->Close();
  }
  
  /////////////////////////////////////////////////////////// 
  
  if( 0 ) {
    cout<<endl<<" ---> fitting with fixed dm2 and t24"<<endl<<endl;
    
    val_dm2_41         = 7.3;
    val_sin2_2theta_14 = 0.1;
    val_sin2_theta_24  = 0.01;

    int    min_status             = 10;
    int    flag_negative          = 0;
    int    idx_dm2                = 0;
    int    idx_t24                = 0;
    double min_chi2               = 0;
    double min_dm2_41_val         = 0;
    double min_sin2_2theta_14_val = 0;
    double min_sin2_theta_24_val  = 0;
    double min_sin2_theta_34_val  = 0;

    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, 0);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready   

    roostr = TString::Format("sub_fit_%02d.root", idm2);
    TFile *subroofile = new TFile(roostr, "recreate");
    TTree *tree = new TTree("tree", "tree");
    
    tree->Branch( "min_status",             &min_status,             "min_status/I" );
    tree->Branch( "flag_negative",          &flag_negative,          "flag_negative/I" );
    tree->Branch( "idx_dm2",                &idx_dm2,                "idx_dm2/I" );
    tree->Branch( "idx_t24",                &idx_t24,                "idx_t24/I" );
    tree->Branch( "min_chi2",               &min_chi2,               "min_chi2/D" );
    tree->Branch( "min_dm2_41_val",         &min_dm2_41_val,         "min_dm2_41_val/D" );
    tree->Branch( "min_sin2_2theta_14_val", &min_sin2_2theta_14_val, "min_sin2_2theta_14_val/D" );
    tree->Branch( "min_sin2_theta_24_val",  &min_sin2_theta_24_val,  "min_sin2_theta_24_val/D" );
    tree->Branch( "min_sin2_theta_34_val",  &min_sin2_theta_34_val,  "min_sin2_theta_34_val/D" );

    /////// X: sin22t14, 1e-4 -> 1   ---> "log10()" ---> -4 -> 0
    /////// Y: m41^2,    1e-1 -> 1e2  ---> "log10()" ---> -1 -> 2
    int bins_theta = 60;
    int bins_dm2   = 60;      
    TH2D *h2_space = new TH2D("h2_space_whole", "h2_space_whole", bins_theta, -4, 0, bins_dm2, -1, 2);
    h2_space->SetBinContent(1,1,1);

    
    for(int itoy=1; itoy<=1; itoy++) {

      cout<<endl<<" ---> proceesing toy "<<itoy<<endl;
      
      osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, 0);  
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT();// meas, CV, COV: all ready
      osc_test->Set_meas2fitdata();
      //osc_test->Set_asimov2fitdata();
      //osc_test->Set_toy_variations( 1 );
      //osc_test->Set_toy2fitdata( 1 );
             
      double dm2_val = 0;
      double t24_val = 0;
	    
      dm2_val = pow(10, h2_space->GetYaxis()->GetBinCenter(idm2) );

      for(int jt24=1; jt24<=bins_theta; jt24++) {
	cout<<" ---> fit1par "<<idm2<<"\t"<<jt24<<endl;
	
	t24_val = pow(10, h2_space->GetXaxis()->GetBinCenter(jt24) );
	
	osc_test->Minimization_OscPars_FullCov(dm2_val, 0.1, t24_val, 0, "dm2_t24");
	
	min_status            = osc_test->minimization_status;
	min_chi2              = osc_test->minimization_chi2;
	min_dm2_41_val        = osc_test->minimization_dm2_41_val;
	min_sin2_2theta_14_val= osc_test->minimization_sin2_2theta_14_val;
	min_sin2_theta_24_val = osc_test->minimization_sin2_theta_24_val;
	min_sin2_theta_34_val = 0;

	idx_dm2 = idm2;
	idx_t24 = jt24;
	
	tree->Fill();
      }
    }

    tree->Write();
    subroofile->Close();
  }  
  
  /////////////////////////////////////////////////////////// 
  
  if( 0 ) {
    cout<<endl<<" ---> Pseudo Expt. and fitting"<<endl<<endl;
    
    val_dm2_41         = 0;
    val_sin2_2theta_14 = 0;
    val_sin2_theta_24  = 0;

    int    min_status             = 10;
    int    flag_negative          = 0;
    double min_chi2               = 0;
    double min_dm2_41_val         = 0;
    double min_sin2_2theta_14_val = 0;
    double min_sin2_theta_24_val  = 0;
    double min_sin2_theta_34_val  = 0;
    double min_dm2_41_err         = 0;
    double min_sin2_2theta_14_err = 0;
    double min_sin2_theta_24_err  = 0;
    double min_sin2_theta_34_err  = 0;
    double true_dm2_41_val         = val_dm2_41;
    double true_sin2_2theta_14_val = val_sin2_2theta_14;
    double true_sin2_theta_24_val  = val_sin2_theta_24;
    double true_sin2_theta_34_val  = 0;
    double chi2_null               = 0;
    
    vector<double> vec_asimov_data;
    vector<double> vec_pseudo_data;
    vector<double> vec_bestfit_pred;
    vector<double> vec_syst_rand;

    
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, 0);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready   
    for(int idx=0; idx<osc_test->matrix_eff_newworld_pred.GetNcols(); idx++) {
      vec_asimov_data.push_back( osc_test->matrix_eff_newworld_pred(0,idx) );
    }
    

    roostr = TString::Format("sub_fit_%06d.root", ifile);
    TFile *subroofile = new TFile(roostr, "recreate");
    TTree *tree = new TTree("tree", "tree");
    
    tree->Branch( "min_status",             &min_status,             "min_status/I" );
    tree->Branch( "flag_negative",          &flag_negative,          "flag_negative/I" );
    tree->Branch( "min_chi2",               &min_chi2,               "min_chi2/D" );
    tree->Branch( "min_dm2_41_val",         &min_dm2_41_val,         "min_dm2_41_val/D" );
    tree->Branch( "min_sin2_2theta_14_val", &min_sin2_2theta_14_val, "min_sin2_2theta_14_val/D" );
    tree->Branch( "min_sin2_theta_24_val",  &min_sin2_theta_24_val,  "min_sin2_theta_24_val/D" );
    tree->Branch( "min_sin2_theta_34_val",  &min_sin2_theta_34_val,  "min_sin2_theta_34_val/D" );
    tree->Branch( "min_dm2_41_err",         &min_dm2_41_err,         "min_dm2_41_err/D" );
    tree->Branch( "min_sin2_2theta_14_err", &min_sin2_2theta_14_err, "min_sin2_2theta_14_err/D" );
    tree->Branch( "min_sin2_theta_24_err",  &min_sin2_theta_24_err,  "min_sin2_theta_24_err/D" );
    tree->Branch( "min_sin2_theta_34_err",  &min_sin2_theta_34_err,  "min_sin2_theta_34_err/D" );
 
    tree->Branch( "true_dm2_41_val",        &true_dm2_41_val,         "true_dm2_41_val/D" );
    tree->Branch( "true_sin2_2theta_14_val",&true_sin2_2theta_14_val, "true_sin2_2theta_14_val/D" );
    tree->Branch( "true_sin2_theta_24_val", &true_sin2_theta_24_val,  "true_sin2_theta_24_val/D" );
    tree->Branch( "true_sin2_theta_34_val", &true_sin2_theta_34_val,  "true_sin2_theta_34_val/D" );
    
    tree->Branch( "chi2_null",               &chi2_null,               "chi2_null/D" );
    
    tree->Branch( "vec_asimov_data",        &vec_asimov_data );
    tree->Branch( "vec_pseudo_data",        &vec_pseudo_data );
    tree->Branch( "vec_bestfit_pred",       &vec_bestfit_pred );
    tree->Branch( "vec_syst_rand",          &vec_syst_rand );

    double dm2_low  = 0;
    double dm2_hgh  = 10;
    double dm2_step = 0.2;
    int num_dm2 = (dm2_hgh-dm2_low)/dm2_step;
    num_dm2 = num_dm2;
    
    double t14_low  = 0;
    double t14_hgh  = 0.5;
    double t14_step = 0.01;
    int num_t14 = (t14_hgh-t14_low)/t14_step;
    num_t14 = num_t14;
    
    double t24_low  = 0;
    double t24_hgh  = 0.5;
    double t24_step = 0.0025;
    int num_t24 = (t24_hgh-t24_low)/t24_step;
    num_t24 = num_t24;

    
    for(int itoy=1; itoy<=1; itoy++) {

      cout<<endl<<" ---> proceesing toy "<<itoy<<endl;
      
      osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, 0);  
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT();// meas, CV, COV: all ready
      //osc_test->Set_asimov2fitdata();	   
      osc_test->Set_toy_variations( 1 );
      osc_test->Set_toy2fitdata( 1 );
      
      flag_negative = 0;
      
      double obj_dm2  = 0;
      double obj_t14  = 0;
      double obj_t24  = 0;
      double obj_chi2 = 1e8;
      
      double pars_4v[4] = {0};      
      chi2_null = osc_test->FCN( pars_4v );;


          
      /////// X: sin22t14, 1e-4 -> 1   ---> "log10()" ---> -4 -> 0
      /////// Y: m41^2,    1e-1 -> 1e2  ---> "log10()" ---> -1 -> 2
      int bins_theta = 60;
      int bins_dm2   = 60;      
      TH2D *h2_space = new TH2D("h2_space_whole", "h2_space_whole", bins_theta, -4, 0, bins_dm2, -1, 2);
      h2_space->SetBinContent(1,1,1);
      
      for(int jdm2=1; jdm2<=bins_theta; jdm2++ ) {
	//for(int jdm2=1; jdm2<=num_dm2; jdm2++ ) {
	cout<<TString::Format(" ---> jdm2 %3d", jdm2)<<endl;
	
        for(int jt14=1; jt14<=bins_theta; jt14++) {
	  //for(int jt14=1; jt14<=num_t14; jt14++) {
          cout<<TString::Format(" ---> jdm2 %3d, jt14 %3d", jdm2, jt14)<<endl;
	  
          for(int jt24=1; jt24<=bins_theta; jt24++) {
	    //for(int jt24=ifile; jt24<=ifile; jt24++) {
	    
            double dm2_val = dm2_low + jdm2*dm2_step;
            double t14_val = t14_low + jt14*t14_step;
            double t24_val = t24_low + jt24*t24_step;

	    dm2_val = pow(10, h2_space->GetYaxis()->GetBinCenter(jdm2));
	    t14_val = pow(10, h2_space->GetXaxis()->GetBinCenter(jt14));
	    t24_val = pow(10, h2_space->GetXaxis()->GetBinCenter(jt24));
	    
            pars_4v[0] = dm2_val;
            pars_4v[1] = t14_val;
            pars_4v[2] = t24_val;

	    //auto time_start = chrono::high_resolution_clock::now();
	    //auto time_stop = chrono::high_resolution_clock::now();
	    //auto time_duration = chrono::duration_cast<chrono::seconds>(time_stop - time_start);
	    //cout<<endl<<" ---> check time duration "<<time_duration.count()<<endl<<endl;
	    //milliseconds, minutes
	    
            double chi2_val = osc_test->FCN( pars_4v );// 1000 times ~ 3 min
	    
	    
            if( chi2_val<0 ) flag_negative++;
            
            if( (obj_chi2>chi2_val) && (chi2_val>=0) ) {
              obj_chi2 = chi2_val;
              obj_dm2 = dm2_val;
              obj_t14 = t14_val;
              obj_t24 = t24_val;
            }
	    
          }// for(int jt24=0; jt24<bins_theta; jt24++)
        }// for(int jt14=0; jt14<bins_theta; jt14++)
      }// for(int jdm2=0; jdm2<bins_dm2; jdm2++ )
      
      //osc_test->Minimization_OscPars_FullCov(obj_dm2, obj_t14, obj_t24, 0, "str_flag_fixpar");

      cout<<" ---> obj "<<obj_chi2<<"\t"<<obj_dm2<<"\t"<<obj_t14<<"\t"<<obj_t24<<endl;
      osc_test->Minimization_OscPars_FullCov(obj_dm2, obj_t14, obj_t24, 0, "str");
      
      //osc_test->Minimization_OscPars_FullCov(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, 0, "t24");
            
      ////////

      min_status            = osc_test->minimization_status;
      min_chi2              = osc_test->minimization_chi2;
      min_dm2_41_val        = osc_test->minimization_dm2_41_val;
      min_sin2_2theta_14_val= osc_test->minimization_sin2_2theta_14_val;
      min_sin2_theta_24_val = osc_test->minimization_sin2_theta_24_val;
      min_sin2_theta_34_val = osc_test->minimization_sin2_theta_34_val;
      min_dm2_41_err        = osc_test->minimization_dm2_41_err;
      min_sin2_2theta_14_err= osc_test->minimization_sin2_2theta_14_err;
      min_sin2_theta_24_err = osc_test->minimization_sin2_theta_24_err;
      min_sin2_theta_34_err = osc_test->minimization_sin2_theta_34_err;
         
      vec_pseudo_data.clear();
      vec_bestfit_pred.clear();
      vec_syst_rand.clear();    
      for(int idx=0; idx<osc_test->matrix_eff_newworld_pred.GetNcols(); idx++) {
	vec_pseudo_data.push_back( osc_test->map_matrix_toy_pred[1](0, idx) );
	vec_bestfit_pred.push_back( osc_test->matrix_eff_newworld_pred(0, idx) );
	//if( (int)(osc_test->vec_rand_syst.size())!=0 ) vec_syst_rand.push_back( osc_test->vec_rand_syst.at(idx) );
      }
        
      tree->Fill();
    }


    tree->Write();
    subroofile->Close();
    cout<<" xyzabc"<<endl;
  }
    
  /////////////////////////////////////////////////////////// Profiling at each grid

  if( 1 ) {
    cout<<endl<<" ---> Profiling at each grid"<<endl;
 
    int bins_theta = 60;
    int bins_dm2   = 60;
  
    /////// X: sin22t14, 1e-4 -> 1   ---> "log10()" ---> -4 -> 0
    /////// Y: m41^2,    1e-1 -> 1e2  ---> "log10()" ---> -1 -> 2
    TH2D *h2_space = new TH2D("h2_space_whole", "h2_space_whole", bins_theta, -4, 0, bins_dm2, -1, 2);

    ///////
    val_dm2_41         = 7.3;
    val_sin2_2theta_14 = 0.1;
    val_sin2_theta_24  = 0.01;  
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready  
    osc_test->Set_meas2fitdata();
    
    double pars_4v[4] = {0};

    ///////
    
    int    min_status             = 10;
    int    usr_idm2               = 0;
    int    usr_it14               = 0;
    double min_chi2               = 0;
    double min_dm2_41_val         = 0;
    double min_sin2_2theta_14_val = 0;
    double min_sin2_theta_24_val  = 0;
    double min_sin2_theta_34_val  = 0;
    double min_dm2_41_err         = 0;
    double min_sin2_2theta_14_err = 0;
    double min_sin2_theta_24_err  = 0;
    double min_sin2_theta_34_err  = 0;
    
    roostr = TString::Format("sub_fit_%02d_%02d.root", idm2, it14);
    TFile *subroofile = new TFile(roostr, "recreate");
    TTree *tree = new TTree("tree", "tree");
    
    tree->Branch( "min_status",             &min_status,             "min_status/I" );
    tree->Branch( "usr_idm2",               &usr_idm2,               "usr_idm2/I" );
    tree->Branch( "usr_it14",               &usr_it14,               "usr_it14/I" );    
    tree->Branch( "min_chi2",               &min_chi2,               "min_chi2/D" );
    tree->Branch( "min_dm2_41_val",         &min_dm2_41_val,         "min_dm2_41_val/D" );
    tree->Branch( "min_sin2_2theta_14_val", &min_sin2_2theta_14_val, "min_sin2_2theta_14_val/D" );
    tree->Branch( "min_sin2_theta_24_val",  &min_sin2_theta_24_val,  "min_sin2_theta_24_val/D" );
    tree->Branch( "min_sin2_theta_34_val",  &min_sin2_theta_34_val,  "min_sin2_theta_34_val/D" );
    tree->Branch( "min_dm2_41_err",         &min_dm2_41_err,         "min_dm2_41_err/D" );
    tree->Branch( "min_sin2_2theta_14_err", &min_sin2_2theta_14_err, "min_sin2_2theta_14_err/D" );
    tree->Branch( "min_sin2_theta_24_err",  &min_sin2_theta_24_err,  "min_sin2_theta_24_err/D" );
    tree->Branch( "min_sin2_theta_34_err",  &min_sin2_theta_34_err,  "min_sin2_theta_34_err/D" );
    
    ///////

    double array_t24[300] = {0};
    for(int ii=0; ii<=100; ii++) {
      array_t24[ii] = ii*0.0001;
    }
    for(int ii=1; ii<200; ii++) {
      array_t24[ii+100] = ii*0.005;
    }
    
    for(int jdm2=idm2; jdm2<=idm2; jdm2++) {      
      for(int jt14=it14; jt14<=it14; jt14++) {
	cout<<" ------> jdm2 "<<jdm2<<" jt14 "<<jt14<<endl;

	double obj_chi2 = 1e10;
	double obj_dm2 = pow(10, h2_space->GetYaxis()->GetBinCenter(jdm2));
	double obj_t14 = pow(10, h2_space->GetXaxis()->GetBinCenter(jt14));
	double obj_t24 = 0;
	
	for(int idx=0; idx<300; idx++) {
	  if(idx%10==0) cout<<"       ---> processing t24idx "<<idx<<endl;
	  
	  //double val_t24 = idx * 0.001;
	  double val_t24 = array_t24[idx];
	  double val_chi2 = 0;
	  
	  if( val_t24<=1 ) {	    
	    pars_4v[0] = obj_dm2;
	    pars_4v[1] = obj_t14;	
	    pars_4v[2] = val_t24;
	    val_chi2 = osc_test->FCN( pars_4v );
	    
	    if( val_chi2 < obj_chi2 ) {
	      obj_chi2 = val_chi2;
	      obj_t24 = val_t24;
	    }
	  }	  
	  
	}// for(int idx=0; idx<100; jdx++)

	if( obj_t24==0 ) obj_t24 = 1e-8;
	cout<<endl<<" ---> objtesta "<<obj_t24<<endl;
	
	osc_test->Minimization_OscPars_FullCov(obj_dm2, obj_t14, obj_t24, 0, "dm2_t14");

	/////////// patch
	/////////// patch
	
	if( osc_test->minimization_status!=0 ) {
	  if( obj_t24<1e-7 ) {
	    cout<<endl<<" ---> pactch_AA_minimization "<<endl<<endl;
	    osc_test->Minimization_OscPars_FullCov(obj_dm2, obj_t14, 1e-7, 0, "dm2_t14");
	  }
	}
	if( osc_test->minimization_status!=0 ) {
	  if( obj_t24<1e-7 ) {
	    cout<<endl<<" ---> pactch_AB_minimization "<<endl<<endl;
	    osc_test->Minimization_OscPars_FullCov(obj_dm2, obj_t14, 1e-6, 0, "dm2_t14");
	  }
	}
	if( osc_test->minimization_status!=0 ) {
	  if( obj_t24<1e-7 ) {
	    cout<<endl<<" ---> pactch_AC_minimization "<<endl<<endl;
	    osc_test->Minimization_OscPars_FullCov(obj_dm2, obj_t14, 1e-5, 0, "dm2_t14");
	  }
	}


	if( osc_test->minimization_status!=0 ) {
	  cout<<endl<<" ---> pactch_BA_minimization "<<endl<<endl;
	  osc_test->Minimization_OscPars_FullCov(obj_dm2, obj_t14, osc_test->minimization_sin2_theta_24_val, 0, "dm2_t14");
	}
	if( osc_test->minimization_status!=0 ) {
	  cout<<endl<<" ---> pactch_BB_minimization "<<endl<<endl;
	  osc_test->Minimization_OscPars_FullCov(obj_dm2, obj_t14, osc_test->minimization_sin2_theta_24_val * 1.1, 0, "dm2_t14");
	}
	
	///////////
	///////////
	
	min_status            = osc_test->minimization_status;
	usr_idm2              = idm2;
	usr_it14              = it14;
	min_chi2              = osc_test->minimization_chi2;
	min_dm2_41_val        = osc_test->minimization_dm2_41_val;
	min_sin2_2theta_14_val= osc_test->minimization_sin2_2theta_14_val;
	min_sin2_theta_24_val = osc_test->minimization_sin2_theta_24_val;
	if(min_sin2_theta_24_val<1e-5) min_sin2_theta_24_val = 0;
	min_sin2_theta_34_val = osc_test->minimization_sin2_theta_34_val;
	min_dm2_41_err        = osc_test->minimization_dm2_41_err;
	min_sin2_2theta_14_err= osc_test->minimization_sin2_2theta_14_err;
	min_sin2_theta_24_err = osc_test->minimization_sin2_theta_24_err;
	min_sin2_theta_34_err = osc_test->minimization_sin2_theta_34_err;
	tree->Fill();
	
      }// for(int jt14=1; jt14<=bins_theta; jt14++)
    }// for(int jdm2=1; jdm2<=bins_dm2; jdm2++)

    tree->Write();
    subroofile->Close();
  }
   
  /////////////////////////////////////////////////////////// FC CLs
  
  if( 0 ) {
    
    cout<<endl;
    cout<<" ---> FC CLs"<<endl;
    
    int bins_theta = 60;
    int bins_dm2   = 60;
  
    /////// X: sin22t14, 1e-4 -> 1   ---> "log10()" ---> -4 -> 0
    /////// Y: m41^2,    1e-1 -> 1e2  ---> "log10()" ---> -1 -> 2
    TH2D *h2_space = new TH2D("h2_space_whole", "h2_space_whole", bins_theta, -4, 0, bins_dm2, -1, 2);

    double dm2_val = pow(10, h2_space->GetYaxis()->GetBinCenter(idm2));
    double t14_val = pow(10, h2_space->GetXaxis()->GetBinCenter(it14));

    double pars_4v[4] = {dm2_val, t14_val, 0, 0};
    double pars_3v[4] = {0};

    ///////////////////////////////////////////////////////////
    
    if( 0 ) {
      cout<<" ---> Asimov/Pseudo CLs "<<endl;

      vector<int> vec_dm2;
      vector<int> vec_t14;      
      vector<double> vec_chi2_4v;
      vector<double> vec_chi2_3v;
      
      roostr = TString::Format("sub_sm_CLs_%06d.root", ifile);
      TFile *subroofile = new TFile(roostr, "recreate");
      TTree *tree = new TTree("tree", "tree");
      
      tree->Branch( "vec_dm2",  &vec_dm2 );
      tree->Branch( "vec_t14",  &vec_t14 );
      tree->Branch( "vec_chi2_4v",  &vec_chi2_4v );
      tree->Branch( "vec_chi2_3v",  &vec_chi2_3v );

      //////// 3v
      osc_test->Set_oscillation_pars(0, 0, 0, 0);  
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT();// meas, CV, COV: all ready
      
      //osc_test->Set_meas2fitdata();
      osc_test->Set_asimov2fitdata();      
      //osc_test->Set_toy_variations( 1 );
      //osc_test->Set_toy2fitdata( 1 );

      for(int jdm2=1; jdm2<=bins_dm2; jdm2++) {
	cout<<" ---> processing jdm2 "<<jdm2<<endl;
	
	for(int jt14=1; jt14<=bins_theta; jt14++) {
	  //cout<<" ---> processing tt "<<jt14<<endl;
	  
	  double dm2 = pow(10, h2_space->GetYaxis()->GetBinCenter(jdm2));
	  double t14 = pow(10, h2_space->GetXaxis()->GetBinCenter(jt14));
	  pars_4v[0] = dm2;
	  pars_4v[1] = t14;

	  double chi2_4v = osc_test->FCN( pars_4v );
	  double chi2_3v = osc_test->FCN( pars_3v );

	  vec_dm2.push_back(jdm2);
	  vec_t14.push_back(jt14);
	  vec_chi2_4v.push_back( chi2_4v );
	  vec_chi2_3v.push_back( chi2_3v );
	}
      }
      
      tree->Fill();

      
      tree->Write();
      subroofile->Close();
	
      return 0;
    }
    
    /////////////////////////////////////////////////////

    
    int ntoys = 2000;

    
    vector<double> vec_chi2_4v_on_4vPseudo;
    vector<double> vec_chi2_3v_on_4vPseudo;      
    vector<double> vec_chi2_4v_on_3vPseudo;
    vector<double> vec_chi2_3v_on_3vPseudo;
      
    roostr = TString::Format("sub_frequentist_CLs_dm2_%02d_tt_%02d.root", idm2, it14);
    TFile *subroofile = new TFile(roostr, "recreate");
    TTree *tree = new TTree("tree", "tree");
    
    tree->Branch( "idm2",          &idm2,     "idm2/I" );
    tree->Branch( "it14",          &it14,     "it14/I" );
    tree->Branch( "vec_chi2_4v_on_4vPseudo",  &vec_chi2_4v_on_4vPseudo );
    tree->Branch( "vec_chi2_3v_on_4vPseudo",  &vec_chi2_3v_on_4vPseudo );
    tree->Branch( "vec_chi2_4v_on_3vPseudo",  &vec_chi2_4v_on_3vPseudo );
    tree->Branch( "vec_chi2_3v_on_3vPseudo",  &vec_chi2_3v_on_3vPseudo );

    //////// 4v
    osc_test->Set_oscillation_pars(dm2_val, t14_val, 0, 0);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    TMatrixD matrix_4v_COV  = osc_test->matrix_eff_newworld_abs_syst_total;
    TMatrixD matrix_4v_pred = osc_test->matrix_eff_newworld_pred;
    
    //////// 3v
    osc_test->Set_oscillation_pars(0, 0, 0, 0);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    TMatrixD matrix_3v_COV  = osc_test->matrix_eff_newworld_abs_syst_total;
    TMatrixD matrix_3v_pred = osc_test->matrix_eff_newworld_pred;
  
    ////////
    
    for(int itoy=1; itoy<=ntoys; itoy++) {

      cout<<" ---> processing itoy "<<itoy<<endl;
      
      /////////////////////// 4v pseudo
    
      // osc_test->Set_oscillation_pars(dm2_val, t14_val, 0, 0);  
      // osc_test->Apply_oscillation();
      // osc_test->Set_apply_POT();// meas, CV, COV: all ready
      osc_test->matrix_eff_newworld_abs_syst_total = matrix_4v_COV;
      osc_test->matrix_eff_newworld_pred = matrix_4v_pred;
      osc_test->Set_toy_variations( 1 );
      osc_test->Set_toy2fitdata( 1 );
      double chi2_4v_on_4vPseudo = osc_test->FCN( pars_4v );
      double chi2_3v_on_4vPseudo = osc_test->FCN( pars_3v );

      vec_chi2_4v_on_4vPseudo.push_back(chi2_4v_on_4vPseudo);
      vec_chi2_3v_on_4vPseudo.push_back(chi2_3v_on_4vPseudo);
      
      /////////////////////// 3v pseudo
    
      // osc_test->Set_oscillation_pars(0, 0, 0, 0);  
      // osc_test->Apply_oscillation();
      // osc_test->Set_apply_POT();// meas, CV, COV: all ready
      osc_test->matrix_eff_newworld_abs_syst_total = matrix_3v_COV;
      osc_test->matrix_eff_newworld_pred = matrix_3v_pred;
      osc_test->Set_toy_variations( 1 );
      osc_test->Set_toy2fitdata( 1 );      
      double chi2_4v_on_3vPseudo = osc_test->FCN( pars_4v );
      double chi2_3v_on_3vPseudo = osc_test->FCN( pars_3v );

      vec_chi2_4v_on_3vPseudo.push_back(chi2_4v_on_3vPseudo);
      vec_chi2_3v_on_3vPseudo.push_back(chi2_3v_on_3vPseudo);

      ////////////////////////

      //cout<<"  ---> "<<chi2_4v_on_3vPseudo - chi2_3v_on_3vPseudo <<endl;
      
    }
    
    tree->Fill();
    tree->Write();
    subroofile->Close();

    delete h2_space;
    
  }
  
  /////////////////////////////////////////////////////////// exclusion

  if( 0 ) {
    
    cout<<endl;
    cout<<" ---> Exclusion processing"<<endl;

    osc_test->Set_oscillation_pars(0, 0, 0, 0);
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();
    osc_test->Set_asimov2noosc();

    ///////
  
    int bins_theta = 60;
    int bins_dm2   = 60;
  
    /////// X: sin22t14, 1e-4 -> 1   ---> "log10()" ---> -4 -> 0
    /////// Y: m41^2,    1e-1 -> 1e2  ---> "log10()" ---> -1 -> 2
    TH2D *h2_space = new TH2D("h2_space_whole", "h2_space_whole", bins_theta, -4, 0, bins_dm2, -1, 2);

    int low_theta = 1;
    int hgh_theta = bins_theta;
    int low_dm2   = 1;
    int hgh_dm2   = bins_dm2;
    int low_t24   = 1;
    int hgh_t24   = bins_theta;
    
    if( it14!=0 ) {
      low_theta = it14; hgh_theta = it14;
    }
    if( idm2!=0 ) {
      low_dm2 = idm2; hgh_dm2 = idm2;
    }
    if( it24!=0) {
      low_t24 = it24; hgh_t24 = it24;
    }
    
    //double user_sin2_theta_24 = 0;// sscan
	
    for(int ibin=low_theta; ibin<=hgh_theta; ibin++) {
      for(int jbin=low_dm2; jbin<=hgh_dm2; jbin++) {
	
	roostr = TString::Format("outfile_theta_%03d_dm2_%03d_t24.txt", ibin, jbin);
	ofstream outfile(roostr, ios::out|ios::trunc);
	
	for(int kbin=low_t24; kbin<=hgh_t24; kbin++) {
	  cout<<TString::Format(" ---> processing theta,dm2,t24: %3d %3d %3d", ibin, jbin, kbin)<<endl;    
	  
	  double xcenter = h2_space->GetXaxis()->GetBinCenter(ibin);// angle 14
	  double ycenter = h2_space->GetYaxis()->GetBinCenter(jbin);// dm2
	  double kcenter = h2_space->GetXaxis()->GetBinCenter(kbin);// angle 24
	  
	  double grid_sin2_2theta_14 = pow( 10, xcenter );
	  double grid_dm2_41         = pow( 10, ycenter );
	  double grid_sin2_theta_24  = pow( 10, kcenter );
	  //grid_sin2_theta_24 = 0.0126;

	  double chi2_4v_on_4vAsimov(0), chi2_3v_on_4vAsimov(0);
	  double chi2_4v_on_3vAsimov(0), chi2_3v_on_3vAsimov(0);
	  double chi2_4v_on_data(0), chi2_3v_on_data(0);

	  double pars_4v[4] = {grid_dm2_41, grid_sin2_2theta_14, grid_sin2_theta_24, 0};
	  double pars_3v[4] = {0};
      
	  /////// 4v Asimov      
	  osc_test->Set_oscillation_pars(pars_4v[0], pars_4v[1], pars_4v[2], pars_4v[3]);
	  osc_test->Apply_oscillation();
	  osc_test->Set_apply_POT();
	  osc_test->Set_asimov2fitdata();
	  chi2_3v_on_4vAsimov = osc_test->FCN( pars_3v );

	  /////// 3v Asimov
	  osc_test->Set_noosc2fitdata();
	  chi2_4v_on_3vAsimov = osc_test->FCN( pars_4v );

	  /////// data
	  osc_test->Set_meas2fitdata();
	  chi2_4v_on_data = osc_test->FCN( pars_4v );
	  chi2_3v_on_data = osc_test->FCN( pars_3v );

	  ///////
	  double dchi2_4vAsimov = chi2_4v_on_4vAsimov - chi2_3v_on_4vAsimov;
	  double dchi2_3vAsimov = chi2_4v_on_3vAsimov - chi2_3v_on_3vAsimov;
	  double dchi2_data     = chi2_4v_on_data - chi2_3v_on_data;
      
	  double delta_4v = dchi2_4vAsimov;
	  double delta_3v = dchi2_3vAsimov;
	  double delta_dd = dchi2_data;
      
	  double data_CL = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_dd) * 100.;
	  double pred_CL = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v) * 100.;
	  double pred_CL_1sigma_plus  = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v-(2*sqrt(fabs(delta_3v))) ) * 100.;
	  double pred_CL_1sigma_minus = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v+(2*sqrt(fabs(delta_3v))) ) * 100.;
	  double pred_CL_2sigma_plus  = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v-2*(2*sqrt(fabs(delta_3v))) ) * 100.;
	  double pred_CL_2sigma_minus = 100 - osc_test->func_CLs(delta_4v, delta_3v, delta_3v+2*(2*sqrt(fabs(delta_3v))) ) * 100.;

	  outfile<<TString::Format("%3d %3d %3d %16.9f %16.9f %16.9f %16.9f %16.9f %16.9f %14.9f %14.9f %14.9f %14.9f %14.9f %14.9f",
				   ibin, jbin, kbin,
				   chi2_4v_on_4vAsimov, chi2_3v_on_4vAsimov,
				   chi2_4v_on_3vAsimov, chi2_3v_on_3vAsimov,
				   chi2_4v_on_data, chi2_3v_on_data,
				   data_CL, pred_CL,
				   pred_CL_1sigma_plus, pred_CL_1sigma_minus,
				   pred_CL_2sigma_plus, pred_CL_2sigma_minus
				   )<<endl;	  

	  //cout<<TString::Format(" ---> pred p-value %10.8f", (100-pred_CL)/100)<<endl;

	}// for(int kbin=low_theta; kbin<=hgh_theta; kbin++)

	outfile.close();
	
      }// for(int jbin=1; jbin<=bins_dm2; jbin++)
    }// for(int ibin=1; ibin<=bins_theta; ibin++)
   
      
  }


  ///////////////////////////////////////////////////////////

  cout<<endl;
  cout<<" ------------------------------ check at the final step ------------------------------"<<endl;
  cout<<" ------------------------------ check at the final step ------------------------------"<<endl;

  cout<<endl;
  cout<<TString::Format(" ---> display(-d) %d, ifile(-f) %d, scaleF_POT_BNB(-pbnb) %6.4f, scaleF_POT_NuMI(-pnumi) %6.4f, theta14(-it14) %d, dm2(-idm2) %d, theta24(-it24) %d",
			display, ifile, osc_test->scaleF_POT_BNB, osc_test->scaleF_POT_NuMI, it14, idm2, it24)<<endl;
  
  cout<<endl;
  cout<<TString::Format(" ---> flag_syst_dirt    %d", osc_test->flag_syst_dirt)<<endl;
  cout<<TString::Format(" ---> flag_syst_mcstat  %d", osc_test->flag_syst_mcstat)<<endl;
  cout<<TString::Format(" ---> flag_syst_flux    %d", osc_test->flag_syst_flux)<<endl;
  cout<<TString::Format(" ---> flag_syst_geant   %d", osc_test->flag_syst_geant)<<endl;
  cout<<TString::Format(" ---> flag_syst_Xs      %d", osc_test->flag_syst_Xs)<<endl;
  cout<<TString::Format(" ---> flag_syst_det     %d", osc_test->flag_syst_det)<<endl;

  cout<<endl;
  cout<<" ---> Finished sucessfully"<<endl;
  
  cout<<endl;
  if( display ) {
    cout<<" Enter Ctrl+c to end the program"<<endl;
    cout<<" Enter Ctrl+c to end the program"<<endl;
    cout<<endl;
    theApp.Run();
  }
  
  return 0;
}
