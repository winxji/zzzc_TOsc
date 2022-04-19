void ec_plot_CLs_GaussCLs_edit()
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
  
  ////////////////////////////////////////////////////////////////////////////////////////

  TString roostr_CLs = "";
  TString roostr_GaussCLs = "";

  roostr_CLs = "outfile_CLs_numuDIS_95.root";
  roostr_GaussCLs = "za_roofile_GaussCLs_numuDIS.root";  
  
  ////////////////

  TString xx_title = "";
  const int ee = 1;
  const int uu = 2;
  const int ue = 3;
  
  int xx = uu;
  
  switch( xx ) {
  case ee:
    xx_title = "sin^{2}2#theta_{ee}"; break;
  case uu:
    xx_title = "sin^{2}2#theta_{#mu#mu}"; break;
  case ue:
    xx_title = "sin^{2}2#theta_{#mue}"; break;
  default:
    cout<<endl<<" ERROR xx_title "<<endl<<endl;
  }
  
  ////////////////

  TFile *file_CLs = new TFile(roostr_CLs, "read");
  TGraph *gh_sensitivity_median = (TGraph*)file_CLs->Get("gh_sensitivity_median");
  TGraph *gh_sensitivity_1s_band = (TGraph*)file_CLs->Get("gh_sensitivity_1s_band");
  TGraph *gh_sensitivity_2s_band = (TGraph*)file_CLs->Get("gh_sensitivity_2s_band");
  TGraph *gh_data_result = (TGraph*)file_CLs->Get("gh_data_result");

  TFile *file_GaussCLs = new TFile(roostr_GaussCLs, "read");
  TGraph *gh_gauss_pred = (TGraph*)file_GaussCLs->Get("gh_CLs_pred_1");
  TGraph *gh_gauss_data = (TGraph*)file_GaussCLs->Get("gh_CLs_data_1");
  TGraph *gh_wilks_pred = (TGraph*)file_GaussCLs->Get("gh_wilk_CL_pred_1");

  TFile *file_CLs_Asimov = new TFile("outfile_CLs_numuDIS_3vAsimovAsMeas.root", "read");
  TGraph *gh_CLs_Asimov = (TGraph*)file_CLs_Asimov->Get("gh_data_result");
  
  ////////////////

  roostr = "canv_results";
  TCanvas *canv_results = new TCanvas(roostr, roostr, 800, 650);
  canv_results->SetBottomMargin(0.15); canv_results->SetTopMargin(0.1); canv_results->SetLeftMargin(0.15); canv_results->SetRightMargin(0.1);
  canv_results->SetLogx();
  canv_results->SetLogy();

  roostr = "h2d_results";
  TH2D *h2d_results = new TH2D(roostr, "", 10, 1e-2, 1, 10, 0.1, 100);
  h2d_results->Draw();
  h2d_results->SetXTitle(xx_title);
  h2d_results->SetYTitle("#Deltam^{2}_{41} (eV^{2})");
  
  h2d_results->GetXaxis()->CenterTitle(); h2d_results->GetXaxis()->SetTitleSize(0.05); h2d_results->GetXaxis()->SetLabelSize(0.05);
  h2d_results->GetYaxis()->CenterTitle(); h2d_results->GetYaxis()->SetTitleSize(0.05); h2d_results->GetYaxis()->SetLabelSize(0.05);  
  h2d_results->GetXaxis()->SetTitleOffset(1.4);
  h2d_results->GetYaxis()->SetTitleOffset(1.2);

  gh_sensitivity_2s_band->Draw("same f");
  gh_sensitivity_1s_band->Draw("same f");
  gh_sensitivity_median->Draw("same l"); gh_sensitivity_median->SetLineColor(kBlue); gh_sensitivity_median->SetLineStyle(7);

  gh_gauss_pred->Draw("same l"); gh_gauss_pred->SetLineColor(kRed); gh_gauss_pred->SetLineWidth(3); gh_gauss_pred->SetLineStyle(7);
  //gh_wilks_pred->Draw("same l"); gh_wilks_pred->SetLineColor(kBlue); gh_wilks_pred->SetLineWidth(3);

  gh_data_result->Draw("same l"); gh_data_result->SetLineColor(kBlack); gh_data_result->SetLineWidth(3);
  gh_gauss_data->Draw("same l"); gh_gauss_data->SetLineColor(kRed); gh_gauss_data->SetLineWidth(3);

  gh_CLs_Asimov->Draw("same l"); gh_CLs_Asimov->SetLineColor(kBlack); gh_CLs_Asimov->SetLineWidth(3); gh_CLs_Asimov->SetLineStyle(7);
  
  TLegend *lg = new TLegend(0.16, 0.20-0.03, 0.44, 0.45-0.03);
  lg->SetTextSize(0.042);
  
  if(xx==ue) {
    lg->SetX1(0.44); lg->SetX2(0.74);
    lg->SetY1(0.6); lg->SetY2(0.85);
    lg->SetTextSize(0.042);
  }
  
  lg->AddEntry( "", "95% CL", "");
  lg->AddEntry( gh_sensitivity_median, "Frequentist CLs: sensitivity", "l");
  //lg->AddEntry( gh_sensitivity_1s_band, "1 sigma", "f");
  //lg->AddEntry( gh_sensitivity_2s_band, "2 sigma", "f");  
  lg->AddEntry( gh_data_result, "Frequentist CLs: data", "l");
  lg->AddEntry( gh_CLs_Asimov, "Frequentist CLs: Asimov", "l");
  lg->AddEntry( gh_gauss_data, "Gaussian CLs: data", "l");
  lg->AddEntry( gh_gauss_pred, "Gaussian CLs + Asimov", "l");  
  //lg->AddEntry( gh_wilks_pred, "Wilks' theorem + Asimov: sensitivity", "l");

  lg->Draw();
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  


  canv_results->SaveAs("canv_results.png");
}
