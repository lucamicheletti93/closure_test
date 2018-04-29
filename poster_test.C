#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>

#include <TROOT.h>
#include <TMinuit.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <string.h>

#include "settings.h"
#include "binning.C"
#endif

double Func_W(double *, double *);
double Func_cost(double *, double *);
double Func_phi(double *, double *);

void generate_polarized_distribution();

const int N_test = 27;
double lambdaTh[N_test] = {-1.,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.,0.2,0.2,0.2,0.2,-0.2,-0.2,-0.2,-0.2,0.6,0.6,0.6,0.6,-0.6,-0.6,-0.6,-0.6};
double lambdaPhi[N_test] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-0.2,-0.1,0.1,0.2,-0.2,-0.1,0.1,0.2,-0.2,-0.1,0.1,0.2,-0.2,-0.1,0.1,0.2};
double lambda_ThPhi[N_test] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

//const int N_test = 11;
//double lambdaTh[N_test] = {-1.,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.};
//double lambdaPhi[N_test] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//double lambda_ThPhi[N_test] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

//const int N_test = 19;
//double lambdaTh[N_test] = {-1.,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.,0.2,0.2,0.2,0.2,-0.2,-0.2,-0.2,-0.2};
//double lambdaPhi[N_test] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-0.2,-0.1,0.1,0.2,-0.2,-0.1,0.1,0.2};
//double lambda_ThPhi[N_test] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

void poster_test(string sample = "LowStat"){

  //============================================================================
  printf("1) Setting main quantities ... \n");
  //============================================================================

  gSystem -> CompileMacro("../settings.h");
  gROOT -> ProcessLine(".x ../binning.C");

  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(1);
  //TGaxis::SetMaxDigits(2);

  char INPUT_FILE_NAME[300];
  if(sample == "LowStat") sprintf(INPUT_FILE_NAME,"/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/POLARIZED_DISTRIBUTIONS/sum_variable_Jpsi_polarization_LowStat.root");
  if(sample == "FullStat") sprintf(INPUT_FILE_NAME,"/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/POLARIZED_DISTRIBUTIONS/sum_variable_Jpsi_polarization_FullStat.root");

  char INPUT_TREE_NAME[100];
  sprintf(INPUT_TREE_NAME,"MCTree");

  double min_fit_range_Cost = -0.5;
  double max_fit_range_Cost = 0.5;
  double min_fit_range_Phi = 0.95;
  double max_fit_range_Phi = 2.2;

  /*double min_fit_range_Cost = -0.6;
  double max_fit_range_Cost = 0.6;
  double min_fit_range_Phi = 0;
  double max_fit_range_Phi = PI;*/

  //============================================================================
  printf("2) Reading the the file.root in ../OUTPUT ... \n");
  //============================================================================

  TFile *input_file = new TFile(INPUT_FILE_NAME,"READ");
  TH2D *hist_CostPhiHE_2pt6_VAR_pol_gen[N_test];
  TH2D *hist_CostPhiHE_2pt6_VAR_pol_rec[N_test], *hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[N_test];
  char hist_name[60];

  for(int i = 0;i < N_test;i++){
    sprintf(hist_name,"hist_gen_polarization%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_gen[i] = (TH2D*) input_file -> Get(hist_name);
    sprintf(hist_name,"hist_rec_polarization%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_rec[i] = (TH2D*) input_file -> Get(hist_name);
    sprintf(hist_name,"hist_rec_polarization_Rebin%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[i] = (TH2D*) input_file -> Get(hist_name);
  }

  //============================================================================
  printf("3) Normalizing the Rebin histo to bin area ... \n"); // AN = Area Normalized
  //============================================================================
  TH2D *hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_AN[N_test];

  //printf("%f +- %f \n",hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[16] -> GetBinContent(2,2),hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[16] -> GetBinError(2,2));

  for(int i = 0;i < N_test;i++){
    sprintf(hist_name,"hist_polarization_area_normalized%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_AN[i] = new TH2D(hist_name,hist_name,N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
    for(int j = 0;j < N_cost_bins_BC;j++){
      for(int k = 0;k < N_phi_bins_BC;k++){
        hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_AN[i] -> SetBinContent(j+1,k+1,(hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[i] -> GetBinContent(j+1,k+1))/bin_area_BC[j][k]);
        hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_AN[i] -> SetBinError(j+1,k+1,(hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[i] -> GetBinError(j+1,k+1))/bin_area_BC[j][k]);
      }
    }
  }

  //printf("%f +- %f \n",hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_AN[16] -> GetBinContent(2,2),hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_AN[16] -> GetBinError(2,2));

  //============================================================================
  printf("4) Correcting REC for Acc x Eff ... \n"); // AC = Acceptance Corrected
  //============================================================================

  sprintf(INPUT_FILE_NAME,"../INPUT/Histo_accxeff_FromOfficialTree_Jpsi_PbPb_Nopol.root");
  TFile *accxeff_file = new TFile(INPUT_FILE_NAME,"READ");

  TH2D *hist_accxeff_HE_2pt6_NO_pol = (TH2D*) accxeff_file -> Get("hist_accxeff_HE_2pt6_phibent");
  TH2D *hist_accxeff_HE_2pt6_NO_pol_Rebin = (TH2D*) accxeff_file -> Get("hist_accxeff_HE_2pt6");

  TH2D *hist_CostPhiHE_2pt6_VAR_pol_rec_AC[N_test], *hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_ANAC[N_test];

  for(int k = 0;k < N_test;k++){
    sprintf(hist_name,"hist_polarization_accxeff_NOpol%i",k);
    hist_CostPhiHE_2pt6_VAR_pol_rec_AC[k] = new TH2D(hist_name,hist_name,100,-1,1,50,0,PI);
    hist_CostPhiHE_2pt6_VAR_pol_rec_AC[k] -> Sumw2();
    hist_CostPhiHE_2pt6_VAR_pol_rec_AC[k] -> Divide(hist_CostPhiHE_2pt6_VAR_pol_rec[k],hist_accxeff_HE_2pt6_NO_pol,1,1);
    //hist_CostPhiHE_2pt6_VAR_pol_rec_AC[k] -> Rebin2D(5,5);

    sprintf(hist_name,"hist_polarization_accxeff_NOpol_Rebin%i",k);
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_ANAC[k] = new TH2D(hist_name,hist_name,N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_ANAC[k] -> Sumw2();
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_ANAC[k] -> Divide(hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_AN[k],hist_accxeff_HE_2pt6_NO_pol_Rebin,1,1);
    //hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_ANAC[k] -> Scale(1/10.);
  }

  //============================================================================
  printf("5) Fit of the REC corrected distributions rebinned with NO pol Acc x Eff correction ... \n");
  //============================================================================
  char func_name[60];
  TGraphErrors *gra_lambdaTh_Phi_Theor = new TGraphErrors(N_test,lambdaTh,lambdaPhi,0,0);
  gra_lambdaTh_Phi_Theor -> SetMarkerColor(kBlack);
  gra_lambdaTh_Phi_Theor -> SetMarkerStyle(21);
  gra_lambdaTh_Phi_Theor -> SetMarkerSize(1.2);
  //============================================================================
  // FIT REC -> NARROW BINNING
  //============================================================================
  TF2 *fit_func_CostPhi2D_HE_VARNO_pol_nrw[N_test]; // REC DISTRIB CORRECTED FOR NO-POL ACCXEFF [NARROW]
  double lambdaTh_rec_nrw[N_test], stat_lambdaTh_rec_nrw[N_test], lambdaPhi_rec_nrw[N_test], stat_lambdaPhi_rec_nrw[N_test];

  for(int i = 0;i < N_test;i++){
    sprintf(func_name,"fit_func_polarization_VAR_pol%i",i);
    fit_func_CostPhi2D_HE_VARNO_pol_nrw[i] = new TF2(func_name,Func_W,min_fit_range_Cost,max_fit_range_Cost,min_fit_range_Phi,max_fit_range_Phi,4);
    fit_func_CostPhi2D_HE_VARNO_pol_nrw[i] -> SetParameters(10000,lambdaTh[i],lambdaPhi[i],0);
    hist_CostPhiHE_2pt6_VAR_pol_rec_AC[i] -> Fit(fit_func_CostPhi2D_HE_VARNO_pol_nrw[i],"RSI0");

    lambdaTh_rec_nrw[i] = fit_func_CostPhi2D_HE_VARNO_pol_nrw[i] -> GetParameter(1);
    stat_lambdaTh_rec_nrw[i] = fit_func_CostPhi2D_HE_VARNO_pol_nrw[i] -> GetParError(1);
    lambdaPhi_rec_nrw[i] = fit_func_CostPhi2D_HE_VARNO_pol_nrw[i] -> GetParameter(2);
    stat_lambdaPhi_rec_nrw[i] = fit_func_CostPhi2D_HE_VARNO_pol_nrw[i] -> GetParError(2);
  }

  TGraphErrors *gra_lambdaTh_Phi_rec_nrw_Filled = new TGraphErrors(N_test,lambdaTh_rec_nrw,lambdaPhi_rec_nrw,stat_lambdaTh_rec_nrw,stat_lambdaPhi_rec_nrw);
  gra_lambdaTh_Phi_rec_nrw_Filled -> SetLineColor(kBlue+1);
  gra_lambdaTh_Phi_rec_nrw_Filled -> SetLineWidth(1);
  gra_lambdaTh_Phi_rec_nrw_Filled -> SetFillColor(kBlue+1);
  gra_lambdaTh_Phi_rec_nrw_Filled -> SetFillStyle(3004);

  TGraphErrors *gra_lambdaTh_Phi_rec_nrw = new TGraphErrors(N_test,lambdaTh_rec_nrw,lambdaPhi_rec_nrw,stat_lambdaTh_rec_nrw,stat_lambdaPhi_rec_nrw);
  gra_lambdaTh_Phi_rec_nrw -> SetMarkerColor(kBlue+1);
  gra_lambdaTh_Phi_rec_nrw -> SetMarkerStyle(20);
  gra_lambdaTh_Phi_rec_nrw -> SetLineColor(kBlue+1);
  gra_lambdaTh_Phi_rec_nrw -> SetLineWidth(1);
  gra_lambdaTh_Phi_rec_nrw -> SetFillColor(kBlue+1);
  gra_lambdaTh_Phi_rec_nrw -> SetFillStyle(3004);

  //============================================================================
  // FIT REC -> LARGE BINNING
  //============================================================================
  TF2 *fit_func_CostPhi2D_HE_VARNO_pol_lrg[N_test]; // REC DISTRIB CORRECTED FOR NO-POL ACCXEFF [NARROW]
  double lambdaTh_rec_lrg[N_test], stat_lambdaTh_rec_lrg[N_test], lambdaPhi_rec_lrg[N_test], stat_lambdaPhi_rec_lrg[N_test];

  for(int i = 0;i < N_test;i++){
    sprintf(func_name,"fit_func_polarization_VAR_pol%i",i);
    fit_func_CostPhi2D_HE_VARNO_pol_lrg[i] = new TF2(func_name,Func_W,min_fit_range_Cost,max_fit_range_Cost,min_fit_range_Phi,max_fit_range_Phi,4);
    fit_func_CostPhi2D_HE_VARNO_pol_lrg[i] -> SetParameters(10000,lambdaTh[i],lambdaPhi[i],0);
    ///fit_func_CostPhi2D_HE_VARNO_pol_lrg[i] -> SetParLimits(1,-1.2,1.2);
    //fit_func_CostPhi2D_HE_VARNO_pol_lrg[i] -> SetParLimits(2,-0.4,0.4);
    fit_func_CostPhi2D_HE_VARNO_pol_lrg[i] -> SetParLimits(1,lambdaTh[i]-0.4,lambdaTh[i]+0.4);
    fit_func_CostPhi2D_HE_VARNO_pol_lrg[i] -> SetParLimits(2,lambdaPhi[i]-0.4,lambdaPhi[i]+0.4);
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_ANAC[i] -> Fit(fit_func_CostPhi2D_HE_VARNO_pol_lrg[i],"RSI0");

    lambdaTh_rec_lrg[i] = fit_func_CostPhi2D_HE_VARNO_pol_lrg[i] -> GetParameter(1);
    stat_lambdaTh_rec_lrg[i] = fit_func_CostPhi2D_HE_VARNO_pol_lrg[i] -> GetParError(1);
    lambdaPhi_rec_lrg[i] = fit_func_CostPhi2D_HE_VARNO_pol_lrg[i] -> GetParameter(2);
    stat_lambdaPhi_rec_lrg[i] = fit_func_CostPhi2D_HE_VARNO_pol_lrg[i] -> GetParError(2);
  }

  TGraphErrors *gra_lambdaTh_Phi_rec_lrg_Filled = new TGraphErrors(N_test,lambdaTh_rec_lrg,lambdaPhi_rec_lrg,stat_lambdaTh_rec_lrg,stat_lambdaPhi_rec_lrg);
  gra_lambdaTh_Phi_rec_lrg_Filled -> SetLineColor(kRed+1);
  gra_lambdaTh_Phi_rec_lrg_Filled -> SetLineWidth(1);
  gra_lambdaTh_Phi_rec_lrg_Filled -> SetFillColor(kRed+1);
  gra_lambdaTh_Phi_rec_lrg_Filled -> SetFillStyle(3005);

  TGraphErrors *gra_lambdaTh_Phi_rec_lrg = new TGraphErrors(N_test,lambdaTh_rec_lrg,lambdaPhi_rec_lrg,stat_lambdaTh_rec_lrg,stat_lambdaPhi_rec_lrg);
  gra_lambdaTh_Phi_rec_lrg -> SetMarkerColor(kRed+1);
  gra_lambdaTh_Phi_rec_lrg -> SetMarkerStyle(20);
  gra_lambdaTh_Phi_rec_lrg -> SetLineColor(kRed+1);
  gra_lambdaTh_Phi_rec_lrg -> SetLineWidth(2);
  //gra_lambdaTh_Phi_rec_lrg -> SetFillColor(kRed+1);
  //gra_lambdaTh_Phi_rec_lrg -> SetFillStyle(0);

  //============================================================================
  printf("6) MIGROS plor ... \n");
  //============================================================================

  TH2D *h_lambdaTh_Phi = new TH2D("h_lambdaTh_Phi","",100,-1.2,1.2,100,-0.4,0.4);
  h_lambdaTh_Phi -> GetXaxis() -> SetNdivisions(12);
  h_lambdaTh_Phi -> GetXaxis() -> SetTitle("#lambda_{#theta}");
  h_lambdaTh_Phi -> GetYaxis() -> SetNdivisions(8);
  h_lambdaTh_Phi -> GetYaxis() -> SetTitle("#lambda_{#phi}");

  TLegend *leg_Th_Phi = new TLegend(0.3,0.80,0.6,0.94,"","brNDC");
  leg_Th_Phi -> SetTextFont(42);
  leg_Th_Phi -> SetTextSize(0.025);
  leg_Th_Phi -> AddEntry(gra_lambdaTh_Phi_Theor,"(#lambda_{#theta},#lambda_{#phi}) #rightarrow input value","P");
  if(sample == "FullStat") leg_Th_Phi -> AddEntry(gra_lambdaTh_Phi_rec_nrw,"(#lambda_{#theta},#lambda_{#phi}) #rightarrow rec value [narrow binning]","PF");
  leg_Th_Phi -> AddEntry(gra_lambdaTh_Phi_rec_lrg,"(#lambda_{#theta},#lambda_{#phi}) #rightarrow rec value [large binning]","PL");

  TCanvas *c_lambdaTh_Phi = new TCanvas("c_lambdaTh_Phi","c_lambdaTh_Phi",4,132,1024,768);
  c_lambdaTh_Phi -> SetGrid();
  h_lambdaTh_Phi -> Draw();
  gra_lambdaTh_Phi_Theor -> Draw("sameP");
  if(sample == "FullStat"){
    //gra_lambdaTh_Phi_rec_nrw_Filled -> Draw("sameP");
    gra_lambdaTh_Phi_rec_nrw -> Draw("samePE5");
  }
  gra_lambdaTh_Phi_rec_lrg_Filled -> Draw("sameP");
  gra_lambdaTh_Phi_rec_lrg -> Draw("sameP");
  leg_Th_Phi -> Draw("same");

  return;

  //============================================================================
  printf("7) Control plots ... \n");
  //============================================================================

  TCanvas *c_grid_plot = new TCanvas("c_grid_plot","c_grid_plot",20,20,600,600);

  TH2D *h_grid_plot = new TH2D("h_grid_plot","",100,-1,1,50,0,PI);
  h_grid_plot -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  h_grid_plot -> GetYaxis() -> SetTitle("#phi_{HE}");
  h_grid_plot -> Draw();
  for(int i = 0;i < N_line_cost_BC;i++){
    if(i < N_line_phi_BC) line_phi_BC[i] -> Draw("same");
    line_cost_BC[i] -> Draw("same");
  }
  TLine *line_cost_min = new TLine(min_fit_range_Cost,min_fit_range_Phi,min_fit_range_Cost,max_fit_range_Phi);
  line_cost_min -> SetLineColor(kRed);
  line_cost_min -> SetLineWidth(2);
  line_cost_min -> Draw("same");
  TLine *line_cost_max = new TLine(max_fit_range_Cost,min_fit_range_Phi,max_fit_range_Cost,max_fit_range_Phi);
  line_cost_max -> SetLineColor(kRed);
  line_cost_max -> SetLineWidth(2);
  line_cost_max -> Draw("same");
  TLine *line_phi_min = new TLine(min_fit_range_Cost,min_fit_range_Phi,max_fit_range_Cost,min_fit_range_Phi);
  line_phi_min -> SetLineColor(kRed);
  line_phi_min -> SetLineWidth(2);
  line_phi_min -> Draw("same");
  TLine *line_phi_max = new TLine(min_fit_range_Cost,max_fit_range_Phi,max_fit_range_Cost,max_fit_range_Phi);
  line_phi_max -> SetLineColor(kRed);
  line_phi_max -> SetLineWidth(2);
  line_phi_max -> Draw("same");

  //============================================================================
  printf("Control plots ... \n");
  //============================================================================

  TH2D *h_lambdaTh = new TH2D("h_lambdaTh","",N_test,0,N_test,100,-1.2,1.2);
  h_lambdaTh -> GetXaxis() -> SetTitle("#lambda_{#theta}");
  TH1D *hist_lambdaTh = new TH1D("hist_lambdaTh","hist_lambdaTh",N_test,0,N_test);
  hist_lambdaTh -> SetLineColor(kBlack);
  TH1D *hist_lambdaTh_rec_nrw = new TH1D("hist_lambdaTh_rec_nrw","hist_lambdaTh_rec_nrw",N_test,0,N_test);
  hist_lambdaTh_rec_nrw -> SetMarkerColor(kBlue+1);
  hist_lambdaTh_rec_nrw -> SetMarkerStyle(20);
  hist_lambdaTh_rec_nrw -> SetLineColor(kBlue+1);
  hist_lambdaTh_rec_nrw -> SetLineWidth(2);
  TH1D *hist_lambdaTh_rec_lrg = new TH1D("hist_lambdaTh_rec_lrg","hist_lambdaTh_rec_lrg",N_test,0,N_test);
  hist_lambdaTh_rec_lrg -> SetMarkerColor(kRed+1);
  hist_lambdaTh_rec_lrg -> SetMarkerStyle(20);
  hist_lambdaTh_rec_lrg -> SetLineColor(kRed+1);
  hist_lambdaTh_rec_lrg -> SetLineWidth(2);

  TH2D *h_lambdaPhi = new TH2D("h_lambdaPhi","",N_test,0,N_test,100,-0.6,0.6);
  h_lambdaPhi -> GetXaxis() -> SetTitle("#lambda_{#phi}");
  TH1D *hist_lambdaPhi = new TH1D("hist_lambdaPhi","hist_lambdaPhi",N_test,0,N_test);
  hist_lambdaPhi -> SetLineColor(kBlack);
  TH1D *hist_lambdaPhi_rec_nrw = new TH1D("hist_lambdaPhi_rec_nrw","hist_lambdaPhi_rec_nrw",N_test,0,N_test);
  hist_lambdaPhi_rec_nrw -> SetMarkerColor(kBlue+1);
  hist_lambdaPhi_rec_nrw -> SetMarkerStyle(20);
  hist_lambdaPhi_rec_nrw -> SetLineColor(kBlue+1);
  hist_lambdaPhi_rec_nrw -> SetLineWidth(2);
  TH1D *hist_lambdaPhi_rec_lrg = new TH1D("hist_lambdaPhi_rec_lrg","hist_lambdaPhi_rec_lrg",N_test,0,N_test);
  hist_lambdaPhi_rec_lrg -> SetMarkerColor(kRed+1);
  hist_lambdaPhi_rec_lrg -> SetMarkerStyle(20);
  hist_lambdaPhi_rec_lrg -> SetLineColor(kRed+1);
  hist_lambdaPhi_rec_lrg -> SetLineWidth(2);

  char bin_label[40];

  for(int i = 0;i < N_test;i++){
    sprintf(bin_label,"%2.1f",lambdaTh[i]);
    h_lambdaTh -> GetXaxis() -> SetBinLabel(i+1,bin_label);
    sprintf(bin_label,"%2.1f",lambdaPhi[i]);
    h_lambdaPhi -> GetXaxis() -> SetBinLabel(i+1,bin_label);

    hist_lambdaTh -> SetBinContent(i+1,lambdaTh[i]);
    hist_lambdaPhi -> SetBinContent(i+1,lambdaPhi[i]);
    hist_lambdaTh_rec_nrw -> SetBinContent(i+1,lambdaTh_rec_nrw[i]); hist_lambdaTh_rec_nrw -> SetBinError(i+1,stat_lambdaTh_rec_nrw[i]);
    hist_lambdaTh_rec_lrg -> SetBinContent(i+1,lambdaTh_rec_lrg[i]); hist_lambdaTh_rec_lrg -> SetBinError(i+1,stat_lambdaTh_rec_lrg[i]);
    hist_lambdaPhi_rec_nrw -> SetBinContent(i+1,lambdaPhi_rec_nrw[i]); hist_lambdaPhi_rec_nrw -> SetBinError(i+1,stat_lambdaPhi_rec_nrw[i]);
    hist_lambdaPhi_rec_lrg -> SetBinContent(i+1,lambdaPhi_rec_lrg[i]); hist_lambdaPhi_rec_lrg -> SetBinError(i+1,stat_lambdaPhi_rec_lrg[i]);
  }

  TCanvas *c_lambdaTh = new TCanvas("c_lambdaTh","c_lambdaTh",4,132,1024,768);
  h_lambdaTh -> Draw();
  hist_lambdaTh -> Draw("same");
  hist_lambdaTh_rec_nrw -> Draw("sameP");
  hist_lambdaTh_rec_lrg -> Draw("sameP");

  TCanvas *c_lambdaPhi = new TCanvas("c_lambdaPhi","c_lambdaPhi",4,132,1024,768);
  h_lambdaPhi -> Draw();
  hist_lambdaPhi -> Draw("same");
  hist_lambdaPhi_rec_nrw -> Draw("sameP");
  hist_lambdaPhi_rec_lrg -> Draw("sameP");

}
//______________________________________________________________________________
void generate_polarized_distribution(){

  //============================================================================
  printf("1) Setting main quantities ... \n");
  //============================================================================

  gSystem -> CompileMacro("../settings.h");
  gROOT -> ProcessLine(".x ../binning.C");

  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(1);

  char INPUT_FILE_NAME[300];
  //sprintf(INPUT_FILE_NAME,"../INPUT/sample_tree_for_closure_test.root");
  sprintf(INPUT_FILE_NAME,"/home/luca/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root");
  char INPUT_TREE_NAME[100];
  sprintf(INPUT_TREE_NAME,"MCTree");

  //============================================================================
  printf("2) Inizializing the tree ... \n");
  //============================================================================

  TFile *input_file = new TFile(INPUT_FILE_NAME,"READ");

  Int_t NDimu_gen;
  Double_t DimuPt_gen[3000], DimuY_gen[3000];
  Double_t CostHE_gen[3000], PhiHE_gen[3000], CostCS_gen[3000], PhiCS_gen[3000];

  TTree *input_tree = (TTree*) input_file -> Get(INPUT_TREE_NAME);
  input_tree -> SetBranchAddress("NDimu_gen",&NDimu_gen);
  input_tree -> SetBranchAddress("DimuPt_gen",DimuPt_gen);
  input_tree -> SetBranchAddress("DimuY_gen",DimuY_gen);
  input_tree -> SetBranchAddress("CostHE_gen",CostHE_gen);
  input_tree -> SetBranchAddress("PhiHE_gen",PhiHE_gen);
  input_tree -> SetBranchAddress("CostCS_gen",CostCS_gen);
  input_tree -> SetBranchAddress("PhiCS_gen",PhiCS_gen);

  Int_t NDimu_rec;
  Double_t DimuPt_rec[3000], DimuY_rec[3000];
  Double_t DimuMass_rec[3000];
  Int_t DimuMatch_rec[3000];
  Double_t CostHE_rec[3000], PhiHE_rec[3000], CostCS_rec[3000], PhiCS_rec[3000];

  input_tree -> SetBranchAddress("NDimu_rec",&NDimu_rec);
  input_tree -> SetBranchAddress("DimuPt_rec",DimuPt_rec);
  input_tree -> SetBranchAddress("DimuY_rec",DimuY_rec);
  input_tree -> SetBranchAddress("DimuMass_rec",DimuMass_rec);
  input_tree -> SetBranchAddress("DimuMatch_rec",DimuMatch_rec);
  input_tree -> SetBranchAddress("CostHE_rec",CostHE_rec);
  input_tree -> SetBranchAddress("PhiHE_rec",PhiHE_rec);
  input_tree -> SetBranchAddress("CostCS_rec",CostCS_rec);
  input_tree -> SetBranchAddress("PhiCS_rec",PhiCS_rec);

  //============================================================================
  printf("3) Filling the histo ... \n"); // INCLUDE FUNCTION TO GENERATE POLARIZED DISTRIBUTIONS
  //============================================================================

  double weight_CostPhi2D = 0;
  int count_Dimu_rec = 0;
  int event_index = 0;
  char func_name[60];
  char hist_name[60];

  TF2 *func_CostPhi2D_HE_VAR_pol[N_test];

  TH2D *hist_CostPhiHE_2pt6_VAR_pol_rec[N_test];
  TH2D *hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[N_test];
  TH2D *hist_CostPhiHE_2pt6_VAR_pol_gen[N_test];
  TH2D *hist_CostPhiHE_2pt6_VAR_pol_gen_Rebin[N_test];

  for(int i = 0;i < N_test;i++){
    //if(i < 11) continue;
    if(i < 19) continue;
    sprintf(func_name,"func_polarization%i",i);
    func_CostPhi2D_HE_VAR_pol[i] = new TF2(func_name,Func_W,-1,1,0,PI,4);
    func_CostPhi2D_HE_VAR_pol[i] -> SetParameters(1000,lambdaTh[i],lambdaPhi[i],lambda_ThPhi[i]);
    sprintf(hist_name,"hist_rec_polarization%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_rec[i] = new TH2D(hist_name,hist_name,100,-1,1,50,0,PI);
    sprintf(hist_name,"hist_rec_polarization_Rebin%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[i] = new TH2D(hist_name,hist_name,N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
    sprintf(hist_name,"hist_gen_polarization%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_gen[i] = new TH2D(hist_name,hist_name,100,-1,1,50,0,PI);
    sprintf(hist_name,"hist_gen_polarization_Rebin%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_gen_Rebin[i] = new TH2D(hist_name,hist_name,N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);


    while(count_Dimu_rec < 1000000){
      input_tree -> GetEntry(event_index);

      for(int k = 0;k < NDimu_gen;k++){
        if(DimuPt_gen[k] > 2 && DimuPt_gen[k] <= 6){
          weight_CostPhi2D = func_CostPhi2D_HE_VAR_pol[i] -> Eval(CostHE_gen[k],TMath::Abs(PhiHE_gen[k]))/func_CostPhi2D_HE_VAR_pol[i] -> GetMaximum();
          hist_CostPhiHE_2pt6_VAR_pol_gen[i] -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k]),weight_CostPhi2D);
          hist_CostPhiHE_2pt6_VAR_pol_gen_Rebin[i] -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k]),weight_CostPhi2D);
        }
      }

      for(int k = 0;k < NDimu_rec;k++){
        if(DimuPt_rec[k] > 2 && DimuPt_rec[k] <= 6){
          if(DimuY_rec[k] > -4. && DimuY_rec[k] < -2.5){
            if(DimuMatch_rec[k] == 2){
              if(DimuMass_rec[k] > 2 && DimuMass_rec[k] < 5){
                //weight_CostPhi2D = func_CostPhi2D_HE_VAR_pol[i] -> Eval(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]))/func_CostPhi2D_HE_VAR_pol[i] -> GetMaximum();
                hist_CostPhiHE_2pt6_VAR_pol_rec[i] -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),weight_CostPhi2D);
                hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[i] -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),weight_CostPhi2D);
                count_Dimu_rec++;
              }
            }
          }
        }
      }
      event_index++;
    }
    printf("POLARIZATION = (%2.1f,%2.1f,%2.1f) \n",lambdaTh[i],lambdaPhi[i],lambda_ThPhi[i]);
    event_index = 0;
    count_Dimu_rec = 0;
  }

  char OUTPUT_FILE_NAME[300];
  sprintf(OUTPUT_FILE_NAME,"/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/POLARIZED_DISTRIBUTIONS/variable_Jpsi_polarization_sample3_FullStat.root");
  TFile *output_file = new TFile(OUTPUT_FILE_NAME,"RECREATE");
  for(int i = 0;i < N_test;i++){
    //if(i < 11) continue;
    if(i < 19) continue;
    hist_CostPhiHE_2pt6_VAR_pol_gen[i] -> Write();
    hist_CostPhiHE_2pt6_VAR_pol_gen_Rebin[i] -> Write();
    hist_CostPhiHE_2pt6_VAR_pol_rec[i] -> Write();
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[i] -> Write();
  }
  output_file -> Close();
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
double Func_W(double *x, double *par){

  double N = par[0];
  double L_th = par[1];
  double L_ph = par[2];
  double L_thph = par[3];

  double costh = x[0];
  double phi = x[1];

  double cosph = TMath::Cos(phi);
  double cos2ph = TMath::Cos(2*phi);

  double W =  (N/(3 + L_th))*(1 + (L_th + L_ph*cos2ph)*costh*costh + 2*L_thph*costh*cosph*TMath::Sqrt(1 - costh*costh) + L_ph*cos2ph);
  return W;
}
//------------------------------------------------------------------------------
double Func_cost(double *x, double *par){

  double N = par[0];
  double L_th = par[1];

  double costh = x[0];

  double W =  (N/(3 + L_th))*(1 + L_th*costh*costh);
  return W;
}
//------------------------------------------------------------------------------
double Func_phi(double *x, double *par){

  double N = par[0];
  double L_th = par[1];
  double L_ph = par[2];

  double phi = x[0];

  double cos2ph = TMath::Cos(2*phi);

  double W =  N*(1 + ((2*L_ph)/(3 + L_th))*cos2ph);
  return W;
}
