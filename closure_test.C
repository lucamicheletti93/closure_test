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
#endif

void closure_test(){

  gSystem -> CompileMacro("../settings.h");
  gROOT -> ProcessLine(".x ../binning.C");

  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(1);
  TGaxis::SetMaxDigits(2);

  //============================================================================
  printf("1) Reading N J/psi histograms ... \n");
  //============================================================================

  TFile *file = new TFile("../OUTPUT/Histo_for_closure_test.root","READ");

  TH2D *hist_CostPhiHE_2pt6_TR_pol_gen = (TH2D*) file -> Get("hist_CostPhiHE_2pt6_TR_pol_gen");
  TH2D *hist_CostPhiHE_2pt6_TR_pol_rec = (TH2D*) file -> Get("hist_CostPhiHE_2pt6_TR_pol_rec");
  TH2D *hist_accxeff_2pt6_TR_pol = (TH2D*) file -> Get("hist_accxeff_2pt6_TR_pol");
  TH2D *hist_CostPhiHE_2pt6_TR_pol_gen_Rebin = (TH2D*) file -> Get("hist_CostPhiHE_2pt6_TR_pol_gen_Rebin");
  TH2D *hist_CostPhiHE_2pt6_TR_pol_rec_Rebin = (TH2D*) file -> Get("hist_CostPhiHE_2pt6_TR_pol_rec_Rebin");
  TH2D *hist_accxeff_2pt6_TR_pol_Rebin = (TH2D*) file -> Get("hist_accxeff_2pt6_TR_pol_Rebin");

  TH2D *hist_CostPhiHE_2pt6_LG_pol_gen = (TH2D*) file -> Get("hist_CostPhiHE_2pt6_LG_pol_gen");
  TH2D *hist_CostPhiHE_2pt6_LG_pol_rec = (TH2D*) file -> Get("hist_CostPhiHE_2pt6_LG_pol_rec");
  TH2D *hist_accxeff_2pt6_LG_pol = (TH2D*) file -> Get("hist_accxeff_2pt6_LG_pol");
  TH2D *hist_CostPhiHE_2pt6_LG_pol_gen_Rebin = (TH2D*) file -> Get("hist_CostPhiHE_2pt6_LG_pol_gen_Rebin");
  TH2D *hist_CostPhiHE_2pt6_LG_pol_rec_Rebin = (TH2D*) file -> Get("hist_CostPhiHE_2pt6_LG_pol_rec_Rebin");
  TH2D *hist_accxeff_2pt6_LG_pol_Rebin = (TH2D*) file -> Get("hist_accxeff_2pt6_LG_pol_Rebin");

  TH2D *hist_CostPhiHE_2pt6_NO_pol_gen = (TH2D*) file -> Get("hist_CostPhiHE_2pt6_NO_pol_gen");
  TH2D *hist_CostPhiHE_2pt6_NO_pol_rec = (TH2D*) file -> Get("hist_CostPhiHE_2pt6_NO_pol_rec");
  TH2D *hist_accxeff_2pt6_NO_pol = (TH2D*) file -> Get("hist_accxeff_2pt6_NO_pol");
  TH2D *hist_CostPhiHE_2pt6_NO_pol_gen_Rebin = (TH2D*) file -> Get("hist_CostPhiHE_2pt6_NO_pol_gen_Rebin");
  TH2D *hist_CostPhiHE_2pt6_NO_pol_rec_Rebin = (TH2D*) file -> Get("hist_CostPhiHE_2pt6_NO_pol_rec_Rebin");
  TH2D *hist_accxeff_2pt6_NO_pol_Rebin = (TH2D*) file -> Get("hist_accxeff_2pt6_NO_pol_Rebin");

  //============================================================================
  printf("2) Normalizing the Rebin histo to bin area ... \n"); // AN = Area Normalized
  //============================================================================

  TH2D *hist_CostPhiHE_2pt6_TR_pol_gen_Rebin_AN = new TH2D("hist_CostPhiHE_2pt6_TR_pol_gen_Rebin_AN","hist_CostPhiHE_2pt6_TR_pol_gen_Rebin_AN",N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
  TH2D *hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_AN = new TH2D("hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_AN","hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_AN",N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
  TH2D *hist_CostPhiHE_2pt6_LG_pol_gen_Rebin_AN = new TH2D("hist_CostPhiHE_2pt6_LG_pol_gen_Rebin_AN","hist_CostPhiHE_2pt6_LG_pol_gen_Rebin_AN",N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
  TH2D *hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_AN = new TH2D("hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_AN","hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_AN",N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
  TH2D *hist_CostPhiHE_2pt6_NO_pol_gen_Rebin_AN = new TH2D("hist_CostPhiHE_2pt6_NO_pol_gen_Rebin_AN","hist_CostPhiHE_2pt6_NO_pol_gen_Rebin_AN",N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
  TH2D *hist_CostPhiHE_2pt6_NO_pol_rec_Rebin_AN = new TH2D("hist_CostPhiHE_2pt6_NO_pol_rec_Rebin_AN","hist_CostPhiHE_2pt6_NO_pol_rec_Rebin_AN",N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);

  for(int i = 0;i< N_cost_bins_BC;i++){
    for(int j = 0;j < N_phi_bins_BC;j++){
      hist_CostPhiHE_2pt6_TR_pol_gen_Rebin_AN -> SetBinContent(i+1,j+1,(hist_CostPhiHE_2pt6_TR_pol_gen_Rebin -> GetBinContent(i+1,j+1))/bin_area_BC[i][j]);
      hist_CostPhiHE_2pt6_TR_pol_gen_Rebin_AN -> SetBinError(i+1,j+1,(hist_CostPhiHE_2pt6_TR_pol_gen_Rebin -> GetBinError(i+1,j+1))/bin_area_BC[i][j]);
      hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_AN -> SetBinContent(i+1,j+1,(hist_CostPhiHE_2pt6_TR_pol_rec_Rebin -> GetBinContent(i+1,j+1))/bin_area_BC[i][j]);
      hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_AN -> SetBinError(i+1,j+1,(hist_CostPhiHE_2pt6_TR_pol_rec_Rebin -> GetBinError(i+1,j+1))/bin_area_BC[i][j]);
      hist_CostPhiHE_2pt6_LG_pol_gen_Rebin_AN -> SetBinContent(i+1,j+1,(hist_CostPhiHE_2pt6_LG_pol_gen_Rebin -> GetBinContent(i+1,j+1))/bin_area_BC[i][j]);
      hist_CostPhiHE_2pt6_LG_pol_gen_Rebin_AN -> SetBinError(i+1,j+1,(hist_CostPhiHE_2pt6_LG_pol_gen_Rebin -> GetBinError(i+1,j+1))/bin_area_BC[i][j]);
      hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_AN -> SetBinContent(i+1,j+1,(hist_CostPhiHE_2pt6_LG_pol_rec_Rebin -> GetBinContent(i+1,j+1))/bin_area_BC[i][j]);
      hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_AN -> SetBinError(i+1,j+1,(hist_CostPhiHE_2pt6_LG_pol_rec_Rebin -> GetBinError(i+1,j+1))/bin_area_BC[i][j]);
      hist_CostPhiHE_2pt6_NO_pol_gen_Rebin_AN -> SetBinContent(i+1,j+1,(hist_CostPhiHE_2pt6_NO_pol_gen_Rebin -> GetBinContent(i+1,j+1))/bin_area_BC[i][j]);
      hist_CostPhiHE_2pt6_NO_pol_gen_Rebin_AN -> SetBinError(i+1,j+1,(hist_CostPhiHE_2pt6_NO_pol_gen_Rebin -> GetBinError(i+1,j+1))/bin_area_BC[i][j]);
      hist_CostPhiHE_2pt6_NO_pol_rec_Rebin_AN -> SetBinContent(i+1,j+1,(hist_CostPhiHE_2pt6_NO_pol_rec_Rebin -> GetBinContent(i+1,j+1))/bin_area_BC[i][j]);
      hist_CostPhiHE_2pt6_NO_pol_rec_Rebin_AN -> SetBinError(i+1,j+1,(hist_CostPhiHE_2pt6_NO_pol_rec_Rebin -> GetBinError(i+1,j+1))/bin_area_BC[i][j]);
    }
  }

  //============================================================================
  printf("3) Correcting REC for Acc x Eff ... \n"); // AC = Acceptance Corrected // NP = No Polarization (Acc x Eff)
  //============================================================================

  TH2D *clone_hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_AN = (TH2D*) hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_AN -> Clone("clone_hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_AN");
  TH2D *clone_hist_accxeff_2pt6_TR_pol_Rebin = (TH2D*) hist_accxeff_2pt6_TR_pol_Rebin -> Clone("clone_hist_accxeff_2pt6_TR_pol_Rebin");
  TH2D *hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_ANAC = new TH2D("hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_ANAC","hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_ANAC",N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
  hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_ANAC -> Divide(clone_hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_AN,clone_hist_accxeff_2pt6_TR_pol_Rebin,1,1);

  TH2D *clone_hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_AN = (TH2D*) hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_AN -> Clone("clone_hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_AN");
  TH2D *clone_hist_accxeff_2pt6_LG_pol_Rebin = (TH2D*) hist_accxeff_2pt6_LG_pol_Rebin -> Clone("clone_hist_accxeff_2pt6_LG_pol_Rebin");
  TH2D *hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_ANAC = new TH2D("hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_ANAC","hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_ANAC",N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
  hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_ANAC -> Divide(clone_hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_AN,clone_hist_accxeff_2pt6_LG_pol_Rebin,1,1);

  TH2D *clone_hist_CostPhiHE_2pt6_NO_pol_rec_Rebin_AN = (TH2D*) hist_CostPhiHE_2pt6_NO_pol_rec_Rebin_AN -> Clone("clone_hist_CostPhiHE_2pt6_NO_pol_rec_Rebin_AN");
  TH2D *clone_hist_accxeff_2pt6_NO_pol_Rebin = (TH2D*) hist_accxeff_2pt6_NO_pol_Rebin -> Clone("clone_hist_accxeff_2pt6_NO_pol_Rebin");
  TH2D *hist_CostPhiHE_2pt6_NO_pol_rec_Rebin_ANAC = new TH2D("hist_CostPhiHE_2pt6_NO_pol_rec_Rebin_ANAC","hist_CostPhiHE_2pt6_NO_pol_rec_Rebin_ANAC",N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
  hist_CostPhiHE_2pt6_NO_pol_rec_Rebin_ANAC -> Divide(clone_hist_CostPhiHE_2pt6_NO_pol_rec_Rebin_AN,clone_hist_accxeff_2pt6_NO_pol_Rebin,1,1);

  // MIXED cases

  TH2D *hist_CostPhiHE_2pt6_TRNO_pol_rec_AC = new TH2D("hist_CostPhiHE_2pt6_TRNO_pol_rec_AC","hist_CostPhiHE_2pt6_TRNO_pol_rec_AC",100,-1,1,50,0,PI);
  hist_CostPhiHE_2pt6_TRNO_pol_rec_AC -> Divide(hist_CostPhiHE_2pt6_TR_pol_rec,hist_accxeff_2pt6_NO_pol,1,1);

  TH2D *hist_CostPhiHE_2pt6_TRNO_pol_rec_Rebin_ANAC = new TH2D("hist_CostPhiHE_2pt6_TRNO_pol_rec_Rebin_ANAC","hist_CostPhiHE_2pt6_TRNO_pol_rec_Rebin_ANAC",N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
  hist_CostPhiHE_2pt6_TRNO_pol_rec_Rebin_ANAC -> Divide(hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_AN,hist_accxeff_2pt6_NO_pol_Rebin,1,1);

  TH2D *hist_CostPhiHE_2pt6_LGNO_pol_rec_AC = new TH2D("hist_CostPhiHE_2pt6_LGNO_pol_rec_AC","hist_CostPhiHE_2pt6_LGNO_pol_rec_AC",100,-1,1,50,0,PI);
  hist_CostPhiHE_2pt6_LGNO_pol_rec_AC -> Divide(hist_CostPhiHE_2pt6_LG_pol_rec,hist_accxeff_2pt6_NO_pol,1,1);

  TH2D *hist_CostPhiHE_2pt6_LGNO_pol_rec_Rebin_ANAC = new TH2D("hist_CostPhiHE_2pt6_LGNO_pol_rec_Rebin_ANAC","hist_CostPhiHE_2pt6_LGNO_pol_rec_Rebin_ANAC",N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
  hist_CostPhiHE_2pt6_LGNO_pol_rec_Rebin_ANAC -> Divide(hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_AN,hist_accxeff_2pt6_NO_pol_Rebin,1,1);

  //============================================================================
  printf("4) TEST 1a : Fit of the GEN distributions rebinned ... \n");
  //============================================================================

  // TRANSVERSE POLARIZATION -> (+1,0,0)

  TF2 *func_W_HE_TR_pol_gen_Rebin_AN = new TF2("func_W_HE_TR_pol_gen_Rebin_AN",Func_W,-1,1,0,PI,4);
  func_W_HE_TR_pol_gen_Rebin_AN -> SetParameters(1000,1,0,0);
  func_W_HE_TR_pol_gen_Rebin_AN -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");
  hist_CostPhiHE_2pt6_TR_pol_gen_Rebin_AN -> Fit(func_W_HE_TR_pol_gen_Rebin_AN,"RLS0");
  printf("> > > (%4.3f,%4.3f,%4.3f) < < < \n",func_W_HE_TR_pol_gen_Rebin_AN -> GetParameter(1),func_W_HE_TR_pol_gen_Rebin_AN -> GetParameter(2),func_W_HE_TR_pol_gen_Rebin_AN -> GetParameter(3));

  // LONGITUDINAL POLARIZATION -> (+1,0,0)

  TF2 *func_W_HE_LG_pol_gen_Rebin_AN = new TF2("func_W_HE_LG_pol_gen_Rebin_AN",Func_W,-1,1,0,PI,4);
  func_W_HE_LG_pol_gen_Rebin_AN -> SetParameters(1000,-1,0,0);
  func_W_HE_LG_pol_gen_Rebin_AN -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");
  hist_CostPhiHE_2pt6_LG_pol_gen_Rebin_AN -> Fit(func_W_HE_LG_pol_gen_Rebin_AN,"RLS0");
  printf("> > > (%4.3f,%4.3f,%4.3f) < < < \n",func_W_HE_LG_pol_gen_Rebin_AN -> GetParameter(1),func_W_HE_LG_pol_gen_Rebin_AN -> GetParameter(2),func_W_HE_LG_pol_gen_Rebin_AN -> GetParameter(3));

  // NO POLARIZATION -> (0,0,0)

  TF2 *func_W_HE_NO_pol_gen_Rebin_AN = new TF2("func_W_HE_NO_pol_gen_Rebin_AN",Func_W,-1,1,0,PI,4);
  func_W_HE_NO_pol_gen_Rebin_AN -> SetParameters(1000,0,0,0);
  func_W_HE_NO_pol_gen_Rebin_AN -> SetParLimits(1,-1,1);
  func_W_HE_NO_pol_gen_Rebin_AN -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");
  hist_CostPhiHE_2pt6_NO_pol_gen_Rebin_AN -> Fit(func_W_HE_NO_pol_gen_Rebin_AN,"RLS0");
  printf("> > > (%4.3f,%4.3f,%4.3f) < < < \n",func_W_HE_NO_pol_gen_Rebin_AN -> GetParameter(1),func_W_HE_NO_pol_gen_Rebin_AN -> GetParameter(2),func_W_HE_NO_pol_gen_Rebin_AN -> GetParameter(3));

  //============================================================================
  printf("4) TEST 1b : Fit of the REC corrected distributions rebinned ... \n");
  //============================================================================

  // TRANSVERSE POLARIZATION -> (+1,0,0)

  TF2 *func_W_HE_TR_pol_rec_Rebin_ANAC = new TF2("func_W_HE_TR_pol_rec_Rebin_ANAC",Func_W,-1,1,0,PI,4);
  func_W_HE_TR_pol_rec_Rebin_ANAC -> SetParameters(1000,1,0,0);
  func_W_HE_TR_pol_rec_Rebin_ANAC -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");
  hist_CostPhiHE_2pt6_TR_pol_rec_Rebin_ANAC -> Fit(func_W_HE_TR_pol_rec_Rebin_ANAC,"RLS0");
  printf("> > > (%4.3f,%4.3f,%4.3f) < < < \n",func_W_HE_TR_pol_rec_Rebin_ANAC -> GetParameter(1),func_W_HE_TR_pol_rec_Rebin_ANAC -> GetParameter(2),func_W_HE_TR_pol_rec_Rebin_ANAC -> GetParameter(3));

  // LONGITUDINAL POLARIZATION -> (+1,0,0)

  TF2 *func_W_HE_LG_pol_rec_Rebin_ANAC = new TF2("func_W_HE_LG_pol_rec_Rebin_ANAC",Func_W,-1,1,0,PI,4);
  func_W_HE_LG_pol_rec_Rebin_ANAC -> SetParameters(1000,-1,0,0);
  func_W_HE_LG_pol_rec_Rebin_ANAC -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");
  hist_CostPhiHE_2pt6_LG_pol_rec_Rebin_ANAC -> Fit(func_W_HE_LG_pol_rec_Rebin_ANAC,"RLS0");
  printf("> > > (%4.3f,%4.3f,%4.3f) < < < \n",func_W_HE_LG_pol_rec_Rebin_ANAC -> GetParameter(1),func_W_HE_LG_pol_rec_Rebin_ANAC -> GetParameter(2),func_W_HE_LG_pol_rec_Rebin_ANAC -> GetParameter(3));

  // NO POLARIZATION -> (0,0,0)

  TF2 *func_W_HE_NO_pol_rec_Rebin_ANAC = new TF2("func_W_HE_NO_pol_rec_Rebin_ANAC",Func_W,-1,1,0,PI,4);
  func_W_HE_NO_pol_rec_Rebin_ANAC -> SetParameters(1000,-1,0,0);
  func_W_HE_NO_pol_rec_Rebin_ANAC -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");
  hist_CostPhiHE_2pt6_NO_pol_rec_Rebin_ANAC -> Fit(func_W_HE_NO_pol_rec_Rebin_ANAC,"RLS0");
  printf("> > > (%4.3f,%4.3f,%4.3f) < < < \n",func_W_HE_NO_pol_rec_Rebin_ANAC -> GetParameter(1),func_W_HE_NO_pol_rec_Rebin_ANAC -> GetParameter(2),func_W_HE_NO_pol_rec_Rebin_ANAC -> GetParameter(3));

  //============================================================================
  printf("5) TEST 2 : Fit of the REC corrected distributions rebinned with NO pol Acc x Eff correction... \n");
  //============================================================================

  // TRANSVERSE POLARIZATION -> (+1,0,0)

  TF2 *func_W_HE_TRNO_pol_rec_Rebin_ANAC = new TF2("func_W_HE_TRNO_pol_rec_Rebin_ANAC",Func_W,-0.6,0.6,0,PI,4);
  func_W_HE_TRNO_pol_rec_Rebin_ANAC -> SetParameters(1000,1,0,0);
  func_W_HE_TRNO_pol_rec_Rebin_ANAC -> SetParLimits(1,-1,1);
  func_W_HE_TRNO_pol_rec_Rebin_ANAC -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");
  hist_CostPhiHE_2pt6_TRNO_pol_rec_Rebin_ANAC -> Fit(func_W_HE_TRNO_pol_rec_Rebin_ANAC,"RLS0");
  printf("> > > (%4.3f,%4.3f,%4.3f) < < < \n",func_W_HE_TRNO_pol_rec_Rebin_ANAC -> GetParameter(1),func_W_HE_TRNO_pol_rec_Rebin_ANAC -> GetParameter(2),func_W_HE_TRNO_pol_rec_Rebin_ANAC -> GetParameter(3));

  TCanvas *c_fit_2pt6_TRNO_pol_rec_Rebin_ANAC = new TCanvas("c_fit_2pt6_TRNO_pol_rec_Rebin_ANAC","c_fit_2pt6_TRNO_pol_rec_Rebin_ANAC",4,132,1024,768);
  hist_CostPhiHE_2pt6_TRNO_pol_rec_Rebin_ANAC -> Draw("SURF1");
  func_W_HE_TRNO_pol_rec_Rebin_ANAC -> Draw("SURFsame");

  // LONGITUDINAL POLARIZATION -> (-1,0,0)

  TF2 *func_W_HE_LGNO_pol_rec_Rebin_ANAC = new TF2("func_W_HE_LGNO_pol_rec_Rebin_ANAC",Func_W,-0.6,0.6,0,PI,4);
  func_W_HE_LGNO_pol_rec_Rebin_ANAC -> SetParameters(1000,-1,0,0);
  func_W_HE_LGNO_pol_rec_Rebin_ANAC -> SetParLimits(-1,-1,1);
  func_W_HE_LGNO_pol_rec_Rebin_ANAC -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");
  hist_CostPhiHE_2pt6_LGNO_pol_rec_Rebin_ANAC -> Fit(func_W_HE_LGNO_pol_rec_Rebin_ANAC,"RLS0");
  printf("> > > (%4.3f,%4.3f,%4.3f) < < < \n",func_W_HE_LGNO_pol_rec_Rebin_ANAC -> GetParameter(1),func_W_HE_LGNO_pol_rec_Rebin_ANAC -> GetParameter(2),func_W_HE_LGNO_pol_rec_Rebin_ANAC -> GetParameter(3));

  TCanvas *c_fit_2pt6_LGNO_pol_rec_Rebin_ANAC = new TCanvas("c_fit_2pt6_LGNO_pol_rec_Rebin_ANAC","c_fit_2pt6_LGNO_pol_rec_Rebin_ANAC",4,132,1024,768);
  hist_CostPhiHE_2pt6_LGNO_pol_rec_Rebin_ANAC -> Draw("SURF1");
  func_W_HE_LGNO_pol_rec_Rebin_ANAC -> Draw("SURFsame");

  //============================================================================
  printf("5) TEST 3 : Fit of the REC corrected distributions NOT rebinned with NO pol Acc x Eff correction... \n");
  //============================================================================

  // TRANSVERSE POLARIZATION -> (+1,0,0)

  TF2 *func_W_HE_TRNO_pol_rec_AC = new TF2("func_W_HE_TRNO_pol_rec_AC",Func_W,-0.7,0.7  ,0,PI,4);
  func_W_HE_TRNO_pol_rec_AC -> SetParameters(1000,1,0,0);
  func_W_HE_TRNO_pol_rec_AC -> SetParLimits(1,0,2);
  func_W_HE_TRNO_pol_rec_AC -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");
  hist_CostPhiHE_2pt6_TRNO_pol_rec_AC -> Fit(func_W_HE_TRNO_pol_rec_AC,"RLS0");
  printf("> > > (%4.3f,%4.3f,%4.3f) < < < \n",func_W_HE_TRNO_pol_rec_AC -> GetParameter(1),func_W_HE_TRNO_pol_rec_AC -> GetParameter(2),func_W_HE_TRNO_pol_rec_AC -> GetParameter(3));

  TCanvas *c_fit_2pt6_TRNO_pol_rec_AC = new TCanvas("c_fit_2pt6_TRNO_pol_rec_AC","c_fit_2pt6_TRNO_pol_rec_AC",4,132,1024,768);
  hist_CostPhiHE_2pt6_TRNO_pol_rec_AC -> Draw("SURF1");
  func_W_HE_TRNO_pol_rec_AC -> Draw("SURFsame");

  // LONGITUDINAL POLARIZATION -> (-1,0,0)

  TF2 *func_W_HE_LGNO_pol_rec_AC = new TF2("func_W_HE_LGNO_pol_rec_AC",Func_W,-1,1,0,PI,4);
  func_W_HE_LGNO_pol_rec_AC -> SetParameters(1000,-1,0,0);
  func_W_HE_LGNO_pol_rec_AC -> SetParLimits(1,-2,0);
  func_W_HE_LGNO_pol_rec_AC -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");
  hist_CostPhiHE_2pt6_LGNO_pol_rec_AC -> Fit(func_W_HE_LGNO_pol_rec_AC,"RLS0");
  printf("> > > (%4.3f,%4.3f,%4.3f) < < < \n",func_W_HE_LGNO_pol_rec_AC -> GetParameter(1),func_W_HE_LGNO_pol_rec_AC -> GetParameter(2),func_W_HE_LGNO_pol_rec_AC -> GetParameter(3));

  TCanvas *c_fit_2pt6_LGNO_pol_rec_AC = new TCanvas("c_fit_2pt6_LGNO_pol_rec_AC","c_fit_2pt6_LGNO_pol_rec_AC",4,132,1024,768);
  hist_CostPhiHE_2pt6_LGNO_pol_rec_AC -> Draw("SURF1");
  func_W_HE_LGNO_pol_rec_AC -> Draw("SURFsame");

}
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
