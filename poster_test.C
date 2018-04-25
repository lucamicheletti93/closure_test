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

const int N_test = 21;
double lambdaTh[N_test] = {-1.,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.,0.2,0.2,0.2,0.2,0.2,-0.2,-0.2,-0.2,-0.2,-0.2};
double lambdaPhi[N_test] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-0.2,-0.1,-0.0,0.1,0.2,-0.2,-0.1,-0.0,0.1,0.2};
double lambda_ThPhi[N_test] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

double min_fit_range_Cost = -0.6;
double max_fit_range_Cost = 0.6;
double min_fit_range_Phi = 0;
double max_fit_range_Phi = PI;

void poster_test(){

  //============================================================================
  printf("1) Setting main quantities ... \n");
  //============================================================================

  gSystem -> CompileMacro("../settings.h");
  gROOT -> ProcessLine(".x ../binning.C");

  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(1);
  //TGaxis::SetMaxDigits(2);

  char INPUT_FILE_NAME[300];
  sprintf(INPUT_FILE_NAME,"/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/POLARIZED_DISTRIBUTIONS/variable_Jpsi_polarization.root");
  char INPUT_TREE_NAME[100];
  sprintf(INPUT_TREE_NAME,"MCTree");

  //============================================================================
  printf("2) Reading the the file.root in ../OUTPUT ... \n");
  //============================================================================

  TFile *input_file = new TFile(INPUT_FILE_NAME,"READ");
  TH2D *hist_CostPhiHE_2pt6_VAR_pol_gen[N_test];
  TH2D *hist_CostPhiHE_2pt6_VAR_pol_rec[N_test], *hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[N_test];
  char hist_name[60];
  char func_name[60];

  for(int i = 0;i < N_test;i++){
    sprintf(hist_name,"hist_gen_polarization%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_gen[i] = (TH2D*) input_file -> Get(hist_name);
    sprintf(hist_name,"hist_rec_polarization%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_rec[i] = (TH2D*) input_file -> Get(hist_name);
    sprintf(hist_name,"hist_rec_polarization_Rebin%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[i] = (TH2D*) input_file -> Get(hist_name);
  }

  //============================================================================
  printf("4) Normalizing the Rebin histo to bin area ... \n"); // AN = Area Normalized
  //============================================================================

  TH2D *hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_AN[N_test];

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

  //============================================================================
  printf("5) Correcting REC for Acc x Eff ... \n"); // AC = Acceptance Corrected
  //============================================================================

  sprintf(INPUT_FILE_NAME,"../INPUT/Histo_accxeff_FromOfficialTree_Jpsi_PbPb_Nopol.root");
  TFile *accxeff_file = new TFile(INPUT_FILE_NAME,"READ");

  TH2D *hist_accxeff_HE_2pt6_NO_pol = (TH2D*) accxeff_file -> Get("hist_accxeff_HE_2pt6_phibent");
  TH2D *hist_accxeff_HE_2pt6_NO_pol_Rebin = (TH2D*) accxeff_file -> Get("hist_accxeff_HE_2pt6");

  TH2D *hist_CostPhiHE_2pt6_VAR_pol_rec_AC[N_test], *hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_ANAC[N_test];

  for(int k = 0;k < N_test;k++){
    sprintf(hist_name,"hist_polarization_accxeff_NOpol%i",k);
    hist_CostPhiHE_2pt6_VAR_pol_rec_AC[k] = new TH2D(hist_name,hist_name,100,-1,1,50,0,PI);
    hist_CostPhiHE_2pt6_VAR_pol_rec_AC[k] -> Divide(hist_CostPhiHE_2pt6_VAR_pol_rec[k],hist_accxeff_HE_2pt6_NO_pol,1,1);

    sprintf(hist_name,"hist_polarization_accxeff_NOpol_Rebin%i",k);
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_ANAC[k] = new TH2D(hist_name,hist_name,N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_ANAC[k] -> Divide(hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_AN[k],hist_accxeff_HE_2pt6_NO_pol_Rebin,1,1);
  }

  //============================================================================
  printf("6) Fit of the REC corrected distributions rebinned with NO pol Acc x Eff correction ... \n");
  //============================================================================

  TF2 *fit_func_CostPhi2D_HE_VARNO_pol[N_test], *fit_func_CostPhi2D_HE_VARNO_pol_Rebin[N_test];
  TF2 *fit_func_CostPhi2D_HE_VAR_pol[N_test]; // fit of generated narrow binning
  double lambdaTh_gen[N_test], lambdaPhi_gen[N_test];
  double lambdaTh_rec[N_test], stat_lambdaTh_rec[N_test], lambdaPhi_rec[N_test], stat_lambdaPhi_rec[N_test];
  double lambdaTh_rec_Rebin[N_test], stat_lambdaTh_rec_Rebin[N_test], lambdaPhi_rec_Rebin[N_test], stat_lambdaPhi_rec_Rebin[N_test];
  double diff_lambdaTh[N_test], diff_lambdaPhi[N_test];
  char bin_label[40];

  TH2D *h_lambdaTh = new TH2D("h_lambdaTh","",N_test,0,N_test,100,-1.2,1.2);
  h_lambdaTh -> GetXaxis() -> SetTitle("#lambda_{#theta}|_{GEN}");
  h_lambdaTh -> GetYaxis() -> SetTitle("#lambda_{#theta}|_{REC}");

  TH1D *hist_lambdaTh_Theor = new TH1D("hist_lambdaTh_Theor","",N_test,0,N_test);
  hist_lambdaTh_Theor -> SetLineColor(kBlue);
  hist_lambdaTh_Theor -> SetFillColor(kBlue);
  hist_lambdaTh_Theor -> SetFillStyle(3005);

  TH1D *hist_lambdaTh = new TH1D("hist_lambdaTh","",N_test,0,N_test);
  hist_lambdaTh -> SetMarkerColor(kBlue);
  hist_lambdaTh -> SetLineColor(kBlue);
  hist_lambdaTh -> SetLineWidth(2);
  hist_lambdaTh -> SetMarkerStyle(20);

  TH1D *hist_lambdaTh_Rebin = new TH1D("hist_lambdaTh_Rebin","",N_test,0,N_test);
  hist_lambdaTh_Rebin -> SetMarkerColor(kRed);
  hist_lambdaTh_Rebin -> SetLineColor(kRed);
  hist_lambdaTh_Rebin -> SetLineWidth(2);
  hist_lambdaTh_Rebin -> SetMarkerStyle(20);

  TH1D *hist_lambdaTh_gen = new TH1D("hist_lambdaTh_gen","",N_test,0,N_test);
  hist_lambdaTh_gen -> SetLineColor(kGreen+1);
  hist_lambdaTh_gen -> SetLineWidth(2);

  TH2D *h_lambdaTh_diff = new TH2D("h_lambdaTh_diff","",N_test,0,N_test,100,-0.2,0.2);
  h_lambdaTh_diff -> GetYaxis() -> SetTitle("#lambda_{#theta}|_{GEN} - #lambda_{#theta}|_{REC}");
  h_lambdaTh_diff -> GetYaxis() -> SetTitleOffset(1.3);

  TH1D *hist_lambdaTh_diff = new TH1D("hist_lambdaTh_diff","",N_test,0,N_test);
  hist_lambdaTh_diff -> SetMarkerColor(kBlue);
  hist_lambdaTh_diff -> SetLineColor(kBlue);
  hist_lambdaTh_diff -> SetLineWidth(2);
  hist_lambdaTh_diff -> SetMarkerStyle(20);

  TH1D *hist_lambdaTh_diff_Rebin = new TH1D("hist_lambdaTh_diff_Rebin","",N_test,0,N_test);
  hist_lambdaTh_diff_Rebin -> SetMarkerColor(kRed);
  hist_lambdaTh_diff_Rebin -> SetLineColor(kRed);
  hist_lambdaTh_diff_Rebin -> SetLineWidth(2);
  hist_lambdaTh_diff_Rebin -> SetMarkerStyle(20);

  TH1D *hist_lambdaTh_gen_diff = new TH1D("hist_lambdaTh_gen_diff","",N_test,0,N_test);
  hist_lambdaTh_gen_diff -> SetLineColor(kGreen+1);
  hist_lambdaTh_gen_diff -> SetLineWidth(2);

  TH2D *h_lambdaPhi = new TH2D("h_lambdaPhi","",N_test,0,N_test,100,-1.2,1.2);
  h_lambdaPhi -> GetXaxis() -> SetTitle("#lambda_{#phi}|_{GEN}");
  h_lambdaPhi -> GetYaxis() -> SetTitle("#lambda_{#phi}|_{REC}");

  TH1D *hist_lambdaPhi_Theor = new TH1D("hist_lambdaPhi_Theor","",N_test,0,N_test);
  hist_lambdaPhi_Theor -> SetMinimum(-1);
  hist_lambdaPhi_Theor -> SetLineColor(kBlue);
  hist_lambdaPhi_Theor -> SetFillColor(kBlue);
  hist_lambdaPhi_Theor -> SetFillStyle(3005);

  TH1D *hist_lambdaPhi = new TH1D("hist_lambdaPhi","",N_test,0,N_test);
  hist_lambdaPhi -> SetMarkerColor(kBlue);
  hist_lambdaPhi -> SetLineColor(kBlue);
  hist_lambdaPhi -> SetLineWidth(2);
  hist_lambdaPhi -> SetMarkerStyle(20);

  TH1D *hist_lambdaPhi_Rebin = new TH1D("hist_lambdaPhi_Rebin","",N_test,0,N_test);
  hist_lambdaPhi_Rebin -> SetMarkerColor(kRed);
  hist_lambdaPhi_Rebin -> SetLineColor(kRed);
  hist_lambdaPhi_Rebin -> SetLineWidth(2);
  hist_lambdaPhi_Rebin -> SetMarkerStyle(20);

  TH1D *hist_lambdaPhi_gen = new TH1D("hist_lambdaPhi_gen","",N_test,0,N_test);
  hist_lambdaPhi_gen -> SetLineColor(kGreen+1);
  hist_lambdaPhi_gen -> SetLineWidth(2);

  TH2D *h_lambdaPhi_diff = new TH2D("h_lambdaPhi_diff","",N_test,0,N_test,100,-0.2,0.2);
  h_lambdaPhi_diff -> GetYaxis() -> SetTitle("#lambda_{#phi}|_{GEN} - #lambda_{#phi}|_{REC}");
  h_lambdaPhi_diff -> GetYaxis() -> SetTitleOffset(1.3);

  TH1D *hist_lambdaPhi_diff = new TH1D("hist_lambdaPhi_diff","",N_test,0,N_test);
  hist_lambdaPhi_diff -> SetMarkerColor(kBlue);
  hist_lambdaPhi_diff -> SetLineColor(kBlue);
  hist_lambdaPhi_diff -> SetLineWidth(2);
  hist_lambdaPhi_diff -> SetMarkerStyle(20);

  TH1D *hist_lambdaPhi_diff_Rebin = new TH1D("hist_lambdaPhi_diff_Rebin","",N_test,0,N_test);
  hist_lambdaPhi_diff_Rebin -> SetMarkerColor(kRed);
  hist_lambdaPhi_diff_Rebin -> SetLineColor(kRed);
  hist_lambdaPhi_diff_Rebin -> SetLineWidth(2);
  hist_lambdaPhi_diff_Rebin -> SetMarkerStyle(20);

  TH1D *hist_lambdaPhi_gen_diff = new TH1D("hist_lambdaPhi_gen_diff","",N_test,0,N_test);
  hist_lambdaPhi_gen_diff -> SetLineColor(kGreen+1);
  hist_lambdaPhi_gen_diff -> SetLineWidth(2);

  /*TH1D *lambdaPhi2_diff = new TH1D("lambdaPhi2_diff","",5,0,5);
  lambdaPhi2_diff -> SetMarkerColor(kBlue+2);
  lambdaPhi2_diff -> SetLineColor(kBlue+2);
  lambdaPhi2_diff -> SetLineWidth(2);
  lambdaPhi2_diff -> SetMarkerStyle(20);

  TH1D *lambdaPhi2_diff_Rebin = new TH1D("lambdaPhi2_diff_Rebin","",5,0,5);
  lambdaPhi2_diff_Rebin -> SetMarkerColor(kRed+2);
  lambdaPhi2_diff_Rebin -> SetLineColor(kRed+2);
  lambdaPhi2_diff_Rebin -> SetLineWidth(2);
  lambdaPhi2_diff_Rebin -> SetMarkerStyle(20);*/


  for(int i = 0;i < N_test;i++){
    lambdaTh_gen[i] = lambdaTh[i];
    lambdaPhi_gen[i] = lambdaPhi[i];
    //if(i < 11){
      sprintf(bin_label,"%2.1f",lambdaTh_gen[i]);
      h_lambdaTh -> GetXaxis() -> SetBinLabel(i+1,bin_label);
      sprintf(bin_label,"#lambda_{#theta} = %2.1f",lambdaTh_gen[i]);
      h_lambdaTh_diff -> GetXaxis() -> SetBinLabel(i+1,bin_label);
    //}
    //if(i >= 11 && i < 16){
      sprintf(bin_label,"%2.1f",lambdaPhi_gen[i]);
      h_lambdaPhi -> GetXaxis() -> SetBinLabel(i+1,bin_label);
      sprintf(bin_label,"#lambda_{#phi} = %2.1f",lambdaPhi_gen[i]);
      h_lambdaPhi_diff -> GetXaxis() -> SetBinLabel(i+1,bin_label);
    //}

    // NARROW BINNING GEN
    sprintf(func_name,"fit_func_polarization_VAR_pol%i",i);
    fit_func_CostPhi2D_HE_VAR_pol[i] = new TF2(func_name,Func_W,min_fit_range_Cost,max_fit_range_Cost,min_fit_range_Phi,max_fit_range_Phi,4);
    fit_func_CostPhi2D_HE_VAR_pol[i] -> SetParameters(1000,lambdaTh_gen[i],lambdaPhi_gen[i],0);
    //fit_func_CostPhi2D_HE_VAR_pol[i] -> FixParameter(3,0.);
    hist_CostPhiHE_2pt6_VAR_pol_gen[i] -> Fit(fit_func_CostPhi2D_HE_VAR_pol[i],"RLS0");

    printf("Narrow binning \n");
    printf("INPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",lambdaTh[i],lambdaPhi[i],lambda_ThPhi[i]);
    printf("OUTPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",fit_func_CostPhi2D_HE_VAR_pol[i] -> GetParameter(1),fit_func_CostPhi2D_HE_VAR_pol[i] -> GetParameter(2),fit_func_CostPhi2D_HE_VAR_pol[i] -> GetParameter(3));
    printf("- - - - - - - - - - - \n");

    lambdaTh_rec[i] = fit_func_CostPhi2D_HE_VAR_pol[i] -> GetParameter(1);
    stat_lambdaTh_rec[i] = fit_func_CostPhi2D_HE_VAR_pol[i] -> GetParError(1);
    lambdaPhi_rec[i] = fit_func_CostPhi2D_HE_VAR_pol[i] -> GetParameter(2);
    stat_lambdaPhi_rec[i] = fit_func_CostPhi2D_HE_VAR_pol[i] -> GetParError(2);

    //if(i < 11){
      diff_lambdaTh[i] = lambdaTh_gen[i] - lambdaTh_rec[i];
      printf("lambdaTh gen - lambdaTh rec = %f \n",diff_lambdaTh[i]);
      hist_lambdaTh_Theor -> SetBinContent(i+1,lambdaTh_gen[i]);
      hist_lambdaTh -> SetBinContent(i+1,lambdaTh_rec[i]);
      hist_lambdaTh -> SetBinError(i+1,fit_func_CostPhi2D_HE_VAR_pol[i] -> GetParError(1));
      hist_lambdaTh_diff -> SetBinContent(i+1,diff_lambdaTh[i]);
      hist_lambdaTh_diff -> SetBinError(i+1,fit_func_CostPhi2D_HE_VAR_pol[i] -> GetParError(1));

      diff_lambdaPhi[i] = lambdaPhi_gen[i] - lambdaPhi_rec[i];
      printf("lambdaPhi gen - lambdaPhi rec = %f \n",diff_lambdaPhi[i]);
      hist_lambdaPhi_Theor -> SetBinContent(i+1,lambdaPhi_gen[i]);
      hist_lambdaPhi -> SetBinContent(i+1,lambdaPhi_rec[i]);
      hist_lambdaPhi -> SetBinError(i+1,fit_func_CostPhi2D_HE_VAR_pol[i] -> GetParError(2));
      hist_lambdaPhi_diff -> SetBinContent(i+1,diff_lambdaPhi[i]);
      hist_lambdaPhi_diff -> SetBinError(i+1,fit_func_CostPhi2D_HE_VAR_pol[i] -> GetParError(2));
    //}

    // NARROW BINNING REC
    sprintf(func_name,"fit_func_polarization_VARNO_pol%i",i);
    fit_func_CostPhi2D_HE_VARNO_pol[i] = new TF2(func_name,Func_W,min_fit_range_Cost,max_fit_range_Cost,min_fit_range_Phi,max_fit_range_Phi,4);
    fit_func_CostPhi2D_HE_VARNO_pol[i] -> SetParameters(1000,lambdaTh_gen[i],lambdaPhi_gen[i],0);
    //fit_func_CostPhi2D_HE_VARNO_pol[i] -> FixParameter(3,0.);
    hist_CostPhiHE_2pt6_VAR_pol_rec_AC[i] -> Fit(fit_func_CostPhi2D_HE_VARNO_pol[i],"RLS0");

    printf("Narrow binning \n");
    printf("INPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",lambdaTh[i],lambdaPhi[i],lambda_ThPhi[i]);
    printf("OUTPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParameter(1),fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParameter(2),fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParameter(3));
    printf("- - - - - - - - - - - \n");

    lambdaTh_rec[i] = fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParameter(1);
    stat_lambdaTh_rec[i] = fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParError(1);
    lambdaPhi_rec[i] = fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParameter(2);
    stat_lambdaPhi_rec[i] = fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParError(2);

    //if(i < 11){
      diff_lambdaTh[i] = lambdaTh_gen[i] - lambdaTh_rec[i];
      printf("lambdaTh gen - lambdaTh rec = %f \n",diff_lambdaTh[i]);
      hist_lambdaTh_gen -> SetBinContent(i+1,lambdaTh_rec[i]);
      hist_lambdaTh_gen -> SetBinError(i+1,fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParError(1));
      hist_lambdaTh_gen_diff -> SetBinContent(i+1,diff_lambdaTh[i]);
      hist_lambdaTh_gen_diff -> SetBinError(i+1,fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParError(1));

      diff_lambdaPhi[i] = lambdaPhi_gen[i] - lambdaPhi_rec[i];
      printf("lambdaPhi gen - lambdaPhi rec = %f \n",diff_lambdaPhi[i]);
      hist_lambdaPhi_gen -> SetBinContent(i+1,lambdaPhi_rec[i]);
      hist_lambdaPhi_gen -> SetBinError(i+1,fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParError(2));
      hist_lambdaPhi_gen_diff -> SetBinContent(i+1,diff_lambdaPhi[i]);
      hist_lambdaPhi_gen_diff -> SetBinError(i+1,fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParError(2));
    //}

    /*if(i >= 11 && i < 16){
      diff_lambdaTh[i] = lambdaTh_gen - lambdaTh_rec;
      diff_lambdaPhi[i] = lambdaPhi_gen - lambdaPhi_rec;
      printf("lambdaPhi gen - lambdaPhi rec = %f \n",diff_lambdaPhi[i]);
      lambdaPhi1_diff -> SetBinContent(i-10,diff_lambdaPhi[i]);
      lambdaPhi1_diff -> SetBinError(i-10,fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParError(2));
    }

    if(i >= 16){
      diff_lambdaTh[i] = lambdaTh_gen - lambdaTh_rec;
      diff_lambdaPhi[i] = lambdaPhi_gen - lambdaPhi_rec;
      printf("lambdaPhi gen - lambdaPhi rec = %f \n",diff_lambdaPhi[i]);
      lambdaPhi2_diff -> SetBinContent(i-15,diff_lambdaPhi[i]);
      lambdaPhi2_diff -> SetBinError(i-15,fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParError(2));
    }*/

    // LARGE BINNING REC
    sprintf(func_name,"fit_func_polarization_VARNO_pol_Rebin%i",i);
    fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] = new TF2(func_name,Func_W,min_fit_range_Cost,max_fit_range_Cost,min_fit_range_Phi,max_fit_range_Phi,4);
    fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> SetParameters(1000,lambdaTh_gen[i],lambdaPhi_gen[i],0);
    fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> FixParameter(3,0.);
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_ANAC[i] -> Fit(fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i],"RLS0");

    printf("Large binning \n");
    printf("INPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",lambdaTh[i],lambdaPhi[i],lambda_ThPhi[i]);
    printf("OUTPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParameter(1),fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParameter(2),fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParameter(3));
    printf("- - - - - - - - - - - \n");

    lambdaTh_rec_Rebin[i] = fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParameter(1);
    lambdaPhi_rec_Rebin[i] = fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParameter(2);

    //if(i < 11){
      diff_lambdaTh[i] = lambdaTh_gen[i] - lambdaTh_rec[i];
      printf("lambdaTh gen - lambdaTh rec = %f \n",diff_lambdaTh[i]);
      hist_lambdaTh_Rebin -> SetBinContent(i+1,lambdaTh_rec[i]);
      hist_lambdaTh_Rebin -> SetBinError(i+1,fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParError(1));
      hist_lambdaTh_diff_Rebin -> SetBinContent(i+1,diff_lambdaTh[i]);
      hist_lambdaTh_diff_Rebin -> SetBinError(i+1,fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParError(1));

      diff_lambdaPhi[i] = lambdaPhi_gen[i] - lambdaPhi_rec[i];
      printf("lambdaPhi gen - lambdaPhi rec = %f \n",diff_lambdaPhi[i]);
      hist_lambdaPhi_Rebin -> SetBinContent(i+1,lambdaPhi_rec[i]);
      hist_lambdaPhi_Rebin -> SetBinError(i+1,fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParError(2));
      hist_lambdaPhi_diff_Rebin -> SetBinContent(i+1,diff_lambdaPhi[i]);
      hist_lambdaPhi_diff_Rebin -> SetBinError(i+1,fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParError(2));
    //}

    /*if(i >= 11 && i < 16){
      diff_lambdaTh[i] = lambdaTh_gen - lambdaTh_rec;
      diff_lambdaPhi[i] = lambdaPhi_gen - lambdaPhi_rec;
      printf("lambdaPhi gen - lambdaPhi rec = %f \n",diff_lambdaPhi[i]);
      lambdaPhi1_diff_Rebin -> SetBinContent(i-10,diff_lambdaPhi[i]);
      lambdaPhi1_diff_Rebin -> SetBinError(i-10,fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParError(2));
    }

    if(i >= 16){
      diff_lambdaTh[i] = lambdaTh_gen - lambdaTh_rec;
      diff_lambdaPhi[i] = lambdaPhi_gen - lambdaPhi_rec;
      printf("lambdaPhi gen - lambdaPhi rec = %f \n",diff_lambdaPhi[i]);
      lambdaPhi2_diff_Rebin -> SetBinContent(i-15,diff_lambdaPhi[i]);
      lambdaPhi2_diff_Rebin -> SetBinError(i-15,fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParError(2));
    }*/
  }

  //============================================================================
  printf("Filling the TGraphErrors ... \n");
  //============================================================================

  TH2D *h_lambdaTh_Phi = new TH2D("h_lambdaTh_Phi","",100,-1.2,1.2,100,-0.3,0.3);
  h_lambdaTh_Phi -> GetXaxis() -> SetNdivisions(12);
  h_lambdaTh_Phi -> GetXaxis() -> SetTitle("#lambda_{#theta}");
  h_lambdaTh_Phi -> GetYaxis() -> SetNdivisions(6);
  h_lambdaTh_Phi -> GetYaxis() -> SetTitle("#lambda_{#phi}");

  TGraphErrors *gra_lambdaTh_Phi_Theor = new TGraphErrors(N_test,lambdaTh,lambdaPhi,0,0);
  gra_lambdaTh_Phi_Theor -> SetMarkerStyle(20);
  gra_lambdaTh_Phi_Theor -> SetMarkerColor(kRed+1);

  TGraphErrors *gra_lambdaTh_Phi_rec = new TGraphErrors(N_test,lambdaTh_rec,lambdaPhi_rec,stat_lambdaTh_rec,stat_lambdaPhi_rec);
  gra_lambdaTh_Phi_rec -> SetMarkerStyle(24);
  gra_lambdaTh_Phi_rec -> SetMarkerColor(kBlue+1);

  TGraphErrors *gra_lambdaTh_Phi_rec_Rebin = new TGraphErrors(N_test,lambdaTh_rec_Rebin,lambdaPhi_rec_Rebin,stat_lambdaTh_rec_Rebin,stat_lambdaPhi_rec_Rebin);
  gra_lambdaTh_Phi_rec_Rebin -> SetMarkerStyle(24);
  gra_lambdaTh_Phi_rec_Rebin -> SetMarkerColor(kRed+1);

  TCanvas *c_lambdaTh_Phi = new TCanvas("c_lambdaTh_Phi","c_lambdaTh_Phi",4,132,1024,768);
  c_lambdaTh_Phi -> SetGrid();
  h_lambdaTh_Phi -> Draw();
  gra_lambdaTh_Phi_Theor -> Draw("sameP");
  gra_lambdaTh_Phi_rec -> Draw("sameP");
  gra_lambdaTh_Phi_rec_Rebin -> Draw("sameP");

  //============================================================================

  TLegend *leg_Th = new TLegend(0.65,0.70,.85,.84,"","brNDC");
  leg_Th -> SetBorderSize(0);
  leg_Th -> SetFillColor(10);
  leg_Th -> SetFillStyle(1);
  leg_Th -> SetLineStyle(0);
  leg_Th -> SetLineColor(0);
  leg_Th -> SetTextFont(42);
  leg_Th -> SetTextSize(0.035);
  leg_Th -> AddEntry(hist_lambdaTh_Theor,"#lambda_{#theta}|_{GEN}","F");
  leg_Th -> AddEntry(hist_lambdaTh_Rebin,"#lambda_{#theta}|_{REC}","PL");

  TCanvas *c_lambdaTh = new TCanvas("c_lambdaTh","c_lambdaTh",4,132,1024,768);
  h_lambdaTh -> Draw();
  hist_lambdaTh_Theor -> Draw("same");
  //hist_lambdaTh_gen -> Draw("samePE");
  hist_lambdaTh -> Draw("samePE");
  hist_lambdaTh_Rebin -> Draw("samePE");
  leg_Th -> Draw("same");

  TLegend *leg_Phi = new TLegend(0.65,0.70,.85,.84,"","brNDC");
  leg_Phi -> SetBorderSize(0);
  leg_Phi -> SetFillColor(10);
  leg_Phi -> SetFillStyle(1);
  leg_Phi -> SetLineStyle(0);
  leg_Phi -> SetLineColor(0);
  leg_Phi -> SetTextFont(42);
  leg_Phi -> SetTextSize(0.035);
  leg_Phi -> AddEntry(hist_lambdaPhi_Theor,"#lambda_{#phi}|_{GEN}","F");
  leg_Phi -> AddEntry(hist_lambdaPhi_Rebin,"#lambda_{#phi}|_{REC}","PL");

  TCanvas *c_lambdaPhi = new TCanvas("c_lambdaPhi","c_lambdaPhi",4,132,1024,768);
  h_lambdaPhi -> Draw();
  hist_lambdaPhi_Theor -> Draw("same");
  //hist_lambdaPhi_gen -> Draw("samePE");
  hist_lambdaPhi -> Draw("samePE");
  hist_lambdaPhi_Rebin -> Draw("samePE");
  leg_Phi -> Draw("same");

  TLine *l_unity_lambdaTh = new TLine(0,0,N_test,0);

  TCanvas *c_lambdaTh_diff = new TCanvas("c_lambdaTh_diff","c_lambdaTh_diff",4,132,1024,768);
  h_lambdaTh_diff -> Draw();
  l_unity_lambdaTh -> Draw("same");
  hist_lambdaTh_gen_diff -> Draw("same");
  hist_lambdaTh_diff -> Draw("samePE");
  hist_lambdaTh_diff_Rebin -> Draw("samePE");

  TLine *l_unity_lambdaPhi = new TLine(0,0,N_test,0);

  TCanvas *c_lambdaPhi_diff = new TCanvas("c_lambdaPhi_diff","c_lambdaPhi_diff",4,132,1024,768);
  h_lambdaPhi_diff -> Draw();
  l_unity_lambdaPhi -> Draw("same");
  hist_lambdaPhi_gen_diff -> Draw("samePE");
  hist_lambdaPhi_diff -> Draw("samePE");
  hist_lambdaPhi_diff_Rebin -> Draw("samePE");
  //for(int i = 0;i < N_test;i++){
    //printf("D_Th = %f ; F_Phi = %f \n",diff_lambdaTh[i],diff_lambdaPhi[i]);
  //}

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

  char INPUT_FILE_NAME[100];
  sprintf(INPUT_FILE_NAME,"../INPUT/sample_tree_for_closure_test.root");
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


    while(count_Dimu_rec < 100000){
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
  sprintf(OUTPUT_FILE_NAME,"/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/POLARIZED_DISTRIBUTIONS/variable_Jpsi_polarization.root");
  TFile *output_file = new TFile(OUTPUT_FILE_NAME,"RECREATE");
  for(int i = 0;i < N_test;i++){
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
