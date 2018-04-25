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
double lambdaTh[N_test] = {-1.,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.,0.2,0.2,0.2,0.2,0.2,0.5,0.5,0.5,0.5,0.5};
double lambdaPhi[N_test] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-0.2,-0.1,-0.0,0.1,0.2,-0.2,-0.1,-0.0,0.1,0.2};
double lambda_ThPhi[N_test] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
//const int N_test = 2;
//double lambdaTh[N_test] = {0.2,0.2};
//double lambdaPhi[N_test] = {-0.5,0.5};
//double lambda_ThPhi[N_test] = {0.,0};

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

  char INPUT_FILE_NAME[100];
  sprintf(INPUT_FILE_NAME,"../OUTPUT/variable_Jpsi_polarization.root");
  char INPUT_TREE_NAME[100];
  sprintf(INPUT_TREE_NAME,"MCTree");

  //============================================================================
  printf("2) Reading the the file.root in ../OUTPUT ... \n");
  //============================================================================

  TFile *input_file = new TFile(INPUT_FILE_NAME,"READ");
  TH2D *hist_CostPhiHE_2pt6_VAR_pol_rec[N_test], *hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[N_test];
  TH2D *hist_CostPhiHE_2pt6_VAR_pol_gen[N_test], *hist_CostPhiHE_2pt6_VAR_pol_gen_Rebin[N_test];
  char hist_name[60];
  char func_name[60];

  for(int i = 0;i < N_test;i++){
    sprintf(hist_name,"hist_rec_polarization%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_rec[i] = (TH2D*) input_file -> Get(hist_name);
    sprintf(hist_name,"hist_rec_polarization_Rebin%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[i] = (TH2D*) input_file -> Get(hist_name);
    sprintf(hist_name,"hist_gen_polarization%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_gen[i] = (TH2D*) input_file -> Get(hist_name);
    sprintf(hist_name,"hist_gen_polarization_Rebin%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_gen_Rebin[i] = (TH2D*) input_file -> Get(hist_name);
  }

  //============================================================================
  printf("4) Normalizing the Rebin histo to bin area ... \n"); // AN = Area Normalized
  //============================================================================

  TH2D *hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_AN[N_test];
  TH2D *hist_CostPhiHE_2pt6_VAR_pol_gen_Rebin_AN[N_test];

  for(int i = 0;i < N_test;i++){
    sprintf(hist_name,"hist_polarization_area_normalized%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_AN[i] = new TH2D(hist_name,hist_name,N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
    sprintf(hist_name,"hist_polarization_area_normalized_gen%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_gen_Rebin_AN[i] = new TH2D(hist_name,hist_name,N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
    for(int j = 0;j < N_cost_bins_BC;j++){
      for(int k = 0;k < N_phi_bins_BC;k++){
        hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_AN[i] -> SetBinContent(j+1,k+1,(hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[i] -> GetBinContent(j+1,k+1))/bin_area_BC[j][k]);
        hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_AN[i] -> SetBinError(j+1,k+1,(hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[i] -> GetBinError(j+1,k+1))/bin_area_BC[j][k]);

        hist_CostPhiHE_2pt6_VAR_pol_gen_Rebin_AN[i] -> SetBinContent(j+1,k+1,(hist_CostPhiHE_2pt6_VAR_pol_gen_Rebin[i] -> GetBinContent(j+1,k+1))/bin_area_BC[j][k]);
        hist_CostPhiHE_2pt6_VAR_pol_gen_Rebin_AN[i] -> SetBinError(j+1,k+1,(hist_CostPhiHE_2pt6_VAR_pol_gen_Rebin[i] -> GetBinError(j+1,k+1))/bin_area_BC[j][k]);
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
  printf("6) Fit of the REC corrected distributions rebinned with NO pol Acc x Eff correction... \n");
  //============================================================================

  TF2 *fit_func_CostPhi2D_HE_VARNO_pol_gen[N_test], *fit_func_CostPhi2D_HE_VARNO_pol_gen_Rebin[N_test];
  TF2 *fit_func_CostPhi2D_HE_VARNO_pol[N_test], *fit_func_CostPhi2D_HE_VARNO_pol_Rebin[N_test];
  double lambdaTh_gen = 0, lambdaPhi_gen = 0;
  double lambdaTh_rec = 0, lambdaPhi_rec = 0;
  double diff_lambdaTh = 0, diff_lambdaPhi = 0;
  char bin_label[40];

  TH2D *h_lambdaTh_gen_rec = new TH2D("h_lambdaTh_gen_rec","",N_test-5,0,N_test-5,100,-0.3,0.3);
  h_lambdaTh_gen_rec -> GetYaxis() -> SetTitle("#lambda_{#theta}|_{GEN} - #lambda_{#theta}|_{REC}");
  h_lambdaTh_gen_rec -> GetYaxis() -> SetTitleOffset(1.3);

  TH1D *lambdaTh_gen_rec = new TH1D("lambdaTh_gen_rec","",N_test-5,0,N_test-5);
  lambdaTh_gen_rec -> SetMarkerColor(kBlue);
  lambdaTh_gen_rec -> SetLineColor(kBlue);
  lambdaTh_gen_rec -> SetLineWidth(2);
  lambdaTh_gen_rec -> SetMarkerStyle(20);

  TH1D *lambdaTh_gen_rec_Rebin = new TH1D("lambdaTh_gen_rec_Rebin","",N_test-5,0,N_test-5);
  lambdaTh_gen_rec_Rebin -> SetMarkerColor(kRed);
  lambdaTh_gen_rec_Rebin -> SetLineColor(kRed);
  lambdaTh_gen_rec_Rebin -> SetLineWidth(2);
  lambdaTh_gen_rec_Rebin -> SetMarkerStyle(20);

  TH1D *lambdaTh_gen_gen = new TH1D("lambdaTh_gen_gen","",N_test-5,0,N_test-5);
  lambdaTh_gen_gen -> SetMarkerColor(kGreen+1);
  lambdaTh_gen_gen -> SetLineColor(kGreen+1);
  lambdaTh_gen_gen -> SetLineWidth(2);
  lambdaTh_gen_gen -> SetMarkerStyle(21);

  TH1D *lambdaTh_gen_gen_Rebin = new TH1D("lambdaTh_gen_gen_Rebin","",5,0,5);
  lambdaTh_gen_gen_Rebin -> SetMarkerColor(kOrange);
  lambdaTh_gen_gen_Rebin -> SetLineColor(kOrange);
  lambdaTh_gen_gen_Rebin -> SetLineWidth(2);
  lambdaTh_gen_gen_Rebin -> SetMarkerStyle(21);

  TH2D *h_lambdaPhi_gen_rec = new TH2D("h_lambdaPhi_gen_rec","",5,0,5,100,-0.3,0.3);
  h_lambdaPhi_gen_rec -> GetYaxis() -> SetTitle("#lambda_{#theta}|_{GEN} - #lambda_{#theta}|_{REC}");
  h_lambdaPhi_gen_rec -> GetYaxis() -> SetTitleOffset(1.3);

  TH1D *lambdaPhi_gen_rec = new TH1D("lambdaPhi_gen_rec","",5,0,5);
  lambdaPhi_gen_rec -> SetMarkerColor(kBlue);
  lambdaPhi_gen_rec -> SetLineColor(kBlue);
  lambdaPhi_gen_rec -> SetLineWidth(2);
  lambdaPhi_gen_rec -> SetMarkerStyle(20);

  TH1D *lambdaPhi_gen_rec_Rebin = new TH1D("lambdaPhi_gen_rec_Rebin","",5,0,5);
  lambdaPhi_gen_rec_Rebin -> SetMarkerColor(kRed);
  lambdaPhi_gen_rec_Rebin -> SetLineColor(kRed);
  lambdaPhi_gen_rec_Rebin -> SetLineWidth(2);
  lambdaPhi_gen_rec_Rebin -> SetMarkerStyle(20);

  TH1D *lambdaPhi_gen_gen = new TH1D("lambdaPhi_gen_gen","",5,0,5);
  lambdaPhi_gen_gen -> SetMarkerColor(kGreen+1);
  lambdaPhi_gen_gen -> SetLineColor(kGreen+1);
  lambdaPhi_gen_gen -> SetLineWidth(2);
  lambdaPhi_gen_gen -> SetMarkerStyle(21);

  TH1D *lambdaPhi_gen_gen_Rebin = new TH1D("lambdaPhi_gen_gen_Rebin","",5,0,5);
  lambdaPhi_gen_gen_Rebin -> SetMarkerColor(kOrange);
  lambdaPhi_gen_gen_Rebin -> SetLineColor(kOrange);
  lambdaPhi_gen_gen_Rebin-> SetLineWidth(2);
  lambdaPhi_gen_gen_Rebin -> SetMarkerStyle(21);

  for(int i = 0;i < N_test;i++){
    lambdaTh_gen = lambdaTh[i];
    lambdaPhi_gen = lambdaPhi[i];
    if(i < 11){
      sprintf(bin_label,"#lambda_{#theta} = %2.1f",lambdaTh_gen);
      h_lambdaTh_gen_rec -> GetXaxis() -> SetBinLabel(i+1,bin_label);
    }
    if(i >= 11 && i < 16){
      sprintf(bin_label,"#lambda_{#phi} = %2.1f",lambdaPhi_gen);
      h_lambdaPhi_gen_rec -> GetXaxis() -> SetBinLabel(i-10,bin_label);
    }
    if(i >= 16){
      sprintf(bin_label,"#lambda_{#phi} = %2.1f",lambdaPhi_gen);
      h_lambdaPhi_gen_rec -> GetXaxis() -> SetBinLabel(i-15,bin_label);
    }

    // NARROW BINNING GEN
    sprintf(func_name,"fit_func_polarization_VARNO_pol_gen%i",i);
    fit_func_CostPhi2D_HE_VARNO_pol_gen[i] = new TF2(func_name,Func_W,min_fit_range_Cost,max_fit_range_Cost,min_fit_range_Phi,max_fit_range_Phi,4);
    fit_func_CostPhi2D_HE_VARNO_pol_gen[i] -> SetParameters(1000,lambdaTh_gen,lambdaPhi_gen,0);
    hist_CostPhiHE_2pt6_VAR_pol_gen[i] -> Fit(fit_func_CostPhi2D_HE_VARNO_pol_gen[i],"RLS0");

    printf("Narrow binning \n");
    printf("INPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",lambdaTh[i],lambdaPhi[i],lambda_ThPhi[i]);
    printf("OUTPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",fit_func_CostPhi2D_HE_VARNO_pol_gen[i] -> GetParameter(1),fit_func_CostPhi2D_HE_VARNO_pol_gen[i] -> GetParameter(2),fit_func_CostPhi2D_HE_VARNO_pol_gen[i] -> GetParameter(3));
    printf("- - - - - - - - - - - \n");

    lambdaTh_rec = fit_func_CostPhi2D_HE_VARNO_pol_gen[i] -> GetParameter(1);
    lambdaPhi_rec = fit_func_CostPhi2D_HE_VARNO_pol_gen[i] -> GetParameter(2);

    if(i < 11){
      diff_lambdaTh = lambdaTh_gen - lambdaTh_rec;
      printf("lambdaTh gen - lambdaTh rec = %f \n",diff_lambdaTh);
      lambdaTh_gen_gen -> SetBinContent(i+1,diff_lambdaTh);
      lambdaTh_gen_gen -> SetBinError(i+1,fit_func_CostPhi2D_HE_VARNO_pol_gen[i] -> GetParError(1));
    }

    if(i >= 11 && i < 16){
      diff_lambdaPhi = lambdaPhi_gen - lambdaPhi_rec;
      printf("lambdaPhi gen - lambdaPhi rec = %f \n",diff_lambdaPhi);
      lambdaPhi_gen_gen -> SetBinContent(i-10,diff_lambdaPhi);
      lambdaPhi_gen_gen -> SetBinError(i-10,fit_func_CostPhi2D_HE_VARNO_pol_gen[i] -> GetParError(2));
    }

    // LARGE BINNING GEN
    sprintf(func_name,"fit_func_polarization_VARNO_pol_Rebin_gen%i",i);
    fit_func_CostPhi2D_HE_VARNO_pol_gen_Rebin[i] = new TF2(func_name,Func_W,min_fit_range_Cost,max_fit_range_Cost,0,PI,4);
    fit_func_CostPhi2D_HE_VARNO_pol_gen_Rebin[i] -> SetParameters(1000,lambdaTh_gen,lambdaPhi_gen,0);
    hist_CostPhiHE_2pt6_VAR_pol_gen_Rebin_AN[i] -> Fit(fit_func_CostPhi2D_HE_VARNO_pol_gen_Rebin[i],"RLS0");

    printf("Large binning \n");
    printf("INPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",lambdaTh[i],lambdaPhi[i],lambda_ThPhi[i]);
    printf("OUTPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",fit_func_CostPhi2D_HE_VARNO_pol_gen_Rebin[i] -> GetParameter(1),fit_func_CostPhi2D_HE_VARNO_pol_gen_Rebin[i] -> GetParameter(2),fit_func_CostPhi2D_HE_VARNO_pol_gen_Rebin[i] -> GetParameter(3));
    printf("- - - - - - - - - - - \n");
    lambdaTh_rec = fit_func_CostPhi2D_HE_VARNO_pol_gen_Rebin[i] -> GetParameter(1);

    if(i < 11){
      diff_lambdaTh = lambdaTh_gen - lambdaTh_rec;
      printf("lambdaTh gen - lambdaTh rec = %f \n",diff_lambdaTh);
      lambdaTh_gen_gen_Rebin -> SetBinContent(i+1,diff_lambdaTh);
      lambdaTh_gen_gen_Rebin -> SetBinError(i+1,fit_func_CostPhi2D_HE_VARNO_pol_gen[i] -> GetParError(1));
    }

    if(i >= 11){
      diff_lambdaPhi = lambdaPhi_gen - lambdaPhi_rec;
      printf("lambdaPhi gen - lambdaPhi rec = %f \n",diff_lambdaPhi);
      lambdaPhi_gen_gen_Rebin -> SetBinContent(i-10,diff_lambdaPhi);
      lambdaPhi_gen_gen_Rebin -> SetBinError(i-10,fit_func_CostPhi2D_HE_VARNO_pol_gen[i] -> GetParError(2));
    }

    //==========================================================================

    // NARROW BINNING REC
    sprintf(func_name,"fit_func_polarization_VARNO_pol%i",i);
    fit_func_CostPhi2D_HE_VARNO_pol[i] = new TF2(func_name,Func_W,min_fit_range_Cost,max_fit_range_Cost,min_fit_range_Phi,max_fit_range_Phi,4);
    fit_func_CostPhi2D_HE_VARNO_pol[i] -> SetParameters(1000,lambdaTh_gen,lambdaPhi_gen,0);
    hist_CostPhiHE_2pt6_VAR_pol_rec_AC[i] -> Fit(fit_func_CostPhi2D_HE_VARNO_pol[i],"RLS0");

    printf("Narrow binning \n");
    printf("INPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",lambdaTh[i],lambdaPhi[i],lambda_ThPhi[i]);
    printf("OUTPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParameter(1),fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParameter(2),fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParameter(3));
    printf("- - - - - - - - - - - \n");

    lambdaTh_rec = fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParameter(1);
    lambdaPhi_rec = fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParameter(2);

    if(i < 11){
      diff_lambdaTh = lambdaTh_gen - lambdaTh_rec;
      printf("lambdaTh gen - lambdaTh rec = %f \n",diff_lambdaTh);
      lambdaTh_gen_rec -> SetBinContent(i+1,diff_lambdaTh);
      lambdaTh_gen_rec -> SetBinError(i+1,fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParError(1));
    }

    if(i >= 11){
      diff_lambdaPhi = lambdaPhi_gen - lambdaPhi_rec;
      printf("lambdaPhi gen - lambdaPhi rec = %f \n",diff_lambdaPhi);
      lambdaPhi_gen_rec -> SetBinContent(i-10,diff_lambdaPhi);
      lambdaPhi_gen_rec -> SetBinError(i-10,fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParError(2));
    }

    // LARGE BINNING REC
    sprintf(func_name,"fit_func_polarization_VARNO_pol_Rebin%i",i);
    fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] = new TF2(func_name,Func_W,min_fit_range_Cost,max_fit_range_Cost,min_fit_range_Phi,max_fit_range_Phi,4);
    fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> SetParameters(1000,lambdaTh_gen,lambdaPhi_gen,0);
    hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin_ANAC[i] -> Fit(fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i],"RLS0");

    printf("Large binning \n");
    printf("INPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",lambdaTh[i],lambdaPhi[i],lambda_ThPhi[i]);
    printf("OUTPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParameter(1),fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParameter(2),fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParameter(3));
    printf("- - - - - - - - - - - \n");

    lambdaTh_rec = fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParameter(1);
    lambdaPhi_rec = fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParameter(2);

    if(i < 11){
      diff_lambdaTh = lambdaTh_gen - lambdaTh_rec;
      printf("lambdaTh gen - lambdaTh rec = %f \n",diff_lambdaTh);
      lambdaTh_gen_rec_Rebin -> SetBinContent(i+1,diff_lambdaTh);
      lambdaTh_gen_rec_Rebin -> SetBinError(i+1,fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParError(1));
    }

    if(i >= 11){
      diff_lambdaPhi = lambdaPhi_gen - lambdaPhi_rec;
      printf("lambdaPhi gen - lambdaPhi rec = %f \n",diff_lambdaPhi);
      lambdaPhi_gen_rec_Rebin -> SetBinContent(i-10,diff_lambdaPhi);
      lambdaPhi_gen_rec_Rebin -> SetBinError(i-10,fit_func_CostPhi2D_HE_VARNO_pol_Rebin[i] -> GetParError(2));
    }
  }

  TLine *l_unity_lambdaTh = new TLine(0,0,11,0);

  TCanvas *c_lambdaTh_gen_rec_Rebin = new TCanvas("c_lambdaTh_gen_rec_Rebin","c_lambdaTh_gen_rec_Rebin",4,132,1024,768);
  h_lambdaTh_gen_rec -> Draw();
  l_unity_lambdaTh -> Draw("same");
  lambdaTh_gen_rec -> Draw("samePE");
  lambdaTh_gen_gen -> Draw("samePE");
  lambdaTh_gen_rec_Rebin -> Draw("samePE");
  lambdaTh_gen_gen_Rebin -> Draw("sameP");

  TLine *l_unity_lambdaPhi = new TLine(0,0,5,0);

  TCanvas *c_lambdaPhi_gen_rec_Rebin = new TCanvas("c_lambdaPhi_gen_rec_Rebin","c_lambdaPhi_gen_rec_Rebin",4,132,1024,768);
  h_lambdaPhi_gen_rec -> Draw();
  l_unity_lambdaPhi -> Draw("same");
  lambdaPhi_gen_rec -> Draw("samePE");
  lambdaPhi_gen_gen -> Draw("samePE");
  lambdaPhi_gen_rec_Rebin -> Draw("samePE");
  lambdaPhi_gen_gen_Rebin -> Draw("sameP");

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
          weight_CostPhi2D = func_CostPhi2D_HE_VAR_pol[i] -> Eval(CostHE_gen[k],PhiHE_gen[k])/func_CostPhi2D_HE_VAR_pol[i] -> GetMaximum();
          hist_CostPhiHE_2pt6_VAR_pol_gen[i] -> Fill(CostHE_gen[k],PhiHE_gen[k],weight_CostPhi2D);
          hist_CostPhiHE_2pt6_VAR_pol_gen_Rebin[i] -> Fill(CostHE_gen[k],PhiHE_gen[k],weight_CostPhi2D);
        }
      }

      for(int k = 0;k < NDimu_rec;k++){
        if(DimuPt_rec[k] > 2 && DimuPt_rec[k] <= 6){
          if(DimuY_rec[k] > -4. && DimuY_rec[k] < -2.5){
            if(DimuMatch_rec[k] == 2){
              if(DimuMass_rec[k] > 2 && DimuMass_rec[k] < 5){
                weight_CostPhi2D = func_CostPhi2D_HE_VAR_pol[i] -> Eval(CostHE_rec[k],PhiHE_rec[k])/func_CostPhi2D_HE_VAR_pol[i] -> GetMaximum();
                hist_CostPhiHE_2pt6_VAR_pol_rec[i] -> Fill(CostHE_rec[k],PhiHE_rec[k],weight_CostPhi2D);
                hist_CostPhiHE_2pt6_VAR_pol_rec_Rebin[i] -> Fill(CostHE_rec[k],PhiHE_rec[k],weight_CostPhi2D);
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

  char OUTPUT_FILE_NAME[50];
  sprintf(OUTPUT_FILE_NAME,"../OUTPUT/variable_Jpsi_polarization.root");
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
