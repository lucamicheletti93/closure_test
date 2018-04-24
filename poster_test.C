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

const int N_test = 11;
double lambdaTh = -1.;
double lambdaPhi = 0.;
double lambda_ThPhi = 0;
double step_width = 0.2;

double min_fit_range_Cost = -0.6;
double max_fit_range_Cost = 0.6;

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

  //TF2 *func_CostPhi2D_HE_VAR_pol = new TF2("func_CostPhi2D_HE_VAR_pol",Func_W,min_Cost,max_Cost,min_Phi,max_Phi,4);
  //func_CostPhi2D_HE_VAR_pol -> SetParameters(pol_par_CostPhi2D_VAR[0],pol_par_CostPhi2D_VAR[1],pol_par_CostPhi2D_VAR[2],pol_par_CostPhi2D_VAR[3]);
  //func_CostPhi2D_HE_VAR_pol -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");

  //============================================================================
  printf("2) Reading the the file.root in ../OUTPUT ... \n");
  //============================================================================

  TFile *input_file = new TFile(INPUT_FILE_NAME,"READ");
  TH2D *hist_CostPhiHE_2pt6_VAR_pol_rec[N_test];
  char hist_name[60];
  char func_name[60];

  for(int i = 0;i < N_test;i++){
    sprintf(hist_name,"hist_polarization%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_rec[i] = (TH2D*) input_file -> Get(hist_name);
  }

  //============================================================================
  printf("4) Normalizing the Rebin histo to bin area ... \n"); // AN = Area Normalized
  //============================================================================

  TH2D *hist_CostPhiHE_2pt6_VAR_pol_rec_AN[N_test];

  for(int i = 0;i < N_test;i++){
    sprintf(hist_name,"hist_polarization_area_normalized%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_rec_AN[i] = new TH2D(hist_name,hist_name,N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
    for(int j = 0;j < N_cost_bins_BC;j++){
      for(int k = 0;k < N_phi_bins_BC;k++){
        hist_CostPhiHE_2pt6_VAR_pol_rec_AN[i] -> SetBinContent(j+1,k+1,(hist_CostPhiHE_2pt6_VAR_pol_rec[i] -> GetBinContent(j+1,k+1))/bin_area_BC[j][k]);
        hist_CostPhiHE_2pt6_VAR_pol_rec_AN[i] -> SetBinError(j+1,k+1,(hist_CostPhiHE_2pt6_VAR_pol_rec[i] -> GetBinError(j+1,k+1))/bin_area_BC[j][k]);
      }
    }
  }

  //TCanvas *c_CostPhiHE_2pt6_VAR_pol_rec_AN = new TCanvas("c_CostPhiHE_2pt6_VAR_pol_rec_AN","c_CostPhiHE_2pt6_VAR_pol_rec_AN",20,20,600,600);
  //c_CostPhiHE_2pt6_VAR_pol_rec_AN -> Divide(2,5);

  //for(int i = 0;i < 1;i++){
    //c_CostPhiHE_2pt6_VAR_pol_rec_AN -> cd(i+1);
    //hist_CostPhiHE_2pt6_VAR_pol_rec_AN[i] -> Draw("COLZ");
  //}

  //============================================================================
  printf("5) Correcting REC for Acc x Eff ... \n"); // AC = Acceptance Corrected
  //============================================================================

  sprintf(INPUT_FILE_NAME,"../INPUT/Histo_accxeff_FromOfficialTree_Jpsi_PbPb_Nopol.root");
  TFile *accxeff_file = new TFile(INPUT_FILE_NAME,"READ");
  TH2D *hist_accxeff_HE_2pt6_NO_pol = (TH2D*) accxeff_file -> Get("hist_accxeff_HE_2pt6");

  //TFile *file = new TFile("../OUTPUT/Histo_for_closure_test.root","READ");
  //TH2D *hist_accxeff_HE_2pt6_NO_pol = (TH2D*) file -> Get("hist_accxeff_2pt6_NO_pol_Rebin");

  TH2D *hist_CostPhiHE_2pt6_VAR_pol_rec_ANAC[N_test];

  for(int k = 0;k < N_test;k++){
    sprintf(hist_name,"hist_accxeff_NOpol%i",k);
    hist_CostPhiHE_2pt6_VAR_pol_rec_ANAC[k] = new TH2D(hist_name,hist_name,N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);
    hist_CostPhiHE_2pt6_VAR_pol_rec_ANAC[k] -> Divide(hist_CostPhiHE_2pt6_VAR_pol_rec_AN[k],hist_accxeff_HE_2pt6_NO_pol,1,1);
  }

  //TCanvas *c_CostPhiHE_2pt6_VAR_pol_rec_ANAC = new TCanvas("c_CostPhiHE_2pt6_VAR_pol_rec_ANAC","c_CostPhiHE_2pt6_VAR_pol_rec_ANAC",20,20,600,600);
  //c_CostPhiHE_2pt6_VAR_pol_rec_ANAC -> Divide(2,5);

  //for(int i = 0;i < 1;i++){
    //c_CostPhiHE_2pt6_VAR_pol_rec_ANAC -> cd(i+1);
    //hist_CostPhiHE_2pt6_VAR_pol_rec_ANAC[i] -> Draw("COLZ");
  //}

  //============================================================================
  printf("6) Fit of the REC corrected distributions rebinned with NO pol Acc x Eff correction... \n");
  //============================================================================

  TF2 *fit_func_CostPhi2D_HE_VARNO_pol[N_test];
  double lambdaTh_gen = 0, lambdaTh_rec = 0;
  double diff_lambdaTh = 0;
  char bin_label[40];

  TH1D *lambdaTh_gen_rec = new TH1D("lambdaTh_gen_rec","",11,0,11);
  lambdaTh_gen_rec -> SetMarkerColor(kRed);
  lambdaTh_gen_rec -> SetMarkerStyle(20);

  TH2D *h_lambdaTh_gen_rec = new TH2D("h_lambdaTh_gen_rec","",11,0,11,100,-0.2,0.2);
  h_lambdaTh_gen_rec -> GetYaxis() -> SetTitle("#lambda_{#theta}|_{GEN} - #lambda_{#theta}|_{REC}");
  h_lambdaTh_gen_rec -> GetYaxis() -> SetTitleOffset(1.3);

  for(int i = 0;i < N_test;i++){
    sprintf(func_name,"fit_func_polarization_VARNOpol_%i",i);
    fit_func_CostPhi2D_HE_VARNO_pol[i] = new TF2(func_name,Func_W,min_fit_range_Cost,max_fit_range_Cost,0,PI,4);
    fit_func_CostPhi2D_HE_VARNO_pol[i] -> SetParameters(1000,-1,0,0);
    hist_CostPhiHE_2pt6_VAR_pol_rec_ANAC[i] -> Fit(fit_func_CostPhi2D_HE_VARNO_pol[i],"RLS0");

    printf("INPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",lambdaTh + step_width*i,lambdaPhi,lambda_ThPhi);
    printf("OUTPUT POLARIZATION = (%4.3f,%4.3f,%4.3f) \n",fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParameter(1),fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParameter(2),fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParameter(3));
    printf("- - - - - - - - - - - \n");
    lambdaTh_gen = lambdaTh + 0.2*i;
    lambdaTh_rec = fit_func_CostPhi2D_HE_VARNO_pol[i] -> GetParameter(1);

    sprintf(bin_label,"#lambda_{#theta} = %2.1f",lambdaTh_gen);
    diff_lambdaTh = lambdaTh_gen - lambdaTh_rec;
    printf("lambdaTh gen - lambdaTh rec = %f \n",diff_lambdaTh);
    lambdaTh_gen_rec -> SetBinContent(i+1,diff_lambdaTh);
    h_lambdaTh_gen_rec -> GetXaxis() -> SetBinLabel(i+1,bin_label);
  }

  TLine *l_unity = new TLine(0,0,11,0);

  TCanvas *c_lambdaTh_gen_rec = new TCanvas("c_lambdaTh_gen_rec","c_lambdaTh_gen_rec",4,132,1024,768);
  h_lambdaTh_gen_rec -> Draw();
  l_unity -> Draw("same");
  lambdaTh_gen_rec -> Draw("sameP");

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

  Int_t NDimu_rec;
  Double_t DimuPt_rec[3000], DimuY_rec[3000];
  Double_t DimuMass_rec[3000];
  Int_t DimuMatch_rec[3000];
  Double_t CostHE_rec[3000], PhiHE_rec[3000], CostCS_rec[3000], PhiCS_rec[3000];

  TTree *input_tree = (TTree*) input_file -> Get(INPUT_TREE_NAME);
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

  TF2 *func_width = new TF2("func_width",Func_W,-1,1,0,PI,4);
  func_width -> SetParameters(1000,lambdaTh,lambdaPhi,lambda_ThPhi);

  for(int i = 0;i < N_test;i++){
    sprintf(func_name,"func_polarization%i",i);
    func_CostPhi2D_HE_VAR_pol[i] = new TF2(func_name,Func_W,-1,1,0,PI,4);
    func_CostPhi2D_HE_VAR_pol[i] -> SetParameters(1000,lambdaTh + step_width*i,lambdaPhi,lambda_ThPhi);
    sprintf(hist_name,"hist_polarization%i",i);
    hist_CostPhiHE_2pt6_VAR_pol_rec[i] = new TH2D(hist_name,hist_name,N_cost_bins_BC,value_cost_BC,N_phi_bins_BC,value_phi_BC);


    while(count_Dimu_rec < 100000){
      input_tree -> GetEntry(event_index);
      for(int k = 0;k < NDimu_rec;k++){
        if(DimuPt_rec[k] > 2 && DimuPt_rec[k] <= 6){


          if(DimuY_rec[k] > -4. && DimuY_rec[k] < -2.5){
            if(DimuMatch_rec[k] == 2){
              if(DimuMass_rec[k] > 2 && DimuMass_rec[k] < 5){
                weight_CostPhi2D = func_CostPhi2D_HE_VAR_pol[i] -> Eval(CostHE_rec[k],PhiHE_rec[k])/func_CostPhi2D_HE_VAR_pol[i] -> GetMaximum();
                hist_CostPhiHE_2pt6_VAR_pol_rec[i] -> Fill(CostHE_rec[k],PhiHE_rec[k],weight_CostPhi2D);
                count_Dimu_rec++;
              }
            }
          }
        }
      }
      event_index++;
    }
    printf("POLARIZATION = (%2.1f,%2.1f,%2.1f) \n",lambdaTh + step_width*i,lambdaPhi,lambda_ThPhi);
    event_index = 0;
    count_Dimu_rec = 0;
  }

  char OUTPUT_FILE_NAME[50];
  sprintf(OUTPUT_FILE_NAME,"../OUTPUT/variable_Jpsi_polarization.root");
  TFile *output_file = new TFile(OUTPUT_FILE_NAME,"RECREATE");
  for(int i = 0;i < N_test;i++){
    hist_CostPhiHE_2pt6_VAR_pol_rec[i] -> Write();
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
