#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>

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

double Func_W(double *, double *);

void generate_polarized_distribution(){

  //============================================================================
  // SET MAIN QUANTITIES
  //============================================================================

  gStyle -> SetOptStat(0);
  double PI = TMath::Pi();
  char INPUT_FILE_NAME[100];
  sprintf(INPUT_FILE_NAME,"../INPUT/sample_tree_for_closure_test.root");
  char INPUT_TREE_NAME[100];
  sprintf(INPUT_TREE_NAME,"MCTree");

  //============================================================================
  printf("Defining the weighting functions ... \n");
  //============================================================================

  double pol_par_TR[4] = {1000,1,0,0}; // TRANSVERSE POLARIZATION -> (+1,0,0)

  TF2 *func_W_HE_TR_pol = new TF2("func_W_HE_TR_pol",Func_W,-1,1,0,PI,4);
  func_W_HE_TR_pol -> SetParameters(pol_par_TR[0],pol_par_TR[1],pol_par_TR[2],pol_par_TR[3]);
  func_W_HE_TR_pol -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");

  double pol_par_LG[4] = {1000,-1,0,0}; // LONGITUDINAL POLARIZATION -> (-1,0,0)

  TF2 *func_W_HE_LG_pol = new TF2("func_W_HE_LG_pol",Func_W,-1,1,0,PI,4);
  func_W_HE_LG_pol -> SetParameters(pol_par_LG[0],pol_par_LG[1],pol_par_LG[2],pol_par_LG[3]);
  func_W_HE_LG_pol -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");

  //============================================================================
  printf("Inizializing the tree ... \n");
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
  Double_t CostHE_rec[3000], PhiHE_rec[3000], CostCS_rec[3000], PhiCS_rec[3000];

  input_tree -> SetBranchAddress("NDimu_rec",&NDimu_rec);
  input_tree -> SetBranchAddress("DimuPt_rec",DimuPt_rec);
  input_tree -> SetBranchAddress("DimuY_rec",DimuY_rec);
  input_tree -> SetBranchAddress("CostHE_rec",CostHE_rec);
  input_tree -> SetBranchAddress("PhiHE_rec",PhiHE_rec);
  input_tree -> SetBranchAddress("CostCS_rec",CostCS_rec);
  input_tree -> SetBranchAddress("PhiCS_rec",PhiCS_rec);

  //============================================================================
  printf("Defining the histos ... \n");
  //============================================================================

  TH2D *hCostPhiHE_2pt6_TR_pol_gen = new TH2D("hCostPhiHE_2pt6_TR_pol_gen","",100,-1,1,100,-PI,PI);
  hCostPhiHE_2pt6_TR_pol_gen -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hCostPhiHE_2pt6_TR_pol_gen -> GetYaxis() -> SetTitle("#phi_{HE}");
  TH2D *hCostPhiHE_2pt6_LG_pol_gen = new TH2D("hCostPhiHE_2pt6_LG_pol_gen","",100,-1,1,100,-PI,PI);
  hCostPhiHE_2pt6_LG_pol_gen -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hCostPhiHE_2pt6_LG_pol_gen -> GetYaxis() -> SetTitle("#phi_{HE}");

  TH2D *hCostPhiHE_2pt6_TR_pol_rec = new TH2D("hCostPhiHE_2pt6_TR_pol_rec","",100,-1,1,100,-PI,PI);
  hCostPhiHE_2pt6_TR_pol_rec -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hCostPhiHE_2pt6_TR_pol_rec -> GetYaxis() -> SetTitle("#phi_{HE}");
  TH2D *hCostPhiHE_2pt6_LG_pol_rec = new TH2D("hCostPhiHE_2pt6_LG_pol_rec","",100,-1,1,100,-PI,PI);
  hCostPhiHE_2pt6_LG_pol_rec -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hCostPhiHE_2pt6_LG_pol_rec -> GetYaxis() -> SetTitle("#phi_{HE}");

  //============================================================================
  printf("Reading the tree & weighting the histos ... \n");
  //============================================================================

  double weight_TR = 0, weight_LG = 0;

  for(int i = 0;i < input_tree -> GetEntries();i++){
    input_tree -> GetEntry(i);

    for(int k = 0;k < NDimu_gen;k++){

      if(DimuPt_gen[k] > 2 && DimuPt_gen[k] <= 6){
        weight_TR = func_W_HE_TR_pol -> Eval(CostHE_gen[k],PhiHE_gen[k])/func_W_HE_TR_pol -> GetMaximum();
        hCostPhiHE_2pt6_TR_pol_gen -> Fill(CostHE_gen[k],PhiHE_gen[k],weight_TR);
        weight_LG = func_W_HE_LG_pol -> Eval(CostHE_gen[k],PhiHE_gen[k])/func_W_HE_LG_pol -> GetMaximum();
        hCostPhiHE_2pt6_LG_pol_gen -> Fill(CostHE_gen[k],PhiHE_gen[k],weight_LG);
      }
    }

    for(int k = 0;k < NDimu_rec;k++){

      if(DimuPt_rec[k] > 2 && DimuPt_rec[k] <= 6){
        hCostPhiHE_2pt6_TR_pol_rec -> Fill(CostHE_rec[k],PhiHE_rec[k],weight_TR);
        hCostPhiHE_2pt6_LG_pol_rec -> Fill(CostHE_rec[k],PhiHE_rec[k],weight_LG);
      }
    }
  }

  //============================================================================
  // PLOTS FOR THE PRESENTATION
  //============================================================================

  char title[100];

  // TRANSVERSE POLARIZATION -> (+1,0,0)

  TPaveText *t_fit_HE_TR = new TPaveText(0.35,0.90,0.65,0.98,"brNDC");
  t_fit_HE_TR -> SetFillColor(kWhite);
  sprintf(title,"#lambda_{#theta} = %2.1f, #lambda_{#phi} = %2.1f, #lambda_{#theta#phi} = %2.1f",pol_par_TR[1],pol_par_TR[2],pol_par_TR[3]);
  t_fit_HE_TR -> AddText(title);

  TCanvas *c_CostPhiHE_2pt6_gen_TR = new TCanvas("c_CostPhiHE_2pt6_gen_TR","c_CostPhiHE_2pt6_gen_TR",4,132,1024,768);
  hCostPhiHE_2pt6_TR_pol_gen -> Draw("COLZ");
  t_fit_HE_TR -> Draw();

  TCanvas *c_CostPhiHE_2pt6_rec_TR = new TCanvas("c_CostPhiHE_2pt6_rec_TR","c_CostPhiHE_2pt6_rec_TR",4,132,1024,768);
  hCostPhiHE_2pt6_TR_pol_rec -> Draw("COLZ");
  t_fit_HE_TR -> Draw();

  // LONGITUDINAL POLARIZATION -> (-1,0,0)

  TPaveText *t_fit_HE_LG = new TPaveText(0.35,0.90,0.65,0.98,"brNDC");
  t_fit_HE_LG -> SetFillColor(kWhite);
  sprintf(title,"#lambda_{#theta} = %2.1f, #lambda_{#phi} = %2.1f, #lambda_{#theta#phi} = %2.1f",pol_par_LG[1],pol_par_LG[2],pol_par_LG[3]);
  t_fit_HE_LG -> AddText(title);

  TCanvas *c_CostPhiHE_2pt6_gen_LG = new TCanvas("c_CostPhiHE_2pt6_gen_LG","c_CostPhiHE_2pt6_gen_LG",4,132,1024,768);
  hCostPhiHE_2pt6_LG_pol_gen -> Draw("COLZ");
  t_fit_HE_LG -> Draw();

  TCanvas *c_CostPhiHE_2pt6_rec_LG = new TCanvas("c_CostPhiHE_2pt6_rec_LG","c_CostPhiHE_2pt6_rec_LG",4,132,1024,768);
  hCostPhiHE_2pt6_LG_pol_rec -> Draw("COLZ");
  t_fit_HE_LG -> Draw();

  //============================================================================
  printf("Writing the resulting histograms ... \n");

  char OUTPUT_FILE_NAME[100];
  sprintf(OUTPUT_FILE_NAME,"../OUTPUT/Histo_for_closure_test.root");

  TFile *output_file = new TFile(OUTPUT_FILE_NAME,"RECREATE");
  hCostPhiHE_2pt6_TR_pol_gen -> Write();
  hCostPhiHE_2pt6_TR_pol_rec -> Write();
  hCostPhiHE_2pt6_LG_pol_gen -> Write();
  hCostPhiHE_2pt6_LG_pol_rec -> Write();
  output_file -> Close();

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
