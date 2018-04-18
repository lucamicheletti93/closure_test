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

double Func_W(double *, double *);

void generate_polarized_distribution(){

  //============================================================================
  printf("1) Setting the binning ... \n");
  //============================================================================

  double diff = 0;

  const int N_cost_bins = 18;
  const int N_TH2_cost_bins = 100;
  const int dim_cost = N_cost_bins + 1;
  double bins_cost[dim_cost] = {0,25,30,35,40,42,44,46,48,50,52,54,56,58,60,65,70,75,100};
  double width_cost[dim_cost];
  double min_cost = -1.;
  double max_cost = 1.;
  double value_cost[dim_cost];
  diff = 0;

  for(int i = 0;i < N_cost_bins;i++){
    width_cost[i+1] = ((max_cost - min_cost)/(double) N_TH2_cost_bins)*(bins_cost[i+1] - bins_cost[i]);
  }

  for(int i = 0;i < dim_cost;i++){
    diff += width_cost[i];
    value_cost[i] = min_cost + diff;
    if(TMath::Abs(value_cost[i]) < 0.01){value_cost[i] = 0.;}
  }


  const int N_phi_bins = 10;
  const int N_TH2_phi_bins = 50;
  const int dim_phi = N_phi_bins + 1;
  double bins_phi[dim_phi] = {0,15,18,20,22,25,28,30,32,35,50};
  double width_phi[dim_phi];
  double min_phi = 0;
  double max_phi = TMath::Pi();
  double value_phi[dim_phi];

  diff = 0;

  for(int i = 0;i < N_phi_bins;i++){
    width_phi[i+1] = ((max_phi - min_phi)/(double) N_TH2_phi_bins)*(bins_phi[i+1] - bins_phi[i]);
  }

  for(int i = 0;i < dim_phi;i++){
    diff += width_phi[i];
    value_phi[i] = min_phi + diff;
    if(TMath::Abs(value_phi[i]) < 0.01){value_phi[i] = 0.;}
  }

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
  printf("2) Defining the weighting functions ... \n");
  //============================================================================

  double pol_par_TR[4] = {1000,1,0,0}; // TRANSVERSE POLARIZATION -> (+1,0,0)

  TF2 *func_W_HE_TR_pol = new TF2("func_W_HE_TR_pol",Func_W,-1,1,0,PI,4);
  func_W_HE_TR_pol -> SetParameters(pol_par_TR[0],pol_par_TR[1],pol_par_TR[2],pol_par_TR[3]);
  func_W_HE_TR_pol -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");

  double pol_par_LG[4] = {1000,-1,0,0}; // LONGITUDINAL POLARIZATION -> (-1,0,0)

  TF2 *func_W_HE_LG_pol = new TF2("func_W_HE_LG_pol",Func_W,-1,1,0,PI,4);
  func_W_HE_LG_pol -> SetParameters(pol_par_LG[0],pol_par_LG[1],pol_par_LG[2],pol_par_LG[3]);
  func_W_HE_LG_pol -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");

  double pol_par_NO[4] = {1000,-1,0,0}; // LONGITUDINAL POLARIZATION -> (-1,0,0)

  TF2 *func_W_HE_NO_pol = new TF2("func_W_HE_NO_pol",Func_W,-1,1,0,PI,4);
  func_W_HE_NO_pol -> SetParameters(pol_par_NO[0],pol_par_NO[1],pol_par_NO[2],pol_par_NO[3]);
  func_W_HE_NO_pol -> SetParNames("N","#lambda_{#theta}","#lambda_{#phi}","#lambda_{#theta#phi}");

  //============================================================================
  printf("3) Inizializing the tree ... \n");
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
  printf("4) Defining the histos ... \n");
  //============================================================================

  TH2D *hist_CostPhiHE_2pt6_TR_pol_gen = new TH2D("hist_CostPhiHE_2pt6_TR_pol_gen","",100,-1,1,50,0,PI);
  hist_CostPhiHE_2pt6_TR_pol_gen -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_CostPhiHE_2pt6_TR_pol_gen -> GetYaxis() -> SetTitle("#phi_{HE}");
  TH2D *hist_CostPhiHE_2pt6_LG_pol_gen = new TH2D("hist_CostPhiHE_2pt6_LG_pol_gen","",100,-1,1,50,0,PI);
  hist_CostPhiHE_2pt6_LG_pol_gen -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_CostPhiHE_2pt6_LG_pol_gen -> GetYaxis() -> SetTitle("#phi_{HE}");
  TH2D *hist_CostPhiHE_2pt6_NO_pol_gen = new TH2D("hist_CostPhiHE_2pt6_NO_pol_gen","",100,-1,1,50,0,PI);
  hist_CostPhiHE_2pt6_NO_pol_gen -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_CostPhiHE_2pt6_NO_pol_gen -> GetYaxis() -> SetTitle("#phi_{HE}");

  TH2D *hist_CostPhiHE_2pt6_TR_pol_gen_Rebin = new TH2D("hist_CostPhiHE_2pt6_TR_pol_gen_Rebin","",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_CostPhiHE_2pt6_TR_pol_gen_Rebin -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_CostPhiHE_2pt6_TR_pol_gen_Rebin -> GetYaxis() -> SetTitle("#phi_{HE}");
  TH2D *hist_CostPhiHE_2pt6_LG_pol_gen_Rebin = new TH2D("hist_CostPhiHE_2pt6_LG_pol_gen_Rebin","",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_CostPhiHE_2pt6_LG_pol_gen_Rebin -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_CostPhiHE_2pt6_LG_pol_gen_Rebin -> GetYaxis() -> SetTitle("#phi_{HE}");
  TH2D *hist_CostPhiHE_2pt6_NO_pol_gen_Rebin = new TH2D("hist_CostPhiHE_2pt6_NO_pol_gen_Rebin","",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_CostPhiHE_2pt6_NO_pol_gen_Rebin -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_CostPhiHE_2pt6_NO_pol_gen_Rebin -> GetYaxis() -> SetTitle("#phi_{HE}");

  TH2D *hist_CostPhiHE_2pt6_TR_pol_rec = new TH2D("hist_CostPhiHE_2pt6_TR_pol_rec","",100,-1,1,50,0,PI);
  hist_CostPhiHE_2pt6_TR_pol_rec -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_CostPhiHE_2pt6_TR_pol_rec -> GetYaxis() -> SetTitle("#phi_{HE}");
  TH2D *hist_CostPhiHE_2pt6_LG_pol_rec = new TH2D("hist_CostPhiHE_2pt6_LG_pol_rec","",100,-1,1,50,0,PI);
  hist_CostPhiHE_2pt6_LG_pol_rec -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_CostPhiHE_2pt6_LG_pol_rec -> GetYaxis() -> SetTitle("#phi_{HE}");
  TH2D *hist_CostPhiHE_2pt6_NO_pol_rec = new TH2D("hist_CostPhiHE_2pt6_NO_pol_rec","",100,-1,1,50,0,PI);
  hist_CostPhiHE_2pt6_NO_pol_rec -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_CostPhiHE_2pt6_NO_pol_rec -> GetYaxis() -> SetTitle("#phi_{HE}");

  TH2D *hist_CostPhiHE_2pt6_TR_pol_rec_Rebin = new TH2D("hist_CostPhiHE_2pt6_TR_pol_rec_Rebin","",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_CostPhiHE_2pt6_TR_pol_rec_Rebin -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_CostPhiHE_2pt6_TR_pol_rec_Rebin -> GetYaxis() -> SetTitle("#phi_{HE}");
  TH2D *hist_CostPhiHE_2pt6_LG_pol_rec_Rebin = new TH2D("hist_CostPhiHE_2pt6_LG_pol_rec_Rebin","",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_CostPhiHE_2pt6_LG_pol_rec_Rebin -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_CostPhiHE_2pt6_LG_pol_rec_Rebin -> GetYaxis() -> SetTitle("#phi_{HE}");
  TH2D *hist_CostPhiHE_2pt6_NO_pol_rec_Rebin = new TH2D("hist_CostPhiHE_2pt6_NO_pol_rec_Rebin","",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_CostPhiHE_2pt6_NO_pol_rec_Rebin -> GetXaxis() -> SetTitle("cos#theta_{HE}");
  hist_CostPhiHE_2pt6_NO_pol_rec_Rebin -> GetYaxis() -> SetTitle("#phi_{HE}");

  //============================================================================
  printf("5) Reading the tree & weighting the histos ... \n");
  //============================================================================

  double weight_TR = 0, weight_LG = 0;

  for(int i = 0;i < input_tree -> GetEntries();i++){
    input_tree -> GetEntry(i);

    for(int k = 0;k < NDimu_gen;k++){

      if(DimuPt_gen[k] > 2 && DimuPt_gen[k] <= 6){
        weight_TR = func_W_HE_TR_pol -> Eval(CostHE_gen[k],TMath::Abs(PhiHE_gen[k]))/func_W_HE_TR_pol -> GetMaximum();
        hist_CostPhiHE_2pt6_TR_pol_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k]),weight_TR);
        hist_CostPhiHE_2pt6_TR_pol_gen_Rebin -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k]),weight_TR);
        weight_LG = func_W_HE_LG_pol -> Eval(CostHE_gen[k],TMath::Abs(PhiHE_gen[k]))/func_W_HE_LG_pol -> GetMaximum();
        hist_CostPhiHE_2pt6_LG_pol_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k]),weight_LG);
        hist_CostPhiHE_2pt6_LG_pol_gen_Rebin -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k]),weight_LG);

        hist_CostPhiHE_2pt6_NO_pol_gen -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k]));
        hist_CostPhiHE_2pt6_NO_pol_gen_Rebin -> Fill(CostHE_gen[k],TMath::Abs(PhiHE_gen[k]));
      }
    }

    for(int k = 0;k < NDimu_rec;k++){

      if(DimuPt_rec[k] > 2 && DimuPt_rec[k] <= 6){
        if(DimuY_rec[k] > -4. && DimuY_rec[k] < -2.5){
          if(DimuMatch_rec[k] == 2){
            if(DimuMass_rec[k] > 2 && DimuMass_rec[k] < 5){
              hist_CostPhiHE_2pt6_TR_pol_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),weight_TR);
              hist_CostPhiHE_2pt6_TR_pol_rec_Rebin -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),weight_TR);
              hist_CostPhiHE_2pt6_LG_pol_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),weight_LG);
              hist_CostPhiHE_2pt6_LG_pol_rec_Rebin -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]),weight_LG);

              hist_CostPhiHE_2pt6_NO_pol_rec -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]));
              hist_CostPhiHE_2pt6_NO_pol_rec_Rebin -> Fill(CostHE_rec[k],TMath::Abs(PhiHE_rec[k]));
            }
          }
        }
      }
    }
  }

  //============================================================================
  printf("6) Computing Acc x Eff ... \n");
  //============================================================================

  TH2D *hist_accxeff_2pt6_TR_pol = new TH2D("hist_accxeff_2pt6_TR_pol","",100,-1,1,50,0,PI);
  hist_accxeff_2pt6_TR_pol -> Divide(hist_CostPhiHE_2pt6_TR_pol_rec,hist_CostPhiHE_2pt6_TR_pol_gen,1,1,"B");

  TH2D *hist_accxeff_2pt6_TR_pol_Rebin = new TH2D("hist_accxeff_2pt6_TR_pol_Rebin","",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_accxeff_2pt6_TR_pol_Rebin -> Divide(hist_CostPhiHE_2pt6_TR_pol_rec_Rebin,hist_CostPhiHE_2pt6_TR_pol_gen_Rebin,1,1,"B");

  TH2D *hist_accxeff_2pt6_LG_pol = new TH2D("hist_accxeff_2pt6_LG_pol","",100,-1,1,50,0,PI);
  hist_accxeff_2pt6_LG_pol -> Divide(hist_CostPhiHE_2pt6_LG_pol_rec,hist_CostPhiHE_2pt6_LG_pol_gen,1,1,"B");

  TH2D *hist_accxeff_2pt6_LG_pol_Rebin = new TH2D("hist_accxeff_2pt6_LG_pol_Rebin","",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_accxeff_2pt6_LG_pol_Rebin -> Divide(hist_CostPhiHE_2pt6_LG_pol_rec_Rebin,hist_CostPhiHE_2pt6_LG_pol_gen_Rebin,1,1,"B");

  TH2D *hist_accxeff_2pt6_NO_pol = new TH2D("hist_accxeff_2pt6_NO_pol","",100,-1,1,50,0,PI);
  hist_accxeff_2pt6_NO_pol -> Divide(hist_CostPhiHE_2pt6_NO_pol_rec,hist_CostPhiHE_2pt6_NO_pol_gen,1,1,"B");

  TH2D *hist_accxeff_2pt6_NO_pol_Rebin = new TH2D("hist_accxeff_2pt6_NO_pol_Rebin","",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_accxeff_2pt6_NO_pol_Rebin -> Divide(hist_CostPhiHE_2pt6_NO_pol_rec_Rebin,hist_CostPhiHE_2pt6_NO_pol_gen_Rebin,1,1,"B");

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
  hist_CostPhiHE_2pt6_TR_pol_gen -> Draw("COLZ");
  t_fit_HE_TR -> Draw();

  TCanvas *c_CostPhiHE_2pt6_rec_TR = new TCanvas("c_CostPhiHE_2pt6_rec_TR","c_CostPhiHE_2pt6_rec_TR",4,132,1024,768);
  hist_CostPhiHE_2pt6_TR_pol_rec -> Draw("COLZ");
  t_fit_HE_TR -> Draw();

  TCanvas *c_accxeff_2pt6_TR_pol = new TCanvas("c_accxeff_2pt6_TR_pol","c_accxeff_2pt6_TR_pol",4,132,1024,768);
  hist_accxeff_2pt6_TR_pol -> Draw("COLZ");
  t_fit_HE_TR -> Draw();

  // LONGITUDINAL POLARIZATION -> (-1,0,0)

  TPaveText *t_fit_HE_LG = new TPaveText(0.35,0.90,0.65,0.98,"brNDC");
  t_fit_HE_LG -> SetFillColor(kWhite);
  sprintf(title,"#lambda_{#theta} = %2.1f, #lambda_{#phi} = %2.1f, #lambda_{#theta#phi} = %2.1f",pol_par_LG[1],pol_par_LG[2],pol_par_LG[3]);
  t_fit_HE_LG -> AddText(title);

  TCanvas *c_CostPhiHE_2pt6_gen_LG = new TCanvas("c_CostPhiHE_2pt6_gen_LG","c_CostPhiHE_2pt6_gen_LG",4,132,1024,768);
  hist_CostPhiHE_2pt6_LG_pol_gen -> Draw("COLZ");
  t_fit_HE_LG -> Draw();

  TCanvas *c_CostPhiHE_2pt6_rec_LG = new TCanvas("c_CostPhiHE_2pt6_rec_LG","c_CostPhiHE_2pt6_rec_LG",4,132,1024,768);
  hist_CostPhiHE_2pt6_LG_pol_rec -> Draw("COLZ");
  t_fit_HE_LG -> Draw();

  TCanvas *c_accxeff_2pt6_LG_pol = new TCanvas("c_accxeff_2pt6_LG_pol","c_accxeff_2pt6_LG_pol",4,132,1024,768);
  hist_accxeff_2pt6_LG_pol -> Draw("COLZ");
  t_fit_HE_LG -> Draw();

  // NO POLARIZATION -> (0,0,0)

  TPaveText *t_fit_HE_NO = new TPaveText(0.35,0.90,0.65,0.98,"brNDC");
  t_fit_HE_NO -> SetFillColor(kWhite);
  sprintf(title,"#lambda_{#theta} = %2.1f, #lambda_{#phi} = %2.1f, #lambda_{#theta#phi} = %2.1f",pol_par_NO[1],pol_par_NO[2],pol_par_NO[3]);
  t_fit_HE_NO -> AddText(title);

  TCanvas *c_CostPhiHE_2pt6_gen_NO = new TCanvas("c_CostPhiHE_2pt6_gen_NO","c_CostPhiHE_2pt6_gen_NO",4,132,1024,768);
  hist_CostPhiHE_2pt6_NO_pol_gen -> Draw("COLZ");
  t_fit_HE_NO -> Draw();

  TCanvas *c_CostPhiHE_2pt6_rec_NO = new TCanvas("c_CostPhiHE_2pt6_rec_NO","c_CostPhiHE_2pt6_rec_NO",4,132,1024,768);
  hist_CostPhiHE_2pt6_NO_pol_rec -> Draw("COLZ");
  t_fit_HE_NO -> Draw();

  TCanvas *c_accxeff_2pt6_NO_pol = new TCanvas("c_accxeff_2pt6_NO_pol","c_accxeff_2pt6_NO_pol",4,132,1024,768);
  hist_accxeff_2pt6_NO_pol -> Draw("COLZ");
  t_fit_HE_NO -> Draw();

  //============================================================================
  // PLOTS FOR THE CLOSURE TEST
  //============================================================================

  // TRANSVERSE POLARIZATION -> (+1,0,0)

  TPaveText *t_fit_HE_TR_Rebin = new TPaveText(0.35,0.90,0.65,0.98,"brNDC");
  t_fit_HE_TR_Rebin -> SetFillColor(kWhite);
  sprintf(title,"#lambda_{#theta} = %2.1f, #lambda_{#phi} = %2.1f, #lambda_{#theta#phi} = %2.1f",pol_par_TR[1],pol_par_TR[2],pol_par_TR[3]);
  t_fit_HE_TR_Rebin -> AddText(title);

  TCanvas *c_CostPhiHE_2pt6_gen_TR_Rebin = new TCanvas("c_CostPhiHE_2pt6_gen_TR_Rebin","c_CostPhiHE_2pt6_gen_TR_Rebin",4,132,1024,768);
  hist_CostPhiHE_2pt6_TR_pol_gen_Rebin -> Draw("COLZtext");
  t_fit_HE_TR_Rebin -> Draw();

  TCanvas *c_CostPhiHE_2pt6_rec_TR_Rebin = new TCanvas("c_CostPhiHE_2pt6_rec_TR_Rebin","c_CostPhiHE_2pt6_rec_TR_Rebin",4,132,1024,768);
  hist_CostPhiHE_2pt6_TR_pol_rec_Rebin -> Draw("COLZtext");
  t_fit_HE_TR_Rebin-> Draw();

  TCanvas *c_accxeff_2pt6_TR_pol_Rebin = new TCanvas("c_accxeff_2pt6_TR_pol_Rebin","c_accxeff_2pt6_TR_pol_Rebin",4,132,1024,768);
  hist_accxeff_2pt6_TR_pol_Rebin -> Draw("COLZtext");
  t_fit_HE_TR_Rebin -> Draw();

  // LONGITUDINAL POLARIZATION -> (-1,0,0)

  TPaveText *t_fit_HE_LG_Rebin = new TPaveText(0.35,0.90,0.65,0.98,"brNDC");
  t_fit_HE_LG_Rebin -> SetFillColor(kWhite);
  sprintf(title,"#lambda_{#theta} = %2.1f, #lambda_{#phi} = %2.1f, #lambda_{#theta#phi} = %2.1f",pol_par_LG[1],pol_par_LG[2],pol_par_LG[3]);
  t_fit_HE_LG_Rebin -> AddText(title);

  TCanvas *c_CostPhiHE_2pt6_gen_LG_Rebin = new TCanvas("c_CostPhiHE_2pt6_gen_LG_Rebin","c_CostPhiHE_2pt6_gen_LG_Rebin",4,132,1024,768);
  hist_CostPhiHE_2pt6_LG_pol_gen_Rebin -> Draw("COLZtext");
  t_fit_HE_LG_Rebin -> Draw();

  TCanvas *c_CostPhiHE_2pt6_rec_LG_Rebin = new TCanvas("c_CostPhiHE_2pt6_rec_LG_Rebin","c_CostPhiHE_2pt6_rec_LG_Rebin",4,132,1024,768);
  hist_CostPhiHE_2pt6_LG_pol_rec_Rebin -> Draw("COLZtext");
  t_fit_HE_LG_Rebin-> Draw();

  TCanvas *c_accxeff_2pt6_LG_pol_Rebin = new TCanvas("c_accxeff_2pt6_LG_pol_Rebin","c_accxeff_2pt6_LG_pol_Rebin",4,132,1024,768);
  hist_accxeff_2pt6_LG_pol_Rebin -> Draw("COLZtext");
  t_fit_HE_LG_Rebin -> Draw();

  // NO POLARIZATION -> (0,0,0)

  TPaveText *t_fit_HE_NO_Rebin = new TPaveText(0.35,0.90,0.65,0.98,"brNDC");
  t_fit_HE_NO_Rebin -> SetFillColor(kWhite);
  sprintf(title,"#lambda_{#theta} = %2.1f, #lambda_{#phi} = %2.1f, #lambda_{#theta#phi} = %2.1f",pol_par_NO[1],pol_par_NO[2],pol_par_NO[3]);
  t_fit_HE_NO_Rebin -> AddText(title);

  TCanvas *c_CostPhiHE_2pt6_gen_NO_Rebin = new TCanvas("c_CostPhiHE_2pt6_gen_NO_Rebin","c_CostPhiHE_2pt6_gen_NO_Rebin",4,132,1024,768);
  hist_CostPhiHE_2pt6_NO_pol_gen_Rebin -> Draw("COLZtext");
  t_fit_HE_NO_Rebin -> Draw();

  TCanvas *c_CostPhiHE_2pt6_rec_NO_Rebin = new TCanvas("c_CostPhiHE_2pt6_rec_NO_Rebin","c_CostPhiHE_2pt6_rec_NO_Rebin",4,132,1024,768);
  hist_CostPhiHE_2pt6_NO_pol_rec_Rebin -> Draw("COLZtext");
  t_fit_HE_NO_Rebin-> Draw();

  TCanvas *c_accxeff_2pt6_NO_pol_Rebin = new TCanvas("c_accxeff_2pt6_NO_pol_Rebin","c_accxeff_2pt6_NO_pol_Rebin",4,132,1024,768);
  hist_accxeff_2pt6_NO_pol_Rebin -> Draw("COLZtext");
  t_fit_HE_NO_Rebin -> Draw();

  //============================================================================
  printf("7) Writing the resulting histograms ... \n");
  //============================================================================

  char OUTPUT_FILE_NAME[100];
  sprintf(OUTPUT_FILE_NAME,"../OUTPUT/Histo_for_closure_test.root");

  TFile *output_file = new TFile(OUTPUT_FILE_NAME,"RECREATE");
  hist_CostPhiHE_2pt6_TR_pol_gen -> Write();
  hist_CostPhiHE_2pt6_TR_pol_rec -> Write();
  hist_accxeff_2pt6_TR_pol -> Write();
  hist_CostPhiHE_2pt6_LG_pol_gen -> Write();
  hist_CostPhiHE_2pt6_LG_pol_rec -> Write();
  hist_accxeff_2pt6_LG_pol -> Write();
  hist_CostPhiHE_2pt6_NO_pol_gen -> Write();
  hist_CostPhiHE_2pt6_NO_pol_rec -> Write();
  hist_accxeff_2pt6_NO_pol -> Write();
  hist_CostPhiHE_2pt6_TR_pol_gen_Rebin -> Write();
  hist_CostPhiHE_2pt6_TR_pol_rec_Rebin -> Write();
  hist_accxeff_2pt6_TR_pol_Rebin -> Write();
  hist_CostPhiHE_2pt6_LG_pol_gen_Rebin -> Write();
  hist_CostPhiHE_2pt6_LG_pol_rec_Rebin -> Write();
  hist_accxeff_2pt6_LG_pol_Rebin -> Write();
  hist_CostPhiHE_2pt6_NO_pol_gen_Rebin -> Write();
  hist_CostPhiHE_2pt6_NO_pol_rec_Rebin -> Write();
  hist_accxeff_2pt6_NO_pol_Rebin -> Write();
  output_file -> Close();

  //============================================================================
  printf("7) Saving the resulting canvas ... \n");
  //============================================================================

  c_CostPhiHE_2pt6_gen_TR -> SaveAs("../OUTPUT/PLOT/c_CostPhiHE_2pt6_gen_TR.png");
  c_CostPhiHE_2pt6_rec_TR -> SaveAs("../OUTPUT/PLOT/c_CostPhiHE_2pt6_rec_TR.png");
  c_accxeff_2pt6_TR_pol -> SaveAs("../OUTPUT/PLOT/c_accxeff_2pt6_TR_pol.png");
  c_CostPhiHE_2pt6_gen_LG -> SaveAs("../OUTPUT/PLOT/c_CostPhiHE_2pt6_gen_LG.png");
  c_CostPhiHE_2pt6_rec_LG -> SaveAs("../OUTPUT/PLOT/c_CostPhiHE_2pt6_rec_LG.png");
  c_accxeff_2pt6_LG_pol -> SaveAs("../OUTPUT/PLOT/c_accxeff_2pt6_LG_pol.png");
  c_CostPhiHE_2pt6_gen_NO -> SaveAs("../OUTPUT/PLOT/c_CostPhiHE_2pt6_gen_NO.png");
  c_CostPhiHE_2pt6_rec_NO -> SaveAs("../OUTPUT/PLOT/c_CostPhiHE_2pt6_rec_NO.png");
  c_accxeff_2pt6_NO_pol -> SaveAs("../OUTPUT/PLOT/c_accxeff_2pt6_NO_pol.png");
  c_CostPhiHE_2pt6_gen_TR_Rebin -> SaveAs("../OUTPUT/PLOT/c_CostPhiHE_2pt6_gen_TR_Rebin.png");
  c_CostPhiHE_2pt6_rec_TR_Rebin -> SaveAs("../OUTPUT/PLOT/c_CostPhiHE_2pt6_rec_TR_Rebin.png");
  c_accxeff_2pt6_TR_pol_Rebin -> SaveAs("../OUTPUT/PLOT/c_accxeff_2pt6_TR_pol_Rebin.png");
  c_CostPhiHE_2pt6_gen_LG_Rebin -> SaveAs("../OUTPUT/PLOT/c_CostPhiHE_2pt6_gen_LG_Rebin.png");
  c_CostPhiHE_2pt6_rec_LG_Rebin -> SaveAs("../OUTPUT/PLOT/c_CostPhiHE_2pt6_rec_LG_Rebin.png");
  c_accxeff_2pt6_LG_pol_Rebin -> SaveAs("../OUTPUT/PLOT/c_CostPhiHE_2pt6_gen_TR.png");
  c_CostPhiHE_2pt6_gen_NO_Rebin -> SaveAs("../OUTPUT/PLOT/c_CostPhiHE_2pt6_gen_NO_Rebin.png");
  c_CostPhiHE_2pt6_rec_NO_Rebin -> SaveAs("../OUTPUT/PLOT/c_CostPhiHE_2pt6_rec_NO_Rebin.png");
  c_accxeff_2pt6_NO_pol_Rebin -> SaveAs("../OUTPUT/PLOT/c_accxeff_2pt6_NO_pol_Rebin.png");

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
