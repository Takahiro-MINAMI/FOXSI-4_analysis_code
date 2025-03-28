//
// Gap loss correction
// _base_cal.root -> _base_cal_gapmod.root
//

#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <fstream>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TSpline.h>
#include <TRandom3.h>
#include <TLine.h>

using namespace std;

int mod_gaploss(const char *inFileName, const char *calFileName, bool save_all = false);

int main(int argc, char **argv)
{

  if (argc == 3)
    mod_gaploss(argv[1], argv[2]);
  else if (argc == 4)
  {
    cout << !strncmp(argv[3], "-all", 4) << endl;
    mod_gaploss(argv[1], argv[2], !strncmp(argv[3], "-all", 4));
  }
  else
  {
    cout << "USAGE: Please input the names of the file to analyze and gap_loss response file" << endl;
    cout << "./path3_mod_gaploss 'input filename' 'Gap-loss response filename'" << endl;
    return -1;
  }
  return 0;
}

int mod_gaploss(const char *inFileName, const char *calFileName, bool save_all)
{
  gStyle->SetOptStat(0);
  gStyle->SetPalette(53);

  cout << "-- Path3: Gap-Loss Correction ------------" << endl;

  // Load Gap-Loss Response function (ROOT, TGraph2D) -------------------------------------
  TFile *fcal = new TFile(calFileName);
  TGraph2D *resp_pt[3];
  TGraph2D *resp_al[3];
  resp_pt[0] = (TGraph2D *)fcal->Get("resp_pt_region0");
  resp_pt[1] = (TGraph2D *)fcal->Get("resp_pt_region1");
  resp_pt[2] = (TGraph2D *)fcal->Get("resp_pt_region2");
  resp_al[0] = (TGraph2D *)fcal->Get("resp_al_region0");
  resp_al[1] = (TGraph2D *)fcal->Get("resp_al_region1");
  resp_al[2] = (TGraph2D *)fcal->Get("resp_al_region2");

  // Input and Set TTree Branch Address -------------------------------------
  TString infilename = TString(inFileName);
  if (!infilename.EndsWith("cal.root"))
  {
    cout << "Error: Input file " << inFileName << " is NOT Cal root file!" << endl;
    return -1;
  }
  double ene_Al_cor = 1.0;
  const char *match = "2022_01";
  if (strstr(calFileName, match) != NULL)
  {
    cout << "Match" << endl;
    ene_Al_cor = 1.034;
  }
  else
  {
    cout << "not Match  " << match << " " << inFileName << endl;
  }

  cout << "  " << endl;
  cout << "Open ROOT file :" << infilename << endl;

  UInt_t ti;
  UInt_t livetime;
  UInt_t integral_livetime;
  UInt_t unixtime;

  Int_t pt_nhit;
  Int_t al_nhit;
  Double_t epi_pt[128] = {0.};
  Double_t pos_pt[128] = {0.0};
  Double_t epi_al[128] = {0.};
  Double_t pos_al[128] = {0.0};

  Double_t merged_epi_pt[128] = {0.};
  Double_t merged_pos_pt[128] = {0.};
  Int_t merged_pt_nhit = 0;
  Int_t merged_pt_nhit_list[128] = {0};
  Double_t merged_epi_al[128] = {0.};
  Double_t merged_pos_al[128] = {0.};
  Int_t merged_al_nhit = 0;
  Int_t merged_al_nhit_list[128] = {0};
  Int_t merged_pt_region[128] = {0};
  Int_t merged_al_region[128] = {0};

  UShort_t flag_pseudo;
  UInt_t pseudo_counter;
  UInt_t survive_pseudo_counter = 0;

  TFile *fin = new TFile(infilename);
  TTree *tin = (TTree *)fin->Get("eventtree");
  TTree *tin_hktree = (TTree *)fin->Get("hktree");
  TTree *tin_threshold_tree = (TTree *)fin->Get("threshold");

  tin->SetBranchAddress("ti", &ti);
  tin->SetBranchAddress("livetime", &livetime);
  tin->SetBranchAddress("integral_livetime", &integral_livetime);
  tin->SetBranchAddress("unixtime", &unixtime);
  tin->SetBranchAddress("flag_pseudo", &flag_pseudo);
  tin->SetBranchAddress("pseudo_counter", &pseudo_counter);
  tin->SetBranchAddress("survive_pseudo_counter", &survive_pseudo_counter);

  tin->SetBranchAddress("pt_nhit", &pt_nhit);
  tin->SetBranchAddress("epi_pt", epi_pt);
  tin->SetBranchAddress("pos_pt", pos_pt);
  tin->SetBranchAddress("merged_pt_nhit", &merged_pt_nhit);
  tin->SetBranchAddress("merged_epi_pt", merged_epi_pt);
  tin->SetBranchAddress("merged_pos_pt", merged_pos_pt);
  tin->SetBranchAddress("merged_pt_nhit_list", merged_pt_nhit_list);

  tin->SetBranchAddress("al_nhit", &al_nhit);
  tin->SetBranchAddress("epi_al", epi_al);
  tin->SetBranchAddress("pos_al", pos_al);
  tin->SetBranchAddress("merged_al_nhit", &merged_al_nhit);
  tin->SetBranchAddress("merged_epi_al", merged_epi_al);
  tin->SetBranchAddress("merged_pos_al", merged_pos_al);
  tin->SetBranchAddress("merged_al_nhit_list", merged_al_nhit_list);

  tin->SetBranchAddress("merged_pt_region", merged_pt_region);
  tin->SetBranchAddress("merged_al_region", merged_al_region);

  // Save org adc information or not -------------------------------------
  if (!save_all)
  {
    tin->SetBranchStatus("cmn*", 0);
    tin->SetBranchStatus("adc0", 0);
    tin->SetBranchStatus("adc1", 0);
    tin->SetBranchStatus("adc2", 0);
    tin->SetBranchStatus("adc3", 0);
    tin->SetBranchStatus("index*", 0);
    tin->SetBranchStatus("hitnum*", 0);
  }

  // Create new ROOT file and TTree  -------------------------------------
  TString outFilename;
  outFilename = inFileName;
  outFilename.ReplaceAll(".root", "_gapmod.root");

  TFile *fout = new TFile(outFilename, "RECREATE");
  TTree *tout = tin->CloneTree(0);
  Double_t gapmod_epi_pt[128] = {0.};
  Double_t gapmod_epi_al[128] = {0.};
  Int_t event_quality_flag = 0;
  tout->Branch("gapmod_epi_pt", gapmod_epi_pt, "gapmod_epi_pt[merged_pt_nhit]/D");
  tout->Branch("gapmod_epi_al", gapmod_epi_al, "gapmod_epi_al[merged_al_nhit]/D");
  tout->Branch("event_quality_flag", &event_quality_flag, "event_quality_flag/I");

  TH1D *histpt_single[3];
  TH1D *histpt_double[3];
  TH1D *histpt_double_mod[3];
  TH2D *hist_diffene_pt_adj[3];
  TH2D *hist_diffene_pt_adj_mod[3];
  for (int j = 0; j < 3; j++)
  {
    histpt_single[j] = new TH1D(Form("histpt_single_region%d", j), Form("Pt, Single, Region%d;Energy [keV];Counts", j), 1100, -10, 100);
    histpt_double[j] = new TH1D(Form("histpt_double_region%d", j), Form("Pt, Double, Region%d;Energy [keV];Counts", j), 1100, -10, 100);
    histpt_double_mod[j] = new TH1D(Form("histpt_double_region%d_mod", j), Form("Pt, Double, Region%d, Gap-Loss Mod.;Energy [keV];Counts", j), 1100, -10, 100);
    hist_diffene_pt_adj[j] = new TH2D(Form("hist_diffene_pt_adj_region%d", j), Form("Pt, Adj. ch, Ediff vs Eave, Region%d;Ediff [keV];Eave [keV]", j), 500, -50, 50, 500, 0, 100);
    hist_diffene_pt_adj_mod[j] = new TH2D(Form("hist_diffene_pt_adj_region%d_mod", j), Form("Pt, Adj. ch, Ediff vs Emod, Region%d, Gap-Loss Mod.;Ediff [keV];Emod [keV]", j), 500, -50, 50, 500, 0, 100);
  }

  TH1D *histal_single[3];
  TH1D *histal_double[3];
  TH1D *histal_double_mod[3];
  TH2D *hist_diffene_al_adj[3];
  TH2D *hist_diffene_al_adj_mod[3];
  for (int j = 0; j < 3; j++)
  {
    histal_single[j] = new TH1D(Form("histal_single_region%d", j), Form("Al, Single, Region%d;Energy[keV];Counts", j), 1100, -10, 100);
    histal_double[j] = new TH1D(Form("histal_double_region%d", j), Form("Al, Double, Region%d;Energy[keV];Counts", j), 1100, -10, 100);
    histal_double_mod[j] = new TH1D(Form("histal_double_region%d_mod", j), Form("Al, Double, Region%d, Gap-Loss Mod.;Energy [keV];Counts", j), 1100, -10, 100);
    hist_diffene_al_adj[j] = new TH2D(Form("hist_diffene_al_adj_region%d", j), Form("Al, Adj. ch, Ediff vs Eave, Region%d;Ediff [keV];Eave [keV]", j), 500, -50, 50, 500, 0, 100);
    hist_diffene_al_adj_mod[j] = new TH2D(Form("hist_diffene_al_adj_region%d_mod", j), Form("Al, Adj. ch, Ediff vs Emod, Region%d, Gap-Loss Mod.;Ediff [keV];Emod [keV]", j), 500, -50, 50, 500, 0, 100);
  }

  TH1D *histpt_single_all = new TH1D("histpt_single_all", "Pt, Single, Allch;Energy [keV];Counts", 1100, -10, 100);
  TH1D *histpt_double_all = new TH1D("histpt_double_all", "Pt, Double, Allch;Energy [keV];Counts", 1100, -10, 100);
  TH1D *histpt_all = new TH1D("histpt_all", "Pt, Single+Double, Allch;Energy [keV];Counts", 1100, -10, 100);
  TH1D *histsingle_onlygood = new TH1D("histsingle_onlygood", "Ediff<0.3 keV, Single;Energy [keV];Counts", 1100, -10, 100);

  TH1D *histal_single_all = new TH1D("histal_single_all", "Al, Single, Allch;Energy [keV];Counts", 1100, -10, 100);
  TH1D *histal_double_all = new TH1D("histal_double_all", "Al, Double, Allch;Energy [keV];Counts", 1100, -10, 100);
  TH1D *histal_all = new TH1D("histal_all", "Al, Single+Double, Allch;Energy [keV];Counts", 1100, -10, 100);

  TH1D *histsum_all = new TH1D("histsum_all", "Allch, Single+Double, Gap-Loss Mod.", 1100, -10, 100);
  TH2D *histall_single = new TH2D("histall_single", "Allch, Single strip (Sum)", 256, -0.5, 255.5, 1000, 0, 100);
  TH2D *histall_double = new TH2D("histall_double", "Allch, Double strip (Sum)", 256, -0.5, 255.5, 1000, 0, 100);

  TH2D *hist_avediff[2][2];
  TH2D *hist_avediff_rot[2][2];
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      hist_avediff[i][j] = new TH2D(Form("hist_avediff_pt%dal%d", i + 1, j + 1), Form("Pt/Al, Eave. vs Ediff, Hit Pt%d Al%d;Eave [keV];Ediff [keV]", i + 1, j + 1), 500 * 2, 0, 100, 500 * 2, -50, 50);
      hist_avediff_rot[i][j] = new TH2D(Form("hist_avediff_rot_pt%dal%d", i + 1, j + 1), Form("Pt/Al, Ediff. vs Eave, Hit Pt%d Al%d;Ediff [keV];Eave [keV]", i + 1, j + 1), 500 * 2, -50, 50, 500 * 2, 0, 100);
    }
  }
  TH1D *hist_event_quality = new TH1D("hist_event_quality", "Event Quality Flag", 10, -0.5, 9.5);

  Long64_t nentries = tin->GetEntriesFast();

  // Main Process  -------------------------------------
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    if (tin->GetEntry(jentry) < 0)
    {
      cerr << "ERROR: error occurs while loading tree" << endl;
      break;
    }
    if (jentry % 100000 == 0)
      cerr << "Proccessed " << jentry << " events(" << jentry * 100. / nentries << "\%)" << endl;

    if (flag_pseudo == 1)
      continue;

    event_quality_flag = 0;

    if (merged_pt_nhit == 1 && merged_al_nhit == 1)
    {
      // Exclude guard-ring events
      if (merged_pos_pt[0] > 120 || merged_pos_pt[0] < 7 || merged_pos_al[0] < 134 - 128 || merged_pos_al[0] > 249 - 128)
      {
        event_quality_flag = -1;
        continue;
      }
      if (merged_pt_region[0] < 0 || merged_al_region[0] < 0)
      {
        event_quality_flag = -1;
        continue;
      }

      // Gap-loss correction // need discussio
      event_quality_flag = 1;
      if (pt_nhit == 1)
        gapmod_epi_pt[0] = merged_epi_pt[0];
      if (pt_nhit == 2)
        gapmod_epi_pt[0] = resp_pt[merged_pt_region[0]]->Interpolate(epi_pt[1] - epi_pt[0], epi_pt[0] + epi_pt[1]);
      if (al_nhit == 1)
        gapmod_epi_al[0] = merged_epi_al[0];
      if (al_nhit == 2)
        gapmod_epi_al[0] = resp_al[merged_al_region[0]]->Interpolate(epi_al[1] - epi_al[0], epi_al[0] + epi_al[1]);

      if (pt_nhit > 2)
      {
        gapmod_epi_pt[0] = merged_epi_pt[0];
        event_quality_flag = 2;
      }
      if (al_nhit > 2)
      {
        gapmod_epi_al[0] = merged_epi_al[0];
        event_quality_flag = 2;
      }
      if (gapmod_epi_al[0] < 1.5 || gapmod_epi_pt[0] < 1.5 || gapmod_epi_al[0] > 32.0 || gapmod_epi_pt[0] > 32.0)
        event_quality_flag = 5;
    }
    else
    {
      int pt_out_index = 0;
      for (int i = 0; i < merged_pt_nhit; i++)
      {
        if (merged_pos_pt[i] > 120 || merged_pos_pt[i] < 7 || merged_pt_region[i] < 0)
        {
          event_quality_flag = -1;
          continue;
        }
        if (merged_pt_nhit_list[i] == 1)
        {
          gapmod_epi_pt[i] = merged_epi_pt[i];
          pt_out_index = pt_out_index + 1;
        }
        if (merged_pt_nhit_list[i] == 2)
        {
          gapmod_epi_pt[i] = resp_pt[merged_pt_region[i]]->Interpolate(epi_pt[pt_out_index + 1] - epi_pt[pt_out_index], epi_pt[pt_out_index] + epi_pt[pt_out_index + 1]);
          pt_out_index = pt_out_index + 2;
        }
        if (merged_pt_nhit_list[i] > 2)
        {
          gapmod_epi_pt[i] = merged_epi_pt[i];
          pt_out_index = pt_out_index + merged_pt_nhit_list[i];
        }
      }
      int al_out_index = 0;
      for (int i = 0; i < merged_al_nhit; i++)
      {
        if (merged_pos_al[i] < 134 - 128 || merged_pos_al[i] > 249 - 128 || merged_al_region[i] < 0)
        {
          event_quality_flag = -1;
          continue;
        }
        if (merged_al_nhit_list[i] == 1)
        {
          gapmod_epi_al[i] = merged_epi_al[i];
          al_out_index = al_out_index + 1;
        }
        if (merged_al_nhit_list[i] == 2)
        {
          gapmod_epi_al[i] = resp_al[merged_al_region[i]]->Interpolate(epi_al[al_out_index + 1] - epi_al[al_out_index], epi_al[al_out_index] + epi_al[al_out_index + 1]);
          al_out_index = al_out_index + 2;
        }
        if (merged_al_nhit_list[i] > 2)
        {
          gapmod_epi_al[i] = merged_epi_al[i];
          al_out_index = al_out_index + merged_al_nhit_list[i];
        }
      }
      if (merged_pt_nhit == merged_al_nhit && event_quality_flag != -1)
        event_quality_flag = 3;
      if (merged_pt_nhit != merged_al_nhit && event_quality_flag != -1)
        event_quality_flag = 4;
    }

    gapmod_epi_al[0] = gapmod_epi_al[0] / ene_Al_cor;

    // for drawing
    if (event_quality_flag == 1)
    {
      if (pt_nhit == 1)
      {
        histpt_single[merged_pt_region[0]]->Fill(epi_pt[0]);
        histpt_single_all->Fill(epi_pt[0]);
        histpt_all->Fill(epi_pt[0]);
        histall_single->Fill(merged_pos_pt[0], epi_pt[0]);
      }
      if (pt_nhit == 2)
      {
        histpt_double[merged_pt_region[0]]->Fill(epi_pt[0] + epi_pt[1]);
        histpt_double_mod[merged_pt_region[0]]->Fill(gapmod_epi_pt[0]);
        histpt_double_all->Fill(gapmod_epi_pt[0]);
        histpt_all->Fill(gapmod_epi_pt[0]);
        hist_diffene_pt_adj[merged_pt_region[0]]->Fill(epi_pt[1] - epi_pt[0], epi_pt[0] + epi_pt[1]);
        hist_diffene_pt_adj_mod[merged_pt_region[0]]->Fill(epi_pt[1] - epi_pt[0], gapmod_epi_pt[0]);
        histall_double->Fill(merged_pos_pt[0], gapmod_epi_pt[0]);
      }
      if (al_nhit == 1)
      {
        histal_single[merged_al_region[0]]->Fill(epi_al[0]);
        histal_single_all->Fill(epi_al[0]);
        histal_all->Fill(epi_al[0]);
        histall_single->Fill(merged_pos_al[0] + 128, epi_al[0]);
      }
      if (al_nhit == 2)
      {
        histal_double[merged_al_region[0]]->Fill(epi_al[0] + epi_al[1]);
        histal_double_mod[merged_al_region[0]]->Fill(gapmod_epi_al[0]);
        histal_double_all->Fill(gapmod_epi_al[0]);
        histal_all->Fill(gapmod_epi_al[0]);
        hist_diffene_al_adj[merged_al_region[0]]->Fill(epi_al[1] - epi_al[0], epi_al[0] + epi_al[1]);
        hist_diffene_al_adj_mod[merged_al_region[0]]->Fill(epi_al[1] - epi_al[0], gapmod_epi_al[0]);
        histall_double->Fill(merged_pos_al[0] + 128, gapmod_epi_al[0]);
      }

      if (pt_nhit == 1 && pt_nhit == 1 && abs(epi_al[0] - epi_pt[0]) < 0.3)
        histsingle_onlygood->Fill((epi_al[0] + epi_pt[0]) / 2.);

      histsum_all->Fill((gapmod_epi_al[0] + gapmod_epi_pt[0]) / 2.0);
      hist_avediff[pt_nhit - 1][al_nhit - 1]->Fill((gapmod_epi_pt[0] + gapmod_epi_al[0]) / 2., (gapmod_epi_pt[0] - gapmod_epi_al[0]) / 2.);
      hist_avediff_rot[pt_nhit - 1][al_nhit - 1]->Fill((gapmod_epi_pt[0] - gapmod_epi_al[0]) / 2., (gapmod_epi_pt[0] + gapmod_epi_al[0]) / 2.);
    }

    hist_event_quality->Fill(event_quality_flag);
    if (event_quality_flag > 0)
      tout->Fill();
  }
  cerr << "Event Loss rate:" << nentries - tout->GetEntries() << " out of " << nentries
       << ":" << 100. - 100. * tout->GetEntries() / nentries << "\%" << endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 1800, 900);
  c1->Divide(2, 1);
  TString pdfname = TString(inFileName).ReplaceAll(".root", "_gapmod.pdf");
  c1->Print(pdfname + "[", "pdf");
  for (int j = 0; j < 3; j++)
  {
    c1->cd(1);
    gPad->SetGrid();
    histpt_single[j]->GetXaxis()->SetRangeUser(0.5, 30);
    histpt_double[j]->GetXaxis()->SetRangeUser(0.5, 30);
    histpt_double_mod[j]->GetXaxis()->SetRangeUser(0.5, 30);
    histpt_single[j]->Scale(1.0 / histpt_single[j]->Integral());
    histpt_double[j]->Scale(1.0 / histpt_double[j]->Integral());
    histpt_double_mod[j]->Scale(1.0 / histpt_double_mod[j]->Integral());
    histpt_double_mod[j]->Draw("HIST");
    histpt_double[j]->Draw("HISTsame");
    histpt_single[j]->SetLineColor(kBlack);
    histpt_double_mod[j]->SetLineColor(kRed);
    histpt_double[j]->SetLineColor(kBlack);

    c1->cd(2);
    gPad->SetGrid();
    histal_single[j]->GetXaxis()->SetRangeUser(0.5, 30);
    histal_double[j]->GetXaxis()->SetRangeUser(0.5, 30);
    histal_double_mod[j]->GetXaxis()->SetRangeUser(0.5, 30);
    histal_single[j]->Scale(1.0 / histal_single[j]->Integral());
    histal_double[j]->Scale(1.0 / histal_double[j]->Integral());
    histal_double_mod[j]->Scale(1.0 / histal_double_mod[j]->Integral());
    histal_double_mod[j]->Draw("HIST");
    histal_double[j]->Draw("HISTsame");
    histal_single[j]->SetLineColor(kBlack);
    histal_double_mod[j]->SetLineColor(kRed);
    histal_double[j]->SetLineColor(kBlack);
    c1->Print(pdfname, "pdf");
  }

  TLine *line3 = new TLine(0, 0, 0, 30);
  TLine *line4 = new TLine(-17.751, 17.751, 17.751, 17.751);
  line3->SetLineStyle(2);
  line4->SetLineStyle(2);
  for (int j = 0; j < 3; j++)
  {
    c1->cd(1);
    gPad->SetGrid(0);
    hist_diffene_pt_adj_mod[j]->Draw("colz");
    hist_diffene_pt_adj_mod[j]->GetXaxis()->SetRangeUser(-30, 30);
    hist_diffene_pt_adj_mod[j]->GetYaxis()->SetRangeUser(0.5, 30);
    line3->Draw("same");
    // line4->Draw("same");
    c1->cd(2);
    gPad->SetGrid(0);
    hist_diffene_al_adj_mod[j]->Draw("colz");
    hist_diffene_al_adj_mod[j]->GetXaxis()->SetRangeUser(-30, 30);
    hist_diffene_al_adj_mod[j]->GetYaxis()->SetRangeUser(0.5, 30);
    line3->Draw("same");
    // line4->Draw("same");
    c1->Print(pdfname, "pdf");
  }

  for (int j = 0; j < 3; j++)
  {
    c1->cd(1);
    gPad->SetGrid(0);
    hist_diffene_pt_adj[j]->Draw("colz");
    hist_diffene_pt_adj[j]->GetXaxis()->SetRangeUser(-30, 30);
    hist_diffene_pt_adj[j]->GetYaxis()->SetRangeUser(0.5, 30);
    line3->Draw("same");
    line4->Draw("same");
    c1->cd(2);
    gPad->SetGrid(0);
    hist_diffene_al_adj[j]->Draw("colz");
    hist_diffene_al_adj[j]->GetXaxis()->SetRangeUser(-30, 30);
    hist_diffene_al_adj[j]->GetYaxis()->SetRangeUser(0.5, 30);
    line3->Draw("same");
    line4->Draw("same");
    c1->Print(pdfname, "pdf");
  }

  c1->cd(1);
  gPad->SetGrid();
  histpt_single_all->GetXaxis()->SetRangeUser(0.5, 30);
  histpt_double_all->GetXaxis()->SetRangeUser(0.5, 30);
  histpt_all->GetXaxis()->SetRangeUser(0.5, 30);
  histpt_all->Draw("HIST");
  histpt_single_all->Draw("HISTsame");
  histpt_double_all->Draw("HISTsame");
  histpt_all->SetLineColor(kBlack);
  histpt_single_all->SetLineColor(kRed);
  histpt_double_all->SetLineColor(kBlue);
  c1->cd(2);
  gPad->SetGrid();
  histal_single_all->GetXaxis()->SetRangeUser(0.5, 30);
  histal_double_all->GetXaxis()->SetRangeUser(0.5, 30);
  histal_all->GetXaxis()->SetRangeUser(0.5, 30);
  histal_all->Draw("HIST");
  histal_single_all->Draw("HISTsame");
  histal_double_all->Draw("HISTsame");
  histal_all->SetLineColor(kBlack);
  histal_single_all->SetLineColor(kRed);
  histal_double_all->SetLineColor(kBlue);
  c1->Print(pdfname, "pdf");

  c1->cd(1);
  gPad->SetGrid();
  histal_all->GetXaxis()->SetRangeUser(0.5, 30);
  histpt_all->GetXaxis()->SetRangeUser(0.5, 30);
  histsum_all->GetXaxis()->SetRangeUser(0.5, 30);
  histsum_all->Draw("HIST");
  histsum_all->SetLineColor(kBlack);
  histpt_all->SetLineColor(kRed);
  histal_all->SetLineColor(kBlue);
  c1->cd(2);
  gPad->SetGrid();
  histsingle_onlygood->GetXaxis()->SetRangeUser(0.5, 30);
  histsingle_onlygood->SetLineColor(kBlack);
  histsingle_onlygood->Draw("HIST");
  c1->Print(pdfname, "pdf");

  c1->cd(1);
  gPad->SetGrid();
  histall_single->Draw("colz");
  c1->cd(2);
  histall_double->Draw("colz");
  c1->Print(pdfname, "pdf");
  c1->Print(pdfname + "]", "pdf");

  // Output File -------------------------------------
  tin_hktree->SetBranchStatus("*", 1);
  TTree *tout_hktree = tin_hktree->CloneTree();
  Long64_t nevent_gapmod = tout->GetEntries();
  tout_hktree->Branch("nevent_gapmod", &nevent_gapmod, "nevent_gapmod/L");
  tout_hktree->GetEntry(0);
  tout_hktree->Fill();

  tin_threshold_tree->SetBranchStatus("*", 1);
  TTree *tout_threshold_tree = tin_threshold_tree->CloneTree();

  TSpline3 *s[256];
  for (Int_t i = 0; i < 256; ++i)
  {
    s[i] = (TSpline3 *)fin->Get(Form("spline%d", i));
    s[i]->Write();
  }

  for (int j = 0; j < 3; j++)
  {
    histpt_single[j]->Write();
    histpt_double[j]->Write();
    histpt_double_mod[j]->Write();
    hist_diffene_pt_adj[j]->Write();
    hist_diffene_pt_adj_mod[j]->Write();
    histal_single[j]->Write();
    histal_double[j]->Write();
    histal_double_mod[j]->Write();
    hist_diffene_al_adj[j]->Write();
    hist_diffene_al_adj_mod[j]->Write();
  }
  histal_all->Write();
  histal_single_all->Write();
  histal_double_all->Write();
  histpt_all->Write();
  histpt_single_all->Write();
  histpt_double_all->Write();
  histsum_all->Write();
  histall_double->Write();
  histall_single->Write();
  histsingle_onlygood->Write();
  hist_event_quality->Write();

  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 2; j++)
    {
      hist_avediff[i][j]->Write();
      hist_avediff_rot[i][j]->Write();
    }
  }

  tout->Write();
  tout_hktree->Write();
  tout_threshold_tree->Write();
  fout->Close();

  cout << "  " << endl;
  cout << "Output ROOT file :" << outFilename << endl;
  cout << "Finish!" << endl;

  return 0;
}
