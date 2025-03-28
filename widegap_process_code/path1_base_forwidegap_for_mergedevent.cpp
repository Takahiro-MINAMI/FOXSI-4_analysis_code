//
// Remap ch, cmn subtraction
// .root -> _base.root
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
#include <TSpline.h>
#include <TRandom3.h>
#include <string.h>
#include <filesystem>

using namespace std;
string path_to_ti2pps = "./cal_func/ti2pps.root";
string path_to_pps2realtime = "./cal_func/pps2realtime.root";

int base(const char *inFileName, Int_t det_number = 0, bool save_all = false);
UInt_t det_ti2timecode_ti(UInt_t det_ti);

int main(int argc, char **argv)
{
  if (argc == 2)
    base(argv[1], NAN);
  else if (argc == 3)
    base(argv[1], Int_t(stoi(argv[2])));
  else if (argc == 4)
    base(argv[1], Int_t(stoi(argv[2])), !strncmp(argv[3], "-all", 4));
  else
  {
    cout << "USAGE:" << endl;
    cout << "./base_forwidegap filename det_number (1-4:flight, 0:Other) flag_all" << endl;
    cout << "ex) ./base_forwidegap ./test.root 1" << endl;
    cout << "ex) ./base_forwidegap ./test.root 1 -all" << endl;
    return -1;
  }

  return 0;
}

UInt_t det_ti2timecode_ti(UInt_t det_ti)
{
  UInt_t shift_ti = det_ti >> 6;
  return shift_ti + 0x4000000;
}

int base(const char *inFileName, Int_t det_number, bool save_all)
{
  gStyle->SetOptStat(0);
  cout << "-- Base Process ------------" << endl;
  if (save_all)
    cout << "All data size mode ...." << endl;
  else
    cout << "Reduce data size mode ...." << endl;

  // Declaration of leaves types -------------------------------------
  Int_t adc_cmn_al[128];
  UChar_t hitnum_al;
  UChar_t index_al[128];
  Int_t adc_cmn_pt[128];
  UChar_t hitnum_pt;
  UChar_t index_pt[128];

  UInt_t ti;
  UInt_t livetime;
  UInt_t integral_livetime;
  UInt_t unixtime;

  UChar_t flag_pseudo;
  UInt_t pseudo_counter;

  // Input and Set TTree Branch Address -------------------------------------
  TString infilename = TString(inFileName);
  if (!infilename.EndsWith(".root"))
  {
    cout << "Error: Input file " << infilename << " is NOT root file!" << endl;
    return -1;
  }
  TFile *fin = new TFile(inFileName);
  TTree *tin = (TTree *)fin->Get("eventtree");

  tin->SetBranchAddress("adc_cmn_al", adc_cmn_al);
  tin->SetBranchAddress("index_al", index_al);
  tin->SetBranchAddress("hitnum_al", &hitnum_al);

  tin->SetBranchAddress("adc_cmn_pt", adc_cmn_pt);
  tin->SetBranchAddress("index_pt", index_pt);
  tin->SetBranchAddress("hitnum_pt", &hitnum_pt);

  tin->SetBranchAddress("ti", &ti);
  tin->SetBranchAddress("livetime", &livetime);
  tin->SetBranchAddress("de_unixtime", &unixtime);
  tin->SetBranchAddress("flag_pseudo", &flag_pseudo);
  tin->SetBranchAddress("pseudo_counter", &pseudo_counter);

  // Save org adc information or not -------------------------------------
  if (!save_all)
  {
    tin->SetBranchStatus("adc_cmn_*", 0);
    tin->SetBranchStatus("index*", 0);
    tin->SetBranchStatus("hitnum*", 0);
  }

  // Load Time cal File for Flight-------------------------------------
  TGraph *g_ti2pps;
  TGraph *g_pps2launchtime;
  bool is_this_flight_data = (det_number >= 1 && det_number <= 4);
  if (is_this_flight_data)
  {
    cout << "Loading ti-to-pps file:  " << path_to_ti2pps << endl;
    if (!std::filesystem::is_regular_file(path_to_ti2pps))
    {
      cout << "ERROR: NO such a ti-to-pps file as " << path_to_ti2pps << " exists!" << endl;
      return -1;
    }
    TFile *ti2pps_calfile = new TFile(TString(path_to_ti2pps));
    TString ti2pps_name;
    cout << "For Flight: Detector Number:   " << det_number << endl;
    ti2pps_name = TString::Format("CdTe%dtivsPPSti", det_number);
    g_ti2pps = (TGraph *)ti2pps_calfile->Get(ti2pps_name);

    cout << "Loading pps-to-realtime file:  " << path_to_pps2realtime << endl;
    if (!std::filesystem::is_regular_file((path_to_pps2realtime)))
    {
      cout << "ERROR: NO such a pps-to-realtime file as " << path_to_pps2realtime << " exists!" << endl;
      return -1;
    }
    TFile *pps2realtime_calfile = new TFile(TString(path_to_pps2realtime));
    g_pps2launchtime = (TGraph *)pps2realtime_calfile->Get("CdTepps2realtime");
  }
  else
    cout << "Not for Flight: Time Sync. not applied " << endl;

  // Initial unixcut for Flight
  Int_t timecut_unix_st[4] = {1685154982, 1685155180, 1685155052, 1685155034};
  Int_t timecut_unix_ed[4] = {1685155005, 1685155189, 1685155070, 1685155056};
  // Int_t timecut_unix_st[4] = {0, 0, 0, 0};
  // Int_t timecut_unix_ed[4] = {2000000000, 2000000000, 2000000000, 2000000000};

  // Create new ROOT file and TTree  -------------------------------------
  TString outFilename = inFileName;
  outFilename.ReplaceAll(".root", "_base.root");
  TFile *fout = new TFile(outFilename, "RECREATE");

  TTree *tout = tin->CloneTree(0);
  tin->SetBranchStatus("cmn*", 1);
  tin->SetBranchStatus("adc*", 1);
  tin->SetBranchStatus("index*", 1);
  tin->SetBranchStatus("hitnum*", 1);

  Int_t pt_adc_num = 0;
  Int_t adc_pt[128] = {};
  Int_t remapch_pt[128] = {0};
  Int_t al_adc_num = 0;
  Int_t adc_al[128] = {};
  Int_t remapch_al[128] = {0};
  UInt_t survive_pseudo_counter = 0;
  Double_t pps_counter = 0.;
  Double_t approx_passed_time_from_launch_sec = 0.;
  ULong64_t integral_livetime_manual;

  Long64_t nevent_org = 0;
  Long64_t nevent_base = 0;

  tout->Branch("pt_adc_num", &pt_adc_num, "pt_adc_num/I");
  tout->Branch("adc_pt", adc_pt, "adc_pt[pt_adc_num]/I");
  tout->Branch("remapch_pt", remapch_pt, "remapch_pt[pt_adc_num]/I");
  tout->Branch("al_adc_num", &al_adc_num, "al_adc_num/I");
  tout->Branch("adc_al", adc_al, "adc_al[al_adc_num]/I");
  tout->Branch("remapch_al", remapch_al, "remapch_al[al_adc_num]/I");
  tout->Branch("survive_pseudo_counter", &survive_pseudo_counter, "survive_pseudo_counter/i");
  tout->Branch("pps_counter", &pps_counter, "pps_counter/D");
  tout->Branch("approx_passed_time_from_launch_sec", &approx_passed_time_from_launch_sec, "approx_passed_time_from_launch_sec/D");
  tout->Branch("integral_livetime_manual", &integral_livetime_manual, "integral_livetime_manual/l");

  TH1D *hist_base[256];
  for (Int_t i = 0; i < 256; i++)
  {
    hist_base[i] = new TH1D(Form("hist_base%d", i), Form("hist_base%d", i), 1020, -20.5, 999.5);
  }
  TH2D *histall_base = new TH2D("histall_base", "histall_base", 256, -0.5, 255.5, 1050, -50.5, 999.5);
  TH1D *histpt_base = new TH1D("histpt_base", "histpt_base", 1020, -20.5, 999.5);
  TH1D *histal_base = new TH1D("histal_base", "histal_base", 1020, -20.5, 999.5);

  Int_t ch;
  Int_t adc;

  ULong64_t integral_livetime_before = 0;
  UInt_t pseudo_counter_before = 0;

  Long64_t nentries = tin->GetEntriesFast();
  nevent_org = nentries;

  cout << "  " << endl;
  cout << "Open ROOT file :" << inFileName << endl;
  cout << "# of input events:  " << nentries << endl;

  // Main Process  -------------------------------------
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    if (tin->GetEntry(jentry) < 0)
    {
      cerr << "ERROR: error occurs while loading tree" << endl;
      break;
    }
    if (is_this_flight_data)
    {
      pps_counter = g_ti2pps->Eval(det_ti2timecode_ti(ti));
      approx_passed_time_from_launch_sec = g_pps2launchtime->Eval(pps_counter) - 887.5;

      if (approx_passed_time_from_launch_sec < 0 && !save_all)
        continue;
    }

    if (pseudo_counter < 0)
    {
      continue;
    }

    pt_adc_num = 0;
    al_adc_num = 0;
    for (Int_t i = 0; i < 128; i++)
    {
      adc_pt[i] = 0;
      adc_al[i] = 0;
      remapch_pt[i] = 0;
      remapch_al[i] = 0;
    }
    for (Int_t i = 0; i < hitnum_pt; i++)
    {
      if (index_pt[i] < 64)
      {
        ch = 63 - index_pt[i];
      }
      else if (index_pt[i] >= 64 && index_pt[i] < 128)
      {
        ch = 127 - (index_pt[i] - 64);
      }
      else
      {
        continue;
      }
      remapch_pt[pt_adc_num] = ch;
      adc = adc_cmn_pt[i];
      adc_pt[pt_adc_num] = adc;
      hist_base[ch]->Fill(adc);
      histall_base->Fill(ch, adc);
      histpt_base->Fill(adc);
      pt_adc_num++;
    }

    for (Int_t i = 0; i < hitnum_al; i++)
    {

      if (index_al[i] < 64)
      {
        ch = 192 + 63 - index_al[i];
      }
      else if (index_al[i] >= 64 && index_al[i] < 128)
      {
        ch = 128 + 63 - (index_al[i] - 64);
      }
      else
      {
        continue;
      }
      remapch_al[al_adc_num] = ch - 128;
      adc = adc_cmn_al[i];
      adc_al[al_adc_num] = adc;
      hist_base[ch]->Fill(adc);
      histall_base->Fill(ch, adc);
      histal_base->Fill(adc);
      al_adc_num++;
    }

    integral_livetime_manual = integral_livetime_before + livetime;
    integral_livetime_before = integral_livetime_manual;

    if (flag_pseudo == 1)
      survive_pseudo_counter = pseudo_counter_before + 1;
    else
      survive_pseudo_counter = pseudo_counter_before;
    pseudo_counter_before = survive_pseudo_counter;

    tout->Fill();

    if (jentry % 100000 == 0)
    {
      cerr << "Proccessed " << jentry << " events(" << jentry * 100. / nentries << "\%)" << endl;
    }
  }

  // HK Information  -------------------------------------
  nevent_base = tout->GetEntriesFast();
  cout << "# of output events:  " << nevent_base << endl;

  TTree *hktree = new TTree("hktree", "hktree");
  hktree->Branch("det_number", &det_number, "det_number/I");
  hktree->Branch("nevent_org", &nevent_org, "nevent_org/L");
  hktree->Branch("nevent_base", &nevent_base, "nevent_base/L");
  hktree->GetEntry(0);
  hktree->Fill();

  // Output File -------------------------------------
  for (Int_t i = 0; i < 256; i++)
    hist_base[i]->Write();
  histall_base->Write();
  histpt_base->Write();
  histal_base->Write();
  tout->Write();
  hktree->Write();
  fout->Close();

  cout << "  " << endl;
  cout << "Output ROOT file :" << outFilename << endl;
  cout << "Finish!" << endl;

  return 0;
}
