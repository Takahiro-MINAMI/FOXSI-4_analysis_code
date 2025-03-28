// pick up single event for calibration file
// this make _base.root -> _base_single.root
// then _base_single.root is used to make en.txt by peaktotxt.C

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
#include <string>
#include <iostream>
#include <vector>
#include <string>
#include <numeric>

using namespace std;

void single(TString inFileName, Int_t *adcThreshIn);

int main(int argc, char **argv)
{

  if (argc == 3)
  {
    ifstream ifs(argv[2]);
    string str;
    int count = 0;
    Int_t adcThreshIn[256] = {};
    // Load Threshold File -------------------------------------
    cout << "Loading threshold file:  " << argv[2] << endl;
    if (ifs.fail())
    {
      cerr << "ERROR: Failed to load threshold file" << endl;
      return -1;
    }
    while (getline(ifs, str))
    {
      if (count > 255)
        break;
      adcThreshIn[count] = stoi(str);
      cerr << count << " " << stoi(str) << endl;
      ++count;
    }
    single(argv[1], adcThreshIn);
  }
  else
  {
    cout << "USAGE:" << endl;
    cout << "./make_single filename threshold_file" << endl;
    cout << "ex) ./make_single ./test_base.root threth_single.txt" << endl;
    return -1;
  }
  return 0;
}

void single(TString inFileName, Int_t *adcThreshIn)
{

  // Input and Set TTree Branch Address -------------------------------------
  if (!inFileName.EndsWith("base.root"))
  {
    cerr << "ERROR:Please input *_base.root file" << endl;
    return;
  }
  TFile *fin = new TFile(inFileName);
  TTree *tin;
  fin->GetObject("eventtree", tin);

  // Declaration of leaves types -------------------------------------
  UShort_t flag_pseudo;
  Int_t pt_num;
  Int_t adc_pt[128];
  Int_t remapch_pt[128];
  Int_t al_num;
  Int_t adc_al[128];
  Int_t remapch_al[128];

  tin->SetBranchAddress("flag_pseudo", &flag_pseudo);
  tin->SetBranchAddress("pt_num", &pt_num);
  tin->SetBranchAddress("adc_pt", adc_pt);
  tin->SetBranchAddress("remapch_pt", remapch_pt);
  tin->SetBranchAddress("al_num", &al_num);
  tin->SetBranchAddress("adc_al", adc_al);
  tin->SetBranchAddress("remapch_al", remapch_al);

  tin->SetBranchStatus("*", 0);
  tin->SetBranchStatus("*pt*", 1);
  tin->SetBranchStatus("*al*", 1);

  // Create new ROOT file and TTree  -------------------------------------
  TString outFilename = inFileName;
  outFilename.ReplaceAll(".root", "_single.root");
  TFile *fout = new TFile(outFilename, "RECREATE");

  TTree *tout = tin->CloneTree(0);
  Long64_t nentries = tin->GetEntries();
  Int_t cnt_al = 0;
  Int_t cnt_pt = 0;
  Int_t chsingle_al = 0;
  Int_t chsingle_pt = 0;
  Int_t ch = 0;
  Int_t adc = 0;
  Int_t index_pt = 0;
  Int_t index_al = 0;
  Long64_t nbytes = 0;
  vector<Double_t> pt_xdata, pt_ydata;
  vector<Double_t> al_xdata, al_ydata;
  Double_t pos_pt, pos_al;
  Double_t ene_pt, ene_al;

  // setting threshold for each ASICs
  Int_t adcThresh[256] = {};
  for (Int_t k = 0; k < 256; ++k)
    adcThresh[k] = adcThreshIn[k];

  tout->Branch("cnt_al", &cnt_al, "cnt_al/I");
  tout->Branch("cnt_pt", &cnt_pt, "cnt_pt/I");
  tout->Branch("chsingle_al", &chsingle_al, "chsingle_al/I");
  tout->Branch("chsingle_pt", &chsingle_pt, "chsingle_pt/I");

  TH1D *hist_single[256];
  for (Int_t i = 0; i < 256; i++)
  {
    hist_single[i] = new TH1D(Form("hist_single%d", i), Form("hist_single%d", i), 1000, -0.5, 999.5);
  }
  TH2D *histall_single = new TH2D("histall_single", "histall_single", 256, -0.5, 255.5, 1000, -0.5, 999.5);
  TH2D *histall_double = new TH2D("histall_double", "histall_double", 256, -0.5, 255.5, 1000, -0.5, 999.5);
  TH2D *histall_tri = new TH2D("histall_tri", "histall_tri", 256, -0.5, 255.5, 1000, -0.5, 999.5);
  TH2D *histall_all = new TH2D("histall_all", "histall_all", 256, -0.5, 255.5, 1000, -0.5, 999.5);
  TH2D *histall_image = new TH2D("histall_image", "histall_image", 128, -0.5, 127.5, 128, -0.5, 127.5);
  TH2D *hist_nhit = new TH2D("hist_nhit", "hist_nhit", 25, 0.5, 25.5, 25, 0.5, 25.5);

  cout << "  " << endl;
  cout << "Open ROOT file :" << inFileName << endl;
  cout << "# of input events:  " << nentries << endl;
  // Main Process  -------------------------------------
  for (Long64_t l = 0; l < nentries; l++)
  {
    nbytes += tin->GetEntry(l);
    if (flag_pseudo == 1)
      continue;
    cnt_al = 0;
    cnt_pt = 0;
    for (Int_t i = 0; i < pt_num; ++i)
    {
      ch = remapch_pt[i];
      if (adc_pt[i] > adcThresh[ch])
      {
        ++cnt_pt;
        pt_xdata.push_back(ch);
        pt_ydata.push_back(adc_pt[i]);
        index_pt = i;
      }
    }
    ene_pt = accumulate(pt_ydata.begin(), pt_ydata.end(), 0.0);
    pos_pt = accumulate(pt_xdata.begin(), pt_xdata.end(), 0.0) / pt_xdata.size();
    if (cnt_pt == 1)
    {
      ch = remapch_pt[index_pt];
      chsingle_pt = ch;
      adc = adc_pt[index_pt];
      hist_single[ch]->Fill(adc);
      histall_single->Fill(ch, adc);
      histall_all->Fill(ch, adc);
    }
    else if (cnt_pt == 2)
    {
      ch = remapch_pt[index_pt];
      chsingle_pt = ch;
      if (index_pt != 0 && (index_pt != 127))
      {
        adc = ene_pt;
        histall_double->Fill(ch, adc);
        histall_all->Fill(ch, adc);
      }
    }
    else if (cnt_pt == 3)
    {
      ch = remapch_pt[index_pt];
      chsingle_pt = ch;
      if (index_pt != 0 && (index_pt != 127))
      {
        adc = ene_pt;
        histall_tri->Fill(ch, adc);
      }
    }
    else
    {
    }
    cnt_al = 0;
    for (Int_t i = 0; i < al_num; ++i)
    {
      ch = remapch_al[i];
      if (adc_al[i] > adcThresh[ch + 128])
      {
        ++cnt_al;
        al_xdata.push_back(ch);
        al_ydata.push_back(adc_al[i]);
        index_al = i;
      }
    }
    ene_al = accumulate(al_ydata.begin(), al_ydata.end(), 0.0);
    pos_al = accumulate(al_xdata.begin(), al_xdata.end(), 0.0) / al_xdata.size();
    if (cnt_al == 1)
    {
      ch = remapch_al[index_al] + 128;
      chsingle_al = ch;
      adc = adc_al[index_al];
      hist_single[ch]->Fill(adc);
      histall_single->Fill(ch, adc);
    }
    else if (cnt_al == 2)
    {
      ch = remapch_al[index_al] + 128;
      ;
      chsingle_al = ch;
      if (index_al != 0 && (index_al != 127))
      {
        adc = ene_al;
        histall_double->Fill(ch, adc);
        histall_all->Fill(ch, adc);
      }
    }
    else if (cnt_al == 3)
    {
      ch = remapch_al[index_al] + 128;
      chsingle_al = ch;
      if (index_al != 0 && (index_al != 127))
      {
        adc = ene_al;
        histall_tri->Fill(ch, adc);
      }
    }
    else
    {
    }
    hist_nhit->Fill(cnt_al, cnt_pt);

    if (cnt_pt == 1 || cnt_al == 1)
      tout->Fill();

    if (l % 100000 == 0)
      cerr << "Proccessed " << l << " events(" << 100. * l / nentries << "\%)" << endl;

    histall_image->Fill(pos_pt, 128 - pos_al);
    pt_xdata.clear();
    pt_ydata.clear();
    al_xdata.clear();
    al_ydata.clear();
  }

  // Output File -------------------------------------
  fout->Write();
  fin->Close();
  fout->Close();

  cout << "  " << endl;
  cout << "Output ROOT file :" << outFilename << endl;
  cout << "Finish!" << endl;
  return;
}
