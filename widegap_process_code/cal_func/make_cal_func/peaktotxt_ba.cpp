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
using namespace std;

void peaktotxt_ba()
{

  string infilename = "./Ba133_single_sum.root";
  string threshname = "./threth_single.txt";

  ofstream outputfile("en.txt", std::ios::app);
  TFile *fin = new TFile(infilename.c_str(), "read");

  Int_t ave_width = 4;
  Int_t npeaks = 4;
  Double_t weight_energy = (34.920 * 5.99 + 34.987 * 11.6 + 35.818 * 3.58) / (5.99 + 11.6 + 3.58);
  Double_t kev_ene[] = {30.973, 7.80, weight_energy, weight_energy - 23.173};

  Int_t adc_start1[] = {300, 60, 380, 100};
  Int_t adc_end1[] = {500, 90, 420, 140};
  Int_t adc_start2[] = {300, 60, 380, 100};
  Int_t adc_end2[] = {500, 90, 420, 140};
  Int_t adc_start3[] = {300, 80, 410, 130};
  Int_t adc_end3[] = {500, 120, 450, 180};
  Int_t adc_start4[] = {300, 80, 400, 130};
  Int_t adc_end4[] = {500, 120, 450, 180};
  Int_t adc_start_ch4[] = {200};
  Int_t adc_end_ch4[] = {300};

  int i = 0;
  TH1D *h1;
  Int_t adcThreshIn[256] = {};
  ifstream ifs(threshname);
  string str;
  int count = 0;
  while (getline(ifs, str))
  {
    if (count > 255)
    {
      break;
    };
    adcThreshIn[count] = stoi(str);
    ++count;
  }
  Double_t tmp_position;
  Double_t tmp_bin;

  // ASIC1 -----------------------------------------------------------------
  npeaks = sizeof(adc_start1) / sizeof(adc_start1[0]);
  int chmin = 5;
  int chmax = 63;
  int j;
  for (j = chmin; j < chmax + 1; ++j)
  {
    if (adcThreshIn[j] > 500)
    {
      continue;
    };
    h1 = (TH1D *)fin->Get(Form("hist_single%d", j));
    for (i = 0; i < npeaks; ++i)
    {
      Int_t xmin = adc_start1[i];
      Int_t xmax = adc_end1[i];
      h1->GetXaxis()->SetRangeUser(xmin, xmax);
      tmp_bin = h1->GetMaximumBin();
      h1->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      tmp_position = h1->GetMean();
      outputfile << j << "    " << kev_ene[i] << "    " << tmp_position << endl;
    }
  }
  // ASIC2----------------------------------------------------------
  npeaks = sizeof(adc_start2) / sizeof(adc_start2[0]);
  chmin = 64;
  chmax = 127;
  for (j = chmin; j < chmax + 1; ++j)
  {
    if (adcThreshIn[j] > 500)
    {
      continue;
    };
    h1 = (TH1D *)fin->Get(Form("hist_single%d", j));

    for (i = 0; i < npeaks; ++i)
    {
      Int_t xmin = adc_start2[i];
      Int_t xmax = adc_end2[i];
      h1->GetXaxis()->SetRangeUser(xmin, xmax);
      tmp_bin = h1->GetMaximumBin();
      h1->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      tmp_position = h1->GetMean();
      outputfile << j << "    " << kev_ene[i] << "    " << tmp_position << endl;
    }
  }
  // ASIC3-----------------------------------------------------------------
  npeaks = sizeof(adc_start3) / sizeof(adc_start3[0]);
  chmin = 128;
  chmax = 191;
  for (j = chmin; j < chmax + 1; ++j)
  {
    if (adcThreshIn[j] > 500)
    {
      continue;
    };
    h1 = (TH1D *)fin->Get(Form("hist_single%d", j));
    for (i = 0; i < npeaks; ++i)
    {
      Int_t xmin = adc_start3[i];
      Int_t xmax = adc_end3[i];
      h1->GetXaxis()->SetRangeUser(xmin, xmax);
      tmp_bin = h1->GetMaximumBin();
      h1->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      tmp_position = h1->GetMean();
      outputfile << j << "    " << kev_ene[i] << "    " << tmp_position << endl;
    }
  }
  // ASIC4-----------------------------------------------------------------
  npeaks = sizeof(adc_start4) / sizeof(adc_start4[0]);
  chmin = 192;
  chmax = 255;
  for (int j = chmin; j < chmax + 1; ++j)
  {
    if (adcThreshIn[j] > 500)
    {
      continue;
    };
    h1 = (TH1D *)fin->Get(Form("hist_single%d", j));
    for (i = 0; i < npeaks; ++i)
    {
      Int_t xmin = adc_start4[i];
      Int_t xmax = adc_end4[i];
      h1->GetXaxis()->SetRangeUser(xmin, xmax);
      tmp_bin = h1->GetMaximumBin();
      h1->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      tmp_position = h1->GetMean();
      outputfile << j << "    " << kev_ene[i] << "    " << tmp_position << endl;
    }
  }

  outputfile.close();
  cout << "OK, Detect peaks! Finish!" << endl;
}
