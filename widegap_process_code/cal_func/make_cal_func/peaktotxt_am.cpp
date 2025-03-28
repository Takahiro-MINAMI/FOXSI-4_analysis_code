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

void peaktotxt_am()
{

  string infilename = "./Am241_single_sum.root";
  string threshname = "./threth_single.txt";

  ofstream outputfile("en.txt", std::ios::app);
  TFile *fin = new TFile(infilename.c_str(), "read");

  Int_t npeaks = 5;
  Int_t ave_width = 4;
  double weight_energy = (20.784 * 1.39 + 21.099 * 0.65 + 21.342 * 0.59 + 21.491 * 0.29) / (1.39 + 0.65 + 0.59 + 0.29);
  Double_t kev_ene[] = {13.946, 17.751, 59.5412, 26.3448, weight_energy};

  Int_t adc_start1[] = {130, 170, 690, 260, 205};
  Int_t adc_end1[] = {155, 210, 740, 300, 240};
  Int_t adc_start2[] = {130, 170, 690, 260, 205};
  Int_t adc_end2[] = {155, 210, 740, 300, 240};
  Int_t adc_start3[] = {150, 205, 600, 300, 255};
  Int_t adc_end3[] = {195, 250, 700, 350, 290};
  Int_t adc_start4[] = {150, 205, 600, 280, 250};
  Int_t adc_end4[] = {190, 250, 700, 350, 280};

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
    // cout<<"Detect peak at ch ="<<j<<endl;
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

  outputfile.close();
  cout << "OK, Detect peaks! Finish!" << endl;
}
