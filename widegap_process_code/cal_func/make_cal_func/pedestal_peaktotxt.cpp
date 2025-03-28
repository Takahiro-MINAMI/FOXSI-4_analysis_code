// To add  pedestal peak to en.txt
// use _base_ped.root

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

void pedestal_peaktotxt()
{
  TFile *fin = new TFile("./Ba133_withped_base_sum.root", "read");
  Int_t npeaks = 1;
  Int_t chmin = 5;
  Int_t chmax = 256;
  Int_t adc_start[] = {-20};
  Int_t adc_end[] = {20};
  Double_t kev_ene[] = {0.0};
  Int_t i = 0;
  TH1D *h1;
  Int_t tmp_bin;

  Int_t adcThreshIn[256] = {};
  ifstream ifs("./threth_single.txt");
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

  ofstream outputfile("en.txt", std::ios::app);
  for (Int_t j = chmin; j < chmax; ++j)
  {
    if (adcThreshIn[j] > 500)
    {
      continue;
    };
    h1 = (TH1D *)fin->Get(Form("hist_base%d", j));
    std::vector<Double_t> position;
    Double_t tmp_position;
    for (i = 0; i < npeaks; ++i)
    {
      Int_t xmin = adc_start[i];
      Int_t xmax = adc_end[i];
      h1->GetXaxis()->SetRangeUser(xmin, xmax);
      tmp_bin = h1->GetMaximumBin();
      h1->GetXaxis()->SetRange(tmp_bin - 2, tmp_bin + 2);
      tmp_position = h1->GetMean();
      outputfile << j << "    " << kev_ene[i] << "    " << tmp_position << endl;
    }
  }
  outputfile.close();
  cout << "OK, Detect peaks! Finish!" << endl;
}
