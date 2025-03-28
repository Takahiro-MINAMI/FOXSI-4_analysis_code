
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
#include <Rtypes.h>
#include <vector>

using namespace std;

void make_doi_func()
{

  gStyle->SetOptStat(0);
  gStyle->SetPalette(53);
  TString inFileName = "test.root";

  TFile *fin = new TFile(inFileName, "read");
  Int_t ave_width = 2;
  Double_t range_peak1[] = {-1.5, 1.3};
  Int_t colorlist[] = {kRed, kOrange, kGreen + 2, kBlue};
  TString pdfname = "test.pdf";
  Double_t kev_ene_low = 13.946;
  Double_t kev_ene_high = 26.3448;
  Double_t x_turn_low[3] = {-0.5, -0.4, -0.5};
  Double_t x_turn_high[3] = {0., 0., -0.0};

  Int_t n_point = 500;
  Int_t niterp = 20;
  Int_t niterp_low = 20;
  Int_t niterp_high = 20;
  Int_t tmp_point = 0;
  Double_t tmpx, tmpy, targety, tmp_peak, ene_ratio;
  Double_t func_xmin = -10.0;
  Double_t func_xmax = 10.0;
  Double_t change_val = 0.15;
  TGraph *g_peak_low = new TGraph();
  TGraph *g_peak_high = new TGraph();

  TF1 *fitfunc_low, *fitfunc2_low, *fitfunc3_low;
  TF1 *fitfunc_high, *fitfunc2_high, *fitfunc3_high;
  TGraph *intp[niterp];
  TGraph *intp_low[niterp_low];
  TGraph *intp_high[niterp_high];
  for (Int_t i = 0; i < niterp; i++)
  {
    intp[i] = new TGraph();
    intp[i]->SetMarkerStyle(7);
    intp[i]->SetLineColor(kRed);
    intp[i]->SetMarkerColor(kRed);
  }
  for (Int_t i = 0; i < niterp_low; i++)
  {
    intp_low[i] = new TGraph();
    intp_low[i]->SetMarkerStyle(7);
    intp_low[i]->SetLineColor(kRed);
    intp_low[i]->SetMarkerColor(kRed);
  }
  for (Int_t i = 0; i < niterp_high; i++)
  {
    intp_high[i] = new TGraph();
    intp_high[i]->SetMarkerStyle(7);
    intp_high[i]->SetLineColor(kRed);
    intp_high[i]->SetMarkerColor(kRed);
  }

  TGraph2D *g_response[2][2];
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      g_response[i][j] = new TGraph2D();

  Double_t offset = 0.0;
  Double_t xmin, xmax;
  Int_t tmp_bin;
  Double_t max, mean, center, kev_ene;
  TH1D *ab;
  TH2D *hist_avediff;

  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
  c1->SetLogz();

  hist_avediff = (TH2D *)fin->Get("hist_avediff_rot_pt1al1");
  xmin = -1.5;
  xmax = 2.0;
  kev_ene = kev_ene_low;
  for (Int_t j = hist_avediff->GetXaxis()->FindBin(xmin); j < hist_avediff->GetXaxis()->FindBin(xmax); j++)
  {
    ab = (TH1D *)hist_avediff->ProjectionY("ab", j, j);
    center = abs(hist_avediff->GetXaxis()->GetBinCenter(j));
    ab->GetXaxis()->SetRangeUser(kev_ene - 2 - center, kev_ene + 1 - center);
    tmp_bin = ab->GetMaximumBin();
    ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
    max = ab->GetMean();
    g_peak_low->SetPoint(g_peak_low->GetN(), hist_avediff->GetXaxis()->GetBinCenter(j), max);
  }
  fitfunc_low = new TF1("fitfunc_low", "[0] +[1]*(x-[2]) ", -100.0, x_turn_low[0]);
  fitfunc2_low = new TF1("fitfunc2_low", "[0]", x_turn_low[0], x_turn_high[0]);
  fitfunc3_low = new TF1("fitfunc3_low", "[0]+[1]*(x-[2]) ", x_turn_high[0], 100.0);
  fitfunc_low->SetNpx(500);
  fitfunc2_low->SetNpx(500);
  fitfunc3_low->SetNpx(500);
  fitfunc_low->SetLineWidth(3);
  fitfunc2_low->SetLineWidth(3);
  fitfunc3_low->SetLineWidth(3);
  g_peak_low->Fit("fitfunc2_low", "", "", x_turn_low[0], x_turn_high[0]);
  fitfunc_low->FixParameter(0, fitfunc2_low->GetParameter(0));
  fitfunc3_low->FixParameter(0, fitfunc2_low->GetParameter(0));
  fitfunc_low->FixParameter(2, x_turn_low[0]);
  fitfunc3_low->FixParameter(2, x_turn_high[0]);
  g_peak_low->Fit("fitfunc_low", "", "", xmin, x_turn_low[0]);
  g_peak_low->Fit("fitfunc3_low", "", "", x_turn_high[0], xmax);

  xmin = -2.0;
  xmax = 1.5;
  kev_ene = kev_ene_high;
  for (Int_t j = hist_avediff->GetXaxis()->FindBin(xmin); j < hist_avediff->GetXaxis()->FindBin(xmax); j++)
  {
    ab = (TH1D *)hist_avediff->ProjectionY("ab", j, j);
    center = abs(hist_avediff->GetXaxis()->GetBinCenter(j));
    ab->GetXaxis()->SetRangeUser(kev_ene - 2 - center, kev_ene + 1 - center);
    tmp_bin = ab->GetMaximumBin();
    ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
    max = ab->GetMean();
    g_peak_high->SetPoint(g_peak_high->GetN(), hist_avediff->GetXaxis()->GetBinCenter(j), max);
  }
  fitfunc_high = new TF1("fitfunc_high", "[0]+[1]*(x-[2]) ", -100.0, x_turn_low[0]);
  fitfunc2_high = new TF1("fitfunc2_high", "[0]", x_turn_low[0], x_turn_high[0]);
  fitfunc3_high = new TF1("fitfunc3_high", "[0]+[1]*(x-[2]) ", x_turn_high[0], 100.0);
  fitfunc_high->SetNpx(500);
  fitfunc2_high->SetNpx(500);
  fitfunc3_high->SetNpx(500);
  fitfunc_high->SetLineWidth(3);
  fitfunc2_high->SetLineWidth(3);
  fitfunc3_high->SetLineWidth(3);
  g_peak_high->Fit("fitfunc2_high", "", "", x_turn_low[0], x_turn_high[0]);
  fitfunc_high->FixParameter(0, fitfunc2_high->GetParameter(0));
  fitfunc3_high->FixParameter(0, fitfunc2_high->GetParameter(0));
  fitfunc_high->FixParameter(2, x_turn_low[0]);
  fitfunc3_high->FixParameter(2, x_turn_high[0]);
  g_peak_high->Fit("fitfunc_high", "", "", xmin, x_turn_low[0]);
  g_peak_high->Fit("fitfunc3_high", "", "", x_turn_high[0], xmax);

  for (int i = 0; i < n_point; i++)
  {
    tmpx = i * (func_xmax - func_xmin) / n_point + func_xmin;
    for (int j = 0; j < niterp; j++)
    {
      targety = kev_ene_high - j * (kev_ene_high - kev_ene_low) / niterp;
      ene_ratio = (targety - kev_ene_low) / (kev_ene_high - kev_ene_low);
      if (tmpx > x_turn_high[0])
        tmp_peak = fitfunc3_low->Eval(tmpx) + (fitfunc3_high->Eval(tmpx) - fitfunc3_low->Eval(tmpx)) * ene_ratio;
      else if (tmpx < x_turn_low[0])
        tmp_peak = fitfunc_low->Eval(tmpx) + (fitfunc_high->Eval(tmpx) - fitfunc_low->Eval(tmpx)) * ene_ratio;
      else
        tmp_peak = fitfunc2_low->Eval(tmpx) + (fitfunc2_high->Eval(tmpx) - fitfunc2_low->Eval(tmpx)) * ene_ratio;
      g_response[0][0]->SetPoint(g_response[0][0]->GetN(), tmpx, tmp_peak, targety);
      intp[j]->SetPoint(i, tmpx, tmp_peak);
    }
    for (int j = 0; j < niterp_low; j++)
    {
      targety = j * (kev_ene_low - (-5.0)) / niterp_low;
      if (tmpx > x_turn_high[0])
        tmp_peak = fitfunc3_low->Eval(tmpx) - targety;
      else if (tmpx < x_turn_low[0])
        tmp_peak = fitfunc_low->Eval(tmpx) - targety;
      else
        tmp_peak = fitfunc2_low->Eval(tmpx) - targety;
      g_response[0][0]->SetPoint(g_response[0][0]->GetN(), tmpx, tmp_peak, kev_ene_low - targety);
      intp_low[j]->SetPoint(i, tmpx, tmp_peak);
    }
    for (int j = 0; j < niterp_high; j++)
    {
      targety = j * (40.0 - kev_ene_high) / niterp_high;
      if (tmpx > x_turn_high[0])
        tmp_peak = fitfunc3_high->Eval(tmpx) + targety;
      else if (tmpx < x_turn_low[0])
        tmp_peak = fitfunc_high->Eval(tmpx) + targety;
      else
        tmp_peak = fitfunc2_high->Eval(tmpx) + targety;
      g_response[0][0]->SetPoint(g_response[0][0]->GetN(), tmpx, tmp_peak, kev_ene_high + targety);
      intp_high[j]->SetPoint(i, tmpx, tmp_peak);
    }
  }

  g_peak_low->Set(0);
  g_peak_high->Set(0);

  hist_avediff = (TH2D *)fin->Get("hist_avediff_rot_pt2al1");
  xmin = -0.7;
  xmax = 2.0;
  kev_ene = kev_ene_low;
  for (Int_t j = hist_avediff->GetXaxis()->FindBin(xmin); j < hist_avediff->GetXaxis()->FindBin(xmax); j++)
  {
    ab = (TH1D *)hist_avediff->ProjectionY("ab", j, j);
    center = abs(hist_avediff->GetXaxis()->GetBinCenter(j));
    ab->GetXaxis()->SetRangeUser(kev_ene - 2 - center, kev_ene + 1 - center);
    tmp_bin = ab->GetMaximumBin();
    ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
    max = ab->GetMean();
    g_peak_low->SetPoint(g_peak_low->GetN(), hist_avediff->GetXaxis()->GetBinCenter(j), max);
  }
  fitfunc_low = new TF1("fitfunc_low", "[0] +[1]*(x-[2]) ", -100.0, x_turn_low[1]);
  fitfunc2_low = new TF1("fitfunc2_low", "[0]", x_turn_low[1], x_turn_high[1]);
  fitfunc3_low = new TF1("fitfunc3_low", "[0]+[1]*(x-[2]) ", x_turn_high[1], 100.0);
  fitfunc_low->SetNpx(500);
  fitfunc2_low->SetNpx(500);
  fitfunc3_low->SetNpx(500);
  fitfunc_low->SetLineWidth(3);
  fitfunc2_low->SetLineWidth(3);
  fitfunc3_low->SetLineWidth(3);
  g_peak_low->Fit("fitfunc2_low", "", "", x_turn_low[1], x_turn_high[1]);
  fitfunc_low->FixParameter(0, fitfunc2_low->GetParameter(0));
  fitfunc3_low->FixParameter(0, fitfunc2_low->GetParameter(0));
  fitfunc_low->FixParameter(2, x_turn_low[1]);
  fitfunc3_low->FixParameter(2, x_turn_high[1]);
  g_peak_low->Fit("fitfunc_low", "", "", xmin, x_turn_low[1]);
  g_peak_low->Fit("fitfunc3_low", "", "", x_turn_high[1], xmax);

  xmin = -1.0;
  xmax = 1.0;
  kev_ene = kev_ene_high;
  for (Int_t j = hist_avediff->GetXaxis()->FindBin(xmin); j < hist_avediff->GetXaxis()->FindBin(xmax); j++)
  {
    ab = (TH1D *)hist_avediff->ProjectionY("ab", j, j);
    center = abs(hist_avediff->GetXaxis()->GetBinCenter(j));
    ab->GetXaxis()->SetRangeUser(kev_ene - 2 - center, kev_ene + 1 - center);
    tmp_bin = ab->GetMaximumBin();
    ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
    max = ab->GetMean();
    g_peak_high->SetPoint(g_peak_high->GetN(), hist_avediff->GetXaxis()->GetBinCenter(j), max);
  }
  fitfunc_high = new TF1("fitfunc_high", "[0]+[1]*(x-[2]) ", -100.0, x_turn_low[1]);
  fitfunc2_high = new TF1("fitfunc2_high", "[0]", x_turn_low[1], x_turn_high[1]);
  fitfunc3_high = new TF1("fitfunc3_high", "[0]+[1]*(x-[2]) ", x_turn_high[1], 100.0);
  fitfunc_high->SetNpx(500);
  fitfunc2_high->SetNpx(500);
  fitfunc3_high->SetNpx(500);
  fitfunc_high->SetLineWidth(3);
  fitfunc2_high->SetLineWidth(3);
  fitfunc3_high->SetLineWidth(3);
  g_peak_high->Fit("fitfunc2_high", "", "", x_turn_low[1], x_turn_high[1]);
  fitfunc_high->FixParameter(0, fitfunc2_high->GetParameter(0));
  fitfunc3_high->FixParameter(0, fitfunc2_high->GetParameter(0));
  fitfunc_high->FixParameter(2, x_turn_low[1]);
  fitfunc3_high->FixParameter(2, x_turn_high[1]);
  g_peak_high->Fit("fitfunc_high", "", "", xmin, x_turn_low[1]);
  g_peak_high->Fit("fitfunc3_high", "", "", x_turn_high[1], xmax);

  for (int i = 0; i < n_point; i++)
  {
    tmpx = i * (func_xmax - func_xmin) / n_point + func_xmin;
    for (int j = 0; j < niterp; j++)
    {
      targety = kev_ene_high - j * (kev_ene_high - kev_ene_low) / niterp;
      ene_ratio = (targety - kev_ene_low) / (kev_ene_high - kev_ene_low);
      if (tmpx > x_turn_high[1])
        tmp_peak = fitfunc3_low->Eval(tmpx) + (fitfunc3_high->Eval(tmpx) - fitfunc3_low->Eval(tmpx)) * ene_ratio;
      else if (tmpx < x_turn_low[1])
        tmp_peak = fitfunc_low->Eval(tmpx) + (fitfunc_high->Eval(tmpx) - fitfunc_low->Eval(tmpx)) * ene_ratio;
      else
        tmp_peak = fitfunc2_low->Eval(tmpx) + (fitfunc2_high->Eval(tmpx) - fitfunc2_low->Eval(tmpx)) * ene_ratio;
      g_response[1][0]->SetPoint(g_response[1][0]->GetN(), tmpx, tmp_peak, targety);
      intp[j]->SetPoint(i, tmpx, tmp_peak);
    }
    for (int j = 0; j < niterp_low; j++)
    {
      targety = j * (kev_ene_low - (-5.0)) / niterp_low;
      if (tmpx > x_turn_high[1])
        tmp_peak = fitfunc3_low->Eval(tmpx) - targety;
      else if (tmpx < x_turn_low[1])
        tmp_peak = fitfunc_low->Eval(tmpx) - targety;
      else
        tmp_peak = fitfunc2_low->Eval(tmpx) - targety;
      g_response[1][0]->SetPoint(g_response[1][0]->GetN(), tmpx, tmp_peak, kev_ene_low - targety);
      intp_low[j]->SetPoint(i, tmpx, tmp_peak);
    }
    for (int j = 0; j < niterp_high; j++)
    {
      targety = j * (40.0 - kev_ene_high) / niterp_high;
      if (tmpx > x_turn_high[1])
        tmp_peak = fitfunc3_high->Eval(tmpx) + targety;
      else if (tmpx < x_turn_low[1])
        tmp_peak = fitfunc_high->Eval(tmpx) + targety;
      else
        tmp_peak = fitfunc2_high->Eval(tmpx) + targety;
      g_response[1][0]->SetPoint(g_response[1][0]->GetN(), tmpx, tmp_peak, kev_ene_high + targety);
      intp_high[j]->SetPoint(i, tmpx, tmp_peak);
    }
  }

  g_peak_low->Set(0);
  g_peak_high->Set(0);

  hist_avediff = (TH2D *)fin->Get("hist_avediff_rot_pt1al2");
  xmin = -1.3;
  xmax = 1.0;
  kev_ene = kev_ene_low;
  ave_width = 4.;
  for (Int_t j = hist_avediff->GetXaxis()->FindBin(xmin); j < hist_avediff->GetXaxis()->FindBin(xmax); j++)
  {
    ab = (TH1D *)hist_avediff->ProjectionY("ab", j, j);
    center = abs(hist_avediff->GetXaxis()->GetBinCenter(j));
    ab->GetXaxis()->SetRangeUser(kev_ene - 2 - center, kev_ene);
    tmp_bin = ab->GetMaximumBin();
    ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
    max = ab->GetMean();
    g_peak_low->SetPoint(g_peak_low->GetN(), hist_avediff->GetXaxis()->GetBinCenter(j), max);
  }
  fitfunc_low = new TF1("fitfunc_low", "[0] +[1]*(x-[2]) ", -100.0, x_turn_low[2]);
  fitfunc2_low = new TF1("fitfunc2_low", "[0]", x_turn_low[2], x_turn_high[2]);
  fitfunc3_low = new TF1("fitfunc3_low", "[0]+[1]*(x-[2]) ", x_turn_high[2], 100.0);
  fitfunc_low->SetNpx(500);
  fitfunc2_low->SetNpx(500);
  fitfunc3_low->SetNpx(500);
  fitfunc_low->SetLineWidth(3);
  fitfunc2_low->SetLineWidth(3);
  fitfunc3_low->SetLineWidth(3);
  g_peak_low->Fit("fitfunc2_low", "", "", x_turn_low[2], x_turn_high[2]);
  fitfunc_low->FixParameter(0, fitfunc2_low->GetParameter(0));
  fitfunc3_low->FixParameter(0, fitfunc2_low->GetParameter(0));
  fitfunc_low->FixParameter(2, x_turn_low[2]);
  fitfunc3_low->FixParameter(2, x_turn_high[2]);
  g_peak_low->Fit("fitfunc_low", "", "", xmin, x_turn_low[2]);
  g_peak_low->Fit("fitfunc3_low", "", "", x_turn_high[2], xmax);

  xmin = -1.5;
  xmax = 0.8;
  kev_ene = kev_ene_high;
  for (Int_t j = hist_avediff->GetXaxis()->FindBin(xmin); j < hist_avediff->GetXaxis()->FindBin(xmax); j++)
  {
    ab = (TH1D *)hist_avediff->ProjectionY("ab", j, j);
    center = abs(hist_avediff->GetXaxis()->GetBinCenter(j));
    ab->GetXaxis()->SetRangeUser(kev_ene - 2 - center, kev_ene + 1 - center);
    tmp_bin = ab->GetMaximumBin();
    ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
    max = ab->GetMean();
    g_peak_high->SetPoint(g_peak_high->GetN(), hist_avediff->GetXaxis()->GetBinCenter(j), max);
  }
  fitfunc_high = new TF1("fitfunc_high", "[0]+[1]*(x-[2]) ", -100.0, x_turn_low[2]);
  fitfunc2_high = new TF1("fitfunc2_high", "[0]", x_turn_low[2], x_turn_high[2]);
  fitfunc3_high = new TF1("fitfunc3_high", "[0]+[1]*(x-[2]) ", x_turn_high[2], 100.0);
  fitfunc_high->SetNpx(500);
  fitfunc2_high->SetNpx(500);
  fitfunc3_high->SetNpx(500);
  fitfunc_high->SetLineWidth(3);
  fitfunc2_high->SetLineWidth(3);
  fitfunc3_high->SetLineWidth(3);
  g_peak_high->Fit("fitfunc2_high", "", "", x_turn_low[2], x_turn_high[2]);
  fitfunc_high->FixParameter(0, fitfunc2_high->GetParameter(0));
  fitfunc3_high->FixParameter(0, fitfunc2_high->GetParameter(0));
  fitfunc_high->FixParameter(2, x_turn_low[2]);
  fitfunc3_high->FixParameter(2, x_turn_high[2]);
  g_peak_high->Fit("fitfunc_high", "", "", xmin, x_turn_low[2]);
  g_peak_high->Fit("fitfunc3_high", "", "", x_turn_high[2], xmax);

  for (int i = 0; i < n_point; i++)
  {
    tmpx = i * (func_xmax - func_xmin) / n_point + func_xmin;
    for (int j = 0; j < niterp; j++)
    {
      targety = kev_ene_high - j * (kev_ene_high - kev_ene_low) / niterp;
      ene_ratio = (targety - kev_ene_low) / (kev_ene_high - kev_ene_low);
      if (tmpx > x_turn_high[2])
        tmp_peak = fitfunc3_low->Eval(tmpx) + (fitfunc3_high->Eval(tmpx) - fitfunc3_low->Eval(tmpx)) * ene_ratio;
      else if (tmpx < x_turn_low[2])
        tmp_peak = fitfunc_low->Eval(tmpx) + (fitfunc_high->Eval(tmpx) - fitfunc_low->Eval(tmpx)) * ene_ratio;
      else
        tmp_peak = fitfunc2_low->Eval(tmpx) + (fitfunc2_high->Eval(tmpx) - fitfunc2_low->Eval(tmpx)) * ene_ratio;
      g_response[0][1]->SetPoint(g_response[0][1]->GetN(), tmpx, tmp_peak, targety);
      intp[j]->SetPoint(i, tmpx, tmp_peak);
    }
    for (int j = 0; j < niterp_low; j++)
    {
      targety = j * (kev_ene_low - (-5.0)) / niterp_low;
      if (tmpx > x_turn_high[2])
        tmp_peak = fitfunc3_low->Eval(tmpx) - targety;
      else if (tmpx < x_turn_low[2])
        tmp_peak = fitfunc_low->Eval(tmpx) - targety;
      else
        tmp_peak = fitfunc2_low->Eval(tmpx) - targety;
      g_response[0][1]->SetPoint(g_response[0][1]->GetN(), tmpx, tmp_peak, kev_ene_low - targety);
      intp_low[j]->SetPoint(i, tmpx, tmp_peak);
    }
    for (int j = 0; j < niterp_high; j++)
    {
      targety = j * (40.0 - kev_ene_high) / niterp_high;
      if (tmpx > x_turn_high[2])
        tmp_peak = fitfunc3_high->Eval(tmpx) + targety;
      else if (tmpx < x_turn_low[2])
        tmp_peak = fitfunc_high->Eval(tmpx) + targety;
      else
        tmp_peak = fitfunc2_high->Eval(tmpx) + targety;
      g_response[0][1]->SetPoint(g_response[0][1]->GetN(), tmpx, tmp_peak, kev_ene_high + targety);
      intp_high[j]->SetPoint(i, tmpx, tmp_peak);
    }
  }

  TFile *fout = new TFile("doi_resp_temp.root", "RECREATE");
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      g_response[i][j]->Write(Form("g_resp_pt%dal%d", i + 1, j + 1));
  fout->Close();
}
