
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

void make_gaploss_func()
{
  TString no_file = "no2022_01";
  gStyle->SetOptStat(0);
  TString am_infilename = "Am241.root";
  TString ba_infilename = "Ba133.root";
  TString fe_infilename = "Fe55.root";
  TString outFilename = "./resp_" + no_file + ".root";
  TString resp_name = "resp_" + no_file;

  Double_t am_pt_range_peak1[] = {9.8, 9.8, 8.8};
  Double_t am_pt_range_peak2[] = {12.0, 12.0, 12.0};
  Double_t am_pt_range_peak3[] = {16.0, 15.8, 15.6};
  Double_t am_pt_range_peak4[] = {20.0, 20.0, 22.0};
  Double_t am_pt_diff1[] = {2.0, 3.0, 3.5};
  Double_t am_pt_diff2[] = {2.5, 3.65, 4.0};
  Double_t am_pt_diff3[] = {2.5, 4.0, 4.5};
  Double_t am_pt_diff4[] = {4.0, 5.5, 6.0};
  Double_t offset[] = {0.0, 0.0, 0.5};

  Double_t ba_pt_range_peak1[] = {26.0, 25.6, 25.6};
  Double_t ba_pt_range_peak2[] = {7.0, 6.5, 5.5};
  Double_t ba_pt_range_peak3[] = {3.5, 3.5, 3.5};
  Double_t fe_pt_range_peak[] = {1.0, 1.0, 1.0};

  Double_t am_al_range_peak1[] = {10.4, 10.4, 9.6};
  Double_t am_al_range_peak2[] = {14.4, 13.5, 12.4};
  Double_t am_al_range_peak3[] = {17.4, 17.4, 15.6};
  Double_t am_al_range_peak4[] = {22.4, 22.4, 20.4};
  Double_t ba_al_range_peak1[] = {27.2, 24.0, 24.2};
  Double_t ba_al_range_peak2[] = {8.5, 8.0, 8.0};
  Double_t ba_al_range_peak3[] = {4.5, 4.0, 4.0};
  Double_t fe_al_range_peak[] = {3.5, 3.5, 3.5};
  Double_t am_al_diff1[] = {0.5, 1.0, 1.5};
  Double_t am_al_diff2[] = {0.8, 1.5, 1.5};
  Double_t am_al_diff3[] = {1.0, 2.0, 2.0};
  Double_t am_al_diff4[] = {2.0, 2.5, 3.0};

  TFile *fin_am = new TFile(am_infilename, "read");
  TFile *fin_ba = new TFile(ba_infilename, "read");
  TFile *fin_fe = new TFile(fe_infilename, "read");

  Int_t ave_width = 0;
  Int_t colorlist[] = {kRed, kOrange, kGreen + 2, kBlue};
  Double_t weight_energy = (34.920 * 5.99 + 34.987 * 11.6 + 35.818 * 3.58) / (5.99 + 11.6 + 3.58);
  Double_t kev_ene[] = {13.946, 17.751, 20.784, 26.3448, 30.973, weight_energy - 23.173, 7.80, weight_energy};
  TString pdfname = outFilename;
  pdfname.ReplaceAll(".root", ".pdf");
  //
  Int_t npeaks = sizeof(kev_ene) / sizeof(kev_ene[0]);

  TH2D *hist_diffene_pt_adj_am;
  TH2D *hist_diffene_pt_adj_ba;
  TH2D *hist_diffene_pt_adj_fe;
  TH2D *hist_diffene_al_adj_am;
  TH2D *hist_diffene_al_adj_ba;
  TH2D *hist_diffene_al_adj_fe;

  TCanvas *c1 = new TCanvas("c1", "c1", 1500, 800);
  c1->Divide(2, 1);
  c1->Print(pdfname + "[", "pdf");
  TGraph2D *g_response_pt[3];
  TGraph2D *g_response_al[3];

  for (int i = 0; i < 3; i++)
  {

    vector<double> x_peak_pt_am1, y_peak_pt_am1;
    vector<double> x_peak_pt_am2, y_peak_pt_am2;
    vector<double> x_peak_pt_am3, y_peak_pt_am3;
    vector<double> x_peak_pt_am4, y_peak_pt_am4;
    vector<double> x_peak_pt_ba1, y_peak_pt_ba1;
    vector<double> x_peak_pt_ba2, y_peak_pt_ba2;
    vector<double> x_peak_pt_ba3, y_peak_pt_ba3;
    vector<double> x_peak_pt_fe, y_peak_pt_fe;

    vector<double> x_peak_al_am1, y_peak_al_am1;
    vector<double> x_peak_al_am2, y_peak_al_am2;
    vector<double> x_peak_al_am3, y_peak_al_am3;
    vector<double> x_peak_al_am4, y_peak_al_am4;
    vector<double> x_peak_al_ba1, y_peak_al_ba1;
    vector<double> x_peak_al_ba2, y_peak_al_ba2;
    vector<double> x_peak_al_ba3, y_peak_al_ba3;
    vector<double> x_peak_al_fe, y_peak_al_fe;

    hist_diffene_pt_adj_am = (TH2D *)fin_am->Get(Form("hist_diffene_pt_adj_region%d", i));
    hist_diffene_pt_adj_ba = (TH2D *)fin_ba->Get(Form("hist_diffene_pt_adj_region%d", i));
    hist_diffene_pt_adj_fe = (TH2D *)fin_fe->Get(Form("hist_diffene_pt_adj_region%d", i));
    hist_diffene_al_adj_am = (TH2D *)fin_am->Get(Form("hist_diffene_al_adj_region%d", i));
    hist_diffene_al_adj_ba = (TH2D *)fin_ba->Get(Form("hist_diffene_al_adj_region%d", i));
    hist_diffene_al_adj_fe = (TH2D *)fin_fe->Get(Form("hist_diffene_al_adj_region%d", i));

    int tmp_bin;
    double max, mean, center, check_range;
    TH1D *ab;
    TH2D *h2_target;
    h2_target = hist_diffene_pt_adj_am;
    check_range = am_pt_range_peak1[i];
    for (int j = h2_target->GetXaxis()->FindBin(-check_range); j < h2_target->GetXaxis()->FindBin(check_range); j++)
    {
      ab = (TH1D *)h2_target->ProjectionY("ab", j, j);
      center = abs(h2_target->GetXaxis()->GetBinCenter(j));
      ab->GetXaxis()->SetRangeUser(9.7, 12.6);
      tmp_bin = ab->GetMaximumBin();
      ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      max = ab->GetMean();
      x_peak_pt_am1.push_back(h2_target->GetXaxis()->GetBinCenter(j));
      y_peak_pt_am1.push_back(max);
    }

    h2_target = hist_diffene_pt_adj_am;
    check_range = am_pt_range_peak2[i];
    for (int j = h2_target->GetXaxis()->FindBin(-check_range); j < h2_target->GetXaxis()->FindBin(check_range); j++)
    {
      ab = (TH1D *)h2_target->ProjectionY("ab", j, j);
      center = abs(h2_target->GetXaxis()->GetBinCenter(j));
      ab->GetXaxis()->SetRangeUser(12.6, 14.0);
      tmp_bin = ab->GetMaximumBin();
      ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      max = ab->GetMean();
      x_peak_pt_am2.push_back(h2_target->GetXaxis()->GetBinCenter(j));
      y_peak_pt_am2.push_back(max);
    }
    h2_target = hist_diffene_pt_adj_am;
    check_range = am_pt_range_peak3[i];
    for (int j = h2_target->GetXaxis()->FindBin(-check_range); j < h2_target->GetXaxis()->FindBin(check_range); j++)
    {
      ab = (TH1D *)h2_target->ProjectionY("ab", j, j);
      center = abs(h2_target->GetXaxis()->GetBinCenter(j));
      ab->GetXaxis()->SetRangeUser(15.0, 18.0);
      tmp_bin = ab->GetMaximumBin();
      ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      max = ab->GetMean();
      x_peak_pt_am3.push_back(h2_target->GetXaxis()->GetBinCenter(j));
      y_peak_pt_am3.push_back(max);
    }
    h2_target = hist_diffene_pt_adj_am;
    check_range = am_pt_range_peak4[i];
    for (int j = h2_target->GetXaxis()->FindBin(-check_range); j < h2_target->GetXaxis()->FindBin(check_range); j++)
    {
      ab = (TH1D *)h2_target->ProjectionY("ab", j, j);
      center = abs(h2_target->GetXaxis()->GetBinCenter(j));
      ab->GetXaxis()->SetRangeUser(18.0, 23.0);
      tmp_bin = ab->GetMaximumBin();
      ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      max = ab->GetMean();
      x_peak_pt_am4.push_back(h2_target->GetXaxis()->GetBinCenter(j));
      y_peak_pt_am4.push_back(max);
    }
    h2_target = hist_diffene_pt_adj_ba;
    check_range = ba_pt_range_peak1[i];
    for (int j = h2_target->GetXaxis()->FindBin(-check_range); j < h2_target->GetXaxis()->FindBin(check_range); j++)
    {
      ab = (TH1D *)h2_target->ProjectionY("ab", j, j);
      center = abs(h2_target->GetXaxis()->GetBinCenter(j));
      ab->GetXaxis()->SetRangeUser(23., 30.5);
      tmp_bin = ab->GetMaximumBin();
      ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      max = ab->GetMean();
      x_peak_pt_ba1.push_back(h2_target->GetXaxis()->GetBinCenter(j));
      y_peak_pt_ba1.push_back(max);
    }
    h2_target = hist_diffene_pt_adj_ba;
    check_range = ba_pt_range_peak2[i];
    for (int j = h2_target->GetXaxis()->FindBin(-check_range); j < h2_target->GetXaxis()->FindBin(check_range); j++)
    {
      ab = (TH1D *)h2_target->ProjectionY("ab", j, j);
      center = abs(h2_target->GetXaxis()->GetBinCenter(j));
      ab->GetXaxis()->SetRangeUser(7.8, 13);
      tmp_bin = ab->GetMaximumBin();
      ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      max = ab->GetMean();
      x_peak_pt_ba2.push_back(h2_target->GetXaxis()->GetBinCenter(j));
      y_peak_pt_ba2.push_back(max);
    }

    h2_target = hist_diffene_al_adj_am;
    check_range = am_al_range_peak1[i];
    for (int j = h2_target->GetXaxis()->FindBin(-check_range); j < h2_target->GetXaxis()->FindBin(check_range); j++)
    {
      ab = (TH1D *)h2_target->ProjectionY("ab", j, j);
      center = abs(h2_target->GetXaxis()->GetBinCenter(j));
      ab->GetXaxis()->SetRangeUser(9.7, 15.0);
      tmp_bin = ab->GetMaximumBin();
      ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      max = ab->GetMean();
      x_peak_al_am1.push_back(h2_target->GetXaxis()->GetBinCenter(j));
      y_peak_al_am1.push_back(max);
    }

    h2_target = hist_diffene_al_adj_am;
    check_range = am_al_range_peak2[i];
    for (int j = h2_target->GetXaxis()->FindBin(-check_range); j < h2_target->GetXaxis()->FindBin(check_range); j++)
    {
      ab = (TH1D *)h2_target->ProjectionY("ab", j, j);
      center = abs(h2_target->GetXaxis()->GetBinCenter(j));
      ab->GetXaxis()->SetRangeUser(16., 19.0);
      tmp_bin = ab->GetMaximumBin();
      ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      max = ab->GetMean();
      x_peak_al_am2.push_back(h2_target->GetXaxis()->GetBinCenter(j));
      y_peak_al_am2.push_back(max);
    }
    h2_target = hist_diffene_al_adj_am;
    check_range = am_al_range_peak3[i];
    for (int j = h2_target->GetXaxis()->FindBin(-check_range); j < h2_target->GetXaxis()->FindBin(check_range); j++)
    {
      ab = (TH1D *)h2_target->ProjectionY("ab", j, j);
      center = abs(h2_target->GetXaxis()->GetBinCenter(j));
      ab->GetXaxis()->SetRangeUser(20, 24.0);
      tmp_bin = ab->GetMaximumBin();
      ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      max = ab->GetMean();
      x_peak_al_am3.push_back(h2_target->GetXaxis()->GetBinCenter(j));
      y_peak_al_am3.push_back(max);
    }
    h2_target = hist_diffene_al_adj_am;
    check_range = am_al_range_peak4[i];
    for (int j = h2_target->GetXaxis()->FindBin(-check_range); j < h2_target->GetXaxis()->FindBin(check_range); j++)
    {
      ab = (TH1D *)h2_target->ProjectionY("ab", j, j);
      center = abs(h2_target->GetXaxis()->GetBinCenter(j));
      ab->GetXaxis()->SetRangeUser(24.5, 30.0);
      tmp_bin = ab->GetMaximumBin();
      ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      max = ab->GetMean();
      x_peak_al_am4.push_back(h2_target->GetXaxis()->GetBinCenter(j));
      y_peak_al_am4.push_back(max);
    }
    h2_target = hist_diffene_al_adj_ba;
    check_range = ba_al_range_peak1[i];
    for (int j = h2_target->GetXaxis()->FindBin(-check_range); j < h2_target->GetXaxis()->FindBin(check_range); j++)
    {
      ab = (TH1D *)h2_target->ProjectionY("ab", j, j);
      center = abs(h2_target->GetXaxis()->GetBinCenter(j));
      ab->GetXaxis()->SetRangeUser(29., 34.0);
      tmp_bin = ab->GetMaximumBin();
      ab->GetXaxis()->SetRange(tmp_bin - ave_width * 2, tmp_bin + ave_width * 2);
      max = ab->GetMean();
      x_peak_al_ba1.push_back(h2_target->GetXaxis()->GetBinCenter(j));
      y_peak_al_ba1.push_back(max);
    }
    h2_target = hist_diffene_al_adj_ba;
    check_range = ba_al_range_peak2[i];
    for (int j = h2_target->GetXaxis()->FindBin(-check_range); j < h2_target->GetXaxis()->FindBin(check_range); j++)
    {
      ab = (TH1D *)h2_target->ProjectionY("ab", j, j);
      center = abs(h2_target->GetXaxis()->GetBinCenter(j));
      ab->GetXaxis()->SetRangeUser(10.4, 14);
      tmp_bin = ab->GetMaximumBin();
      ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      max = ab->GetMean();
      x_peak_al_ba2.push_back(h2_target->GetXaxis()->GetBinCenter(j));
      y_peak_al_ba2.push_back(max);
    }
    h2_target = hist_diffene_al_adj_ba;
    check_range = ba_al_range_peak3[i];
    for (int j = h2_target->GetXaxis()->FindBin(-check_range); j < h2_target->GetXaxis()->FindBin(check_range); j++)
    {
      ab = (TH1D *)h2_target->ProjectionY("ab", j, j);
      center = abs(h2_target->GetXaxis()->GetBinCenter(j));
      ab->GetXaxis()->SetRangeUser(6.8, 10);
      tmp_bin = ab->GetMaximumBin();
      ab->GetXaxis()->SetRange(tmp_bin - ave_width, tmp_bin + ave_width);
      max = ab->GetMean();
      x_peak_al_ba3.push_back(h2_target->GetXaxis()->GetBinCenter(j));
      y_peak_al_ba3.push_back(max);
    }

    TGraph *g_pt_am1 = new TGraph(x_peak_pt_am1.size(), &x_peak_pt_am1[0], &y_peak_pt_am1[0]);
    TGraph *g_pt_am2 = new TGraph(x_peak_pt_am2.size(), &x_peak_pt_am2[0], &y_peak_pt_am2[0]);
    TGraph *g_pt_am3 = new TGraph(x_peak_pt_am3.size(), &x_peak_pt_am3[0], &y_peak_pt_am3[0]);
    TGraph *g_pt_am4 = new TGraph(x_peak_pt_am4.size(), &x_peak_pt_am4[0], &y_peak_pt_am4[0]);
    TGraph *g_pt_ba1 = new TGraph(x_peak_pt_ba1.size(), &x_peak_pt_ba1[0], &y_peak_pt_ba1[0]);
    TGraph *g_pt_ba2 = new TGraph(x_peak_pt_ba2.size(), &x_peak_pt_ba2[0], &y_peak_pt_ba2[0]);

    TGraph *g_al_am1 = new TGraph(x_peak_al_am1.size(), &x_peak_al_am1[0], &y_peak_al_am1[0]);
    TGraph *g_al_am2 = new TGraph(x_peak_al_am2.size(), &x_peak_al_am2[0], &y_peak_al_am2[0]);
    TGraph *g_al_am3 = new TGraph(x_peak_al_am3.size(), &x_peak_al_am3[0], &y_peak_al_am3[0]);
    TGraph *g_al_am4 = new TGraph(x_peak_al_am4.size(), &x_peak_al_am4[0], &y_peak_al_am4[0]);
    TGraph *g_al_ba1 = new TGraph(x_peak_al_ba1.size(), &x_peak_al_ba1[0], &y_peak_al_ba1[0]);
    TGraph *g_al_ba2 = new TGraph(x_peak_al_ba2.size(), &x_peak_al_ba2[0], &y_peak_al_ba2[0]);
    TGraph *g_al_ba3 = new TGraph(x_peak_al_ba3.size(), &x_peak_al_ba3[0], &y_peak_al_ba3[0]);

    TF1 *fitfunc_pt_am1 = new TF1("fitfunc_pt_am1", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4]) - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4])  +[3]*x*x*x*x", -kev_ene[0], kev_ene[0]);
    fitfunc_pt_am1->FixParameter(0, kev_ene[0]);
    fitfunc_pt_am1->SetParameter(1, 0);
    fitfunc_pt_am1->SetParameter(3, 0);
    fitfunc_pt_am1->SetParameter(4, 10);
    fitfunc_pt_am1->SetParLimits(4, 0, 1000);
    fitfunc_pt_am1->SetNpx(500);
    g_pt_am1->Fit("fitfunc_pt_am1", "EMR", "", -am_pt_range_peak1[i], am_pt_range_peak1[i]);
    g_pt_am1->Fit("fitfunc_pt_am1", "EMR", "", -am_pt_range_peak1[i], am_pt_range_peak1[i]);

    TF1 *fitfunc_pt_am2 = new TF1("fitfunc_pt_am2", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4])  +[3]*x*x*x*x", -kev_ene[1], kev_ene[1]);
    fitfunc_pt_am2->FixParameter(0, kev_ene[1]);
    fitfunc_pt_am2->SetParameter(1, 0);
    fitfunc_pt_am2->SetParameter(3, 0);
    fitfunc_pt_am2->SetParameter(4, 12);
    fitfunc_pt_am2->SetParLimits(4, 0, 1000);
    fitfunc_pt_am2->SetNpx(500);
    g_pt_am2->Fit("fitfunc_pt_am2", "EMR", "", -am_pt_range_peak2[i], am_pt_range_peak2[i]);
    g_pt_am2->Fit("fitfunc_pt_am2", "EMR", "", -am_pt_range_peak2[i], am_pt_range_peak2[i]);

    TF1 *fitfunc_pt_am3 = new TF1("fitfunc_pt_am3", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4])  +[3]*x*x*x*x", -kev_ene[2], kev_ene[2]);
    fitfunc_pt_am3->FixParameter(0, kev_ene[2]);
    fitfunc_pt_am3->SetParameter(1, 0);
    fitfunc_pt_am3->SetParameter(3, 0);
    fitfunc_pt_am3->SetParameter(4, 17);
    fitfunc_pt_am3->SetParLimits(4, 0, 1000);
    fitfunc_pt_am3->SetNpx(500);
    g_pt_am3->Fit("fitfunc_pt_am3", "EMR", "", -am_pt_range_peak3[i], am_pt_range_peak3[i]);
    g_pt_am3->Fit("fitfunc_pt_am3", "EMR", "", -am_pt_range_peak3[i], am_pt_range_peak3[i]);

    TF1 *fitfunc_pt_am4 = new TF1("fitfunc_pt_am4", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4])  +[3]*x*x*x*x", -kev_ene[3], kev_ene[3]);
    fitfunc_pt_am4->FixParameter(0, kev_ene[3]);
    fitfunc_pt_am4->SetParameter(1, 0);
    fitfunc_pt_am4->SetParameter(3, 0);
    fitfunc_pt_am4->SetParameter(4, 22);
    fitfunc_pt_am4->SetParLimits(4, 0, 1000);
    fitfunc_pt_am4->SetNpx(500);
    g_pt_am4->Fit("fitfunc_pt_am4", "EMR", "", -am_pt_range_peak4[i], am_pt_range_peak4[i]);
    g_pt_am4->Fit("fitfunc_pt_am4", "EMR", "", -am_pt_range_peak4[i], am_pt_range_peak4[i]);

    TF1 *fitfunc_pt_ba1 = new TF1("fitfunc_pt_ba1", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4]) + +[3]*x*x*x*x", -kev_ene[4], kev_ene[4]);
    fitfunc_pt_ba1->FixParameter(0, kev_ene[4]);
    fitfunc_pt_ba1->SetParameter(1, 0);
    fitfunc_pt_ba1->SetParameter(3, 0);
    fitfunc_pt_ba1->SetParameter(4, 25);
    fitfunc_pt_ba1->SetParLimits(4, 0, 1000);
    fitfunc_pt_ba1->SetNpx(500);
    g_pt_ba1->Fit("fitfunc_pt_ba1", "EMR", "", -ba_pt_range_peak1[i], ba_pt_range_peak1[i]);
    g_pt_ba1->Fit("fitfunc_pt_ba1", "EMR", "", -ba_pt_range_peak1[i], ba_pt_range_peak1[i]);

    TF1 *fitfunc_pt_ba2 = new TF1("fitfunc_pt_ba2", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4]) + +[3]*x*x*x*x", -kev_ene[5], kev_ene[5]);
    fitfunc_pt_ba2->FixParameter(0, kev_ene[5]);
    fitfunc_pt_ba2->SetParameter(1, 0);
    fitfunc_pt_ba2->FixParameter(3, 0);
    fitfunc_pt_ba2->SetParameter(4, 25);
    fitfunc_pt_ba2->SetParLimits(4, 0, 1000);
    fitfunc_pt_ba2->SetNpx(500);
    g_pt_ba2->Fit("fitfunc_pt_ba2", "EMR", "", -ba_pt_range_peak2[i], ba_pt_range_peak2[i]);
    g_pt_ba2->Fit("fitfunc_pt_ba2", "EMR", "", -ba_pt_range_peak2[i], ba_pt_range_peak2[i]);

    TF1 *fitfunc_pt_bound_max = new TF1("fitfunc_pt_bound_max", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4]) + +[3]*x*x*x*x + [5]", -30, 30);
    fitfunc_pt_bound_max->FixParameter(3, fitfunc_pt_ba1->GetParameter(3));
    fitfunc_pt_bound_max->FixParameter(2, fitfunc_pt_ba1->GetParameter(2));
    fitfunc_pt_bound_max->FixParameter(1, fitfunc_pt_ba1->GetParameter(1));
    fitfunc_pt_bound_max->FixParameter(0, fitfunc_pt_ba1->GetParameter(0));
    fitfunc_pt_bound_max->FixParameter(4, fitfunc_pt_ba1->GetParameter(4));
    fitfunc_pt_bound_max->FixParameter(5, 10.0);
    fitfunc_pt_bound_max->SetNpx(500);

    TF1 *fitfunc_pt_bound_min = new TF1("fitfunc_pt_bound_min", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4]) + +[3]*x*x*x*x + [5]", -kev_ene[5], kev_ene[5]);
    fitfunc_pt_bound_min->FixParameter(3, fitfunc_pt_ba2->GetParameter(3));
    fitfunc_pt_bound_min->FixParameter(2, fitfunc_pt_ba2->GetParameter(2));
    fitfunc_pt_bound_min->FixParameter(1, fitfunc_pt_ba2->GetParameter(1));
    fitfunc_pt_bound_min->FixParameter(0, fitfunc_pt_ba2->GetParameter(0));
    fitfunc_pt_bound_min->FixParameter(4, fitfunc_pt_ba2->GetParameter(4));
    fitfunc_pt_bound_min->FixParameter(5, -7.0);
    fitfunc_pt_bound_min->SetNpx(500);

    TF1 *fitfunc_al_am1 = new TF1("fitfunc_al_am1", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4]) - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4])  +[3]*x*x*x*x", -kev_ene[0], kev_ene[0]);
    fitfunc_al_am1->SetParameter(0, kev_ene[0]);
    fitfunc_al_am1->SetParameter(1, 0);
    fitfunc_al_am1->SetParameter(3, 0);
    fitfunc_al_am1->SetParameter(4, 10);
    fitfunc_al_am1->SetParLimits(4, 0, 1000);
    fitfunc_al_am1->SetNpx(500);
    g_al_am1->Fit("fitfunc_al_am1", "EMR", "", -am_al_range_peak1[i], am_al_range_peak1[i]);
    g_al_am1->Fit("fitfunc_al_am1", "EMR", "", -am_al_range_peak1[i], am_al_range_peak1[i]);

    TF1 *fitfunc_al_am2 = new TF1("fitfunc_al_am2", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4])  +[3]*x*x*x*x", -kev_ene[1], kev_ene[1]);
    fitfunc_al_am2->SetParameter(0, kev_ene[1]);
    fitfunc_al_am2->SetParameter(1, 0);
    fitfunc_al_am2->SetParameter(3, 0);
    fitfunc_al_am2->SetParameter(4, 12);
    fitfunc_al_am2->SetParLimits(4, 0, 1000);
    fitfunc_al_am2->SetNpx(500);
    g_al_am2->Fit("fitfunc_al_am2", "EMR", "", -am_al_range_peak2[i], am_al_range_peak2[i]);
    g_al_am2->Fit("fitfunc_al_am2", "EMR", "", -am_al_range_peak2[i], am_al_range_peak2[i]);

    TF1 *fitfunc_al_am3 = new TF1("fitfunc_al_am3", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4])  +[3]*x*x*x*x", -kev_ene[2], kev_ene[2]);
    fitfunc_al_am3->SetParameter(0, kev_ene[2]);
    fitfunc_al_am3->SetParameter(1, 0);
    fitfunc_al_am3->SetParameter(3, 0);
    fitfunc_al_am3->SetParameter(4, 17);
    fitfunc_al_am3->SetParLimits(4, 0, 1000);
    fitfunc_al_am3->SetNpx(500);
    g_al_am3->Fit("fitfunc_al_am3", "EMR", "", -am_al_range_peak3[i], am_al_range_peak3[i]);
    g_al_am3->Fit("fitfunc_al_am3", "EMR", "", -am_al_range_peak3[i], am_al_range_peak3[i]);

    TF1 *fitfunc_al_am4 = new TF1("fitfunc_al_am4", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4])  +[3]*x*x*x*x", -kev_ene[3], kev_ene[3]);
    fitfunc_al_am4->SetParameter(0, kev_ene[3]);
    fitfunc_al_am4->SetParameter(1, 0);
    fitfunc_al_am4->SetParameter(3, 0);
    fitfunc_al_am4->SetParameter(4, 22);
    fitfunc_al_am4->SetParLimits(4, 0, 1000);
    fitfunc_al_am4->SetNpx(500);
    g_al_am4->Fit("fitfunc_al_am4", "EMR", "", -am_al_range_peak4[i], am_al_range_peak4[i]);
    g_al_am4->Fit("fitfunc_al_am4", "EMR", "", -am_al_range_peak4[i], am_al_range_peak4[i]);

    TF1 *fitfunc_al_ba1 = new TF1("fitfunc_al_ba1", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4]) + +[3]*x*x*x*x", -kev_ene[4], kev_ene[4]);
    fitfunc_al_ba1->SetParameter(0, kev_ene[4]);
    fitfunc_al_ba1->SetParameter(1, 0);
    fitfunc_al_ba1->SetParameter(3, 0);
    fitfunc_al_ba1->SetParameter(4, 25);
    fitfunc_al_ba1->SetParLimits(4, 0, 1000);
    fitfunc_al_ba1->SetNpx(500);
    g_al_ba1->Fit("fitfunc_al_ba1", "EMR", "", -ba_al_range_peak1[i], ba_al_range_peak1[i]);
    g_al_ba1->Fit("fitfunc_al_ba1", "EMR", "", -ba_al_range_peak1[i], ba_al_range_peak1[i]);

    TF1 *fitfunc_al_ba2 = new TF1("fitfunc_al_ba2", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4])  +[3]*x*x*x*x", -kev_ene[5], kev_ene[5]);
    fitfunc_al_ba2->SetParameter(0, kev_ene[5]);
    fitfunc_al_ba2->FixParameter(2, 0.1);
    fitfunc_al_ba2->FixParameter(1, 0.01);
    fitfunc_al_ba2->FixParameter(3, 0);
    fitfunc_al_ba2->FixParameter(4, 0);
    fitfunc_al_ba2->SetParLimits(4, 0, 10000);
    fitfunc_al_ba2->SetNpx(500);
    g_al_ba2->Fit("fitfunc_al_ba2", "", "", -ba_al_range_peak2[i], ba_al_range_peak2[i]);
    fitfunc_al_ba2->ReleaseParameter(2);
    fitfunc_al_ba2->ReleaseParameter(4);
    fitfunc_al_ba2->SetParLimits(4, 0, 10000);
    g_al_ba2->Fit("fitfunc_al_ba2", "EMR", "", -ba_al_range_peak2[i], ba_al_range_peak2[i]);
    g_al_ba2->Fit("fitfunc_al_ba2", "EMR", "", -ba_al_range_peak2[i], ba_al_range_peak2[i]);

    TF1 *fitfunc_al_ba3 = new TF1("fitfunc_al_ba3", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4]) +[3]*x*x*x*x", -kev_ene[6], kev_ene[6]);
    if (i > 0)
      fitfunc_al_ba3->FixParameter(0, kev_ene[6]);
    fitfunc_al_ba3->SetParameter(0, kev_ene[6]);
    fitfunc_al_ba3->SetParameter(1, 5.0);
    fitfunc_al_ba3->SetParameter(2, 1.0);
    fitfunc_al_ba3->FixParameter(3, 0);
    fitfunc_al_ba3->FixParameter(2, 0);
    fitfunc_al_ba3->SetParameter(4, 8);
    fitfunc_al_ba3->SetParLimits(4, 0, 100);
    fitfunc_al_ba3->SetNpx(500);
    g_al_ba3->Fit("fitfunc_al_ba3", "EMR", "", -ba_al_range_peak3[i], ba_al_range_peak3[i]);
    g_al_ba3->Fit("fitfunc_al_ba3", "EMR", "", -ba_al_range_peak3[i], ba_al_range_peak3[i]);

    TF1 *fitfunc_al_bound_max = new TF1("fitfunc_al_bound_max", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4]) + +[3]*x*x*x*x + [5]", -30, 30);
    fitfunc_al_bound_max->FixParameter(3, fitfunc_al_ba1->GetParameter(3));
    fitfunc_al_bound_max->FixParameter(2, fitfunc_al_ba1->GetParameter(2));
    fitfunc_al_bound_max->FixParameter(1, fitfunc_al_ba1->GetParameter(1));
    fitfunc_al_bound_max->FixParameter(0, fitfunc_al_ba1->GetParameter(0));
    fitfunc_al_bound_max->FixParameter(4, fitfunc_al_ba1->GetParameter(4));
    fitfunc_al_bound_max->FixParameter(5, 10.0);
    fitfunc_al_bound_max->SetNpx(500);

    TF1 *fitfunc_al_bound_min = new TF1("fitfunc_al_bound_min", "-[1]*[0]*[0] - [2]*sqrt([0]*[0]+[4])   - [3]*[0]*[0]*[0]*[0]+ [0] +[1]*x*x +[2]*sqrt(x*x+[4]) + +[3]*x*x*x*x + [5]", -kev_ene[5], kev_ene[5]);
    fitfunc_al_bound_min->FixParameter(3, fitfunc_al_ba2->GetParameter(3));
    fitfunc_al_bound_min->FixParameter(2, fitfunc_al_ba2->GetParameter(2));
    fitfunc_al_bound_min->FixParameter(1, fitfunc_al_ba2->GetParameter(1));
    fitfunc_al_bound_min->FixParameter(0, fitfunc_al_ba2->GetParameter(0));
    fitfunc_al_bound_min->FixParameter(4, fitfunc_al_ba2->GetParameter(4));
    fitfunc_al_bound_min->FixParameter(5, -7.0);
    fitfunc_al_bound_min->SetNpx(500);

    g_response_pt[i] = new TGraph2D();
    g_response_al[i] = new TGraph2D();

    Int_t n_point = 500;
    Int_t tmp_point = 0;
    Double_t tmpx = 0.0;

    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[5] / n_point - kev_ene[5];
      g_response_pt[i]->SetPoint(tmp_point, tmpx, fitfunc_pt_bound_min->Eval(tmpx), kev_ene[5] - 7.0);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[0] / n_point - kev_ene[0];
      g_response_pt[i]->SetPoint(tmp_point, tmpx, fitfunc_pt_am1->Eval(tmpx), kev_ene[0]);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[1] / n_point - kev_ene[1];
      g_response_pt[i]->SetPoint(tmp_point, tmpx, fitfunc_pt_am2->Eval(tmpx), kev_ene[1]);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[2] / n_point - kev_ene[2];
      g_response_pt[i]->SetPoint(tmp_point, tmpx, fitfunc_pt_am3->Eval(tmpx), kev_ene[2]);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[3] / n_point - kev_ene[3];
      g_response_pt[i]->SetPoint(tmp_point, tmpx, fitfunc_pt_am4->Eval(tmpx), kev_ene[3]);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[4] / n_point - kev_ene[4];
      g_response_pt[i]->SetPoint(tmp_point, tmpx, fitfunc_pt_ba1->Eval(tmpx), kev_ene[4]);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[5] / n_point - kev_ene[5];
      g_response_pt[i]->SetPoint(tmp_point, tmpx, fitfunc_pt_ba2->Eval(tmpx), kev_ene[5]);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * (kev_ene[4] + 10.0) / n_point - (kev_ene[4] + 10.0);
      g_response_pt[i]->SetPoint(tmp_point, tmpx, fitfunc_pt_bound_max->Eval(tmpx), kev_ene[4] + 10.0);
      tmp_point = tmp_point + 1;
    }

    tmp_point = 0.0;
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[5] / n_point - kev_ene[5];
      g_response_al[i]->SetPoint(tmp_point, tmpx, fitfunc_al_bound_min->Eval(tmpx), kev_ene[5] - 7.0);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[0] / n_point - kev_ene[0];
      g_response_al[i]->SetPoint(tmp_point, tmpx, fitfunc_al_am1->Eval(tmpx), kev_ene[0]);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[1] / n_point - kev_ene[1];
      g_response_al[i]->SetPoint(tmp_point, tmpx, fitfunc_al_am2->Eval(tmpx), kev_ene[1]);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[2] / n_point - kev_ene[2];
      g_response_al[i]->SetPoint(tmp_point, tmpx, fitfunc_al_am3->Eval(tmpx), kev_ene[2]);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[3] / n_point - kev_ene[3];
      g_response_al[i]->SetPoint(tmp_point, tmpx, fitfunc_al_am4->Eval(tmpx), kev_ene[3]);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[4] / n_point - kev_ene[4];
      g_response_al[i]->SetPoint(tmp_point, tmpx, fitfunc_al_ba1->Eval(tmpx), kev_ene[4]);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[5] / n_point - kev_ene[5];
      g_response_al[i]->SetPoint(tmp_point, tmpx, fitfunc_al_ba2->Eval(tmpx), kev_ene[5]);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * kev_ene[6] / n_point - kev_ene[6];
      g_response_al[i]->SetPoint(tmp_point, tmpx, fitfunc_al_ba3->Eval(tmpx), kev_ene[6]);
      tmp_point = tmp_point + 1;
    }
    for (int k = 0; k < 2 * n_point + 1; k++)
    {
      tmpx = k * (kev_ene[4] + 10.0) / n_point - (kev_ene[4] + 10.0);
      g_response_al[i]->SetPoint(tmp_point, tmpx, fitfunc_al_bound_max->Eval(tmpx), kev_ene[4] + 10.0);
      tmp_point = tmp_point + 1;
    }
  }

  cout << "Output file ...." << endl;
  c1->Print(pdfname + "]", "pdf");

  TFile *fout = new TFile(outFilename, "RECREATE");
  for (int i = 0; i < 3; i++)
    g_response_pt[i]->Write(Form("resp_pt_region%d", i));
  for (int i = 0; i < 3; i++)
    g_response_al[i]->Write(Form("resp_al_region%d", i));
  fout->Close();
}
