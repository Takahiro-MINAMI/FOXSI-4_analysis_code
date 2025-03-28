//
// DoI loss correction
// _base_cal_gapmod.root -> _base_cal_gapmod_doimod.root
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
#include <string.h>

using namespace std;

TGraph2D *g_response_pt1al1;
TGraph2D *g_response_pt1al2;
TGraph2D *g_response_pt2al1;
TGraph2D *g_response_pt2al2;

int mod_doiloss(const char *inFileName, const char *calFileName);
double reconstruct_energypt1al1(double ave, double diff);
double reconstruct_energypt2al1(double ave, double diff);
double reconstruct_energypt1al2(double ave, double diff);
double reconstruct_energypt2al2(double ave, double diff);

int main(int argc, char **argv)
{

	if (argc == 3)
	{
		mod_doiloss(argv[1], argv[2]);
	}
	else
	{
		cout << "USAGE: Please input the names of the file to analyze and doi_loss response file" << endl;
		cout << "./path4_mod_doiloss 'input filename' 'DoI response filename'" << endl;
		return -1;
	}
	return 0;
}

double reconstruct_energypt1al1(double ave, double diff)
{
	return g_response_pt1al1->Interpolate(diff, ave);
}
double reconstruct_energypt1al2(double ave, double diff)
{
	return g_response_pt1al2->Interpolate(diff, ave);
}

double reconstruct_energypt2al1(double ave, double diff)
{
	return g_response_pt2al1->Interpolate(diff, ave);
}

double reconstruct_energypt2al2(double ave, double diff)
{
	return g_response_pt2al2->Interpolate(diff, ave);
}

int mod_doiloss(const char *inFileName, const char *calFileName)
{
	gStyle->SetOptStat(0);
	gStyle->SetPalette(53);

	// Load DoI-Loss Response function (ROOT, TGraph2D) -------------------------------------
	TFile *fcal = new TFile(calFileName);
	g_response_pt1al1 = (TGraph2D *)fcal->Get("g_resp_pt1al1");
	g_response_pt1al2 = (TGraph2D *)fcal->Get("g_resp_pt1al2");
	g_response_pt2al1 = (TGraph2D *)fcal->Get("g_resp_pt2al1");
	g_response_pt2al2 = (TGraph2D *)fcal->Get("g_resp_pt2al2");

	// Input and Set TTree Branch Address -------------------------------------
	TString infilename = TString(inFileName);
	if (!infilename.EndsWith("gapmod.root"))
	{
		cout << "Error: Input file " << inFileName << " is NOT Gap-Mod root file!" << endl;
		return -1;
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

	Double_t gapmod_epi_pt[128] = {0.};
	Double_t gapmod_epi_al[128] = {0.};
	Int_t event_quality_flag = 0;

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

	tin->SetBranchAddress("gapmod_epi_pt", gapmod_epi_pt);
	tin->SetBranchAddress("gapmod_epi_al", gapmod_epi_al);
	tin->SetBranchAddress("event_quality_flag", &event_quality_flag);

	// Create new ROOT file and TTree  -------------------------------------
	TString outFilename;
	outFilename = inFileName;
	outFilename.ReplaceAll(".root", "_doimod.root");

	TFile *fout = new TFile(outFilename, "RECREATE");
	TTree *tout = tin->CloneTree(0);
	Double_t reconst_epi = 0.0;
	tout->Branch("reconst_epi", &reconst_epi, "reconst_epi/D");

	TH1D *hist_moddoi[2][2];
	TH1D *hist_ave[2][2];
	TH1D *hist_pt[2][2];
	TH1D *hist_al[2][2];
	TH2D *hist_avediff[2][2];
	TH2D *hist_reconst[2][2];

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			hist_moddoi[i][j] = new TH1D(Form("hist_moddoi_pt%dal%d", i + 1, j + 1), Form("Erec, DoI Mod., Pt%dAl%d;Energy [keV];Counts", i + 1, j + 1), 1100, -10, 100);
			hist_ave[i][j] = new TH1D(Form("hist_ave_pt%dal%d", i + 1, j + 1), Form("Eave, Ave., Pt%dAl%d;Energy [keV];Counts", i + 1, j + 1), 1100, -10, 100);
			hist_pt[i][j] = new TH1D(Form("hist_pt_pt%dal%d", i + 1, j + 1), Form("Pt, Pt%dAl%d;Energy [keV];Counts", i + 1, j + 1), 1100, -10, 100);
			hist_al[i][j] = new TH1D(Form("hist_al_pt%dal%d", i + 1, j + 1), Form("Al, Pt%dAl%d;Energy [keV];Counts", i + 1, j + 1), 1100, -10, 100);
			hist_avediff[i][j] = new TH2D(Form("hist_avediff_pt%dal%d", i + 1, j + 1), Form("Eave vs Ediff, Pt%dAl%d;Eave [keV];Ediff [keV]", i + 1, j + 1), 500, 0, 100, 1000, -100, 100);
			hist_reconst[i][j] = new TH2D(Form("hist_reconst_pt%dal%d", i + 1, j + 1), Form("Erec vs Ediff, Pt%dAl%d;Erec [keV];Ediff [keV]", i + 1, j + 1), 500, 0, 100, 1000, -100, 100);
		}
	}
	TH1D *histsum_reconst = new TH1D("histsum_reconst", "Erec, DoI Mod.;Energy [keV];Counts", 1100 * 2, -10, 100);
	TH1D *histsum_ave = new TH1D("histsum_ave", "Eave, Ave.;Energy [keV];Counts", 1100 * 2, -10, 100);
	TH1D *histsum_pt = new TH1D("histsum_pt", "Pt;Energy [keV];Counts", 1100 * 2, -10, 100);
	TH1D *histsum_al = new TH1D("histsum_al", "Al;Energy [keV];Counts", 1100 * 2, -10, 100);
	TH1D *histsingle_onlygood = new TH1D("histsingle_onlygood", "Ediff<0.2 keV, Single;;Energy [keV];Counts", 1100 * 2, -10, 100);
	TH2D *hist_avediff_all = new TH2D("hist_avediff_all", "Eave vs Ediff;Eave [keV];Ediff [keV]", 500, 0, 100, 1000, -100, 100);
	TH2D *hist_reconst_all = new TH2D("hist_reconst_all", "Erec vs Ediff;Erec [keV];Ediff [keV]", 500, 0, 100, 1000, -100, 100);

	TH2D *hist_each_all = new TH2D("hist_each_all", "Erec, All ch;ch;Erec [keV]", 256, -0.5, 255.5, 550 * 3, -10, 100);

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

		if (event_quality_flag != 1 && event_quality_flag != 2)
			continue;

		// DoI-loss correction
		if (pt_nhit == 1 && al_nhit == 1)
			reconst_epi = reconstruct_energypt1al1((gapmod_epi_pt[0] + gapmod_epi_al[0]) / 2.0, (gapmod_epi_pt[0] - gapmod_epi_al[0]) / 2.0);

		if (pt_nhit == 1 && al_nhit == 2)
			reconst_epi = reconstruct_energypt1al2((gapmod_epi_pt[0] + gapmod_epi_al[0]) / 2.0, (gapmod_epi_pt[0] - gapmod_epi_al[0]) / 2.0);

		if (pt_nhit == 2 && al_nhit == 1)
			reconst_epi = reconstruct_energypt2al1((gapmod_epi_pt[0] + gapmod_epi_al[0]) / 2.0, (gapmod_epi_pt[0] - gapmod_epi_al[0]) / 2.0);

		if (pt_nhit == 2 && al_nhit == 2)
			reconst_epi = reconstruct_energypt2al2((gapmod_epi_pt[0] + gapmod_epi_al[0]) / 2.0, (gapmod_epi_pt[0] - gapmod_epi_al[0]) / 2.0);

		if (al_nhit > 2 || pt_nhit > 2)
			reconst_epi = (gapmod_epi_pt[0] + gapmod_epi_al[0]) / 2.0;

		// for drawing
		if (pt_nhit < 3 && al_nhit < 3)
		{
			hist_moddoi[pt_nhit - 1][al_nhit - 1]->Fill(reconst_epi);
			hist_ave[pt_nhit - 1][al_nhit - 1]->Fill((gapmod_epi_pt[0] + gapmod_epi_al[0]) / 2.0);
			hist_pt[pt_nhit - 1][al_nhit - 1]->Fill(gapmod_epi_pt[0]);
			hist_al[pt_nhit - 1][al_nhit - 1]->Fill(gapmod_epi_al[0]);
			hist_reconst[pt_nhit - 1][al_nhit - 1]->Fill(reconst_epi, (gapmod_epi_pt[0] - gapmod_epi_al[0]) / 2.0);
			hist_avediff[pt_nhit - 1][al_nhit - 1]->Fill((merged_epi_al[0] + merged_epi_pt[0]) / 2.0, (gapmod_epi_pt[0] - gapmod_epi_al[0]) / 2.0);
		}

		if (pt_nhit == 1 && al_nhit == 1 && (epi_al[0] - epi_pt[0]) < 0.2 && (epi_al[0] - epi_pt[0]) > -0.2)
			histsingle_onlygood->Fill((epi_al[0] + epi_pt[0]) / 2.0);

		histsum_ave->Fill((merged_epi_al[0] + merged_epi_pt[0]) / 2.0);
		hist_reconst_all->Fill(reconst_epi, (gapmod_epi_pt[0] - gapmod_epi_al[0]) / 2.0);
		hist_avediff_all->Fill((gapmod_epi_pt[0] + gapmod_epi_al[0]) / 2.0, (gapmod_epi_pt[0] - gapmod_epi_al[0]) / 2.0);
		histsum_reconst->Fill(reconst_epi);
		histsum_pt->Fill(gapmod_epi_pt[0]);
		histsum_al->Fill(gapmod_epi_al[0]);

		hist_each_all->Fill(merged_pos_pt[0], gapmod_epi_pt[0]);
		hist_each_all->Fill(merged_pos_pt[0] + 128, gapmod_epi_al[0]);
		tout->Fill();
	}
	cerr << "Event Loss rate:" << nentries - tout->GetEntries() << " out of " << nentries
		 << ":" << 100. - 100. * tout->GetEntries() / nentries << "\%" << endl;

	TCanvas *c1 = new TCanvas("c1", "c1", 2400, 900);
	c1->Divide(3, 1);
	TString pdfname = TString(outFilename).ReplaceAll(".root", ".pdf");
	c1->Print(pdfname + "[", "pdf");

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			c1->cd(1);
			gPad->SetGrid();
			hist_moddoi[i][j]->GetXaxis()->SetRangeUser(0.5, 30);
			hist_moddoi[i][j]->Draw("HIST");
			hist_pt[i][j]->Draw("HISTsame");
			hist_al[i][j]->Draw("HISTsame");
			hist_reconst[i][j]->SetLineColor(kBlack);
			hist_pt[i][j]->SetLineColor(kRed);
			hist_al[i][j]->SetLineColor(kBlue);

			c1->cd(2);
			gPad->SetLogz(1);
			gPad->SetGrid();
			hist_reconst[i][j]->GetXaxis()->SetRangeUser(0.5, 30);
			hist_reconst[i][j]->GetYaxis()->SetRangeUser(-10, 10);
			hist_reconst[i][j]->Draw("colz");

			c1->cd(3);
			gPad->SetLogz(1);
			gPad->SetGrid();
			hist_avediff[i][j]->GetXaxis()->SetRangeUser(0.5, 30);
			hist_avediff[i][j]->GetYaxis()->SetRangeUser(-10, 10);
			hist_avediff[i][j]->Draw("colz");
			c1->Print(pdfname, "pdf");
		}
	}
	c1->cd(1);
	gPad->SetGrid();
	gPad->SetLogz(0);
	histsum_reconst->GetXaxis()->SetRangeUser(0.5, 30);
	histsum_reconst->SetLineColor(kRed);
	histsum_ave->SetLineColor(kBlack);
	histsum_reconst->Draw("HIST");
	histsum_ave->Draw("HISTsame");

	c1->cd(2);
	gPad->SetGrid();
	gPad->SetLogz(1);
	hist_avediff_all->GetXaxis()->SetRangeUser(0.5, 30);
	hist_avediff_all->GetYaxis()->SetRangeUser(-10, 10);
	hist_avediff_all->Draw("colz");

	c1->cd(3);
	gPad->SetLogz(1);
	gPad->SetGrid();
	hist_reconst_all->GetXaxis()->SetRangeUser(0.5, 30);
	hist_reconst_all->GetYaxis()->SetRangeUser(-10, 10);
	hist_reconst_all->Draw("colz");
	c1->Print(pdfname, "pdf");

	c1->cd(1);
	gPad->SetTickx();
	gPad->SetTicky();
	histsingle_onlygood->GetYaxis()->SetMaxDigits(3);
	histsingle_onlygood->GetXaxis()->SetRangeUser(0.5, 30);
	histsingle_onlygood->SetLineColor(kBlack);
	histsingle_onlygood->Draw("HIST");
	c1->Print(pdfname, "pdf");

	c1->Clear();
	hist_each_all->Draw("colz");
	hist_each_all->GetYaxis()->SetRangeUser(4, 32);
	c1->Print(pdfname, "pdf");

	c1->Print(pdfname + "]", "pdf");

	// Output File -------------------------------------
	tin_hktree->SetBranchStatus("*", 1);
	TTree *tout_hktree = tin_hktree->CloneTree();
	Long64_t nevent_doimod = tout->GetEntries();
	tout_hktree->Branch("nevent_doimod", &nevent_doimod, "nevent_doimod/L");
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

	histsum_reconst->Write();
	histsum_pt->Write();
	histsum_al->Write();
	hist_avediff_all->Write();
	hist_reconst_all->Write();

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			hist_moddoi[i][j]->Write();
			hist_ave[i][j]->Write();
			hist_pt[i][j]->Write();
			hist_al[i][j]->Write();
			hist_avediff[i][j]->Write();
			hist_reconst[i][j]->Write();
		}
	}
	histsum_ave->Write();
	histsingle_onlygood->Write();
	hist_each_all->Write();

	tout->Write();
	tout_hktree->Write();
	tout_threshold_tree->Write();
	fout->Close();

	cout << "  " << endl;
	cout << "Output ROOT file :" << outFilename << endl;
	cout << "Finish!" << endl;

	return 0;
}
