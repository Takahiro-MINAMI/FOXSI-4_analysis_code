//
// Position Calibration (det um -> arcsec)
// _base_cal_gapmod_doimod.root -> _base_cal_gapmod_doimod_poscal.root
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

int path5_poscal(const char *inFileName, const char *calFileName, Int_t det_number = 0);

int main(int argc, char **argv)
{

	if (argc == 4)
	{
		path5_poscal(argv[1], argv[2], Int_t(stoi(argv[3])));
	}
	else
	{
		cout << "USAGE: Please input the names of the file to analyze and Position file" << endl;
		cout << "./path5_poscal 'input filename' 'Position calibration filename' det_number (1-4:flight)" << endl;
		return -1;
	}
	return 0;
}

int path5_poscal(const char *inFileName, const char *calFileName, Int_t det_number)
{
	gStyle->SetOptStat(0);
	gStyle->SetPalette(53);

	// Load Poscal Response function (ROOT, TH1D) -------------------------------------
	TFile *fcal = new TFile(calFileName);
	const Double_t PI = TMath::Pi();
	const Double_t um2arcsec = 0.103;
	TString cal_name;
	Double_t angle_default = 0.0;
	Double_t center_cor_x = 0.0;
	Double_t center_cor_y = 0.0;
	Double_t center_cor_x_shift = -373.7907494362859;
	Double_t center_cor_y_shift = -56.4252368257264 + 10;

	TString pos_calfilenams = "T1_posxposy_CdTe1to4.root";
	TFile *pos_calfile = new TFile("T1_posxposy_CdTe1to4.root");
	TH1D *hist_pos_cal = (TH1D *)pos_calfile->Get("T1_posxposy_CdTe1to4");

	Double_t det_rot_x_um_shift = hist_pos_cal->GetBinContent((det_number - 1) * 2 + 1);
	Double_t det_rot_y_um_shift = hist_pos_cal->GetBinContent((det_number - 1) * 2 + 2);

	cout << "Road: Position File: " << calFileName << "  for det#" << det_number << endl;
	if (det_number == 1)
	{
		cal_name = "pos_cdte0";
		angle_default = -180. * PI / 180.;
	}
	else if (det_number == 2)
	{
		cal_name = "pos_cdte1";
		angle_default = -60 * PI / 180.;
	}
	else if (det_number == 3)
	{
		cal_name = "pos_cdte2";
		angle_default = -120. * PI / 180.;
	}
	else if (det_number == 4)
	{
		cal_name = "pos_cdte3";
		angle_default = 0. * PI / 180.;
	}
	else
	{
		cout << "ERROR: Specify the detector number" << endl;
		return -1;
	}

	TH2D *h_poscal = (TH2D *)fcal->Get(cal_name);
	Double_t center_x = h_poscal->GetBinContent(0);
	Double_t center_y = h_poscal->GetBinContent(1);
	// Double_t angle_x_from_x_ax = h_poscal->GetBinContent(2);
	// Double_t angle_y_from_x_ax = h_poscal->GetBinContent(3);

	TString substrip_calfilenams[3] = {"./poscal_resp/substrip_pos_60um.root", "./poscal_resp/substrip_pos_80um.root", "./poscal_resp/substrip_pos_100um.root"};

	TFile *fsubstrip_calfiles[3];
	TF1 *pt_single_norm[3];
	TH2D *pt_double_norm[3];
	TF1 *al_single_norm[3];
	TH2D *al_double_norm[3];
	for (int i = 0; i < 3; i++)
	{
		fsubstrip_calfiles[i] = new TFile(substrip_calfilenams[i]);
		pt_single_norm[i] = (TF1 *)fsubstrip_calfiles[i]->Get("pt_single_norm");
		pt_double_norm[i] = (TH2D *)fsubstrip_calfiles[i]->Get("pt_double_norm");
		al_single_norm[i] = (TF1 *)fsubstrip_calfiles[i]->Get("al_single_norm");
		al_double_norm[i] = (TH2D *)fsubstrip_calfiles[i]->Get("al_double_norm");
	}

	// Input and Set TTree Branch Address -------------------------------------
	TString infilename = TString(inFileName);
	if (!infilename.EndsWith("doimod.root"))
	{
		cout << "Error: Input file " << inFileName << " is NOT DoI-Mod root file!" << endl;
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
	Double_t approx_passed_time_from_launch_sec = 0.;
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
	tin->SetBranchAddress("approx_passed_time_from_launch_sec", &approx_passed_time_from_launch_sec);

	tin->SetBranchAddress("event_quality_flag", &event_quality_flag);

	// Create new ROOT file and TTree  -------------------------------------
	TString outFilename;
	outFilename = inFileName;
	outFilename.ReplaceAll(".root", "_poscal.root");

	TFile *fout = new TFile(outFilename, "RECREATE");
	TTree *tout = tin->CloneTree(0);

	Double_t det_x_um = 0.;
	Double_t det_y_um = 0.;
	Double_t det_rot_x_um = 0.;
	Double_t det_rot_y_um = 0.;
	Double_t det_rot_x_arcsec = 0.;
	Double_t det_rot_y_arcsec = 0.;
	Double_t helio_x_arcsec = 0.;
	Double_t helio_y_arcsec = 0.;
	Double_t helio_x_arcsec_shift = 0.;
	Double_t helio_y_arcsec_shift = 0.;
	Double_t livetime_ratio = 1.;
	Double_t pre_pseudo_counter = 0.;
	Double_t pre_survive_pseudo_counter = 0.;

	tout->Branch("det_x_um", &det_x_um, "det_x_um/D");
	tout->Branch("det_y_um", &det_y_um, "det_y_um/D");
	tout->Branch("det_rot_x_um", &det_rot_x_um, "det_rot_x_um/D");
	tout->Branch("det_rot_y_um", &det_rot_y_um, "det_rot_y_um/D");
	tout->Branch("det_rot_x_arcsec", &det_rot_x_arcsec, "det_rot_x_arcsec/D");
	tout->Branch("det_rot_y_arcsec", &det_rot_y_arcsec, "det_rot_y_arcsec/D");
	tout->Branch("helio_x_arcsec", &helio_x_arcsec, "helio_x_arcsec/D");
	tout->Branch("helio_y_arcsec", &helio_y_arcsec, "helio_y_arcsec/D");
	tout->Branch("helio_x_arcsec_shift", &helio_x_arcsec_shift, "helio_x_arcsec_shift/D");
	tout->Branch("helio_y_arcsec_shift", &helio_y_arcsec_shift, "helio_y_arcsec_shift/D");
	tout->Branch("livetime_ratio", &livetime_ratio, "livetime_ratio/D");

	// For Drawing  -------------------------------------
	Double_t pos_strip_center[128] = {0.};
	Double_t gap_width[128] = {0.};
	pos_strip_center[64] = 30;
	pos_strip_center[63] = -30;
	gap_width[64] = 60;
	gap_width[63] = 60;

	Double_t strip_gap = 0;
	for (int i = 62; i >= 0; i--)
	{
		if (i >= 49)
			strip_gap = 60;
		else if (i < 49 && i >= 29)
			strip_gap = 80;
		else if (i < 29 && i >= 0)
			strip_gap = 100;
		pos_strip_center[i] = pos_strip_center[i + 1] - strip_gap;
		pos_strip_center[64 - i + 63] = pos_strip_center[63 - i + 63] + strip_gap;
		gap_width[i] = strip_gap;
		gap_width[64 - i + 63] = strip_gap;
	}

	const Int_t bins = Int_t((pos_strip_center[127] - pos_strip_center[0] + 100) / 64 * 2);
	const Double_t start = pos_strip_center[0];
	const Double_t end = pos_strip_center[127];
	const Double_t start_arcsec = um2arcsec * start;
	const Double_t end_arcsec = um2arcsec * end;
	const Double_t helio_para_x_center = 0;
	const Double_t helio_para_y_center = 0;
	const Double_t helio_rot_x_center = 0.;
	const Double_t helio_rot_y_center = 0.0;
	const Double_t helio_x_center = um2arcsec * start + helio_para_x_center;
	const Double_t helio_y_center = um2arcsec * start + helio_para_y_center;

	TH2D *image_det_um = new TH2D("image_det_um", "image_det_um", bins, start, end, bins, start, end);
	TH2D *image_det_para_um = new TH2D("image_det_para_um", "image_det_para_um", bins, start, end, bins, start, end);
	TH2D *image_det_rot_um = new TH2D("image_det_rot_um", "image_det_rot_um", bins, start, end, bins, start, end);
	TH2D *image_det_rot_arcsec = new TH2D("image_det_rot_arcsec", "image_det_rot_arcsec", bins, start_arcsec, end_arcsec, bins, start_arcsec, end_arcsec);
	TH2D *image_helio_arcsec = new TH2D("helio_arcsec", "helio_arcsec", bins, helio_x_center, end_arcsec + helio_para_x_center, bins, helio_y_center, end_arcsec + helio_para_y_center);
	TH2D *image_helio_arcsec_shift = new TH2D("helio_arcsec_shift", "helio_arcsec_shift", bins, helio_x_center, end_arcsec + helio_para_x_center, bins, helio_y_center, end_arcsec + helio_para_y_center);

	TH1D *al_forfind = (TH1D *)al_double_norm[0]->ProjectionY("ab", 1, 1);
	TH1D *pt_forfind = (TH1D *)pt_double_norm[0]->ProjectionY("abc", 1, 1);
	TH1D *hist_tmp = new TH1D("", "", 80, -40, 40);
	TH1D *p_prob;

	Long64_t nentries = tin->GetEntriesFast();
	Double_t energy_ratio_x = 0.;
	Double_t energy_ratio_y = 0.;
	Double_t det_pitch_um[3] = {60., 80., 100.};

	// Main Process  -------------------------------------
	TRandom3 rndm;
	Double_t dx = 0.0;
	Double_t dy = 0.0;
	Double_t tmp_pos = 0.0;
	Int_t tmp_ene_bin = 0;

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
		if (al_nhit > 3 || pt_nhit > 3)
			continue;

		// Position correction
		if (pt_nhit == 1 && al_nhit == 1)
		{
			energy_ratio_x = 1.;
			tmp_pos = pt_single_norm[merged_pt_region[0]]->GetRandom();
			det_x_um = -(pos_strip_center[int(pos_pt[0])] + tmp_pos);

			energy_ratio_y = 1.;
			tmp_pos = al_single_norm[merged_al_region[0]]->GetRandom();
			det_y_um = -(pos_strip_center[int(pos_al[0])] + tmp_pos);
		}

		if (pt_nhit == 1 && al_nhit == 2)
		{
			energy_ratio_x = 1.;
			tmp_pos = pt_single_norm[merged_pt_region[0]]->GetRandom();
			det_x_um = -(pos_strip_center[int(pos_pt[0])] + tmp_pos);

			energy_ratio_y = epi_al[1] / (epi_al[0] + epi_al[1]);
			tmp_ene_bin = al_forfind->FindBin(energy_ratio_y);
			p_prob = (TH1D *)al_double_norm[merged_al_region[0]]->ProjectionX("abcz", tmp_ene_bin, tmp_ene_bin);
			hist_tmp->FillRandom(p_prob, 1);
			tmp_pos = hist_tmp->GetMean() + det_pitch_um[merged_al_region[0]] / 2.0;
			det_y_um = -(pos_strip_center[int(pos_al[0])] + tmp_pos);
			hist_tmp->Reset();
		}

		if (pt_nhit == 2 && al_nhit == 1)
		{
			energy_ratio_x = epi_pt[1] / (epi_pt[0] + epi_pt[1]);
			tmp_ene_bin = pt_forfind->FindBin(energy_ratio_x);
			p_prob = (TH1D *)pt_double_norm[merged_pt_region[0]]->ProjectionX("abcz", tmp_ene_bin, tmp_ene_bin);
			hist_tmp->FillRandom(p_prob, 1);
			tmp_pos = (hist_tmp->GetMean() + det_pitch_um[merged_pt_region[0]] / 2.0);
			det_x_um = -(pos_strip_center[int(pos_pt[0])] + tmp_pos);
			hist_tmp->Reset();

			energy_ratio_y = 1.;
			tmp_pos = al_single_norm[merged_al_region[0]]->GetRandom();
			det_y_um = -(pos_strip_center[int(pos_al[0])] + tmp_pos);
		}

		if (pt_nhit == 2 && al_nhit == 2)
		{
			energy_ratio_x = epi_pt[1] / (epi_pt[0] + epi_pt[1]);
			tmp_ene_bin = pt_forfind->FindBin(energy_ratio_x);
			p_prob = (TH1D *)pt_double_norm[merged_pt_region[0]]->ProjectionX("abcz", tmp_ene_bin, tmp_ene_bin);
			hist_tmp->FillRandom(p_prob, 1);
			tmp_pos = (hist_tmp->GetMean() + det_pitch_um[merged_pt_region[0]] / 2.0);
			det_x_um = -(pos_strip_center[int(pos_pt[0])] + tmp_pos);
			hist_tmp->Reset();

			energy_ratio_y = epi_al[1] / (epi_al[0] + epi_al[1]);
			tmp_ene_bin = al_forfind->FindBin(energy_ratio_y);
			p_prob = (TH1D *)al_double_norm[merged_al_region[0]]->ProjectionX("abcz", tmp_ene_bin, tmp_ene_bin);
			hist_tmp->FillRandom(p_prob, 1);
			tmp_pos = hist_tmp->GetMean() + det_pitch_um[merged_al_region[0]] / 2.0;
			det_y_um = -(pos_strip_center[int(pos_al[0])] + tmp_pos);
			hist_tmp->Reset();
		}

		center_x = 0;
		center_y = 0;
		det_rot_x_um = TMath::Cos(angle_default) * (det_x_um + center_x) + TMath::Sin(angle_default) * (det_y_um - center_y);
		det_rot_y_um = -TMath::Sin(angle_default) * (det_x_um + center_x) + TMath::Cos(angle_default) * (det_y_um - center_y);
		det_rot_x_arcsec = (det_rot_x_um - det_rot_x_um_shift) * um2arcsec;
		det_rot_y_arcsec = (det_rot_y_um - det_rot_y_um_shift) * um2arcsec;
		helio_x_arcsec = TMath::Cos(-PI / 3.0) * (det_rot_x_arcsec - helio_rot_x_center) + TMath::Sin(-PI / 3.0) * (det_rot_y_arcsec - helio_rot_y_center) + center_cor_x_shift;
		helio_y_arcsec = -TMath::Sin(-PI / 3.0) * (det_rot_x_arcsec - helio_rot_x_center) + TMath::Cos(-PI / 3.0) * (det_rot_y_arcsec - helio_rot_y_center) + center_cor_y_shift;

		// Fine Position Calibration

		if (approx_passed_time_from_launch_sec >= 100.0 && approx_passed_time_from_launch_sec <= 320.0)
		{
			dx = (-0.828961 * approx_passed_time_from_launch_sec + 90.8787) / (0.97 * 11.0) * 11.0 * um2arcsec;
			dy = (-0.0110396 * dx * dx + 0.381066 * dx + 0.0259491);
			helio_x_arcsec_shift = helio_x_arcsec - dx * 11.0 * um2arcsec;
			helio_y_arcsec_shift = helio_y_arcsec - dy * 11.0 * um2arcsec;
		}
		else if (approx_passed_time_from_launch_sec > 320.0 && approx_passed_time_from_launch_sec <= 375.0)
		{
			helio_x_arcsec_shift = helio_x_arcsec + (-324.21 + 311.21);
			helio_y_arcsec_shift = helio_y_arcsec + (420.15 + 154.85);
			cout << approx_passed_time_from_launch_sec << " " << helio_x_arcsec << " " << helio_y_arcsec << " " << helio_x_arcsec_shift << " " << helio_y_arcsec_shift << endl;
		}
		else if (approx_passed_time_from_launch_sec > 395 && approx_passed_time_from_launch_sec <= 430)
		{
			helio_x_arcsec_shift = helio_x_arcsec + (-324.21 + 311.21);
			helio_y_arcsec_shift = helio_y_arcsec; //(-137.85 + 154.85);
												   // helio_x_arcsec_shift = helio_x_arcsec + (-324.21 + 311.21) + 29;
												   // helio_y_arcsec_shift = helio_y_arcsec + 31; //(-137.85 + 154.85);
		}
		// livetime_ratio = 1.;
		// pre_pseudo_counter = 0.;
		// pre_survive_pseudo_counter = 0.;

		// if (jentry == 0)
		// {
		// 	pre_pseudo_counter = pseudo_counter;
		// 	pre_survive_pseudo_counter = survive_pseudo_counter;
		// 	livetime_ratio = 0.;
		// }
		// else
		// {
		// 	if ((pseudo_counter - pre_pseudo_counter) < 1)
		// 		livetime_ratio = 0.;
		// 	else
		// 	{
		// 		livetime_ratio = (survive_pseudo_counter - pre_survive_pseudo_counter) / (pseudo_counter - pre_pseudo_counter);
		// 		pre_pseudo_counter = pseudo_counter;
		// 		pre_survive_pseudo_counter = survive_pseudo_counter;
		// 	}
		// }

		// for drawing
		image_det_um->Fill(det_x_um, det_y_um);
		image_det_para_um->Fill(det_x_um + center_x, det_y_um - center_y);
		image_det_rot_um->Fill(det_rot_x_um, det_rot_y_um);
		image_det_rot_arcsec->Fill(det_rot_x_arcsec, det_rot_y_arcsec);
		image_helio_arcsec->Fill(helio_x_arcsec, helio_y_arcsec);
		image_helio_arcsec_shift->Fill(helio_x_arcsec_shift, helio_y_arcsec_shift);
		tout->Fill();
	}

	cerr << "Event Loss rate:" << nentries - tout->GetEntries() << " out of " << nentries
		 << ":" << 100. - 100. * tout->GetEntries() / nentries << "\%" << endl;

	TCanvas *c1 = new TCanvas("c1", "c1", 1800, 900);
	c1->Divide(2, 1);
	TString pdfname = TString(outFilename).ReplaceAll(".root", ".pdf");
	c1->Print(pdfname + "[", "pdf");

	c1->cd(1);
	gPad->SetGrid();
	image_det_um->Draw("colz");
	c1->cd(2);
	gPad->SetGrid();
	image_det_para_um->Draw("colz");
	c1->Print(pdfname, "pdf");

	c1->cd(1);
	gPad->SetGrid();
	image_det_rot_um->Draw("colz");
	c1->cd(2);
	gPad->SetGrid();
	image_det_rot_arcsec->Draw("colz");
	c1->Print(pdfname, "pdf");

	c1->cd(1);
	gPad->SetGrid();
	image_helio_arcsec->Draw("colz");
	c1->cd(2);
	gPad->SetGrid();
	image_helio_arcsec_shift->Draw("colz");
	c1->Print(pdfname, "pdf");

	c1->Print(pdfname + "]", "pdf");

	// Output File -------------------------------------
	tin_hktree->SetBranchStatus("*", 1);
	TTree *tout_hktree = tin_hktree->CloneTree();
	Long64_t nevent_poscal = tout->GetEntries();
	tout_hktree->Branch("nevent_poscal", &nevent_poscal, "nevent_poscal/L");
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

	image_det_um->Write();
	image_det_para_um->Write();
	image_det_rot_um->Write();
	image_det_rot_arcsec->Write();
	image_helio_arcsec->Write();
	image_helio_arcsec_shift->Write();

	tout->Write();
	tout_hktree->Write();
	tout_threshold_tree->Write();
	fout->Close();

	cout << "  " << endl;
	cout << "Output ROOT file :" << outFilename << endl;
	cout << "Finish!" << endl;

	return 0;
}
