//
// calibration (adc->Energy)
// _base.root -> _base_cal.root
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
#include <Rtypes.h>
#include <sstream>
#include <TPaletteAxis.h>
#include <filesystem>

using namespace std;

int pdfout = 1;

int cal(const char *inFileName, const char *calfilename, const char *threshfilename = "", bool save_all = false);

int main(int argc, char **argv)
{
	if (argc == 4)
		cal(argv[1], argv[2], argv[3]);
	else if (argc == 5)
	{
		cout << !strncmp(argv[4], "-all", 4) << endl;
		cal(argv[1], argv[2], argv[3], !strncmp(argv[4], "-all", 4));
	}
	else
	{
		cout << "USAGE: Please input the names of the file to analyze and cal file" << endl;
		cout << "./path2_cal_and_merge_eachside 'input filename' 'cal filename' 'threshold filename' -all" << endl;
		return -1;
	}

	return 0;
}

int cal(const char *inFileName, const char *calfilename, const char *threshfilename, bool save_all)
{
	gStyle->SetOptStat(0);
	gStyle->SetPalette(53);

	cout << "-- Path2: Each ch Calibration ------------" << endl;
	if (save_all)
		cout << "All data size mode ...." << endl;
	else
		cout << "Reduce data size mode ...." << endl;

	// Load Calibration function (ROOT, Spline3) -------------------------------------
	cout << "Loading Calibration file:  " << calfilename << endl;
	TFile *fCal = new TFile(calfilename);
	if (!fCal->IsOpen())
	{
		cerr << "ERROR: NO such calibration file as " << calfilename << endl;
		return -1;
	}
	TSpline3 *s[256];

	for (Int_t i = 0; i < 256; ++i)
		s[i] = (TSpline3 *)fCal->Get(Form("spline%d", i));
	fCal->Close();

	// Load Threshold information -------------------------------------
	cout << "Loading Threshold file:  " << threshfilename << endl;
	if (!std::filesystem::is_regular_file(threshfilename))
	{
		cout << "ERROR: NO such a threshold file as " << threshfilename << " exists!" << endl;
		return -1;
	}
	ifstream ifs(threshfilename);
	Double_t adcThreshIn[256], adcThreshIn_peak[256], thresh_each, thresh_peak;
	string str;
	Int_t count = 0;
	TTree *threshold_tree = new TTree("threshold", "threshold");
	threshold_tree->Branch("thresh", &thresh_each, "thresh_each/D");
	threshold_tree->Branch("thresh_peak", &thresh_peak, "thresh_peak/D");

	while (getline(ifs, str))
	{
		std::istringstream ss(str);
		if (ss >> thresh_each >> thresh_peak)
		{
			adcThreshIn[count] = thresh_each;
			adcThreshIn_peak[count] = thresh_peak;
			cout << "ch " << count << "  " << adcThreshIn[count] << "  " << adcThreshIn_peak[count] << endl;
		}
		threshold_tree->Fill();
		if (count > 255)
			break;
		++count;
	}
	if (count < 255)
	{
		cout << "ERROR: Threshold file must contain 256 or more thresholds!" << endl;
		return -1;
	}

	// Input and Set TTree Branch Address -------------------------------------
	TString infilename = TString(inFileName);
	if (!infilename.EndsWith("base.root"))
	{
		cout << "Error: Input file " << inFileName << " is NOT Base root file!" << endl;
		return -1;
	}

	UInt_t ti;
	UInt_t livetime;
	UInt_t integral_livetime;
	UInt_t unixtime;
	Int_t pt_adc_num;
	Int_t adc_pt[128];
	Int_t remapch_pt[128];
	Int_t al_adc_num;
	Int_t adc_al[128];
	Int_t remapch_al[128];

	UShort_t flag_pseudo;
	UInt_t pseudo_counter;
	UInt_t survive_pseudo_counter = 0;

	TFile *fin = new TFile(inFileName);
	TTree *tin = (TTree *)fin->Get("eventtree");
	TTree *tin_hktree = (TTree *)fin->Get("hktree");

	tin->SetBranchAddress("ti", &ti);
	tin->SetBranchAddress("livetime", &livetime);
	tin->SetBranchAddress("integral_livetime", &integral_livetime);
	tin->SetBranchAddress("unixtime", &unixtime);
	tin->SetBranchAddress("pt_adc_num", &pt_adc_num);
	tin->SetBranchAddress("adc_pt", adc_pt);
	tin->SetBranchAddress("remapch_pt", remapch_pt);
	tin->SetBranchAddress("al_adc_num", &al_adc_num);
	tin->SetBranchAddress("adc_al", adc_al);
	tin->SetBranchAddress("remapch_al", remapch_al);
	tin->SetBranchAddress("flag_pseudo", &flag_pseudo);
	tin->SetBranchAddress("pseudo_counter", &pseudo_counter);
	tin->SetBranchAddress("survive_pseudo_counter", &survive_pseudo_counter);

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
	TString outFilename = inFileName;
	outFilename.ReplaceAll(".root", "_cal.root");

	TFile *fout = new TFile(outFilename, "RECREATE");
	TTree *tout = tin->CloneTree(0);
	Int_t pt_nhit = 0;
	Double_t epi_pt[128] = {0.};
	Double_t pos_pt[128] = {0.0};
	Int_t merged_pt_nhit = 0;
	Double_t merged_epi_pt[128] = {0.};
	Double_t merged_pos_pt[128] = {0.};
	Int_t merged_pt_nhit_list[128] = {0};

	Int_t al_nhit = 0;
	Double_t epi_al[128] = {0.};
	Double_t pos_al[128] = {0.0};
	Int_t merged_al_nhit = 0;
	Double_t merged_epi_al[128] = {0.};
	Double_t merged_pos_al[128] = {0.};
	Int_t merged_al_nhit_list[128] = {0};
	Int_t merged_pt_region[128] = {0};
	Int_t merged_al_region[128] = {0};

	Double_t pt_dx[128] = {0.};
	Double_t al_dx[128] = {0.};

	tout->Branch("pt_nhit", &pt_nhit, "pt_nhit/I");
	tout->Branch("epi_pt", epi_pt, "epi_pt[pt_nhit]/D");
	tout->Branch("pos_pt", pos_pt, "pos_pt[pt_nhit]/D");
	tout->Branch("merged_pt_nhit", &merged_pt_nhit, "merged_pt_nhit/I");
	tout->Branch("merged_epi_pt", merged_epi_pt, "merged_epi_pt[merged_pt_nhit]/D");
	tout->Branch("merged_pos_pt", merged_pos_pt, "merged_pos_pt[merged_pt_nhit]/D");
	tout->Branch("merged_pt_nhit_list", merged_pt_nhit_list, "merged_pt_nhit_list[merged_pt_nhit]/I");
	tout->Branch("merged_pt_region", merged_pt_region, "merged_pt_region[merged_pt_nhit]/I");
	tout->Branch("pt_dx", pt_dx, "pt_dx[merged_pt_nhit]/D");

	tout->Branch("al_nhit", &al_nhit, "al_nhit/I");
	tout->Branch("epi_al", epi_al, "epi_al[al_nhit]/D");
	tout->Branch("pos_al", pos_al, "pos_al[al_nhit]/D");
	tout->Branch("merged_al_nhit", &merged_al_nhit, "merged_al_nhit/I");
	tout->Branch("merged_epi_al", merged_epi_al, "merged_epi_al[merged_al_nhit]/D");
	tout->Branch("merged_pos_al", merged_pos_al, "merged_pos_al[merged_al_nhit]/D");
	tout->Branch("merged_al_nhit_list", merged_al_nhit_list, "merged_al_nhit_list[merged_al_nhit]/I");
	tout->Branch("merged_al_region", merged_al_region, "merged_al_region[merged_al_nhit]/I");
	tout->Branch("al_dx", al_dx, "al_dx[merged_al_nhit]/D");

	Double_t adc, enEval;
	Int_t ch;

	TH1D *hist_cal[256];
	TH1D *hist_merged[256];
	for (Int_t i = 0; i < 256; i++)
	{
		hist_cal[i] = new TH1D(Form("hist_cal%d", i), Form("hist_cal%d", i), 1100, -10, 100);
		hist_merged[i] = new TH1D(Form("hist_merged%d", i), Form("hist_merged%d", i), 1010, -10, 100);
	}
	TH2D *histall_cal = new TH2D("histall_cal", "Calibration spectrum, Allch", 256, -0.5, 255.5, 1000, 0, 100);
	TH1D *histpt_cal = new TH1D("histpt_cal", "Calibration spectrum, Pt side sum;Energy [keV];Counts", 1100, -10, 100);
	TH1D *histal_cal = new TH1D("histal_cal", "Calibration spectrum, Al side sum;Energy [keV];Counts", 1100, -10, 100);
	TH1D *histpt_cal_single = new TH1D("histpt_cal_single", "Calibration spectrum, Pt side sum, Single strip;Energy [keV];Counts", 1100, -10, 100);
	TH1D *histal_cal_single = new TH1D("histal_cal_single", "Calibration spectrum, Al side sum, Single strip;Energy [keV];Counts", 1100, -10, 100);
	TH1D *histpt_cal_double = new TH1D("histpt_cal_double", "Calibration spectrum, Pt side sum, Double strip (Sum);Energy [keV];Counts", 1100, -10, 100);
	TH1D *histal_cal_double = new TH1D("histal_cal_double", "Calibration spectrum, Al side sum, Double strip (Sum);Energy [keV];Counts", 1100, -10, 100);

	TH2D *histall_merged = new TH2D("histall_merged", "Calibration spectrum, Allch, All event (Sum);ch;Energy [keV]", 256, -0.5, 255.5, 1000, 0, 100);
	TH2D *histall_single = new TH2D("histall_single", "Calibration spectrum, Allch, Single strip;ch;Energy [keV]", 256, -0.5, 255.5, 1000, 0, 100);
	TH2D *histall_double = new TH2D("histall_double", "Calibration spectrum, Allch, Double strip (Sum);ch;Energy [keV]", 256, -0.5, 255.5, 1000, 0, 100);

	TH2D *hist_image_ch = new TH2D("hist_image_ch", "Image (ch scale);Pt ch;Al ch", 128, -0.5, 127.5, 128, -0.5, 127.5);
	TH2D *hist_image_ch_single = new TH2D("hist_image_ch_single", "Image (ch scale, Single strip);Pt ch;Al ch", 128, -0.5, 127.5, 128, -0.5, 127.5);

	TH1D *hist_nhit_pt = new TH1D("hist_nhit_pt", "# of hit, Pt side;hit ch", 26, -0.5, 25.5);
	TH1D *hist_nhit_al = new TH1D("hist_nhit_al", "# of hit, Al side;hit ch", 26, -0.5, 25.5);
	TH2D *hist_merged_nhit = new TH2D("hist_merged_nhit", "# of Multiplicity;Pt ch;Al ch", 26, -0.5, 25.5, 26, -0.5, 25.5);

	TH1D *histpt_cal_single_each[3];
	TH1D *histpt_cal_double_each[3];
	TH1D *histal_cal_single_each[3];
	TH1D *histal_cal_double_each[3];
	TH2D *hist_vsene_pt_adj[3];
	TH2D *hist_vsene_al_adj[3];
	TH2D *hist_diffene_pt_adj[3];
	TH2D *hist_diffene_al_adj[3];

	for (Int_t i = 0; i < 3; i++)
	{
		histpt_cal_single_each[i] = new TH1D(Form("histpt_cal_single_each_region%d", i), Form("histpt_cal_single_each_region%d", i), 500, 0, 100);
		histpt_cal_double_each[i] = new TH1D(Form("histpt_cal_double_each_region%d", i), Form("histpt_cal_double_each_region%d", i), 500, 0, 100);
		histal_cal_single_each[i] = new TH1D(Form("histal_cal_single_each_region%d", i), Form("histal_cal_single_each_region%d", i), 500, 0, 100);
		histal_cal_double_each[i] = new TH1D(Form("histal_cal_double_each_region%d", i), Form("histal_cal_double_each_region%d", i), 500, 0, 100);
		hist_vsene_pt_adj[i] = new TH2D(Form("hist_vsene_pt_adj_region%d", i), Form("hist_vsene_pt_adj_region%d", i), 1000, 0, 100, 1000, 0, 100);
		hist_vsene_al_adj[i] = new TH2D(Form("hist_vsene_al_adj_region%d", i), Form("hist_vsene_al_adj_region%d", i), 1000, 0, 100, 1000, 0, 100);
		hist_diffene_pt_adj[i] = new TH2D(Form("hist_diffene_pt_adj_region%d", i), Form("hist_diffene_pt_adj_region%d", i), 1000, -50, 50, 1000, 0, 100);
		hist_diffene_al_adj[i] = new TH2D(Form("hist_diffene_al_adj_region%d", i), Form("hist_diffene_al_adj_region%d", i), 1000, -50, 50, 1000, 0, 100);
	}

	Double_t xbins[129] = {0.0};
	Double_t xbins_center[128] = {0.0};
	Double_t tmp_x = -5.380;
	xbins[0] = tmp_x;
	for (int ch = 0; ch < 128; ++ch)
	{
		if (ch >= 0 && ch <= 3)
			tmp_x = -4930. - 100. * (4. - ch);
		if (ch >= 4 && ch <= 27)
			tmp_x = -2630. - 100. * (27. - ch);
		if (ch >= 28 && ch <= 47)
			tmp_x = -1010. - 80. * (47. - ch);
		if (ch >= 48 && ch <= 63)
			tmp_x = -30. - 60. * (63. - ch);
		if (ch >= 64 && ch <= 79)
			tmp_x = 30. + 60. * (ch - 64);
		if (ch >= 80 && ch <= 99)
			tmp_x = 1010. + 80. * (ch - 80);
		if (ch >= 100 && ch <= 123)
			tmp_x = 2630. + 100. * (ch - 100);
		if (ch >= 124 && ch <= 128)
			tmp_x = 4930. + 100. * (ch - 123);
		xbins_center[ch] = tmp_x / 1000.; // convert to mm
	}
	for (int ch = 1; ch < 128; ++ch)
		xbins[ch] = (xbins_center[ch - 1] + xbins_center[ch]) / 2.0;

	xbins[0] = -5350. - 50.0;
	xbins[128] = 5350. + 50.0;

	TH2D *hist_image_dxdy = new TH2D("hist_image_dxdy", "Image (mm scale)", 128, xbins, 128, xbins);
	TH2D *hist_image_dxdy_scale = new TH2D("hist_image_dxdy_scale", "Image (mm scale, Area scaled)", 128, xbins, 128, xbins);

	TH1D *hist_slice_dx_scale = new TH1D("hist_slice_dx_scale", "Image (mm scale, Area scaled), Xslice", 128, xbins);
	TH1D *hist_slice_dy_scale = new TH1D("hist_slice_dy_scale", "Image (mm scale, Area scaled), Yslice", 128, xbins);

	TH1D *hist_slice_dx = new TH1D("hist_slice_dx", "Image (mm scale), Xslice", 128, -0.5, 127.5);
	TH1D *hist_slice_dy = new TH1D("hist_slice_dy", "Image (mm scale), Yslice", 128, -0.5, 127.5);

	TH1D *dist_hit = new TH1D("dist_hit", "dist_hit", 128, -0.5, 127.5);

	TRandom3 rndm;
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

		Int_t hit_al = 0;
		Int_t hit_pt = 0;
		merged_al_nhit = 0;
		merged_pt_nhit = 0;
		Int_t lowbin, upbin, maxbin;

		if (flag_pseudo == 1)
			continue;

		// Pt Side Process  -------------------------------------
		for (Int_t i = 0; i < pt_adc_num; ++i)
		{
			ch = remapch_pt[i];
			adc = adc_pt[i] + rndm.Uniform(-0.5, 0.5);
			enEval = s[ch]->Eval(adc);
			if (adcThreshIn[ch] > enEval)
				continue;
			pos_pt[hit_pt] = 1.0 * remapch_pt[i];
			epi_pt[hit_pt] = enEval;
			++hit_pt;
			hist_cal[ch]->Fill(enEval);
			histall_cal->Fill(ch, enEval);
			histpt_cal->Fill(enEval);
		}

		// Merge  -------------------------------------
		dist_hit->Reset();
		for (Int_t j = 0; j < hit_pt; j++)
			dist_hit->Fill(pos_pt[j], epi_pt[j]);

		hit_pt = 0;
		while (dist_hit->GetMaximum() > adcThreshIn_peak[(dist_hit->GetMaximumBin() - 1)])
		{
			maxbin = dist_hit->GetMaximumBin();
			for (lowbin = maxbin - 1; lowbin >= 1; lowbin--)
			{
				if (dist_hit->GetBinContent(lowbin) > 0.0)
				{
				}
				else
				{
					lowbin++;
					break;
				}
			}
			for (upbin = maxbin + 1; upbin <= 128; upbin++)
			{
				if (dist_hit->GetBinContent(upbin) > 0.0)
				{
				}
				else
				{
					upbin--;
					break;
				}
			}
			dist_hit->GetXaxis()->SetRange(lowbin, upbin);
			merged_epi_pt[merged_pt_nhit] = dist_hit->Integral(lowbin, upbin); // energy sum
			merged_pos_pt[merged_pt_nhit] = dist_hit->GetMean(1);							 // energy weighted mean
			pt_dx[merged_pt_nhit] = xbins_center[(int)merged_pos_pt[merged_pt_nhit]];

			// Numbering by the strip-pitch region (-1: Guard ring, 0: 60um, 1: 80um, 2: 100um)
			Double_t tmp_pos_pt = merged_pos_pt[merged_pt_nhit];
			if (tmp_pos_pt < 4)
				merged_pt_region[merged_pt_nhit] = -1;
			if (tmp_pos_pt >= 4 && tmp_pos_pt < 28)
				merged_pt_region[merged_pt_nhit] = 2;
			if (tmp_pos_pt >= 28 && tmp_pos_pt < 48)
				merged_pt_region[merged_pt_nhit] = 1;
			if (tmp_pos_pt >= 48 && tmp_pos_pt < 79)
				merged_pt_region[merged_pt_nhit] = 0;
			if (tmp_pos_pt >= 79 && tmp_pos_pt < 99)
				merged_pt_region[merged_pt_nhit] = 1;
			if (tmp_pos_pt >= 99 && tmp_pos_pt < 124)
				merged_pt_region[merged_pt_nhit] = 2;
			if (tmp_pos_pt >= 123)
				merged_pt_region[merged_pt_nhit] = -1;

			hist_merged[(int)merged_pos_pt[merged_pt_nhit]]->Fill(merged_epi_pt[merged_pt_nhit]);
			histall_merged->Fill(merged_pos_pt[merged_pt_nhit], merged_epi_pt[merged_pt_nhit]);
			merged_pt_nhit_list[merged_pt_nhit] = (upbin - lowbin) + 1;

			for (Int_t j = lowbin; j <= upbin; j++)
			{
				pos_pt[hit_pt] = dist_hit->GetBinCenter(j);
				epi_pt[hit_pt] = dist_hit->GetBinContent(j);
				dist_hit->SetBinContent(j, 0.0);
				hit_pt = hit_pt + 1;
			}
			dist_hit->GetXaxis()->SetRange(1, 128);
			merged_pt_nhit++;
		}
		pt_nhit = hit_pt;

		// Al Side Process  -------------------------------------
		for (Int_t i = 0; i < al_adc_num; ++i)
		{
			ch = remapch_al[i] + 128;
			adc = adc_al[i] + rndm.Uniform(-0.5, 0.5);
			enEval = s[ch]->Eval(adc);
			if (adcThreshIn[ch] > enEval)
				continue;
			pos_al[hit_al] = 1.0 * remapch_al[i];
			epi_al[hit_al] = enEval;
			++hit_al;
			hist_cal[ch]->Fill(enEval);
			histall_cal->Fill(ch, enEval);
			histal_cal->Fill(enEval);
		}

		// Merge  -------------------------------------
		dist_hit->Reset();
		for (Int_t j = 0; j < hit_al; j++)
			dist_hit->Fill(pos_al[j], epi_al[j]);
		hit_al = 0;

		while (dist_hit->GetMaximum() > adcThreshIn_peak[(dist_hit->GetMaximumBin() - 1) + 128])
		{
			maxbin = dist_hit->GetMaximumBin();
			for (lowbin = maxbin - 1; lowbin >= 1; lowbin--)
			{
				if (dist_hit->GetBinContent(lowbin) > 0.0)
				{
				}
				else
				{
					lowbin++;
					break;
				}
			}
			for (upbin = maxbin + 1; upbin <= 128; upbin++)
			{
				if (dist_hit->GetBinContent(upbin) > 0.0)
				{
				}
				else
				{
					upbin--;
					break;
				}
			}
			dist_hit->GetXaxis()->SetRange(lowbin, upbin);
			merged_epi_al[merged_al_nhit] = dist_hit->Integral(lowbin, upbin);
			merged_pos_al[merged_al_nhit] = dist_hit->GetMean(1);
			al_dx[merged_al_nhit] = xbins_center[(int)merged_pos_al[merged_al_nhit]];

			Double_t tmp_pos_al = merged_pos_al[merged_al_nhit];
			if (tmp_pos_al < 4)
				merged_al_region[merged_al_nhit] = -1;
			if (tmp_pos_al >= 4 && tmp_pos_al < 28)
				merged_al_region[merged_al_nhit] = 2;
			if (tmp_pos_al >= 28 && tmp_pos_al < 48)
				merged_al_region[merged_al_nhit] = 1;
			if (tmp_pos_al >= 48 && tmp_pos_al < 79)
				merged_al_region[merged_al_nhit] = 0;
			if (tmp_pos_al >= 79 && tmp_pos_al < 99)
				merged_al_region[merged_al_nhit] = 1;
			if (tmp_pos_al >= 99 && tmp_pos_al < 124)
				merged_al_region[merged_al_nhit] = 2;
			if (tmp_pos_al >= 123)
				merged_al_region[merged_al_nhit] = -1;

			hist_merged[(int)merged_pos_al[merged_al_nhit] + 128]->Fill(merged_epi_al[merged_al_nhit]);
			histall_merged->Fill(merged_pos_al[merged_al_nhit] + 128, merged_epi_al[merged_al_nhit]);
			merged_al_nhit_list[merged_al_nhit] = (upbin - lowbin) + 1;

			for (Int_t j = lowbin; j <= upbin; j++)
			{
				pos_al[hit_al] = dist_hit->GetBinCenter(j);
				epi_al[hit_al] = dist_hit->GetBinContent(j);
				dist_hit->SetBinContent(j, 0.0);
				hit_al = hit_al + 1;
			}
			dist_hit->GetXaxis()->SetRange(1, 128);
			merged_al_nhit++;
		}
		al_nhit = hit_al;

		hist_merged_nhit->Fill(merged_pt_nhit, merged_al_nhit);
		for (Int_t i = 0; i < merged_pt_nhit; i++)
			hist_nhit_pt->Fill(merged_pt_nhit_list[i]);

		for (Int_t i = 0; i < merged_al_nhit; i++)
			hist_nhit_al->Fill(merged_al_nhit_list[i]);

		if (!save_all && (merged_pt_nhit < 1 || merged_al_nhit < 1))
			continue;

		if (merged_pt_nhit == 1 && merged_al_nhit == 1)
		{
			if (pt_nhit == 1 && al_nhit > 0)
			{
				histpt_cal_single->Fill(merged_epi_pt[0]);
				if (merged_pt_region[0] != -1)
					histpt_cal_single_each[merged_pt_region[0]]->Fill(merged_epi_pt[0]);
				histall_single->Fill(merged_pos_pt[0], merged_epi_pt[0]);
			}
			if (al_nhit == 1 && pt_nhit > 0)
			{
				histal_cal_single->Fill(merged_epi_al[0]);
				if (merged_al_region[0] != -1)
					histal_cal_single_each[merged_al_region[0]]->Fill(merged_epi_al[0]);
				histall_single->Fill(merged_pos_al[0] + 128, merged_epi_al[0]);
			}
			if (pt_nhit == 2 && al_nhit > 0)
			{
				histpt_cal_double->Fill(merged_epi_pt[0]);
				histall_double->Fill(merged_pos_pt[0], merged_epi_pt[0]);
				if (merged_pt_region[0] != -1)
				{
					histpt_cal_double_each[merged_pt_region[0]]->Fill(merged_epi_pt[0]);
					hist_vsene_pt_adj[merged_pt_region[0]]->Fill(epi_pt[0], epi_pt[1]);
					hist_diffene_pt_adj[merged_pt_region[0]]->Fill(epi_pt[0] - epi_pt[1], epi_pt[0] + epi_pt[1]);
				}
			}
			if (al_nhit == 2 && pt_nhit > 0)
			{
				histal_cal_double->Fill(merged_epi_al[0]);
				histall_double->Fill(merged_pos_al[0] + 128, merged_epi_al[0]);
				if (merged_al_region[0] != -1)
				{
					histal_cal_double_each[merged_al_region[0]]->Fill(merged_epi_al[0]);
					hist_vsene_al_adj[merged_al_region[0]]->Fill(epi_al[0], epi_al[1]);
					hist_diffene_al_adj[merged_al_region[0]]->Fill(epi_al[0] - epi_al[1], epi_al[0] + epi_al[1]);
				}
			}

			hist_image_ch->Fill(merged_pos_pt[0], merged_pos_al[0]);
			if (al_nhit == 1 && pt_nhit == 1)
				hist_image_ch_single->Fill(merged_pos_pt[0], merged_pos_al[0]);

			hist_image_dxdy->Fill(xbins_center[(int)merged_pos_pt[0]], xbins_center[(int)merged_pos_al[0]]);
			hist_image_dxdy_scale->Fill(xbins_center[(int)merged_pos_pt[0]], xbins_center[(int)merged_pos_al[0]], 1.0 / ((-xbins[(int)merged_pos_pt[0]] + xbins[(int)merged_pos_pt[0] + 1]) * (-xbins[(int)merged_pos_al[0]] + xbins[(int)merged_pos_al[0] + 1])));
			hist_slice_dx->Fill(merged_pos_pt[0]);
			hist_slice_dy->Fill(merged_pos_al[0]);
			hist_slice_dx_scale->Fill(xbins_center[(int)merged_pos_pt[0]], 1.0 / (-xbins[(int)merged_pos_pt[0]] + xbins[(int)merged_pos_pt[0] + 1]));
			hist_slice_dy_scale->Fill(xbins_center[(int)merged_pos_al[0]], 1.0 / (-xbins[(int)merged_pos_al[0]] + xbins[(int)merged_pos_al[0] + 1]));
		}

		tout->Fill();
	}

	cerr << "Event Loss rate:" << nentries - tout->GetEntries() << " out of " << nentries
			 << ":" << 100. - 100. * tout->GetEntries() / nentries << "\%" << endl;

	hist_image_ch->GetXaxis()->SetTitle("Pt ch");
	hist_image_ch->GetYaxis()->SetTitle("Al ch");

	hist_merged_nhit->GetXaxis()->SetTitle("Pt ch");
	hist_merged_nhit->GetYaxis()->SetTitle("Al ch");

	Int_t hist_color[3] = {kRed, kGreen + 2, kBlue};
	if (pdfout == 1)
	{
		TCanvas *c1 = new TCanvas("c1", "c1", 1500, 800);
		TString pdfname(inFileName);
		pdfname.ReplaceAll(".root", "_cal.pdf");
		pdfname.ReplaceAll("base/test_no", "fig/test_no");
		c1->Print(pdfname + "[");

		c1->SetLogz();
		histall_cal->Draw("colz");
		c1->Print(pdfname, "pdf");
		histall_single->Draw("colz");
		c1->Print(pdfname, "pdf");
		c1->Clear();

		c1->Divide(2, 1);
		for (Int_t i = 2; i > -1; i--)
		{
			c1->cd(1);
			gPad->SetGrid();
			histpt_cal_single_each[i]->GetXaxis()->SetRangeUser(0, 70);
			histpt_cal_single_each[i]->SetLineColor(hist_color[i]);
			histpt_cal_single_each[i]->Scale(1. / histpt_cal_single_each[i]->Integral());
			if (i == 2)
				histpt_cal_single_each[i]->Draw("HIST");
			else
				histpt_cal_single_each[i]->Draw("HISTsame");
			c1->cd(2);
			gPad->SetGrid();
			histal_cal_single_each[i]->GetXaxis()->SetRangeUser(0, 70);
			histal_cal_single_each[i]->SetLineColor(hist_color[i]);
			histal_cal_single_each[i]->Scale(1. / histal_cal_single_each[i]->Integral());
			if (i == 2)
				histal_cal_single_each[i]->Draw("HIST");
			else
				histal_cal_single_each[i]->Draw("HISTsame");
		}
		c1->Print(pdfname, "pdf");
		for (Int_t i = 2; i > -1; i--)
		{
			c1->cd(1);
			gPad->SetGrid();
			histpt_cal_double_each[i]->GetXaxis()->SetRangeUser(0, 70);
			histpt_cal_double_each[i]->SetLineColor(hist_color[i]);
			histpt_cal_double_each[i]->Scale(1. / histpt_cal_single_each[i]->Integral());
			if (i == 2)
				histpt_cal_double_each[i]->Draw("HIST");
			else
				histpt_cal_double_each[i]->Draw("HISTsame");
			c1->cd(2);
			gPad->SetGrid();
			histal_cal_double_each[i]->GetXaxis()->SetRangeUser(0, 70);
			histal_cal_double_each[i]->SetLineColor(hist_color[i]);
			histal_cal_double_each[i]->Scale(1. / histal_cal_single_each[i]->Integral());
			if (i == 2)
				histal_cal_double_each[i]->Draw("HIST");
			else
				histal_cal_double_each[i]->Draw("HISTsame");
		}
		c1->Print(pdfname, "pdf");
		for (Int_t i = 0; i < 3; i++)
		{
			c1->cd(1);
			gPad->SetGrid();
			hist_diffene_pt_adj[i]->GetXaxis()->SetRangeUser(-20, 20);
			hist_diffene_pt_adj[i]->GetYaxis()->SetRangeUser(0, 40);
			hist_diffene_pt_adj[i]->Draw("colz");
			c1->cd(2);
			gPad->SetGrid();
			hist_diffene_al_adj[i]->GetXaxis()->SetRangeUser(-20, 20);
			hist_diffene_al_adj[i]->GetYaxis()->SetRangeUser(0, 40);
			hist_diffene_al_adj[i]->Draw("colz");
			c1->Print(pdfname, "pdf");
		}
		c1->Print(pdfname + "]");
	}

	// Output File -------------------------------------
	tin_hktree->SetBranchStatus("*", 1);
	TTree *tout_hktree = tin_hktree->CloneTree();
	Long64_t nevent_cal = tout->GetEntries();
	tout_hktree->Branch("nevent_cal", &nevent_cal, "nevent_cal/L");
	tout_hktree->GetEntry(0);
	tout_hktree->Fill();

	for (Int_t i = 0; i < 256; ++i)
		s[i]->Write();
	histall_cal->Write();
	histpt_cal->Write();
	histal_cal->Write();
	histpt_cal_single->Write();
	histal_cal_single->Write();
	histpt_cal_double->Write();
	histal_cal_double->Write();

	histall_merged->Write();
	histall_single->Write();
	histall_double->Write();

	hist_image_ch->Write();
	hist_image_ch_single->Write();

	hist_image_dxdy->Write();
	hist_image_dxdy_scale->Write();

	hist_slice_dx->Write();
	hist_slice_dy->Write();
	hist_slice_dx_scale->Write();
	hist_slice_dy_scale->Write();

	hist_merged_nhit->Write();
	hist_nhit_pt->Write();
	hist_nhit_al->Write();

	for (Int_t i = 0; i < 3; i++)
	{
		histpt_cal_single_each[i]->Write();
		histpt_cal_double_each[i]->Write();
		histal_cal_single_each[i]->Write();
		histal_cal_double_each[i]->Write();
		hist_vsene_pt_adj[i]->Write();
		hist_vsene_al_adj[i]->Write();
		hist_diffene_pt_adj[i]->Write();
		hist_diffene_al_adj[i]->Write();
	}

	tout->Write();
	tout_hktree->Write();
	threshold_tree->Write();
	fout->Close();

	cout << "  " << endl;
	cout << "Output ROOT file :" << outFilename << endl;
	cout << "Finish!" << endl;

	return 0;
}
