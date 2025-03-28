// generate genCalDataTree.root for caliburation based on en.txt.
// first use co_peaktotxt.C & pedestal_peaktotxt.C to make en.txt
// en.txt -> energy vs adc

#include <iostream>
#include <utility>
#include <vector>
#include <array>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TString.h>

void sort_multi(std::vector<Double_t> v1, std::vector<Double_t> v2){
    int n = v1.size();
    for(int i=0; i<n-1; ++i){
        for(int j=i+1; j<n; ++j){
            if(v1[i]>v1[j]){
                std::swap(v1[i],v1[j]);
                std::swap(v2[i],v2[j]);
            }
        }
    }
}

template <typename T>
void putToLast(T* array,Int_t init, Int_t end){

  for (int i = init; i < end; ++i){
    swap(array[i],array[i+1]);
  }
}

Int_t genCalTree(TString ifn1 = "en.txt", TString ofn = "enCalDataTree.root", int useReverse = 0){

  std::cerr << "making enCalDataTree.root" << std::endl;
  // STEP1: Set input filename
	// TString ifn1 = "enCal.dat";

  // STEP2: Create TTree
  TTree *tree1 = new TTree("enCal", "tree using ReadFile()");

  // STEP3: Read data using TTree::ReadFile(...) method
  // STEP4: Create TFile to save TTree
  //  TString ofn = "enCalDataTree.root";
  TFile *fout = new TFile(ofn, "recreate");


  // make tree1 from text file
  tree1->ReadFile(ifn1.Data(), "ch/I:keV/D:adc/D");

  std::array< std::vector<Double_t>, 256> adc;
  std::array< std::vector<Double_t>, 256> keV;


  // read out tree1
  Int_t ch_temp,remapped_ch;
  Double_t adc_temp,keV_temp;
  tree1->SetBranchAddress("ch", &ch_temp);
  tree1->SetBranchAddress("adc", &adc_temp);
  tree1->SetBranchAddress("keV", &keV_temp);

  // set tree1 data (ch, line energy, adc value) into vectors
  Int_t nentries = tree1->GetEntries();
  //  std::cout << nentries << std::endl;
  for (int index = 0; index < nentries; ++index) {
    tree1->GetEntry(index);
    // if(adc_temp<0){continue;}

    if(useReverse){  // pt side channel non-reverse mode
      //      std::cout << "use reverse mode" << std::endl;
      if(ch_temp < 128){
        remapped_ch = ch_temp;
      }else if(ch_temp < 192){
        remapped_ch = ch_temp+64;
      }else if(ch_temp < 256){
        remapped_ch = ch_temp-64;
      }
    }else{            // pt side channel non-reverse mode
      remapped_ch = ch_temp;
    }
    adc[remapped_ch].push_back(adc_temp);
    keV[remapped_ch].push_back(keV_temp);
  }

  for (int i = 0; i < 256; ++i) {
  	if(adc[i].size() <= 1){
  	  adc[i].push_back(100);
  	  keV[i].push_back(0);
  	  adc[i].push_back(10);
  	  keV[i].push_back(0);
  	}
  }

  // sort adc in energy-incremental way
  for (int index = 0; index < 256; ++index) {
    // sort_multi(keV[index],adc[index]);
    std::sort(keV[index].begin(),keV[index].end() );
    std::sort(adc[index].begin(),adc[index].end() );

    // extrapolate linearly from two last points
    int n = keV[index].size();
    double x1 = keV[index][n-2];
    double x2 = keV[index][n-1];
    keV[index].push_back(2.*x2-x1);

    int m = adc[index].size();
    double y1 = adc[index][m-2];
    double y2 = adc[index][m-1];
    adc[index].push_back(2.*y2-y1);
  }


  // make TGraph and spline function
  TSpline3 *s[256];
  TGraph* gr[256];

  for (int i = 0; i < 256; ++i) {
    gr[i] = new TGraph(adc[i].size(),&adc[i][0],&keV[i][0]);
    s[i] = new TSpline3(Form("grs%d",i),gr[i],"e2b2");
    gr[i]->SetTitle(Form("gr%d",i));
    gr[i]->SetName(Form("gr%d",i));
    s[i]->SetTitle(Form("spline%d",i));
    s[i]->SetName(Form("spline%d",i));
  }


  TCanvas *c1 = new TCanvas("c1","title",1000,600);
  for(Int_t i=0;i<256;i++){
    gr[i]->Write();
    s[i]->Write();
    // gr[i]->Draw("same");
  }

  const Int_t panelNum = 4;
  TCanvas *c[panelNum];
  for(Int_t i=0;i<panelNum;i++){
    c[i] = new  TCanvas(Form("c%d",i),"Spline curves",1200,1500);
    c[i]->Divide(8,256/8/panelNum);
    c[i]->SetGrid();
  }

  for(Int_t k=0;k<panelNum;k++){
    for(Int_t i=0;i<256/panelNum;i++){
      c[k]->cd(i+1)->SetGrid();
      gr[i+k*256/panelNum]->SetMarkerColor(2);
      gr[i+k*256/panelNum]->Draw("AP*");
      s[i+k*256/panelNum]->Draw("same");
    }
  }
  ofn.ReplaceAll(".root",".pdf");
  TString pdfname = ofn;
  c[0]->SaveAs(pdfname+"(");
  for(Int_t k=1;k<panelNum-1;k++){
    c[k]->SaveAs(pdfname);
  }
  c[panelNum-1]->SaveAs(pdfname+")");

  // STEP5: Write TTree to TFile
  tree1->Write();

  // STEP6: Close TFile
  fout->Close();

	return 0;
}

int main(int argc, char const *argv[]){
  if(argc==2){
    genCalTree(argv[1]);
  }else if(argc==3){
    genCalTree(argv[1],argv[2]);
  }else if(argc==4){
    int opt=0;
    if(strncmp(argv[3],"-reverse",8)==0){
      opt = 1;
    }
    genCalTree(argv[1],argv[2],opt);
  }else{
    std::cout << "USAGE: one input file is required" << std::endl;
  }
  return 0;
}
