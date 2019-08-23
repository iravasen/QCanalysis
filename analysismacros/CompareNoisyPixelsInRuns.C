#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <TH2.h>
#include <TFile.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TKey.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TText.h>
#include <TSystem.h>

using namespace std;

//Functions
std::array<long int,4> CompareTwoRuns(TH2 *href, TH2 *h2);
void SetStyle(TGraphErrors *ge, Color_t col);
void DoAnalysis(string filepath, const int nChips, bool isIB, long int refrun, int stavenum);

void CompareNoisyPixelsInRuns(){
  string fpath;
  int nchips=9;
  cout<<"\n\nAvailable file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data -Art | tail -n 500");
  cout<<"\nCopy file name: ";
  cin>>fpath;
  cout<<endl;

  string layernum = fpath.substr(fpath.find("Layer")+5, 1);
  string stavenum = fpath.substr(fpath.find("Stave")+5, 1);

  if(stoi(layernum)>=0 && stoi(layernum)<=2) nchips = 9;
  else if (stoi(layernum)==3 && stoi(layernum)==4) nchips = 54*2;
  else nchips = 98*2;

  bool isIB;
  if(nchips==9) isIB=kTRUE;
  else isIB=kFALSE;

  cout<<"Available (non-empty) runs in your file:\n"<<endl;
  TFile *infile=new TFile(Form("../Data/%s",fpath.c_str()));
  TList *list = infile->GetListOfKeys();
  TIter next(list);
  TKey *key;
  TObject *obj;
  TH2 *h2;
  while((key = ((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1"))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    h2 = (TH2*)obj->Clone(obj->GetName());
    if(!h2->GetEntries()) continue;
    string objname = (string)obj->GetName();
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_")+1, 13) : objname.substr(objname.find("_",3)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("_",1)+1, objname.find("_",3)-3);
    cout<<runnum<<" (time: "<<timestamp<<")"<<endl;
  }


  long int runnum;
  cout<<"\n\n=>Insert a run you want to use as a reference for the comparison with all the others: \n"<<endl;
  cin>>runnum;

  DoAnalysis("../Data/"+fpath, nchips, isIB, runnum, std::stoi(stavenum));
}

//
// Analysis
//
void DoAnalysis(string filepath, const int nChips, bool isIB, long int refrun, int stavenum){

  gStyle->SetOptStat(0000);

  std::vector<TH2*> hmaps;
  std::vector<string> timestamps, runnumbers;
  std::vector<array<long int,4>> noisypix;
  int posrefrun=0;

  //Read the file and the list of plots with entries
  TFile *infile=new TFile(filepath.c_str());
  TList *list = infile->GetListOfKeys();
  list->ls();
  TIter next(list);
  TKey *key;
  TObject *obj;
  TH2 *h2;
  while((key = ((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1"))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    h2 = (TH2*)obj->Clone(obj->GetName());
    if(!h2->GetEntries()) continue;
    cout<<"... Reading "<<obj->GetName()<<endl;
    hmaps.push_back(h2);
    string objname = (string)obj->GetName();
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_")+1, 13) : objname.substr(objname.find("_",3)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("_",1)+1, objname.find("_",3)-3);
    if(timestamp==std::to_string(refrun)) posrefrun=(int)hmaps.size()-1;
    if(runnum.find(std::to_string(refrun))!=string::npos) posrefrun=(int)hmaps.size()-1;
    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);
  }

  //Compare all the runs (non-empty ones) with the reference run chosen by the user
  vector<string> timelabel, runlabel;
  for(int irun=0; irun<(int)hmaps.size(); irun++){
    if(timestamps[irun]==std::to_string(refrun)) continue; // do not compare refrun with itself
    else if(runnumbers[irun].find(std::to_string(refrun))!=string::npos) continue;
    if(!hmaps[irun]->GetEntries()) continue; // do not compare ref run with empty run (= empty maps)
    noisypix.push_back(CompareTwoRuns(hmaps[posrefrun], hmaps[irun]));
    if(filepath.find("run")==string::npos) noisypix[noisypix.size()-1][0] = std::stol(timestamps[irun]);//save timestamp on run2 (the ref one is known)
    else {
      string realrunnum = runnumbers[irun];
      noisypix[noisypix.size()-1][0] = std::stol(realrunnum.erase(0,3)); //else save runnumber of run2 --> erase removes the word "run"
    }
    //cout<<noisypix[noisypix.size()-1][0]<<"  "<<noisypix[noisypix.size()-1][1]<<"  "<<noisypix[noisypix.size()-1][2]<<"  "<<noisypix[noisypix.size()-1][3]<<endl;
    timelabel.push_back(timestamps[irun]);
    runlabel.push_back(runnumbers[irun]);
  }

  //Print a warning if the run is not compared to any other run
  if(!timelabel.size())
    cout<<"\n\n WARNING: no available runs to compare with your reference run "<<refrun<<"\n\n"<<endl;

  //Make plot
  TGraphErrors *ge_nref = new TGraphErrors();
  TGraphErrors *ge_n2 = new TGraphErrors();
  TGraphErrors *ge_ncom1 = new TGraphErrors();
  TGraphErrors *ge_ncom2 = new TGraphErrors();
  double xshift = 3.;
  double max = -1;
  double min = 1e35;

  for(int icomp=(int)noisypix.size()-1; icomp>=0; icomp--){//first the older data and last the most recent
    //first couple of bar on the left
    int ipoint = (int)noisypix.size()-icomp-1;
    if(!ipoint) xshift=1.;
    else xshift = 3.;
    ge_nref->SetPoint(ipoint, ipoint*xshift, (double)noisypix[icomp][3]/2.+(double)noisypix[icomp][1]/2.);
    ge_nref->SetPointError(ipoint, 0.5, (double)noisypix[icomp][1]/2.);
    ge_ncom1->SetPoint(ipoint, ipoint*xshift, 0.);
    ge_ncom1->SetPointError(ipoint, 0.5, (double)noisypix[icomp][3]/2.);
    if((double)noisypix[icomp][3]/2.+(double)noisypix[icomp][1]/2. > max) max = (double)noisypix[icomp][3]/2.+(double)noisypix[icomp][1]/2.;

    //second couple of bar on the right
    ge_n2->SetPoint(ipoint, ipoint*xshift+1, -(double)noisypix[icomp][3]/2.-(double)noisypix[icomp][2]/2.);
    ge_n2->SetPointError(ipoint, 0.5, (double)noisypix[icomp][2]/2.);
    ge_ncom2->SetPoint(ipoint, ipoint*xshift+1, 0.);
    ge_ncom2->SetPointError(ipoint, 0.5, (double)noisypix[icomp][3]/2.);
    if(-(double)noisypix[icomp][3]/2.-(double)noisypix[icomp][2]/2. < min) min = -(double)noisypix[icomp][3]/2.-(double)noisypix[icomp][2]/2.;
  }

  //Style
  SetStyle(ge_nref, kBlue);
  SetStyle(ge_ncom1, kBlack);
  SetStyle(ge_ncom2, kBlack);
  SetStyle(ge_n2, kRed+2);

  //Draw
  TCanvas *canvas = new TCanvas("mycanvas", "mycanvas", 1300, 800);
  canvas->SetMargin(0.0092, 0.1271, 0.1759, 0.0996);
  canvas->cd();

  //fake histo (just for the axes)
  double x2,y2;
  ge_ncom2->GetPoint(ge_ncom2->GetN()-1, x2,y2);
  TH1F *hfake = new TH1F("hfake","hfake", (int)x2+6, -3, x2+3);
  hfake->Draw();
  //canvas->SetLogy();
  hfake->SetTitle(Form("%s%06ld compared to all",filepath.find("run")==string::npos? "":"run",refrun));
  ge_nref->Draw("P E2 same");
  ge_ncom1->Draw("E2 same");
  ge_ncom2->Draw("E2 same");
  ge_n2->Draw("E2 same");
  hfake->GetYaxis()->SetRangeUser(min+0.1*min, max+0.1*max);
  hfake->GetYaxis()->SetLabelColor(kWhite);
  hfake->GetYaxis()->SetTickLength(0.005);
  TLine *lineref = new TLine(-0.5, 0, x2+0.5, 0);
  lineref->SetLineColor(kGray-1);
  lineref->SetLineStyle(2);
  lineref->Draw("same");

  //draw labels on x axis
  //TAxis *ax = ge_nref->GetHistogram()->GetXaxis();
  int counter = (int)timelabel.size()-1;
  for(Int_t k=4;k<=hfake->GetNbinsX()-3;k+=3){
    hfake->GetXaxis()->SetBinLabel(k,filepath.find("run")==string::npos ? timelabel[counter].c_str():runlabel[counter].c_str());
    counter--;
  }

  //draw legend
  TLegend *leg = new TLegend(0.876,0.176, 0.994, 0.902);
  leg->SetLineColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(ge_nref, "#splitline{#noisy pix}{ref. run only}", "f");
  leg->AddEntry(ge_n2, "#splitline{#noisy pix}{2nd run only}", "f");
  leg->AddEntry(ge_ncom1, "#splitline{#noisy pix}{both}");
  leg->Draw("same");

  canvas->SaveAs(Form("../Plots/%sStave%d_NoisyPixComparison_%s%ld_compared_to_run_%s.pdf", isIB?"IB":"OB",stavenum,filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  canvas->SaveAs(Form("../Plots/%sStave%d_NoisyPixComparison_%s%ld_compared_to_run_%s.root", isIB?"IB":"OB",stavenum,filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

}

//
// Function to compare two hitmaps --> returns an arrays with timestamp of run2, noisyPixInRefRun, noisyPixInRun2, noisyPixInCommon
//
std::array<long int,4> CompareTwoRuns(TH2 *href, TH2 *h2){

  std::array<long int,4> noisypix = {0, 0, 0, 0};
  //number of noisy pix in refrun_only and in common
  for(int ixbin=1; ixbin<=href->GetXaxis()->GetNbins(); ixbin++){
    for(int iybin=1; iybin<=href->GetYaxis()->GetNbins(); iybin++){
      if(href->GetBinContent(ixbin, iybin)>1 && h2->GetBinContent(ixbin, iybin)>1){//noisy in both runs
        noisypix[3]++;
      }
      else if(href->GetBinContent(ixbin, iybin)>1 && h2->GetBinContent(ixbin, iybin)==0){//noisy only in ref run
        noisypix[1]++;
      }
      else if(href->GetBinContent(ixbin, iybin)==0 && h2->GetBinContent(ixbin, iybin)>1){//noisy only in second run
        noisypix[2]++;
      }
      else continue;
    }
  }

  return noisypix;
}

//
// Style
//
void SetStyle(TGraphErrors *ge, Color_t col){
  ge->SetMarkerStyle(0);
  ge->SetMarkerColor(col);
  ge->SetFillColor(col);
  ge->SetLineColor(col);
}
