#include <sstream>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <string>
#include <iostream>
#include <vector>
#include <ctime>
#include <TH2.h>
#include <TFile.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TKey.h>
#include <TLine.h>

//
// Function to compare two hitmaps --> returns an arrays with timestamp of run2, noisyPixInRefRun, noisyPixInRun2, noisyPixInCommon
//
std::array<long int,5> CompareTwoRuns(TH2 *href, TH2 *h2){

  std::array<long int,5> noisypix = {0, 0, 0, 0, 0};
  //number of noisy pix in refrun_only and in common
  for(int ixbin=1; ixbin<=href->GetXaxis()->GetNbins(); ixbin++){
    for(int iybin=1; iybin<=href->GetYaxis()->GetNbins(); iybin++){
      if(href->GetBinContent(ixbin, iybin)>16 && h2->GetBinContent(ixbin, iybin)>16){//noisy in both runs
        noisypix[2]++;
      }
      else if(href->GetBinContent(ixbin, iybin)>16 && h2->GetBinContent(ixbin, iybin)==0){//noisy only in ref run
        noisypix[0]++;
      }
      else if(href->GetBinContent(ixbin, iybin)==0 && h2->GetBinContent(ixbin, iybin)>16){//noisy only in second run
        noisypix[1]++;
      }
      else continue;
    }
  }

  /*cout<<"Stave/run: "<<h2->GetName()<<endl;
  cout<<"Ref run:  "<<noisypix[0]<<endl;
  cout<<"Sec run:  "<<noisypix[1]<<endl;
  cout<<"Both run: "<<noisypix[2]<<endl;*/

  return noisypix;
}

//
//Function to return the number of active (i.e. enabled) chips in a stave
//
int GetNchipsActive(TH2 *hmap, int maxchip){
  int ix1=1;
  int ix2=256;
  int activechips = maxchip;
  while(ix2<=hmap->GetNbinsX()){
    if((hmap->Integral(ix1,ix2,1,128)<1e-20))
      activechips--;
    ix1=ix2+1;
    ix2+=256;
  }
  return activechips;
}

//
//Set Style of a TGraph
//
void SetStyle(TGraph *h, Int_t col, Style_t mkr){
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.4);
  h->SetMarkerColor(col);
  //h->SetFillStyle(0);
  //h->SetFillColorAlpha(col,0.8);
}

//
//Set Style of a TH1
//
void SetStyle(TH1 *h, Int_t col, Style_t mkr){
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.4);
  h->SetMarkerColor(col);
  //h->SetFillStyle(0);
  //h->SetFillColorAlpha(col,0.8);
}

//
// Set Style of TGraphErrors
//
void SetStyle(TGraphErrors *ge, Color_t col){
  ge->SetMarkerStyle(0);
  ge->SetMarkerColor(col);
  ge->SetFillColor(col);
  ge->SetLineColor(col);
}


//
//Function to return the number of runs in which a stave has no hits
//
int GetNrunsWOhits(TH2 *hFhrStv){
  int runswohits = 0;
  for(int iy=1; iy<=hFhrStv->GetNbinsY(); iy++){
    if(hFhrStv->GetBinContent(1,iy)<1e-15)
      runswohits++;
  }
  return runswohits;
}

//
// Get current date and time
//
string GetCurrentDateTime(bool useLocalTime) {
    stringstream currentDateTime;
    // current date/time based on current system
    time_t ttNow = time(0);
    tm * ptmNow;

    if (useLocalTime)
        ptmNow = localtime(&ttNow);
    else
        ptmNow = gmtime(&ttNow);

    currentDateTime << 1900 + ptmNow->tm_year;

    //month
    if (ptmNow->tm_mon < 9)
        //Fill in the leading 0 if less than 10
        currentDateTime << "0" << 1 + ptmNow->tm_mon;
    else
        currentDateTime << (1 + ptmNow->tm_mon);

    //day
    if (ptmNow->tm_mday < 10)
        currentDateTime << "0" << ptmNow->tm_mday << "_";
    else
        currentDateTime <<  ptmNow->tm_mday << "_";

    //hour
    if (ptmNow->tm_hour < 10)
        currentDateTime << "0" << ptmNow->tm_hour;
    else
        currentDateTime << ptmNow->tm_hour;

    //min
    if (ptmNow->tm_min < 10)
        currentDateTime << "0" << ptmNow->tm_min;
    else
        currentDateTime << ptmNow->tm_min;

    //sec
    if (ptmNow->tm_sec < 10)
        currentDateTime << "0" << ptmNow->tm_sec;
    else
        currentDateTime << ptmNow->tm_sec;


    return currentDateTime.str();
}

//
// Function to compare two hitmaps --> returns an arrays with timestamp of run2, noisyPixInRefRun, noisyPixInRun2, noisyPixInCommon
//
std::array<float,nMasked+1> GetFHRwithMasking(TH2 *hmap, const int nchips, double ntrig, TH2 *hhotmap){

  array<float,nMasked+1> fhrstave;
  TH2 *hmapclone = (TH2*)hmap->Clone(Form("hmapclone_%s",hmap->GetName()));

  for(int iter=0; iter<nMasked+1; iter++){
    long int totalhits = hmapclone->Integral();
    float fhr = (float)totalhits / (512.*1024.*nchips*ntrig);
    if(!nchips) fhr=0.;
    //cout<<fhr<<endl;
    fhrstave[iter] = fhr;
    int binmax = hmapclone->GetMaximumBin();
    int binmax_x, binmax_y, binmax_z;
    hmapclone->GetBinXYZ(binmax, binmax_x, binmax_y, binmax_z);
    hmapclone->SetBinContent(binmax_x, binmax_y, 0);
    if(totalhits!=0) hhotmap->SetBinContent(binmax_x, binmax_y, 1);
    //cout<<iter<<" FHR: "<<fhr<<endl;
    //cout<<"Masking binx: "<<binmax_x<<"  biny: "<<binmax_y<<endl;
  }

  delete hmapclone;

  return fhrstave;

}
