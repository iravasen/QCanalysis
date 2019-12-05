#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TText.h>
#include <TSystem.h>
#include <TKey.h>

using namespace std;

const int nMasked = 10;
//Functions
std::array<float,nMasked+1> GetFHRwithMasking(TH2 *hmap, const int nchips, double ntrig);
void SetStyle(TH1 *h, Int_t col, Style_t mkr);
void DoAnalysis(string filepath_hit, const int nChips, bool isIB);


void MaskNoisyPixelsInRuns(){
  string fpath;

  int nchips=9;
  cout<<"\n\nAvailable file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*FHRMAPS_HITMAPS* -Art | tail -n 500");
  cout<<"\nCopy file name: ";
  cin>>fpath;
  cout<<endl;

  bool isIB;
  if(fpath.find("IB")!=string::npos){
    isIB = kTRUE;
  }
  else{
    string layernum = fpath.substr(fpath.find("Layer")+5, 1);
    if(stoi(layernum)>=0 && stoi(layernum)<=2) nchips = 9;
    else if (stoi(layernum)==3 && stoi(layernum)==4) nchips = 54*2;
    else nchips = 98*2;
    if(nchips==9) isIB=kTRUE;
    else isIB=kFALSE;
  }

  DoAnalysis(fpath, nchips, isIB);
}

//
// Analysis
//
void DoAnalysis(string filepath_hit, const int nChips, bool isIB){

  gStyle->SetOptStat(0000);

  Int_t col[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#0033a1")};

  std::vector<TH2*> hmaps;
  std::vector<string> timestamps, runnumbers, stavenums, laynums;
  vector<int> posrefrun;
  int nLayers=1, nRuns=1;
  vector<int> nStavesInLay;
  vector<string> runlabel;
  int cstv=0;

  //Read the file and the list of plots with entries (hitmaps!)
  TFile *infile=new TFile(filepath_hit.c_str());
  TList *list = (TList*)infile->GetListOfKeys();
  TIter next(list);
  TKey *key;
  TObject *obj;
  while((key=(TKey*)next())){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1"))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("Stv")==string::npos) continue;
    TH2 *h2 = (TH2*)obj;
    //if(!h2->GetEntries()) continue;
    cout<<"... Reading "<<obj->GetName()<<endl;
    hmaps.push_back(h2);
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);
    string stvnum = objname.substr(objname.find("Stv")+3,2);
    if(stvnum.find("_")!=string::npos){
      stvnum = objname.substr(objname.find("Stv")+3,1);
    }

    if((int)stavenums.size()>1 && stvnum!=stavenums[stavenums.size()-1])
      cstv++;

    if((int)laynums.size()>1 && laynum!=laynums[laynums.size()-1]){
      nLayers++;
      nStavesInLay.push_back(cstv);
      cstv=0;
    }

    if((int)stavenums.size()>1 && stvnum==stavenums[stavenums.size()-1])
      nRuns++;
    else nRuns=1;

    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);
    laynums.push_back(laynum);
    stavenums.push_back(stvnum);

    if((int)stavenums.size()>1 && stvnum!=stavenums[stavenums.size()-2])
      continue;
    runlabel.push_back(runnum);

  }
  nStavesInLay.push_back(cstv+1);//in case of 1 layer or for last layer

  //Read file with fhr maps for each layer --> ONLY TO EXTRACT NUMBER OF TRIGGERS!!
  TList *list2 = (TList*)infile->GetListOfKeys();
  TIter next2(list2);
  TH2 *h2_2[nRuns];
  TKey *key2;
  TObject *obj2;
  int cntrun=0;
  while((key2=((TKey*)next2()))){
    obj2 = key2->ReadObj();
    if ((strcmp(obj2->IsA()->GetName(),"TProfile")!=0)
         && (!obj2->InheritsFrom("TH2"))
	       && (!obj2->InheritsFrom("TH1"))
       ) {
            cout<<"<W> Object "<<obj2->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    string objname2 = (string)obj2->GetName();
    if(objname2.find("Stv")!=string::npos) break;
    h2_2[cntrun] = (TH2*)obj2;
    cntrun++;
    if(cntrun==nRuns) break;
    //if(!h2->GetEntries()) continue;
    cout<<"... Reading "<<obj2->GetName()<<endl;
  }

  //Calculate number of triggers for each run
  cout<<endl;
  cout<<"... Extract the number of triggers (events) from for each run"<<endl;
  vector<double> ntrig;
  for(int ir=0; ir<nRuns; ir++){
    double fhr_run = h2_2[ir]->GetBinContent(1,1);
    double hits_chip = hmaps[ir]->Integral(1,256,1,512);
    ntrig.push_back(hits_chip/(512.*1024.*fhr_run));
    if(!ir) cout<<"~"<<ntrig[ntrig.size()-1]<<" triggers"<<endl;
  }

  //Start masking hottest pixels for each stave in each run
  cout<<endl;
  cout<<"... Analysing FHR with (hot) pixel masking"<<endl;
  vector<array<float,nMasked+1>> fhrall;
  int irun = nRuns-1;
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){ //start from the bottom in order to start with the oldest run
    fhrall.push_back(GetFHRwithMasking(hmaps[ihist],nChips,ntrig[irun]));
    irun--;
    if(ihist>0){
      if(stavenums[ihist-1]!=stavenums[ihist]){
        irun=nRuns-1;
      }
    }
  }

  TH2F *hFhrStv[nLayers][100];
  for(int ilay=0; ilay<nLayers; ilay++)
    for(int is=0; is<nStavesInLay[ilay]; is++)
      hFhrStv[ilay][is] = new TH2F(Form("h2FhrStv_%s_%d", laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),is), Form("Layer-%s - Stave-%d; # Hot Pixel Masked;Run", laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),is),nMasked+1, -0.5, nMasked+0.5, nRuns, 0.5, nRuns+0.5);

  //Fill histogram
  int ilayer = nLayers-1;
  int istave = nStavesInLay[ilayer]-1;
  irun=0;
  for(int i=0; i<(int)fhrall.size(); i++){

    for(int ifhr=0; ifhr<(int)fhrall[i].size(); ifhr++){
      hFhrStv[ilayer][istave]->SetBinContent(ifhr+1, irun+1, fhrall[i][ifhr]);
    }
    irun++;
    if(i<(int)fhrall.size()-1){// in case ref run is the first into the list of runs
      if(stavenums[fhrall.size()-i-2]!=stavenums[fhrall.size()-i-1]){
        istave--;
        irun=0;
      }
      if(laynums[fhrall.size()-i-2]!=laynums[fhrall.size()-i-1]){
        ilayer--;
        istave = nStavesInLay[ilayer]-1;
      }
    }
  }

  //Make FHR (averaged on all runs) vs #masked pix for all staves in a layer
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas cnv(Form("cnv_%d",ilay), Form("cnv_%d",ilay));
    cnv.cd();
    cnv.SetLogy();
    cnv.SetTickx();
    cnv.SetTicky();
    cnv.SetMargin(0.0988,0.1,0.1,0.0993);
    TH1F *hframe = cnv.DrawFrame(-0.5,1e-14,10.5,1e-3,Form("Layer %s - Average FHR %s; # Hot Pixel masked ; FHR (/event/pixel)",laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str()));
    //legend
    TLegend leg(0.904, 0.127,0.997,0.898);
    for(int is=0; is<nStavesInLay[ilay];is++){
      TH1F *proj = (TH1F*)hFhrStv[ilay][is]->ProjectionX(Form("proj_%d%d",ilay,is));
      proj->Scale(1./nRuns); //Divide by the number of runs
      SetStyle(proj, col[is<nStavesInLay[ilay]/2 ? is : is-nStavesInLay[ilay]/2],is<nStavesInLay[ilay]/2 ? 24:26);
      proj->Draw("PL same");
      leg.AddEntry(proj, Form("Stv%d",is),"pl");
    }
    leg.Draw("same");
    cnv.SaveAs(Form("../Plots/Layer%s_FHRpixmask_%s.pdf", laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str()));
    cnv.SaveAs(Form("../Plots/Layer%s_FHRpixmask_%s.root", laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str()));
  }

}

//
// Function to compare two hitmaps --> returns an arrays with timestamp of run2, noisyPixInRefRun, noisyPixInRun2, noisyPixInCommon
//
std::array<float,nMasked+1> GetFHRwithMasking(TH2 *hmap, const int nchips, double ntrig){

  array<float,nMasked+1> fhrstave;

  for(int iter=0; iter<nMasked+1; iter++){
    long int totalhits = hmap->Integral();
    float fhr = (float)totalhits / (512.*1024.*nchips*ntrig);
    fhrstave[iter] = fhr;
    int binmax = hmap->GetMaximumBin();
    int binmax_x, binmax_y, binmax_z;
    hmap->GetBinXYZ(binmax, binmax_x, binmax_y, binmax_z);
    hmap->SetBinContent(binmax_x, binmax_y, 0);
    //cout<<iter<<" FHR: "<<fhr<<endl;
    //cout<<"Masking binx: "<<binmax_x<<"  biny: "<<binmax_y<<endl;
  }

  return fhrstave;

}
//
//Set Style
//
void SetStyle(TH1 *h, Int_t col, Style_t mkr){
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.4);
  h->SetMarkerColor(col);
  //h->SetFillStyle(0);
  //h->SetFillColorAlpha(col,0.8);
}
