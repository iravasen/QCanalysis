#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <TH2.h>
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
#include <THnSparse.h>

#include "inc/constants.h"

using namespace std;

//Functions
std::array<long int,5> CompareTwoRuns(THnSparse *href, THnSparse *h2);
void SetStyle(TGraphErrors *ge, Color_t col);
void DoAnalysis(string filepath, const int nChips, bool isIB, string skipruns, long int refrun);

void CompareDeadPixelsInRuns(){
  string fpath;
  int nchips=9;
  cout<<"\n\nAvailable file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*THRMAPS_DEADPIXMAPS* -Art | tail -n 500");
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

  //Choose whether to skip runs
  string skipans, skipruns;
  cout<<endl;
  cout<<"Would you like to skip some run(s)? [y/n] ";
  cin>>skipans;
  if(skipans=="y" || skipans=="Y"){
    cout<<endl;
    cout<<"Specify run number(s) separated by comma (no white spaces!):";
    cin>>skipruns;
    cout<<endl;
  }
  else
    skipruns=" ";

  cout<<"Available runs in your file:\n"<<endl;
  TFile *infile=new TFile(fpath.c_str());
  TList *list = (TList*)infile->GetListOfKeys();
  TIter next(list);
  TObject *obj;
  TKey *key;
  vector<string> stavenums;
  while((key = ((TKey*)next()))){
    obj=key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1")
         && (!obj->InheritsFrom("THnSparse")))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D, 2D, Sparse histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("Stv")==string::npos) continue;

    string runnum =  objname.substr(objname.find("run")+3, 6);
    string stvnum = objname.substr(objname.find("Stv")+3,2);
    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user
    if(stvnum.find("_")!=string::npos)
      stvnum = objname.substr(objname.find("Stv")+3,1);
    stavenums.push_back(stvnum);
    //cout<<"run: "<<runnum<<"   timestamp: "<<timestamp<<"    laynum: "<<laynum<<endl;
    if((int)stavenums.size()>1 && stvnum!=stavenums[stavenums.size()-2])
      break;

    cout<<runnum<<endl;
  }

  long int refrun;
  cout<<"\n\n=>Insert a run you want to use as a reference for the comparison with all the others: \n"<<endl;
  cin>>refrun;

  DoAnalysis(fpath, nchips, isIB, skipruns, refrun);
}

//
// Analysis
//
void DoAnalysis(string filepath, const int nChips, bool isIB, string skipruns, long int refrun){

  gStyle->SetOptStat(0000);

  std::vector<THnSparse*> hmaps;
  std::vector<string> timestamps, runnumbers, stavenums, laynums;
  int nLayers=1, nRuns=1;
  vector<string> runlabel;

  //Read the file and the list of plots with entries
  TFile *infile=new TFile(filepath.c_str());
  TList *list = (TList*)infile->GetListOfKeys();
  TIter next(list);
  TObject *obj;
  TKey *key;
  THnSparse *hsparse;
  while((key=((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1")
         && (!obj->InheritsFrom("THnSparse")))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D, 2D, sparse histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("Stv")==string::npos) continue;
    hsparse = (THnSparse*)obj;
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);
    string stvnum = objname.substr(objname.find("Stv")+3,2);
    if(stvnum.find("_")!=string::npos){
      stvnum = objname.substr(objname.find("Stv")+3,1);
    }

    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user

    cout<<"... Reading "<<obj->GetName()<<endl;
    hmaps.push_back(hsparse);

    if((int)laynums.size()>1 && laynum!=laynums[laynums.size()-1]){
      nLayers++;
    }

    if((int)stavenums.size()>1 && stvnum==stavenums[stavenums.size()-1]){
      nRuns++;
    }
    else {
      nRuns=1;
    }

    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);
    laynums.push_back(laynum);
    stavenums.push_back(stvnum);


  }

  //set runlabels
  for(int i=0; i<nRuns; i++){
    runlabel.push_back(runnumbers[i]);
  }

  //Compare all the runs (non-empty ones) with the reference run chosen by the user
  long int first[nLayers][100], second[nLayers][100], both[nLayers][100];
  bool filled[nLayers][100];
  vector<array<long int,5>> noisypix;
  for(int ilay=0; ilay<nLayers; ilay++)
    for(int i=0; i<100; i++){
      first[ilay][i]=0; second[ilay][i]=0; both[ilay][i]=0;
      filled[ilay][i] = false;
    }
  int posrefrun[nLayers>1 ? 48:nStavesInLay[stoi(laynums[0])]];
  int stop = (nLayers>1) ? 48:nStavesInLay[stoi(laynums[0])];
  for(int is=0; is < stop; is++)
    posrefrun[is] = -1;
  //find ref run positions
  for(int ihist=0; ihist<(int)hmaps.size(); ihist++){
    string hname = hmaps[ihist]->GetName();
    string runn =  hname.substr(hname.find("run")+3, 6);
    string stvnum = hname.substr(hname.find("Stv")+3,2);
    if(stvnum.find("_")!=string::npos){
      stvnum = hname.substr(hname.find("Stv")+3,1);
    }
    int snum = stoi(stvnum);
    int lnum = stoi(hname.substr(hname.find("L")+1,1));
    int staveindex = snum;
    if(nLayers>1){
      staveindex = lnum==0 ? snum: lnum==1? snum+nStavesInLay[0]:snum+nStavesInLay[0]+nStavesInLay[1];
    }
    if(stol(runn)==refrun){ //position of refence run
      posrefrun[staveindex] = ihist;
    }
  }

  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){ //start from the bottom in order to start with the oldest run

    string hname = hmaps[ihist]->GetName();
    string runn =  hname.substr(hname.find("run")+3, 6);
    string lnum =  hname.substr(hname.find("L")+1,1);
    string stvnum = hname.substr(hname.find("Stv")+3,2);
    if(stvnum.find("_")!=string::npos){
      stvnum = hname.substr(hname.find("Stv")+3,1);
    }
    int snum = stoi(stvnum);
    auto itf = find(runlabel.begin(), runlabel.end(), runn);
    int irun = distance(runlabel.begin(), itf);
    cout<<irun<<endl;
    int ilayer = nLayers>1 ? stoi(lnum) : 0;
    int istave = snum;
    if(nLayers>1)
      istave = (lnum=="0") ? snum : (lnum=="1")? snum+nStavesInLay[0]:snum+nStavesInLay[0]+nStavesInLay[1];

    if(hmaps[ihist]->GetEntries()>1e4){
      cout<<"L"<<lnum<<"_"<<snum<<" - run"<<runn<<" skipped because has more than 10000 entries (probably a bad run)."<<endl;
      continue;
    }

    if(runnumbers[ihist].find(std::to_string(refrun))!=string::npos){
      continue;
    }
    if(posrefrun[istave]==-1){
      continue; //skip comparison if no run is found
    }

    noisypix.push_back(CompareTwoRuns(hmaps[posrefrun[istave]], hmaps[ihist]));
    noisypix[noisypix.size()-1][3] = stol(stavenums[ihist]);
    noisypix[noisypix.size()-1][4] = stol(laynums[ihist]);
    first[ilayer][irun]+=noisypix[noisypix.size()-1][0];
    second[ilayer][irun]+=noisypix[noisypix.size()-1][1];
    both[ilayer][irun]+=noisypix[noisypix.size()-1][2];
    filled[ilayer][irun] = true;

    //cout<<"L"<<ilayer<<"_"<<istave<<" "<<noisypix[noisypix.size()-1][0]<<"  "<<noisypix[noisypix.size()-1][1]<<"  "<<noisypix[noisypix.size()-1][2]<<"  "<<noisypix[noisypix.size()-1][3]<<"  "<<noisypix[noisypix.size()-1][4]<<endl;

  }//end loop on histograms

  //Make plot for each layer and for each stave in the root file
  TGraphErrors *ge_nref[nLayers];
  TGraphErrors *ge_n2[nLayers];
  TGraphErrors *ge_ncom1[nLayers];
  TGraphErrors *ge_ncom2[nLayers];
  TGraphErrors *ge_nref_stave[nLayers][100];//100 is an abitrary large number of staves
  TGraphErrors *ge_n2_stave[nLayers][100];
  TGraphErrors *ge_ncom1_stave[nLayers][100];
  TGraphErrors *ge_ncom2_stave[nLayers][100];
  double xshift = 3.;
  double max[nLayers];
  double min[nLayers];
  for(int ilay=0; ilay<nLayers; ilay++){
    ge_nref[ilay] = new TGraphErrors();
    ge_n2[ilay] = new TGraphErrors();
    ge_ncom1[ilay] = new TGraphErrors();
    ge_ncom2[ilay] = new TGraphErrors();

    max[ilay] = -1.;
    min[ilay] = 1e35;

    int ipoint = 0;
    for(int ir=0; ir<100; ir++){//first the older data and last the most recent
      if(!filled[ilay][ir]) continue; //skip if not filled
      if(!ipoint) xshift=1.;
      else xshift = 3.;
      ge_nref[ilay]->SetPoint(ipoint, ipoint*xshift, (double)both[ilay][ir]/2.+(double)first[ilay][ir]/2.);
      ge_nref[ilay]->SetPointError(ipoint, 0.5, (double)first[ilay][ir]/2.);
      ge_ncom1[ilay]->SetPoint(ipoint, ipoint*xshift, 0.);
      ge_ncom1[ilay]->SetPointError(ipoint, 0.5, (double)both[ilay][ir]/2.);
      if((double)both[ilay][ir]/2.+(double)first[ilay][ir] > max[ilay]) max[ilay] = (double)both[ilay][ir]/2.+(double)first[ilay][ir];

      //second couple of bar on the right
      ge_n2[ilay]->SetPoint(ipoint, ipoint*xshift+1, -(double)both[ilay][ir]/2.-(double)second[ilay][ir]/2.);
      ge_n2[ilay]->SetPointError(ipoint, 0.5, (double)second[ilay][ir]/2.);
      ge_ncom2[ilay]->SetPoint(ipoint, ipoint*xshift+1, 0.);
      ge_ncom2[ilay]->SetPointError(ipoint, 0.5, (double)both[ilay][ir]/2.);
      if(-(double)both[ilay][ir]/2.-(double)second[ilay][ir] < min[ilay]) min[ilay] = -(double)both[ilay][ir]/2.-(double)second[ilay][ir];
      ipoint++;
    }//end first loop on runs

    //Style
    SetStyle(ge_nref[ilay], kBlue);
    SetStyle(ge_ncom1[ilay], kBlack);
    SetStyle(ge_ncom2[ilay], kBlack);
    SetStyle(ge_n2[ilay], kRed+2);
  }//end loop on layers

  //FILL PLOTS FOR EACH STAVE
  //fill the plots for each stave
  /*int cnt = 0;
  int nstaves = posrefrun.size();
  double maxs[nLayers][100];
  double mins[nLayers][100];
  while(cnt<(int)noisypix.size()){
    int isidx = noisypix[cnt][3];
    int ilay  = noisypix[cnt][4];
    maxs[ilay][isidx] = -1.;
    mins[ilay][isidx] = 1e35;
    ge_nref_stave[ilay][isidx] = new TGraphErrors();
    ge_ncom1_stave[ilay][isidx] = new TGraphErrors();
    ge_n2_stave[ilay][isidx] = new TGraphErrors();
    ge_ncom2_stave[ilay][isidx] = new TGraphErrors();
    for(int ir=0; ir<nRuns-1; ir++){//first the older data and last the most recent
      if(!ir) xshift=1.;
      else xshift = 3.;
      ge_nref_stave[ilay][isidx]->SetPoint(ir, ir*xshift, (double)noisypix[cnt][2]/2.+(double)noisypix[cnt][0]/2.);
      ge_nref_stave[ilay][isidx]->SetPointError(ir, 0.5, (double)noisypix[cnt][0]/2.);
      ge_ncom1_stave[ilay][isidx]->SetPoint(ir, ir*xshift, 0.);
      ge_ncom1_stave[ilay][isidx]->SetPointError(ir, 0.5, (double)noisypix[cnt][2]/2.);
      if((double)noisypix[cnt][2]/2.+(double)noisypix[cnt][0] > maxs[ilay][isidx]) maxs[ilay][isidx] = (double)noisypix[cnt][2]/2.+(double)noisypix[cnt][0];

      //second couple of bar on the right
      ge_n2_stave[ilay][isidx]->SetPoint(ir, ir*xshift+1, -(double)noisypix[cnt][2]/2.-(double)noisypix[cnt][1]/2.);
      ge_n2_stave[ilay][isidx]->SetPointError(ir, 0.5, (double)noisypix[cnt][1]/2.);
      ge_ncom2_stave[ilay][isidx]->SetPoint(ir, ir*xshift+1, 0.);
      ge_ncom2_stave[ilay][isidx]->SetPointError(ir, 0.5, (double)noisypix[cnt][2]/2.);
      if(-(double)noisypix[cnt][2]/2.-(double)noisypix[cnt][1] < mins[ilay][isidx]) mins[ilay][isidx] = -(double)noisypix[cnt][2]/2.-(double)noisypix[cnt][1];
      cnt++;
    }
    SetStyle(ge_nref_stave[ilay][isidx], kBlue);
    SetStyle(ge_ncom1_stave[ilay][isidx], kBlack);
    SetStyle(ge_ncom2_stave[ilay][isidx], kBlack);
    SetStyle(ge_n2_stave[ilay][isidx], kRed+2);
  }*/

  //Legend
  TLegend *leg = new TLegend(0.876,0.176, 0.994, 0.902);
  leg->SetLineColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(ge_nref[0], "#splitline{#dead pix}{ref. run only}", "f");
  leg->AddEntry(ge_n2[0], "#splitline{#dead pix}{2nd run only}", "f");
  leg->AddEntry(ge_ncom1[0], "#splitline{#dead pix}{both}");

  //Draw plot for each layer
  for(int ilay=0; ilay<nLayers; ilay++){

    TCanvas *canvas = new TCanvas(Form("mycanvas_%d",ilay), Form("mycanvas_%d",ilay), 1300, 800);
    canvas->SetMargin(0.08, 0.1271, 0.1759, 0.0996);
    canvas->cd();

    //fake histo (just for the axes)
    double x2,y2;
    ge_ncom2[ilay]->GetPoint(ge_ncom2[ilay]->GetN()-1, x2,y2);
    TH1F *hfake = new TH1F("hfake","hfake", (int)x2+6, -3, x2+3);
    //draw labels on x axis
    int counter = runlabel.size()-1;
    for(Int_t k=4;k<=hfake->GetNbinsX()-3;k+=3){
      if(stol(runlabel[counter])==refrun){
        k-=3;
        counter--;
        continue;
      }
      hfake->GetXaxis()->SetBinLabel(k, Form("run%s",runlabel[counter].c_str()));
      counter--;
    }

    hfake->Draw();
    //canvas->SetLogy();
    hfake->SetTitle(Form("Layer-%s - %s%06ld compared to all",nLayers>1 ? to_string(ilay).c_str():laynums[0].c_str(), filepath.find("run")==string::npos? "":"run",refrun));
    ge_nref[ilay]->Draw("P E2 same");
    ge_ncom1[ilay]->Draw("E2 same");
    ge_ncom2[ilay]->Draw("E2 same");
    ge_n2[ilay]->Draw("E2 same");
    hfake->GetYaxis()->SetRangeUser(min[ilay]+0.1*min[ilay], max[ilay]+0.1*max[ilay]);
    //hfake->GetYaxis()->SetLabelColor(kWhite);
    hfake->GetYaxis()->SetTickLength(0.005);
    hfake->GetYaxis()->SetMaxDigits(4);
    TLine *lineref = new TLine(-0.5, 0, x2+0.5, 0);
    lineref->SetLineColor(kGray-1);
    lineref->SetLineStyle(2);
    lineref->Draw("same");

    //draw legend
    leg->Draw("same");

    canvas->SaveAs(Form("../Plots/Layer%s_DeadPixComparison_%s%ld_compared_to_run_%s.pdf", nLayers>1 ? to_string(ilay).c_str() : laynums[0].c_str(),filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    canvas->SaveAs(Form("../Plots/Layer%s_DeadPixComparison_%s%ld_compared_to_run_%s.root", nLayers>1 ? to_string(ilay).c_str() : laynums[0].c_str(),filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

    delete canvas;
    delete hfake;
    delete lineref;
  }//end loop on layers

  //Draw plot for each stave
  //Draw
  /*cnt=0;
  cout<<"noisypix size: "<<noisypix.size()<<endl;
  while(cnt<(int)noisypix.size()){
  //for(int ilay=0; ilay<nLayers; ilay++){
    //for(int is=0; is<nstaves; is++){
      int isidx = (int)noisypix[cnt][3];//stave number
      int ilay = (int)noisypix[cnt][4];//layer number
      cout<<"cnt: "<<cnt<<"  ilay: "<<ilay<<"  istave: "<<isidx<<endl;
      TCanvas *canvas = new TCanvas(Form("mycanvas_%d_%d",ilay,isidx), Form("mycanvas_%d_%d",ilay,isidx), 1300, 800);
      canvas->SetMargin(0.08, 0.1271, 0.1759, 0.0996);
      canvas->cd();

      //fake histo (just for the axes)
      double x2,y2;
      ge_ncom2_stave[ilay][isidx]->GetPoint(ge_ncom2_stave[ilay][isidx]->GetN()-1, x2,y2);
      TH1F *hfake = new TH1F("hfake","hfake", (int)x2+6, -3, x2+3);
      //draw labels on x axis
      int counter = 0;
      for(Int_t k=4;k<=hfake->GetNbinsX()-3;k+=3){
        hfake->GetXaxis()->SetBinLabel(k, Form("run%s",runlabel[counter].c_str()));
        counter++;
      }
      hfake->Draw();
      //canvas->SetLogy();
      hfake->SetTitle(Form("Layer-%s - Stave-%d - %s%06ld compared to all",nLayers>1? to_string(ilay).c_str():laynums[0].c_str(), isidx, filepath.find("run")==string::npos? "":"run",refrun));
      ge_nref_stave[ilay][isidx]->Draw("P E2 same");
      ge_ncom1_stave[ilay][isidx]->Draw("E2 same");
      ge_ncom2_stave[ilay][isidx]->Draw("E2 same");
      ge_n2_stave[ilay][isidx]->Draw("E2 same");
      hfake->GetYaxis()->SetRangeUser(mins[ilay][isidx]+0.1*mins[ilay][isidx], maxs[ilay][isidx]+0.1*maxs[ilay][isidx]);
      //hfake->GetYaxis()->SetLabelColor(kWhite);
      hfake->GetYaxis()->SetTickLength(0.005);
      hfake->GetYaxis()->SetMaxDigits(4);
      TLine *lineref = new TLine(-0.5, 0, x2+0.5, 0);
      lineref->SetLineColor(kGray-1);
      lineref->SetLineStyle(2);
      lineref->Draw("same");

      leg->Draw("same");

      if(!cnt) canvas->SaveAs(Form("../Plots/%sAllStaves_DeadPixComparison_%s%ld_compared_to_run_%s.pdf[",nLayers==1 ? Form("Layer%s_",laynums[0].c_str()) : "AllLayers_",filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      canvas->SaveAs(Form("../Plots/%sAllStaves_DeadPixComparison_%s%ld_compared_to_run_%s.pdf",nLayers==1 ? Form("Layer%s_",laynums[0].c_str()) : "AllLayers_",filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      if(cnt==((int)noisypix.size())-(nRuns-1)) canvas->SaveAs(Form("../Plots/%sAllStaves_DeadPixComparison_%s%ld_compared_to_run_%s.pdf]",nLayers==1 ? Form("Layer%s_",laynums[0].c_str()) : "AllLayers_",filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

      delete canvas;
      delete hfake;
      delete lineref;

      cnt+=(nRuns-1);
    //}
  }*/

}

//
// Function to compare two hitmaps --> returns an arrays with timestamp of run2, noisyPixInRefRun, noisyPixInRun2, noisyPixInCommon
//
std::array<long int,5> CompareTwoRuns(THnSparse *href, THnSparse *h2){

  std::array<long int,5> noisypix = {0, 0, 0, 0, 0};
  //number of noisy pix in refrun_only and in common
  for(int iybin=1; iybin<=512; iybin++){
    href->GetAxis(1)->SetRange(iybin,iybin); //select a row
    h2->GetAxis(1)->SetRange(iybin,iybin);//select a row
    TH1D *hrefproj = (TH1D*)href->Projection(0);
    TH1D *h2proj = (TH1D*)h2->Projection(0);
    for(int ixbin=1; ixbin<=9216; ixbin++){

      if(hrefproj->GetBinContent(ixbin)==1 && h2proj->GetBinContent(ixbin)==1){//dead in both runs
        noisypix[2]++;
      }
      else if(hrefproj->GetBinContent(ixbin)==1 && h2proj->GetBinContent(ixbin)==0){//dead only in ref run
        noisypix[0]++;
      }
      else if(hrefproj->GetBinContent(ixbin)==0 && h2proj->GetBinContent(ixbin)==1){//dead only in second run
        noisypix[1]++;
      }
      else continue;
    }
    delete hrefproj;
    delete h2proj;
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
