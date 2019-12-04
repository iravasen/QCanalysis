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

using namespace std;

//Functions
std::array<long int,5> CompareTwoRuns(TH2 *href, TH2 *h2);
void SetStyle(TGraphErrors *ge, Color_t col);
void DoAnalysis(string filepath, const int nChips, bool isIB, long int refrun);

void CompareNoisyPixelsInRuns(){
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

  cout<<"Available runs in your file:\n"<<endl;
  TFile *infile=new TFile(fpath.c_str());
  TList *list = (TList*)infile->GetListOfKeys();
  TIter next(list);
  TObject *obj;
  TKey *key;
  TH2 *h2;
  vector<string> stavenums;
  while((key = ((TKey*)next()))){
    obj=key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1"))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("Stv")==string::npos) continue;
    h2 = (TH2*)obj->Clone(obj->GetName());
    //if(!h2->GetEntries()) continue;
    string runnum =  objname.substr(objname.find("run")+3, 6);
    string stvnum = objname.substr(objname.find("Stv")+3,2);
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

  DoAnalysis(fpath, nchips, isIB, refrun);
}

//
// Analysis
//
void DoAnalysis(string filepath, const int nChips, bool isIB, long int refrun){

  gStyle->SetOptStat(0000);

  std::vector<TH2*> hmaps;
  std::vector<string> timestamps, runnumbers, stavenums, laynums;
  vector<int> posrefrun;
  int nLayers=1, nRuns=1;
  vector<int> nStavesInLay;
  int cstv=0;

  //Read the file and the list of plots with entries
  TFile *infile=new TFile(filepath.c_str());
  TList *list = (TList*)infile->GetListOfKeys();
  TIter next(list);
  TObject *obj;
  TKey *key;
  TH2 *h2;
  while((key=((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1"))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("Stv")==string::npos) continue;
    h2 = (TH2*)obj->Clone(obj->GetName());
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
    if(stol(runnum)==refrun) //position of refence run
      posrefrun.push_back((int)hmaps.size()-1);

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
  }
  nStavesInLay.push_back(cstv+1);//in case of 1 layer or for last layer

  /*cout<<"Staves in Lay: "<<nStavesInLay[0]<<endl;
  cout<<"nRuns: "<<nRuns<<endl;
  cout<<"nLayers: "<<nLayers<<endl;
  cout<<"SizeOf posrefrun: "<<posrefrun.size()<<endl;*/

  //Compare all the runs (non-empty ones) with the reference run chosen by the user
  vector<string> runlabel;
  int istave = 0;
  for(int ilay=0; ilay<(int)nStavesInLay.size(); ilay++)
    istave+=nStavesInLay[ilay];
  istave--;
  int ilayer = nLayers-1;
  long int first[nLayers][nRuns], second[nLayers][nRuns], both[nLayers][nRuns];
  int irun=0;
  vector<array<long int,5>> noisypix;
  for(int ilay=0; ilay<nLayers; ilay++)
    for(int i=0; i<nRuns; i++){
      first[ilay][i]=0; second[ilay][i]=0; both[ilay][i]=0;
    }
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){ //start from the bottom in order to start with the oldest run
    if(runnumbers[ihist].find(std::to_string(refrun))!=string::npos){
      if(ihist>0){// in case ref run is the first into the list of runs
        if(stavenums[ihist-1]!=stavenums[ihist]){
          istave--;
          irun=0;
        }
        if(laynums[ihist-1]!=laynums[ihist]){
          ilayer--;
        }
      }
      continue;
    }
    //if(!hmaps[ihist]->GetEntries()) continue; // do not compare ref run with empty run (= empty maps)

    noisypix.push_back(CompareTwoRuns(hmaps[posrefrun[istave]], hmaps[ihist]));
    noisypix[noisypix.size()-1][3] = stol(stavenums[ihist]);
    noisypix[noisypix.size()-1][4] = stol(laynums[ihist]);
    first[ilayer][irun]+=noisypix[noisypix.size()-1][0];
    second[ilayer][irun]+=noisypix[noisypix.size()-1][1];
    both[ilayer][irun]+=noisypix[noisypix.size()-1][2];
    irun++;

    //cout<<noisypix[noisypix.size()-1][0]<<"  "<<noisypix[noisypix.size()-1][1]<<"  "<<noisypix[noisypix.size()-1][2]<<"  "<<noisypix[noisypix.size()-1][3]<<endl;
    if(ilayer==nLayers-1 && stavenums[ihist]=="0")
      runlabel.push_back(runnumbers[ihist]);

    if(ihist>0){
      if(stavenums[ihist-1]!=stavenums[ihist]){
        istave--;
        irun=0;
      }
      if(laynums[ihist-1]!=laynums[ihist]){
        ilayer--;
      }
    }

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
  double max = -1;
  double min = 1e35;
  for(int ilay=0; ilay<nLayers; ilay++){
    ge_nref[ilay] = new TGraphErrors();
    ge_n2[ilay] = new TGraphErrors();
    ge_ncom1[ilay] = new TGraphErrors();
    ge_ncom2[ilay] = new TGraphErrors();

    for(int ir=0; ir<nRuns-1; ir++){//first the older data and last the most recent
      //first couple of bar on the left
      //int ipoint = (int)noisypix.size()-icomp-1;
      if(!ir) xshift=1.;
      else xshift = 3.;
      ge_nref[ilay]->SetPoint(ir, ir*xshift, (double)both[ilay][ir]/2.+(double)first[ilay][ir]/2.);
      ge_nref[ilay]->SetPointError(ir, 0.5, (double)first[ilay][ir]/2.);
      ge_ncom1[ilay]->SetPoint(ir, ir*xshift, 0.);
      ge_ncom1[ilay]->SetPointError(ir, 0.5, (double)both[ilay][ir]/2.);
      if((double)both[ilay][ir]/2.+(double)first[ilay][ir] > max) max = (double)both[ilay][ir]/2.+(double)first[ilay][ir];

      //second couple of bar on the right
      ge_n2[ilay]->SetPoint(ir, ir*xshift+1, -(double)both[ilay][ir]/2.-(double)second[ilay][ir]/2.);
      ge_n2[ilay]->SetPointError(ir, 0.5, (double)second[ilay][ir]/2.);
      ge_ncom2[ilay]->SetPoint(ir, ir*xshift+1, 0.);
      ge_ncom2[ilay]->SetPointError(ir, 0.5, (double)both[ilay][ir]/2.);
      if(-(double)both[ilay][ir]/2.-(double)second[ilay][ir] < min) min = -(double)both[ilay][ir]/2.-(double)second[ilay][ir];
    }//end first loop on runs

    //Style
    SetStyle(ge_nref[ilay], kBlue);
    SetStyle(ge_ncom1[ilay], kBlack);
    SetStyle(ge_ncom2[ilay], kBlack);
    SetStyle(ge_n2[ilay], kRed+2);
  }//end loop on layers

  //FILL PLOTS FOR EACH STAVE
  //fill the plots for each stave
  int cnt = 0;
  double maxs[nLayers][100];
  double mins[nLayers][100];
  for(int ilay=nLayers-1; ilay>=0; ilay--){
    for(int is=nStavesInLay[ilay]-1; is>=0; is--){
      maxs[ilay][is] = -1.;
      mins[ilay][is] = 1e35;
      ge_nref_stave[ilay][is] = new TGraphErrors();
      ge_ncom1_stave[ilay][is] = new TGraphErrors();
      ge_n2_stave[ilay][is] = new TGraphErrors();
      ge_ncom2_stave[ilay][is] = new TGraphErrors();
      for(int ir=0; ir<nRuns-1; ir++){//first the older data and last the most recent
        if(!ir) xshift=1.;
        else xshift = 3.;
        ge_nref_stave[ilay][is]->SetPoint(ir, ir*xshift, (double)noisypix[cnt][2]/2.+(double)noisypix[cnt][0]/2.);
        ge_nref_stave[ilay][is]->SetPointError(ir, 0.5, (double)noisypix[cnt][0]/2.);
        ge_ncom1_stave[ilay][is]->SetPoint(ir, ir*xshift, 0.);
        ge_ncom1_stave[ilay][is]->SetPointError(ir, 0.5, (double)noisypix[cnt][2]/2.);
        if((double)noisypix[cnt][2]/2.+(double)noisypix[cnt][0] > maxs[ilay][is]) maxs[ilay][is] = (double)noisypix[cnt][2]/2.+(double)noisypix[cnt][0];

        //second couple of bar on the right
        ge_n2_stave[ilay][is]->SetPoint(ir, ir*xshift+1, -(double)noisypix[cnt][2]/2.-(double)noisypix[cnt][1]/2.);
        ge_n2_stave[ilay][is]->SetPointError(ir, 0.5, (double)noisypix[cnt][1]/2.);
        ge_ncom2_stave[ilay][is]->SetPoint(ir, ir*xshift+1, 0.);
        ge_ncom2_stave[ilay][is]->SetPointError(ir, 0.5, (double)noisypix[cnt][2]/2.);
        if(-(double)noisypix[cnt][2]/2.-(double)noisypix[cnt][1] < mins[ilay][is]) mins[ilay][is] = -(double)noisypix[cnt][2]/2.-(double)noisypix[cnt][1];
        cnt++;
      }
      SetStyle(ge_nref_stave[ilay][is], kBlue);
      SetStyle(ge_ncom1_stave[ilay][is], kBlack);
      SetStyle(ge_ncom2_stave[ilay][is], kBlack);
      SetStyle(ge_n2_stave[ilay][is], kRed+2);
    }
  }

  //Legend
  TLegend *leg = new TLegend(0.876,0.176, 0.994, 0.902);
  leg->SetLineColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(ge_nref[0], "#splitline{#noisy pix}{ref. run only}", "f");
  leg->AddEntry(ge_n2[0], "#splitline{#noisy pix}{2nd run only}", "f");
  leg->AddEntry(ge_ncom1[0], "#splitline{#noisy pix}{both}");

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
    int counter = 0;
    for(Int_t k=4;k<=hfake->GetNbinsX()-3;k+=3){
      hfake->GetXaxis()->SetBinLabel(k, Form("run%s",runlabel[counter].c_str()));
      counter++;
    }
    hfake->Draw();
    //canvas->SetLogy();
    hfake->SetTitle(Form("Layer-%s - %s%06ld compared to all",laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(), filepath.find("run")==string::npos? "":"run",refrun));
    ge_nref[ilay]->Draw("P E2 same");
    ge_ncom1[ilay]->Draw("E2 same");
    ge_ncom2[ilay]->Draw("E2 same");
    ge_n2[ilay]->Draw("E2 same");
    hfake->GetYaxis()->SetRangeUser(min+0.1*min, max+0.1*max);
    //hfake->GetYaxis()->SetLabelColor(kWhite);
    hfake->GetYaxis()->SetTickLength(0.005);
    hfake->GetYaxis()->SetMaxDigits(4);
    TLine *lineref = new TLine(-0.5, 0, x2+0.5, 0);
    lineref->SetLineColor(kGray-1);
    lineref->SetLineStyle(2);
    lineref->Draw("same");

    //draw legend
    leg->Draw("same");

    canvas->SaveAs(Form("../Plots/Layer%s_NoisyPixComparison_%s%ld_compared_to_run_%s.pdf", laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    canvas->SaveAs(Form("../Plots/Layer%s_NoisyPixComparison_%s%ld_compared_to_run_%s.root", laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

    delete canvas;
    delete hfake;
    delete lineref;
  }//end loop on layers

  //Draw plot for each stave
  //Draw
  for(int ilay=0; ilay<nLayers; ilay++){
    for(int is=0; is<nStavesInLay[ilay]; is++){
      TCanvas *canvas = new TCanvas(Form("mycanvas_%d_%d",ilay,is), Form("mycanvas_%d_%d",ilay,is), 1300, 800);
      canvas->SetMargin(0.08, 0.1271, 0.1759, 0.0996);
      canvas->cd();

      //fake histo (just for the axes)
      double x2,y2;
      ge_ncom2_stave[ilay][is]->GetPoint(ge_ncom2_stave[ilay][is]->GetN()-1, x2,y2);
      TH1F *hfake = new TH1F("hfake","hfake", (int)x2+6, -3, x2+3);
      //draw labels on x axis
      int counter = 0;
      for(Int_t k=4;k<=hfake->GetNbinsX()-3;k+=3){
        hfake->GetXaxis()->SetBinLabel(k, Form("run%s",runlabel[counter].c_str()));
        counter++;
      }
      hfake->Draw();
      //canvas->SetLogy();
      hfake->SetTitle(Form("Layer-%s - Stave-%d - %s%06ld compared to all",laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(), is, filepath.find("run")==string::npos? "":"run",refrun));
      ge_nref_stave[ilay][is]->Draw("P E2 same");
      ge_ncom1_stave[ilay][is]->Draw("E2 same");
      ge_ncom2_stave[ilay][is]->Draw("E2 same");
      ge_n2_stave[ilay][is]->Draw("E2 same");
      hfake->GetYaxis()->SetRangeUser(mins[ilay][is]+0.1*mins[ilay][is], maxs[ilay][is]+0.1*maxs[ilay][is]);
      //hfake->GetYaxis()->SetLabelColor(kWhite);
      hfake->GetYaxis()->SetTickLength(0.005);
      hfake->GetYaxis()->SetMaxDigits(4);
      TLine *lineref = new TLine(-0.5, 0, x2+0.5, 0);
      lineref->SetLineColor(kGray-1);
      lineref->SetLineStyle(2);
      lineref->Draw("same");

      leg->Draw("same");

      if(!ilay && !is) canvas->SaveAs(Form("../Plots/%sAllStaves_NoisyPixComparison_%s%ld_compared_to_run_%s.pdf[",nLayers==1 ? Form("Layer%s_",laynums[ilay*nRuns*nStavesInLay[ilay]].c_str()) : "AllLayers_",filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      canvas->SaveAs(Form("../Plots/%sAllStaves_NoisyPixComparison_%s%ld_compared_to_run_%s.pdf",nLayers==1 ? Form("Layer%s_",laynums[ilay*nRuns*nStavesInLay[ilay]].c_str()) : "AllLayers_",filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      if(ilay==nLayers-1 && is==nStavesInLay[ilay]-1) canvas->SaveAs(Form("../Plots/%sAllStaves_NoisyPixComparison_%s%ld_compared_to_run_%s.pdf]",nLayers==1 ? Form("Layer%s_",laynums[ilay*nRuns*nStavesInLay[ilay]].c_str()) : "AllLayers_",filepath.find("run")==string::npos? "":"run",refrun, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

      delete canvas;
      delete hfake;
      delete lineref;
    }
  }

}

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
// Style
//
void SetStyle(TGraphErrors *ge, Color_t col){
  ge->SetMarkerStyle(0);
  ge->SetMarkerColor(col);
  ge->SetFillColor(col);
  ge->SetLineColor(col);
}
