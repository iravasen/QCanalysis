#include <string>
#include <iostream>
#include <vector>
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

using namespace std;

void SetStyle(TGraph *h, Int_t col, Style_t mkr);
void DoAnalysis(string filepath, const int nChips, string skipruns, bool isIB);

//
// MAIN
//
void AnalyzeLayerOccupancy(){
  string fpath;
  int nchips=9;
  cout<<"\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*FHRMAPS_HITMAPS* -Art | tail -n 500");
  cout<<"\nCopy file name: ";
  cin>>fpath;
  cout<<endl;

  bool isIB;

  if(fpath.find("IB")!=string::npos){
    isIB = kTRUE;
  }
  else if (fpath.find("OB")!=string::npos){
    isIB = kFALSE;
  }
  else{
    string layernum = fpath.substr(fpath.find("Layer")+5, 1);
    if(stoi(layernum)>=0 && stoi(layernum)<=2) nchips = 9; 
    else if (stoi(layernum)==3 || stoi(layernum)==4) nchips = 8; //number of HICs 
    else nchips = 14; //number of HICs
    if(nchips==9) isIB=kTRUE;
    else isIB=kFALSE;
  }

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


  //Call
  DoAnalysis(fpath, nchips, skipruns, isIB);

}

//
//Set Style
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
// Analyse data
//
void DoAnalysis(string filepath, int nChips, string skipruns, bool isIB){

  gStyle->SetOptStat(0000);

  std::vector<TH2*> hmaps;
  std::vector<string> timestamps, runnumbers, laynums;
  int nTimes=0, nRuns=1;
  Int_t col[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#0033a1"), TColor::GetColor("#708090"), TColor::GetColor("#9400D3")};

  //Read the file and the list of plots with entries
  TFile *infile=new TFile(filepath.c_str());
  TList *list = infile->GetListOfKeys();
  TKey *key;
  TObject *obj;
  TIter next(list);
  TH2 *h2;
  while((key = ((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
	&& (!obj->InheritsFrom("TH2"))
	&& (!obj->InheritsFrom("TH1"))
	) {
      cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
    }
    string objname = (string)obj->GetName();
    if(objname.find("Stv")!=string::npos) break;
    h2 = (TH2*)obj;
    if(!h2->GetEntries()) continue;
    cout<<"... Reading "<<obj->GetName()<<endl;
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);
    if(skipruns.find(runnum)!=string::npos) continue;// eventually skip runs if specified
    hmaps.push_back(h2);
    //cout<<runnum<<"  "<<timestamp<<endl;
    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);

    laynums.push_back(laynum);
    //cout<<"run: "<<runnum<<"   timestamp: "<<timestamp<<"    laynum: "<<laynum<<endl;
    nTimes++;
    if(nTimes>1 && laynum==laynums[laynums.size()-2])
      nRuns++;
    else nRuns=1;
  }

  const int nLayersIB = (int)hmaps.size()==nRuns ? 1 : stoi(laynums[laynums.size()-1])+1;
  const int nLayersOB = (int)hmaps.size()==nRuns ? 1 : stoi(laynums[laynums.size()-1])+1-3;
  int nLayers;
  if (isIB) nLayers = nLayersIB;
  else nLayers=nLayersOB;

  //const int nRuns = (int)runnumbers.size() / nLayers;
  //  cout<<"Lay: "<<nLayers<<"  Runs: "<<nRuns<<endl;
  TGraph *trend[nLayers][100];
  int ilayer=nLayers-1;
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){
    for(int ibiny=1; ibiny<=hmaps[ihist]->GetNbinsY(); ibiny++){
      trend[ilayer][ibiny-1] = new TGraph();
    }
    if(ihist>0)
      if(laynums[ihist-1]!=laynums[ihist]){
        ilayer--;
      }
  }
  ilayer=nLayers-1;
  TH1F *hproj = new TH1F();
  string histname = hmaps[0]->GetName();
  int irun=0;

  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){// start from the last in order to have the runs from the oldest to the newest
    if (!isIB){ //nChips depend on the layer for the OB
      if (stoi(laynums[ihist]) ==3 || stoi(laynums[ihist]) ==4)    nChips=8;
      else if (stoi(laynums[ihist]) ==5 || stoi(laynums[ihist]) ==6)    nChips=14;
    }

    for(int ibiny=1; ibiny<=hmaps[ihist]->GetNbinsY(); ibiny++){//loop on y bins (stave)s
      TH1D *hproj = hmaps[ihist]->ProjectionX("proj",ibiny,ibiny); //single stave
      trend[ilayer][ibiny-1]->SetName(Form("gr_L%s_stave%d",laynums[ihist].c_str(),ibiny-1));
      int deadchips = 0;
      for(int ibinx=1; ibinx<=hmaps[ihist]->GetNbinsX(); ibinx++){//evaluate the number of disabled chips
        if(hmaps[ihist]->GetBinContent(ibinx,ibiny)<1e-15)
          deadchips++;
      }
      if(deadchips>0){
	if (isIB) {
	  cout<<"\nLayer "<<laynums[ihist]<<" Stave "<<ibiny-1<<" Run: "<<runnumbers[ihist]<<" Chips number: " << nChips<<" --> Chips active:"<<nChips-deadchips<<endl;
	} else {
	  cout<<"\nLayer "<<laynums[ihist]<<" Stave "<<ibiny-1<<" Run: "<<runnumbers[ihist]<<" Hic number: " << nChips << " --> HICs active:"<<nChips-deadchips<<endl;
	}
      }
      if(deadchips!=nChips)
        trend[ilayer][ibiny-1]->SetPoint(irun, irun, hproj->Integral()/(nChips-deadchips));
      else
        trend[ilayer][ibiny-1]->SetPoint(irun, irun, 0.);

      double x=0;
      double y=0;

      if(isIB){
	if((ibiny-1)<hmaps[ihist]->GetNbinsY()/2)
	  SetStyle(trend[ilayer][ibiny-1], col[ibiny-1], 24);
	else
	  SetStyle(trend[ilayer][ibiny-1], col[ibiny-1-hmaps[ihist]->GetNbinsY()/2], 26);
      }
      else{
	if (stoi(laynums[ihist]) ==3 || stoi(laynums[ihist])==4){
	  if((ibiny-1)<hmaps[ihist]->GetNbinsY()/3){
	    SetStyle(trend[ilayer][ibiny-1], col[ibiny-1], 24);
	  }
	  else if ((ibiny-1)<hmaps[ihist]->GetNbinsY()*2/3){
	    SetStyle(trend[ilayer][ibiny-1], col[ibiny-1-hmaps[ihist]->GetNbinsY()/3], 26);
	  }
	  else{
	    SetStyle(trend[ilayer][ibiny-1], col[ibiny-1-hmaps[ihist]->GetNbinsY()*2/3], 25);	    
	  }
	}
	else if (stoi(laynums[ihist]) ==5 || stoi(laynums[ihist])==6){
	  if((ibiny-1)<int(hmaps[ihist]->GetNbinsY()/4))
	    SetStyle(trend[ilayer][ibiny-1], col[ibiny-1], 24);
	  else if ((ibiny-1)<2*int(hmaps[ihist]->GetNbinsY()/4))
	    SetStyle(trend[ilayer][ibiny-1], col[ibiny-1-int(hmaps[ihist]->GetNbinsY()/4)], 26);
	  else if ((ibiny-1)<3*hmaps[ihist]->GetNbinsY()/4)
	    SetStyle(trend[ilayer][ibiny-1], col[ibiny-1-2*int(hmaps[ihist]->GetNbinsY()/4)], 25);
	  else
	    SetStyle(trend[ilayer][ibiny-1], col[ibiny-1-int(hmaps[ihist]->GetNbinsY()*3/4)], 30);
	}

      }
    }
    irun++;
    if(ihist>0)
      if(laynums[ihist-1]!=laynums[ihist]){
        irun=0;
        ilayer--;
      }
  }

  int npoints = trend[0][0]->GetN();
  TH1F *hfake = new TH1F("hfake", "; Run; Fake-hit Rate (/event/pixel)", npoints, -0.5, (double)npoints-0.5);

  for(int ir=0; ir<(int)runnumbers.size()/nLayers; ir++)
    hfake->GetXaxis()->SetBinLabel(ir+1, Form("run%06d", stoi(runnumbers[(int)runnumbers.size()/nLayers-1-ir])));


  //Draw
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas *canvas = new TCanvas();
    canvas->cd();
    canvas->SetLogy();
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetMargin(0.0988,0.15,0.194,0.0993);
    TLegend *leg = new TLegend(0.857, 0.197,0.997,0.898);

    //setting a two column legend for L5 and L6
    if (!isIB && nLayers==4 && (ilay==1 || ilay==2 || ilay==3)) leg->SetNColumns(2);
    else  if (stoi(laynums[0]) ==4 || stoi(laynums[0]) ==5 || stoi(laynums[0]) ==6) leg->SetNColumns(2);
    
    for(int istave=0; istave<hmaps[ilay*nRuns]->GetNbinsY(); istave++)
      leg->AddEntry(trend[ilay][istave], Form("Stv%d",istave), "p");
    hfake->GetYaxis()->SetRangeUser(1e-14, 1e-3);
    hfake->GetXaxis()->SetTitleOffset(2.8);
    hfake->SetTitle(Form("Layer-%s, %s",laynums[ilay*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    hfake->Draw();
    for(int istave=0; istave<hmaps[ilay*nRuns]->GetNbinsY(); istave++)
      trend[ilay][istave]->Draw("P same");
    leg->Draw("same");
    canvas->SaveAs(Form("../Plots/Layer%s_fakehitrate_%s.pdf", laynums[ilay*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    canvas->SaveAs(Form("../Plots/Layer%s_fakehitrate_%s.root", laynums[ilay*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    delete canvas;
    delete leg;
  }

  //Make GIF with TH2 for each run and for each layer
  irun=1;
  ilayer=nLayers-1;
  gStyle->SetPalette(1);
  for(int ilay=0; ilay<nLayers; ilay++)//remove images if they exist already
    gSystem->Unlink(Form("../Plots/Layer%s_fakehitratemap_%s.gif", laynums[ilay*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){// start from the last in order to have the runs from the oldest to the newest
    TCanvas *canvas = new TCanvas();
    canvas->cd();
    canvas->SetLogz();
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.15);
    hmaps[ihist]->Draw("colz");
    hmaps[ihist]->SetMinimum(1e-14);
    hmaps[ihist]->SetMaximum(1e-3);
    hmaps[ihist]->GetZaxis()->SetTitle("Fake-hit Rate (/event/pixel)");
    hmaps[ihist]->SetTitle(Form("Layer-%s, Run %06d (%d/%d)", laynums[ilayer*nRuns].c_str(), stoi(runnumbers[ihist]), irun,nRuns));
    canvas->Print(Form("../Plots/Layer%s_fakehitratemap_%s.gif+40", laynums[ilayer*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    irun++;
    if(!ihist && nLayers==1){
      canvas->Print(Form("../Plots/Layer%s_fakehitratemap_%s.gif++40++", laynums[ilayer*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      break;
    }
    if(nLayers>1){
      if(ihist>0){
        if(laynums[ihist-1]!=laynums[ihist]){
          canvas->Print(Form("../Plots/Layer%s_fakehitratemap_%s.gif++40++", laynums[ilayer*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
          ilayer--;
          irun=1;
        }
      }
      else if(!ihist && !ilayer){
        canvas->Print(Form("../Plots/Layer%s_fakehitratemap_%s.gif++40++", laynums[ilayer*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      }
    }

    delete canvas;
  }

}
