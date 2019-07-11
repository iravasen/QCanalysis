#include <string>
#include <iostream>
#include <vector>
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
#include <TSystem.h>

using namespace std;

void SetStyle(TH1 *h, Int_t col);
void DoAnalysis(string filepath, const int nChips, bool isIB, int layernum, int stavenum);

//
// MAIN
//
void AnalyzeStaveHitmaps(){
  string fpath;
  int nchips=9;
  cout<<"\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n"<<endl;
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

  //Call
  DoAnalysis("../Data/"+fpath, nchips, isIB, stoi(layernum), stoi(stavenum));

}

//
//Set Style
//
void SetStyle(TH1 *h, Int_t col){
  h->SetLineColor(col);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(1.4);
  h->SetMarkerColor(col);
  //h->SetFillStyle(0);
  //h->SetFillColorAlpha(col,0.8);
}

//
// Analyse data
//
void DoAnalysis(string filepath, const int nChips, bool isIB, int layernum, int stavenum){

  gStyle->SetOptStat(0000);

  std::vector<TH2S*> hmaps;
  std::vector<string> timestamps, runnumbers;
  int nTimes=0;
  Int_t col[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0")};

  //Read the file and the list of plots with entries
  TFile *infile=new TFile(filepath.c_str());
  TList *list = infile->GetListOfKeys();
  list->ls();
  TIter next(list);
  TKey *key;
  TObject *obj;
  TH2S *h2;
  while((key = ((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1"))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    h2 = (TH2S*)obj->Clone(obj->GetName());
    if(!h2->GetEntries()) continue;
    cout<<"... Reading "<<obj->GetName()<<endl;
    hmaps.push_back(h2);
    string objname = (string)obj->GetName();
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_")+1, 13) : objname.substr(objname.find("_",3)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("_",1)+1, objname.find("_",3)-3);
    //cout<<runnum<<"  "<<timestamp<<endl;
    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);
    nTimes++;
  }

  TH1F *hnoiseocc[nChips];
  TH1F *hproj = new TH1F();
  string histname = hmaps[0]->GetName();
  for(int itime=0; itime<(int)hmaps.size(); itime++){
    for(int ichip=0; ichip<nChips; ichip++){
      hproj = (TH1F*)hmaps[itime]->ProjectionY(Form("hproj_%d_%d",ichip,itime), 1+256*ichip, 256+256*ichip);
      double ntrigger = std::stol(timestamps[itime])>1560426360000 ? 300.*50000. : 600.*50000.; //after 13 June 2019 13.46, the trigger is 50kHz for 5 min, before was 50kHz for 10 min
      //cout<<ntrigger<<endl;
      double noiseocc = (double)(hproj->GetEntries()) / (double)(512.*1024.*ntrigger);
      //cout<<"entries: "<<(double)(hproj->GetEntries())<<"  noise occ: "<<noiseocc<<endl;
      if(!itime) hnoiseocc[ichip] = new TH1F(Form("hnoiseocc_chip%d",ichip), Form("hnoiseocc_chip%d; %s; Fake Hit Rate (/pixel/event)",ichip,histname.find("run")==string::npos? "Time (timestamp)":"Run"), nTimes, 0.5, (double)nTimes+0.5);
      hnoiseocc[ichip]->SetBinContent(nTimes-itime, noiseocc);
      double errrel_num = TMath::Sqrt(hproj->GetEntries()) / hproj->GetEntries();//poissonian error
      double errrel_den = 50000. / ntrigger; // because error or number of triggers is +- 1s*50kHz
      hnoiseocc[ichip]->SetBinError(nTimes-itime, noiseocc*TMath::Sqrt((errrel_num*errrel_num+errrel_den*errrel_den)));//nTimes-itime because the oldest first
      SetStyle(hnoiseocc[ichip], col[ichip]);
      hproj->Reset();
    }
  }

  //Set labels
  for(int ichip=0; ichip<nChips; ichip++)
    for(int ibin=1; ibin<=hnoiseocc[ichip]->GetNbinsX(); ibin++)
      hnoiseocc[ichip]->GetXaxis()->SetBinLabel(nTimes+1-ibin, histname.find("run")==string::npos? timestamps[ibin-1].c_str():runnumbers[ibin-1].c_str());


  //Draw
  TCanvas *canvas = new TCanvas();
  canvas->cd();
  canvas->SetLogy();
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->SetMargin(0.0988,0.1,0.194,0.0993);
  TLegend *leg = new TLegend(0.904, 0.197,0.997,0.898);
  for(int ichip=0; ichip<nChips; ichip++)
    leg->AddEntry(hnoiseocc[ichip], Form("Ch%d",ichip));
  hnoiseocc[0]->GetYaxis()->SetRangeUser(1e-12, 1e-3);
  hnoiseocc[0]->GetXaxis()->SetTitleOffset(2.8);
  hnoiseocc[0]->SetTitle(Form("%s Layer-%d Stave-%d, %s",isIB?"IB":"OB",layernum,stavenum, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  hnoiseocc[0]->Draw("P E1 X0");
  hnoiseocc[0]->Draw("same Lhist");
  for(int ichip=1; ichip<nChips; ichip++){
    hnoiseocc[ichip]->Draw("P E1 X0 same");
    hnoiseocc[ichip]->Draw("same Lhist");
  }
  leg->Draw();
  canvas->SaveAs(Form("../Plots/%s_layer%d_stave%d_fakehitrate_%s.pdf", isIB?"IB":"OB", layernum, stavenum, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  canvas->SaveAs(Form("../Plots/%s_layer%d_stave%d_fakehitrate_%s.root", isIB?"IB":"OB", layernum, stavenum, filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));

}
