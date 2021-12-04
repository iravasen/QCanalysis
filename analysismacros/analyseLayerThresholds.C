#include "itsAnalysis.hh"

void SetStyle(TGraph *h, Int_t col, Style_t mkr){
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.4);
  h->SetMarkerColor(col);
}

Int_t col[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#0033a1")};

// Main function
void analyseLayerThresholds(){
  itsAnalysis myAnalysis("Threshold");

  auto nLayers      = myAnalysis.nLayers();
  auto laynums      = myAnalysis.Layers();
  auto runNumbers   = myAnalysis.Runs();
  auto hmaps        = myAnalysis.loadedHists();

  TGraph *trend[nLayers][100]; //list of trends
  for (int ilayer = 0; ilayer < nLayers; ++ilayer){ //loop over layers
    for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){ //loop over number of histograms
      for(int ibiny=1; ibiny<=hmaps[ihist]->GetNbinsY(); ibiny++){ // Loop over y-bins (staves)
        trend[ilayer][ibiny-1] = new TGraph();
      }
    }
  }
  

  TH1F *hproj = new TH1F();
  for (string layer : laynums){ // loop over layers
    int irun = myAnalysis.nRuns()-1; // Count back to have right order in histogram
    for(auto hist : myAnalysis.loadLayer(stoi(layer))){ //loop hist for layer
      //cout<<myAnalysis.loadLayer(stoi(layer)).size()<<endl;

      for(int ibiny=1; ibiny<=hist->GetNbinsY(); ibiny++){//loop on y bins (staves)
        TH1D *hproj = hist->ProjectionX("proj",ibiny,ibiny); //single stave
        trend[stoi(layer)][ibiny-1]->SetName(Form("gr_L%s_stave%d",layer.c_str(),ibiny-1));

        int nChips = myAnalysis.nChips(stoi(layer));
        int deadchips = 0;
        for(int ibinx=1; ibinx<=hist->GetNbinsX(); ibinx++){//evaluate the number of disabled chips
          if(hist->GetBinContent(ibinx,ibiny)<1e-15) deadchips++;
        }
        if(deadchips>0)
          cout<<"Layer: "<<layer<<" Stave: "<<ibiny-1<<" Run: "<<myAnalysis.getRunNumber(hist)<<" --> Chips active:"<<nChips-deadchips<<endl;
        if(deadchips!=nChips)
          trend[stoi(layer)][ibiny-1]->SetPoint(irun, irun, hproj->Integral()/(nChips-deadchips)); // average per chip
        else
          trend[stoi(layer)][ibiny-1]->SetPoint(irun, irun, 0.);
        if((ibiny-1)<hist->GetNbinsY()/2)
          SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1], 24);
        else
          SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1-hist->GetNbinsY()/2], 26);
        }
      irun--;
      }
    }
  
  
  int npoints = myAnalysis.nRuns();
  TH1F *hfake = new TH1F("hfake", "; Run; Avg. Threshold (DAC)", npoints, -0.5, (double)npoints-0.5);
  for(int ir=0; ir<npoints; ir++) // Set bin labels to run numbers
      hfake->GetXaxis()->SetBinLabel(npoints-(ir), Form("run%06d",stoi(myAnalysis.Runs()[ir])));
  
  for (string layer : laynums){ // loop over layers
    TCanvas *canvas = new TCanvas();
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetMargin(0.0988,0.1,0.194,0.0993);
    TLegend *leg = new TLegend(0.904, 0.197,0.997,0.898);
    hfake->GetYaxis()->SetRangeUser(8.5, 12.5);
    hfake->GetXaxis()->SetTitleOffset(2.8);
    hfake->SetStats(0);
    hfake->SetTitle(Form("Layer-%s, from run%s to run%s",layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
    //hfake->SetTitle(Form("Layer-%s",layer.c_str()));
    hfake->Draw();
    for(int istave=0; istave<myAnalysis.stavesInLayer(stoi(layer)); istave++){
      leg->AddEntry(trend[stoi(layer)][istave], Form("Stv%d",istave), "p");
      trend[stoi(layer)][istave]->Draw("P same");
    }
    leg->Draw("same");
    canvas->SaveAs(Form("../Plots/Layer%s_thresholds_run%s-run%s.pdf", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
    canvas->SaveAs(Form("../Plots/Layer%s_thresholds_run%s-run%s.root", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
    delete canvas;
    delete leg;

  }

} // end of analyseLayerThresholds()