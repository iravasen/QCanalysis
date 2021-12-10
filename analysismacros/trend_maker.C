
void trend_maker(){
gROOT->SetBatch(kTRUE);
string fpath;
cout<<"Available file(s) for the analysis of track task:"<<endl;
gSystem->Exec("ls ../Data/*TrackTask* -Art | tail -n 500");
cout<<"\nCopy the file name to analyse:"<<endl;
cin>>fpath;
cout<<endl;


TFile *file = TFile::Open(fpath.c_str());

   TList * keyslist=(TList *) file->GetListOfKeys();
   TIter next(keyslist);
   TKey *key;
   vector<string>runs_clusage, runs_eta, runs_phi, runs_ncl, runs_occ;
   vector<double>fcluster, etacounts, relative_eta, phicounts, phicounts1, phicounts2, relative_phi1, relative_phi2, ncluster, rms_ncl, first_bin_edge_ncl, last_bin_edge_ncl, occ;
   vector<double>cycle_clusage, cycle_eta, cycle_phi, cycle_ncl, cycle_occ;
   double i,j,k,m,n=0;

   cout<<"List of runs:"<<endl;
   for(auto&& keyAsObj : *file->GetListOfKeys()){
       auto keysub = (TKey*) keyAsObj;
       string name= keysub->GetName();
       if(name.find("ClusterUsage")!=string::npos){
        string run = name.substr(name.find("run")+3, 6);
        cout<<run<<endl;  
       }
    }

   string skipans, skipruns;
   cout<<"Would you like to skip some run(s)? [y/n]"<<endl;
   cin>>skipans;
   if(skipans=="y" || skipans=="Y"){
     cout<<endl;
     cout<<"Specify run number(s) to skip separated by comma (no white spaces!):"<<endl;
     cin>>skipruns;
     cout<<endl;
   }
   else
     skipruns=" ";

   while ((key=(TKey*)next())){

      string keyname = (string)key->GetName();
      string runnum = keyname.substr(keyname.find("run")+3, 6);
      if(skipruns.find(runnum)!=string::npos) continue;

      if (keyname.find("ClusterUsage")!=string::npos){
          i++;
         // key->Print();
          cout<<keyname<<endl;
          cout<<"Run number: "<<runnum<<endl;

          runs_clusage.push_back(runnum);
          cycle_clusage.push_back(i);

          TH1D *h= (TH1D *) key->ReadObj();

          //Extract the mean value in x axis
          double x_mean = h->GetMean(1);
          //Find the bin number for a value in x axis
          TAxis* xaxis = h->GetXaxis();
          int binx = xaxis->FindBin(x_mean);

          //Find the value in y axis
          double bin = h->GetBinContent(binx);
          //std::cout<<"Fraction of clusters= "<< bin << std::endl;
          fcluster.push_back(bin);
     }

//////////////////////////////////////////////////////////////////////////////////////////////

       if (keyname.find("EtaDistribution")!=string::npos){
          j++;
          //key->Print();
          cout<<keyname<<endl;
          cout<<"Run number: "<<runnum<<endl;

          runs_eta.push_back(runnum);
          cycle_eta.push_back(j);

          TH1D *h= (TH1D *) key->ReadObj();

          //Extract the mean value in x axis
          double x_mean = h->GetMean(1);
          
          //Integral of the distribution 
          double integral = h->Integral();
          
          //Find the bin number for a value in x axis
          TAxis* xaxis = h->GetXaxis();
          int binx = xaxis->FindBin(x_mean);

          //Find the value in y axis
          double bin = h->GetBinContent(binx);
          //std::cout<<"etaCounts_y_value= "<< bin << std::endl;
          etacounts.push_back(bin);
          
          //Relative values of eta counts
          double relative;
          if (integral != 0){
            relative = (bin/integral);
          }
          else {relative = 0;}
          
          relative_eta.push_back(relative);
     }

//////////////////////////////////////////////////////////////////////////////////////////////

       if (keyname.find("PhiDistribution")!=string::npos){
          k++;
          //key->Print();
          cout<<keyname<<endl;
          cout<<"Run number: "<<runnum<<endl;

          runs_phi.push_back(runnum);
          cycle_phi.push_back(k);

          TH1D *h= (TH1D *) key->ReadObj();
          
          //Integral in range 0-2pi
          double integral = h->Integral();
                    
          //Extract the mean value in x axis
          double x_mean = h->GetMean(1);
          //Find the bin number for a value in x axis
          TAxis* xaxis = h->GetXaxis();
          int binx = xaxis->FindBin(x_mean);
          //Find the value in y axis
          double bin = h->GetBinContent(binx);
          //cout<<"phi_y_value= "<<bin<<endl;
          phicounts.push_back(bin);

          TH1D *hc1 = (TH1D *)h->Clone();
          hc1->GetXaxis()->SetRangeUser(0, 3);
          
          //Integral in range 0-pi
          double integral1 = hc1->Integral();
          
          //Extract the mean value in x axis
          double x_mean1 = hc1->GetMean(1);
          //Find the bin number for a value in x axis
          TAxis* xaxis1 = hc1->GetXaxis();
          int binx1 = xaxis1->FindBin(x_mean1);
          //Find the value in y axis
          double bin1 = hc1->GetBinContent(binx1);
          //cout<<"y value clone 1= "<<bin1<<endl;
          phicounts1.push_back(bin1);

          TH1D *hc2 = (TH1D *)h->Clone();
          hc2->GetXaxis()->SetRangeUser(3, 6);
          
          //Integral in range pi-2pi
          double integral2 = hc2->Integral();
          
          //Extract the mean value in x axis
          double x_mean2 = hc2->GetMean(1);
          //Find the bin number for a value in x axis
          TAxis* xaxis2 = hc2->GetXaxis();
          int binx2 = xaxis2->FindBin(x_mean2);
          //Find the value in y axis
          double bin2 = hc2->GetBinContent(binx2);
          //cout<<"y value clone 2= "<<bin2<<endl;
          phicounts2.push_back(bin2);
          
          //Relative values
          double rel_phi1;
          if (integral != 0){
            rel_phi1 = (integral1/integral);
          }
          else {rel_phi1 = 0;}
          
          relative_phi1.push_back(rel_phi1);
          
          double rel_phi2;
          if (integral != 0){
            rel_phi2 = (integral2/integral);
          }
          else {rel_phi2 = 0;}
          
          relative_phi2.push_back(rel_phi2);
      }

//////////////////////////////////////////////////////////////////////////////////////////////

       if (keyname.find("NClusters")!=string::npos){
          m++;
          //key->Print();
          cout<<keyname<<endl;
          cout<<"Run number: "<<runnum<<endl;

          runs_ncl.push_back(runnum);
          cycle_ncl.push_back(m);

          TH1D *h= (TH1D *) key->ReadObj();         

         //Extract the mean value in x axis
         double x_mean = h->GetMean(1);
         ncluster.push_back(x_mean);
         
         //Extract RMS of the distribution
         double rms = h->GetRMS();
         rms_ncl.push_back(rms);
         
         //bin low edge of the leftmost bin filled
         double first_edge;
         if ((h->FindFirstBinAbove())>0){
           first_edge = h->GetXaxis()->GetBinLowEdge(h->FindFirstBinAbove());         
         }
         else{ 
           first_edge =0;
         }           
         first_bin_edge_ncl.push_back(first_edge);       
         
         //bin low edge of the rightmost bin filled
         double last_edge;
         if ((h->FindLastBinAbove())>0){
           last_edge = h->GetXaxis()->GetBinLowEdge(h->FindLastBinAbove());
         }
         else{
           last_edge = 0;
         }         
         last_bin_edge_ncl.push_back(last_edge);

     }

//////////////////////////////////////////////////////////////////////////////////////////////

       if (keyname.find("OccupancyROF")!=string::npos){
          n++;
          //key->Print();
          cout<<keyname<<endl;
          cout<<"Run number: "<<runnum<<endl;

          runs_occ.push_back(runnum);
          cycle_occ.push_back(n);

          TH1D *h= (TH1D *) key->ReadObj();

          //Extract the mean value in x axis
          double x_mean = h->GetMean(1);

          //Find the bin number for a value in x axis
          TAxis* xaxis = h->GetXaxis();
          int binx = xaxis->FindBin(x_mean);

          //Find the value in y axis
          double bin = h->GetBinContent(binx);
          //std::cout<<"OccupancyROF= "<< bin << std::endl;
          occ.push_back(bin);
     }
   }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::reverse(runs_clusage.begin(),runs_clusage.end());
  std::reverse(fcluster.begin(),fcluster.end());

  auto c1 = new TCanvas("Fraction of clusters used in tracking");
  c1->SetGrid(0,1);
  auto gr = new TGraph(cycle_clusage.size(),&cycle_clusage[0],&fcluster[0]);
  double b=0;
  gr->GetXaxis()->SetNdivisions(cycle_clusage.size());
  gr->GetXaxis()->Set(1+cycle_clusage.size(),gr->GetXaxis()->GetXmin(),1+cycle_clusage.size());

  for (auto itr : runs_clusage){
     b++;
     int binIndex=gr->GetXaxis()->FindBin(b);
      gr->GetXaxis()->SetBinLabel(binIndex,itr.data());
      gr->GetXaxis()->ChangeLabel(binIndex,60,-1,39,-1,-1);
  }
  gr->GetXaxis()->SetLabelSize(0.027);
  gr->SetNameTitle("Fraction of clusters used in tracking", "Fraction of clusters used in tracking");
  gr->GetXaxis()->SetTitle("Run");
  gr->GetYaxis()->SetTitle("nCluster in track/total cluster");
  gr->SetMarkerStyle(20);
  gr->Draw("APL");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::reverse(runs_eta.begin(),runs_eta.end());
  std::reverse(relative_eta.begin(),relative_eta.end());

  auto c2 = new TCanvas("Relative number of tracks in Eta distributions");
  c2->SetGrid(0,1);
  auto gr2 = new TGraph(cycle_eta.size(),&cycle_eta[0],&relative_eta[0]);
  double b2=0;
  gr2->GetXaxis()->SetNdivisions(cycle_eta.size());
  gr2->GetXaxis()->Set(1+cycle_eta.size(),gr->GetXaxis()->GetXmin(),1+cycle_eta.size());

  for (auto itr : runs_eta){
     b2++;
     int binIndex=gr2->GetXaxis()->FindBin(b2);
      gr2->GetXaxis()->SetBinLabel(binIndex,itr.data());
      gr2->GetXaxis()->ChangeLabel(binIndex,60,-1,39,-1,-1);
  }
  gr2->GetXaxis()->SetLabelSize(0.027);
  gr2->SetNameTitle("Relative number of tracks in Eta distributions", "Relative number of tracks in Eta distributions");
  gr2->GetXaxis()->SetTitle("Run");
  gr2->GetYaxis()->SetTitle("mean counts / total counts");
  gr2->SetMarkerStyle(20);
  gr2->Draw("APL");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::reverse(runs_phi.begin(),runs_phi.end());
  std::reverse(phicounts.begin(),phicounts.end());

  auto c3 = new TCanvas("Mean number of tracks in Phi distributions (0 < Phi < 2#pi)");
  c3->SetGrid(0,1);

  auto gr3 = new TGraph(cycle_phi.size(),&cycle_phi[0],&phicounts[0]);
  double b3=0;

  gr3->GetXaxis()->SetNdivisions(cycle_phi.size());
  gr3->GetXaxis()->Set(1+cycle_phi.size(),gr3->GetXaxis()->GetXmin(),1+cycle_phi.size());

  for (auto itr : runs_phi){
     b3++;
     int binIndex=gr3->GetXaxis()->FindBin(b3);
      gr3->GetXaxis()->SetBinLabel(binIndex,itr.data());
      gr3->GetXaxis()->ChangeLabel(binIndex,60,-1,39,-1,-1);
  }

  gr3->GetXaxis()->SetLabelSize(0.027);
  gr3->SetNameTitle("Mean number of tracks in Phi distributions (0 < Phi < 2#pi)", "Mean number of tracks in Phi distributions (0 < Phi < 2#pi)");
  gr3->GetXaxis()->SetTitle("Run");
  gr3->GetYaxis()->SetTitle("<counts>");
  gr3->SetMarkerStyle(20);
  gr3->Draw("APL");

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::reverse(phicounts1.begin(),phicounts1.end());
  std::reverse(phicounts2.begin(),phicounts2.end());
  std::reverse(relative_phi1.begin(),relative_phi1.end());
  std::reverse(relative_phi2.begin(),relative_phi2.end());
  
  auto c31 = new TCanvas("Relative number of tracks in Phi distributions");
  c31->SetGrid(0,1);

  auto gr31 = new TGraph(cycle_phi.size(),&cycle_phi[0],&relative_phi1[0]);
  auto gr32 = new TGraph(cycle_phi.size(),&cycle_phi[0],&relative_phi2[0]);
  double b31=0;

  gr31->GetXaxis()->SetNdivisions(cycle_phi.size());
  gr31->GetXaxis()->Set(1+cycle_phi.size(),gr31->GetXaxis()->GetXmin(),1+cycle_phi.size());

  for (auto itr : runs_phi){
     b31++;
     int binIndex=gr31->GetXaxis()->FindBin(b31);
      gr31->GetXaxis()->SetBinLabel(binIndex,itr.data());
      gr31->GetXaxis()->ChangeLabel(binIndex,60,-1,39,-1,-1);
  }

  gr31->GetXaxis()->SetLabelSize(0.027);
  gr31->SetNameTitle("Relative number of tracks in Phi distributions", "Relative number of tracks in Phi distributions");
  gr31->GetXaxis()->SetTitle("Run");
  gr31->GetYaxis()->SetTitle("Number of counts in subrange / Total counts 0-2#pi");
  gr31->SetMarkerStyle(20);
  gr31->SetMarkerColor(kRed);
  gr31->SetLineColor(kRed);
  gr31->Draw("APL");
  gr32->SetMarkerStyle(20);
  gr32->SetMarkerColor(kBlue);
  gr32->SetLineColor(kBlue);
  gr32->Draw("PL SAME");

  auto legend=new TLegend(0.769784,0.842456,0.91207,0.927515);
  legend->SetFillColor(0);
  legend->SetHeader("Phi subrange","");
  legend->AddEntry(gr31,"0-#pi","lp");
  legend->AddEntry(gr32,"#pi-2#pi","lp");
  legend->Draw();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::reverse(runs_ncl.begin(),runs_ncl.end());
  std::reverse(ncluster.begin(),ncluster.end());

  auto c4 = new TCanvas("hNClusters");
  c4->SetGrid(0,1);
  auto gr4 = new TGraph(cycle_ncl.size(),&cycle_ncl[0],&ncluster[0]);
  double b4=0;
  gr4->GetXaxis()->SetNdivisions(cycle_ncl.size());
  gr4->GetXaxis()->Set(1+cycle_ncl.size(),gr4->GetXaxis()->GetXmin(),1+cycle_ncl.size());

  for (auto itr : runs_ncl){
     b4++;
     int binIndex=gr4->GetXaxis()->FindBin(b4);
      gr4->GetXaxis()->SetBinLabel(binIndex,itr.data());
      gr4->GetXaxis()->ChangeLabel(binIndex,60,-1,39,-1,-1);
  }
  gr4->GetXaxis()->SetLabelSize(0.027);
  gr4->SetNameTitle("hNClusters", "hNClusters");
  gr4->GetXaxis()->SetTitle("Run");
  gr4->GetYaxis()->SetTitle("<Number of clusters per track>");
  gr4->SetMarkerStyle(20);
  gr4->Draw("APL");
  
//////////////////////////////////////////NClusters RMS//////////////////////////////////////////////////////////////////////////////////////////

  std::reverse(rms_ncl.begin(),rms_ncl.end());

  auto c41 = new TCanvas("RMS_NClustersDistribution");
  c41->SetGrid(0,1);
  auto gr41 = new TGraph(cycle_ncl.size(),&cycle_ncl[0],&rms_ncl[0]);
  double b41=0;
  gr41->GetXaxis()->SetNdivisions(cycle_ncl.size());
  gr41->GetXaxis()->Set(1+cycle_ncl.size(),gr41->GetXaxis()->GetXmin(),1+cycle_ncl.size());

  for (auto itr : runs_ncl){
     b41++;
     int binIndex=gr41->GetXaxis()->FindBin(b41);
      gr41->GetXaxis()->SetBinLabel(binIndex,itr.data());
      gr41->GetXaxis()->ChangeLabel(binIndex,60,-1,39,-1,-1);
  }
  gr41->GetXaxis()->SetLabelSize(0.027);
  gr41->SetNameTitle("RMS_NClustersDistribution", "RMS_NClustersDistribution");
  gr41->GetXaxis()->SetTitle("Run");
  gr41->GetYaxis()->SetTitle("<RMS>");
  gr41->SetMarkerStyle(20);
  gr41->Draw("APL");
  

//////////////////////////////////////////NClusters. Bin low edge of the first and last bins filled//////////////////////////////////////

  std::reverse(first_bin_edge_ncl.begin(),first_bin_edge_ncl.end());
  std::reverse(last_bin_edge_ncl.begin(),last_bin_edge_ncl.end());
  
  auto c43 = new TCanvas("Low edges of first and last bins filled in NClustersDistribution");
  c43->SetGrid(0,1);
  auto gr43 = new TGraph(cycle_ncl.size(),&cycle_ncl[0],&last_bin_edge_ncl[0]);
  auto gr42 = new TGraph(cycle_ncl.size(),&cycle_ncl[0],&first_bin_edge_ncl[0]);
  double b43=0;
  gr43->GetXaxis()->SetNdivisions(cycle_ncl.size());
  gr43->GetXaxis()->Set(1+cycle_ncl.size(),gr43->GetXaxis()->GetXmin(),1+cycle_ncl.size());

  for (auto itr : runs_ncl){
     b43++;
     int binIndex=gr43->GetXaxis()->FindBin(b43);
      gr43->GetXaxis()->SetBinLabel(binIndex,itr.data());
      gr43->GetXaxis()->ChangeLabel(binIndex,60,-1,39,-1,-1);
  }
  gr43->GetXaxis()->SetLabelSize(0.027);
  gr43->SetNameTitle("Low edges of first and last bins filled in NClustersDistribution", "Low edges of first and last bins filled in NClustersDistribution");
  gr43->GetXaxis()->SetTitle("Run");
  gr43->GetYaxis()->SetTitle("Value in x axis");
  gr43->SetMarkerStyle(20);
  gr43->SetMarkerColor(kRed);
  gr43->SetLineColor(kRed);
  gr43->GetHistogram()->SetMinimum(0.); //Set minimum y value
  gr43->Draw("APL");
  gr42->SetMarkerStyle(20);
  gr42->SetMarkerColor(kBlue);
  gr42->SetLineColor(kBlue);
  gr42->Draw("PL SAME"); 
  auto legend1=new TLegend(0.809353,0.830621,0.951639,0.914941);
  legend1->SetFillColor(0);
  //legend1->SetHeader("Phi range","");
  legend1->AddEntry(gr42,"leftmost bin filled","lp");
  legend1->AddEntry(gr43,"rightmost bin filled","lp");
  legend1->Draw();
 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::reverse(runs_occ.begin(),runs_occ.end());
  std::reverse(occ.begin(),occ.end());

  auto c5 = new TCanvas("Track occupancy in ROF");
  c5->SetGrid(0,1);
  auto gr5 = new TGraph(cycle_occ.size(),&cycle_occ[0],&occ[0]);
  double b5=0;
  gr5->GetXaxis()->SetNdivisions(cycle_occ.size());
  gr5->GetXaxis()->Set(1+cycle_occ.size(),gr5->GetXaxis()->GetXmin(),1+cycle_occ.size());

  for (auto itr : runs_occ){
     b5++;
     int binIndex=gr5->GetXaxis()->FindBin(b5);
      gr5->GetXaxis()->SetBinLabel(binIndex,itr.data());
      gr5->GetXaxis()->ChangeLabel(binIndex,60,-1,39,-1,-1);
  }
  gr5->GetXaxis()->SetLabelSize(0.027);
  gr5->SetNameTitle("Track occupancy in ROF", "Track occupancy in ROF");
  gr5->GetXaxis()->SetTitle("Run");
  gr5->GetYaxis()->SetTitle("nTracks/ROF");
  gr5->SetMarkerStyle(20);
  gr5->Draw("APL");
  
TFile *f = new TFile("../Plots/TrackTask_Result_trends.root", "RECREATE");
f->cd();
c1->Write();
c2->Write();
c3->Write();
c31->Write();
c4->Write();
c41->Write();
c43->Write();
c5->Write();
f->Close();

}



