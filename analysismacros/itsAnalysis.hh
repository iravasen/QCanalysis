#ifndef itsAnalysis_hh
#define itsAnalysis_hh

#include <string>
#include <iostream>
#include <vector>
#include <TFile.h>
#include <TList.h>

//using namespace std;

//---------------------------------------------------------------------------
// Description:
// Class for loading ROOT files to analyse in ALICE ITS2 Threshold QC
// Author: B. Hofman
//---------------------------------------------------------------------------

class itsAnalysis {

private :
	vector<string> runNumbers_;
	vector<string> layerNumbers_;

	vector<string>  runNumbers_skips;
	std::vector<TH2*> hmaps_skips;

	string skip_ans,skip_runs,skip_layers;

	std::vector<TH2*> hmaps;

public :
	itsAnalysis(string histType) { // Function to load in a ROOT file
		// We will ask for the file to analyse and load that
		cout<<"\n\n=> Available file(s) for the analysis: \n"<<endl;
		gSystem->Exec("ls ../Data/*THRMAPS_DEADPIXMAPS* -Art | tail -n 500");
		cout<<"\nCopy file name: ";
		
		string filepath;
		cin>>filepath;
		cout<<endl;

		TFile *file = TFile::Open(filepath.c_str());
		if (!file) {std::cout<<"Could not open root file"<<std::endl;}

		TObject *obj; TKey *key;
		TIter next( file->GetListOfKeys());
		while ((key = (TKey *) next())) {
		    obj = file->Get(key->GetName()); // copy object to memory
		    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
		         && (!obj->InheritsFrom("TH2"))
			       && (!obj->InheritsFrom("TH1"))
		         && (!obj->InheritsFrom("THnSparse"))
		       ) { cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;}   
			string objname = (string)obj->GetName();
			string objtitle = (string)obj->GetTitle();
			if(objname.find("Stv")!=string::npos) break;
			if(objtitle.find(histType)==string::npos) continue;

			//Load in histograms into vector<TH2*>
			hmaps.push_back((TH2*)obj);

			//Find runNumber if unique
			string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
			if (std::find(runNumbers_.begin(), runNumbers_.end(), runnum) == runNumbers_.end()) {
			  runNumbers_.push_back(runnum);
			}
			//Find layerNumbers_ if unique
			//string laynum = objname.substr(objname.find("L")+1,1);
			string laynum =  objname.find("L")==string::npos ? "noLayer":objname.substr(objname.find("L")+1, 1);
			if (std::find(layerNumbers_.begin(), layerNumbers_.end(), laynum) == layerNumbers_.end()) {
			  layerNumbers_.push_back(laynum);
			}

		} // end of loop over objects

		cout << "Found run numbers: ";
		for (auto i : runNumbers_){ cout << i <<", "; }
		cout << endl << "Found layers: ";
		for (auto i : layerNumbers_){ cout << i <<", "; }
		cout << endl<<endl;

		cout<<"Would you like to skip some run(s)? [y/n]";
		cin>>skip_ans;
		if(skip_ans=="y" || skip_ans=="Y"){
		  cout<<"Specify run number(s) separated by comma (no white spaces!):";
		  cin>>skip_runs;
		  cout<<endl;
		}

		for (auto hist: hmaps)	{
			string objname = (string)hist->GetName();
			string run     = objname.substr(objname.find("run")+3, 6);
			if (skip_runs.find(run) == std::string::npos)
				hmaps_skips.push_back(hist);
		}

		for (auto run : runNumbers_) {
			if (skip_runs.find(run) == std::string::npos)
				runNumbers_skips.push_back(run);
		}


	} // end of loadFile()s

	int nChips(int layer) { // Function to return number of chips in a its Layer
		int nChips_ = -999;
		if (layer == 0 || layer == 1 || layer == 2) nChips_ = 9;
		if (layer == 3 || layer == 4) nChips_ = 8;	// Modules instead of chips
		if (layer == 5 || layer == 6) nChips_ = 14; // Modules instead of chips
		return nChips_;
	}

	int stavesInLayer(int layer)	{
		int staves;
		if (layer == 0) staves = 11;
		if (layer == 1) staves = 16;
		if (layer == 2) staves = 20;
		if (layer == 3 || layer == 4) staves = 54;
		if (layer == 5 || layer == 6) staves = 90;
		return staves;
	}

	vector<string> Layers() {
		return layerNumbers_;
	}

	vector<string> Runs() {
		return runNumbers_skips;
	}

	int nLayers() {
		return layerNumbers_.size();
	}

	int nRuns() {
		return runNumbers_skips.size();
	}

	std::vector<TH2*> loadedHists() {
		return hmaps_skips;
	}

	std::vector<TH2*> loadLayer(int layer) {
		std::vector<TH2*> hmaps_per_layer; // empty list, histograms of layer
		for (auto hist: hmaps_skips)	{
			string objname = (string)hist->GetName();
			if (stoi(objname.substr(objname.find("L")+1,1))==layer){ // Check title for layer
				hmaps_per_layer.push_back(hist);
			}
		}
		return hmaps_per_layer;
	}

	string getRunNumber(TH2* hist) {
		string objname = hist->GetName();
		string runN = objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
		return runN;
	}

}; // end of class def

#endif