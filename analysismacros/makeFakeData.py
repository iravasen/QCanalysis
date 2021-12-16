import ROOT
import numpy as np
import random
ROOT.EnableImplicitMT() # Tell ROOT you want to go parallel

fakeDataFile = ROOT.TFile.Open("xx","RECREATE")

layer = [0,1,2,3,4,5,6]
chips = [9,9,9,8,8,14,14]
stave = [12,16,20,48,60,84,96]

for run in range(0,10):
	for i in range (0,7):
		hist = ROOT.TH2F("hist", "hist", chips[layer[i]], 0, chips[layer[i]],stave[layer[i]],0,stave[layer[i]])
		hist.SetStats(0)
		hist.SetTitle("DeadPixel")
		#hist.SetTitle("Threshold")
		hist.SetName('h2_L{}_run{}'.format(layer[i],run+100000))

		for j in range (0,chips[layer[i]]):
			for k in range (0,stave[layer[i]]):
				hist.SetBinContent(j+1,k+1,random.random()*4+8)

		hist.Write()