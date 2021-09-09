import yaml
import argparse
import ROOT
import array as a
import numpy as np

def main():

    coltochipid0 = {0:16, 1:17, 2:18, 3:19, 4:20, 5:21, 6:22,
                    7:32, 8:33, 9:34, 10:35, 11:36, 12:37, 13:38,
                    14:48, 15:49, 16:50, 17:51, 18:52, 19:53, 20:54,
                    21:64, 22:65, 23:66, 24:67, 25:68, 26:69, 27:70,
                    28:80, 29:81, 30:82, 31:83, 32:84, 33:85, 34:86,
                    35:96, 36:97, 37:98, 38:99, 39:100, 40:101, 41:102,
                    42:112, 43:113, 44:114, 45:115, 46:116, 47:117, 48:118}

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True, help="Input file to be analysed")
    args = parser.parse_args()
    print(f"Analysing file: {args.file}")

    #open file
    infl = ROOT.TFile.Open(args.file, "READ")
    ntriggers = GetTriggers(infl)
    print(f"Number of triggers: {ntriggers}")
    #Loop over all THnSparse and prepare yaml file with noisy pixels
    NOISECUT = 1e-6
    for key in infl.GetListOfKeys():
        obj=key.ReadObj()
        if obj.InheritsFrom("THnSparse"):
            if obj.GetEntries()-1 < 0: #skip empty THnSparse
                continue
            ## Loop of pixels which fired
            layer = obj.GetName()[9:10]
            dict = {}
            npix = 0
            for ibin in range(obj.GetNbins()):
                coord = np.array([0,0], dtype=np.int32)
                pixelhits = obj.GetBinContent(ibin, np.asarray(coord))
                fhr = pixelhits/ntriggers #fhr of the pixel
                if fhr < NOISECUT: ##noise cut
                    continue
                npix = npix+1
                chipid = 0
                if int(layer)<3: #ib
                    chipid = int((coord[0]-1)/1024)
                    if chipid not in dict:
                        dict.update({chipid:[[int(coord[0]-1-chipid*1024),int(coord[1]-1),fhr]]})
                    else:
                        dict[chipid].append([int(coord[0]-1-chipid*1024),int(coord[1]-1),fhr])
                else: #ob - to be added
                    continue
            ##save yaml
            stavenum = obj.GetName()[14:15] if obj.GetName()[15:16] == "_" else obj.GetName()[14:16]
            print(f"L{layer}_{int(stavenum):02d}: {npix} hot pixels above cut")
            with open(f"../yaml/noise_masks/L{layer}_{int(stavenum):02d}.yml", 'w') as f:
                yaml.dump(dict, f)


## Function to get number of triggers
def GetTriggers(infl):
    fhr_chip_ib = 0.
    nhits_chip_ib = 0.
    for key in infl.GetListOfKeys():
        obj=key.ReadObj()
        name = obj.GetName()
        if obj.InheritsFrom("TH2") and ("L0" in name): ##just take L0
            fhr_chip_ib = obj.GetBinContent(1,1)
            if fhr_chip_ib<1e-15:
                print("Please change chip, this bin is empty")
        if obj.InheritsFrom("THnSparse") and ("L0_Stv0" in name): ##just take L0_00
            h2 = obj.Projection(1,0)
            nhits_chip_ib = h2.Integral(1,1024,1,512)
            del h2
            break
    return nhits_chip_ib / (512.*1024.*fhr_chip_ib)

if __name__ == '__main__':
    main()
