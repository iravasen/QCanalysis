# QC analysis instructions
In the following some instructions about how to use the scripts are presented. 

## Purpose of this repository
The purpose of the code developed in this repository is to analyse data coming from the commissioning runs on the available ITS layers. The data are downloaded from the CCDB database or read directly from the computer (FLP) on which data are saved after each run. 

## Prepare the environment
In order to be able to run the scripts you will need to build the QualityControl (origin/master) and O2 (origin/dev) software. Be sure to **build the mentioned versions**. The installation is based on *aliBuild*. For installing aliBuild, see the pre-requisites [here](https://alice-doc.github.io/alice-analysis-tutorial/building/custom.html#prerequisites) and follow [these](https://alice-doc.github.io/alice-analysis-tutorial/building/custom.html#get-or-upgrade-alibuild) instructions for the installation. Then, you need to build the packages. To do so you can follow the standard guide [here](https://alice-doc.github.io/alice-analysis-tutorial/building/build.html#%F0%9F%9B%A0-build-the-packages). The main things to consider are:
1. [Preparation](https://alice-doc.github.io/alice-analysis-tutorial/building/build.html#prepare-your-source-code) of the souce code. Consider the command for O2. For QualityControl you will need to run the command: 
```bash
aliBuild init QualityControl@master --defaults o2
```
2. Check the [pre-requisites](https://alice-doc.github.io/alice-analysis-tutorial/building/build.html#check-your-prerequisites-skip-if-using-alidock) with aliDoctor following the command for O2. 
3. [Build](https://alice-doc.github.io/alice-analysis-tutorial/building/build.html#build-and-rebuild) the O2 software following the specific command line for O2. For QualityControl software use the following command:
```bash
aliBuild build QualityControl --defaults o2
```
The script *startanalysis.sh* will load the QualityControl environment with:
```bash
alienv load QualityControl/latest
```
**Be sure that this command works on your machine**. Then, you need to have permissions to execute the script *startanalysis.sh*. To do so, just do on your terminal (in the folder where the script is):
```bash
chmod +x startanalysis.sh
```
Inside the folder there is also a MakeFile needed to compile *getObject.cxx*. As a **temporary solution**, you need to edit the lines 46,47,48,49,55,56 where you see a path to the "sw/" folder. Edit them by writing the correct path of the sw folder on your computer. 

## Start the analysis
To start the download (optional) and analysis of the data you simply need to run the script *startanalysis.sh*:
```bash
./startanalysis.sh
```
Then you need to follow the instructions that appear on the terminal window. The workflow is described in the following. **Remember that you can exit the program with Ctrl+c at any moment**. 

### Software update
Everytime you run the script *startanalysis.sh*, a pull of this repository is done automatically. In case you performed local modifications (as for the Makefile), they will be kept (git stash).

### Choose what to do (main menu)
When starting the script, it will automatically load the QualityControl modules (after the git repository update). After this, a first menu will appear on the terminal allowing the choice of the following options:

1. Download and analyse data from CCDB
3. Analyse data on flp

You just need to enter the option number (1 or 2). Any other character that will be inserted will be recognized as invalid and the option will have to be retyped. In general, at present, options 1 is meant to analyse data coming from fake-hit rate scans saved into the CCDB while option 3 is for threshold scan and fake-hit scan runs on FLP. 

#### Option 1 - Download and analyse data from CCDB
In case option 1 is chosen, two options are presented
1. Download data from CCDB and perform the analysis
2. Analyse data only (I have already a data sample)
##### Option 1 - Download data from CCDB and perform the analysis
If option 1 (**Download data from CCDB and perform the analysis**) is chosen, the script compiles automatically the code for database access (getObject.cxx) provided that the MakeFile has been properly edited. First, a list of all the available *tasks* in the database is shown (the database is accessible also from [here](http://ccdb-test.cern.ch:8080/browse/)). You need to type a *task name* choosing among the one proposed (type *quit* to exit). The tasks used during shifts are **ITSRAWDS** and **ITSQCTrhesholdTask**. They might be the ones you want to type. 
After this, a list of all the *objects* inside the chosen task is shown. You need to copy a single *object name* from where to download the data (type *quit* to exit). Have a look to the end of this paragraph to know which are the correct object names to select depending on the analysis you want to perform. 

After selecting an object, a menu with the following options will appear:

1. Enter run numbers
2. Enter timestamps

You need to choose one of the options by typing the corresponding number (any other character will be recognized as invalid and the correct option will have to be retyped).
If option 1 (**Enter run numbers**) is chosen, you will need to type an interval of runs from which to download data (i.e. if the interval is 100 - 200, all the data related to runs having a number between 100 and 200 will be downloaded). **Note that the run number is available in the database from July 3rd, 2019 only!**. For older data you will need to choose option 2 (**Enter timestamps**).  
Instead, if option 2 (**Enter timestamps**) is chosen, you will need to type a time interval from which to download data: the time interval will include a date (DD MM YY) and a time in H24 format (HH MM SS). 
In both cases, the final result will be a *.root* file saved in the directory *Data/* containing the histrograms of each run or timestamp. This is the file to use for the analysis (next step). 
After the downloading of the data, a menu with the available analyses is shown (**!IMPORTANT! Other analyses will be added soon**):

1. Fake-hit rate run by run
2. Compare noisy pixels between runs
3. Option 1 and 2 together

If option 1 (**Fake-hit rate run by run**) is chosen, the software will show you the fake-hit rate as a function of the run number (or timestamp) for all the chips in a Stave. The ROOT macro is run automatically and a list of (good) files for the analysis is shown. Copy and paste the file name in the input line that will appear. **A good file for the analysis must be downloaded from objects named as *ITSQC/Occupancy/LayerXX/LayerXXStaveYYHITMAP* where XX is the layer number, YY the stave number**. The final plot is automatically saved in *pdf* and *root* format in the repository *Plots/*. 
If option 2 (**Compare noisy pixels between runs**) is chosen, the software will compare the number of noisy pixels run by run. The ROOT macro is run automatically and a list of files for the analysis is shown. Copy and paste the file name in the input line that will appear. Then, a reference run number (**choose a run among the ones in the input file!**) has to be typed. All the runs in the input file will be compared with the chosen reference run. **A good file for the analysis must be downloaded from objects named as *ITSQC/Occupancy/LayerXX/LayerXXStaveYYHITMAP* where XX is the layer number, YY the stave number**. The plot is automatically saved in *pdf* and *root* format in the repository *Plots/*. 
Finally, if option 3 is chosen, the two previous options will be performed at the same time since the data sample can be the same. 

##### Option 2 - Analyse data only (I have already a data sample)
For this option you will need to have already a file containing a set of histograms to analyse (downloaded from CCDB!). If this is the case, choosing this option, you will access directly the analysis menu. Choose a single analysis as explained in the previous paragraph (after database downloading). 

#### Option 2 - Analyse data on flp 
This option allows you to analyse data directly on FLP skipping the download from the database. This option can be chosen both if you are working on your computer and on FLP. First, you need to choose on which FLP the data you want to analyse are. **Type only the number of the FLP**. Depending on your choice, if you are not working already on the selected FLP, the script will ask you to type your CERN username. This is used to make a connection via ssh to the selected flp (through lxplus). You just need to insert first your CERN account password and then the one of the "its" account on the selected flp. 
Then, the script asks to type the layer number. After this, an update of the repository is automatically done keeping any local (on flp) modifications. 
Later, the QualityControl environment is loaded on the flp and the script asks to type the path in which data are saved (tipycally: */data/shifts/* or */data/L0_shifts/*). When the path is selected, the script asks you to choose an analysis type:

1. Average stave thresholds run by run
2. Compare dead pixels between runs
3. Option 1 and 2 together
4. Fake-hit rate run by run
5. Compare noisy pixels between runs
6. Option 4 and 5 together

In general, once you have selected one of the options above the script asks you whether you want to prepare and analyse a data sample or you want to analyse an already existing sample. You just need to enter the option you want to perform. 
When you choose to prepare and analyse a data sample, for all the options 1-6 above, the script will ask you to type a run interval. If you choose "starting run = 100" and "final run = 200", all the runs between 100 and 200 will be included in the sample and hence in the analysis. Then, for options 1, 3, 4, 6 you need to type also the number of the first stave in the layer (for example "6" if you are analysing the layer 0 currently under commissioning -> stave numbered from 6 to 11). Now, the details of each option (1-6) will be described. 
##### Option 1 - Average stave thresholds run by run
If this option is chosen, the software will show you the average threshold (stave by stave on the layer) as a function of the run number. The ROOT macro is run automatically and a list of files for the analysis is shown. Copy and paste the file name in the input line that will appear. The plot is automatically saved in *pdf* and *root* format in the repository *Plots/* **of the flp selected**.
##### Option 2 - Compare dead pixels between runs
If this option is chosen, the software will compare the number of dead pixels (i.e. having null threshold) run by run. The ROOT macro is run automatically and a list of files for the analysis is shown. Copy and paste the file name in the input line that will appear. Then, a reference run number (**choose a run among the ones in the input file!**) has to be typed. All the runs in the input file will be compared with the reference run. The plot is automatically saved in *pdf* and *root* format in the repository *Plots/* **of the flp selected**.
##### Option 3 - Option 1 and 2 together 
This option will simply do the option 1 and option 2 at the same time since the data sample for the analysis can be the same. 
##### Option 4 - Fake-hit rate run by run
If this option is chosen, the software will show you the fake-hit rate as a function of the run number for all the staves composing the layer. The ROOT macro is run automatically and a list of files for the analysis is shown. Copy and paste the file name in the input line that will appear. The final plot is automatically saved in *pdf* and *root* format in the repository *Plots/* **of the flp selected**.
##### Option 5 - Compare noisy pixels between runs
If this option is chosen, the software will compare the number of noisy pixels run by run on all the layer (no stave-by-stave analysis at present). The ROOT macro is run automatically and a list of files for the analysis is shown. Copy and paste the file name in the input line that will appear. Then, a reference run number (**choose a run among the ones in the input file!**) has to be typed. All the runs in the input file will be compared with the reference run. The plot is automatically saved in *pdf* and *root* format in the repository *Plots/* **of the flp selected**.
##### Option 6 - Option 4 and 5 together
This option will simply do the option 4 and option 5 at the same time since the data sample for the analysis can be the same.
## Contacts
For any issue, please contact <ivan.ravasenga@cern.ch>.
