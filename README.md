# QC analysis instructions
In the following some instructions about how to use the scripts are presented. 

## Purpose of this repository
The purpose of the code developed in this repository is to analyse data coming from the commissioning runs on the available ITS layers. The data are downloaded from the CCDB database or read directly from the compurer on which data are saved after each run. 

## Prepare the environment
In order to be able to run the scripts you will need to build the QualityControl and O2 software. To this aim, you can follow the instruction provided [here](https://github.com/MYOMAO/QualityControl/blob/master/README.md#installing-qc-with-alibuild) for the QC software (installation based on aliBuild). The script *startanalysis.sh* will load the QualityControl environment with:
```bash
alienv load QualityControl/latest
```
**Be sure that this command works on your machine**. If it works only in the installation directory you might want to add the following lines to your *.bashrc*:
```bash
export ALIBUILD_WORK_DIR="$HOME/alice/sw"
eval "`alienv shell-helper`"
```
If this is not the case you will have to copy the folder QCanalysis in the folder where the *alienv* command works. 
Then, you need to have permissions to execute the script *startanalysis.sh*. To do so, just do on your terminal (in the folder where the script is):
```bash
chmod +x startanalysis.sh
```

## Start the analysis
To start the download (optional) and analysis of the data you simply need to run the script *startanalysis.sh*:
```bash
./startanalysis.sh
```
Then you need to follow the instructions that appear on the terminal window. The workflow is described in the following. **Remember that you can exit the program with Ctrl+c at any moment**. 

### Choose what to do (main menu)
When starting the script, it will load automatically the QualityControl modules. After this, a first menu will appear on the terminal allowing the choice of the following options:
1. Download and analyse data
2. Analyse data only

You just need to enter the option number (1 or 2). Any other character that will be inserted will be recognized as invalid and the option will have to be retyped. 
#### Option 1 - Download and analyse data
In case option 1 from the menu above is chosen, the script compiles automatically the code for database access. First, a list of all the available *tasks* in the database is shown (the database is accessible also from [here](http://ccdb-test.cern.ch:8080/browse/)). You need to type a *task name* choosing among the one proposed (type *quit* to exit). The tasks used during shifts are **ITSRAWDS** and **ITSQCTrhesholdTask**. They might be the ones you want to type. 
After this, a list of all the *objects* inside the chosen task is shown. You need to copy a single *object name* from where to download the data (type *quit* to exit). Have a look to the end of this paragraph to know which are the correct object names depending on the analysis you want to perform. 

After selecting an object, a menu with the following options will appear:
1. Enter run numbers
2. Enter timestamps

You need to choose one of the options by typing the corresponding number (any other character inserted will be recognized as invalid and the correct option will have to be retyped).
If **option 1** is chosen, you need to type an interval of runs from which to download data (i.e. if the interval is run100 - run200, all the data related to runs having a number between 100 and 200 will be downloaded. **Note the run number is available in the database from July 3rd, 2019 only!**. For older data you will need to choose option 2.  
Instead, if **option 2** is chosen, you will need to type a time interval from which to download data: the time interval will include a date (DD MM YY) and a time in H24 format (HH MM SS). 
The final result will be a *.root* saved in the directory *Data/* containing the histrograms of each run or timestamp. This is the file to use for the analysis of the data. 
After the downloading of the data, a menu with the **available analyses** is shown (**!IMPORTANT! Other analyses will be added soon**):
1. Fake-hit rate run by run
2. Compare noisy pixels between runs

If **option 1** is chosen, the software will show you the fake-hit rate as a function of the run number for all the chips in a Stave. The ROOT macro is run automatically and a list of (good) files for the analysis is shown. Copy and paste the file name in the input line that will appear. **A good file for the analysis must be downloaded from objects named as *ITSQC/Occupancy/LayerXX/LayerXXStaveYYHITMAP* where XX is the layer number, YY the stave number**. The plot is automatically saved in *pdf* and *root* format in the repository *Plots/*. 
If **option 2** is chosen, the software will compare the number of noisy pixels run by run. The ROOT macro is run automatically and a list of (good) files for the analysis is shown. Copy and paste the file name in the input line that will appear. Then, a reference run number (**choose a run among the ones in the input file!**) has to be typed. All the runs in the input file will be compared with the chosen reference run. **A good file for the analysis must be downloaded from objects named as *ITSQC/Occupancy/LayerXX/LayerXXStaveYYHITMAP* where XX is the layer number, YY the stave number**. The plot is automatically saved in *pdf* and *root* format in the repository *Plots/*. 

#### Option 2 - Analyse data only
For this option you will need to have already a file containing a set of histograms to analyse. If this is the case, choosing this option, you will access directly the analysis menu. Choose a single analysis as explained in the paragraph above

## Contacts
For any issue, please contact <ivan.ravasenga@cern.ch>.
