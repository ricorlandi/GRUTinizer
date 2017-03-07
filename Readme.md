GRUTinizer
==========

A generic unpacker and analysis package for gamma-ray spectroscopy.

The doxygen documentation can be found [here](https://pcbend.github.io/GRUTinizer).


--
The most important thing to know is that in order to use GRUTinizer, you must *first* run the ```grutsh``` command in your shell session. (e.g. after you have logged in to cagragr@miho-1 or aino-1, you type ```grutsh```. This will source a number of need scripts and directories to correctly setup your environment to compile and run GRUTinizer. 

Next, copy the ```~/ana/template_campaign_december/``` directory to your own analysis space. Compile it with 

``` make -j22 ```

Note, you do not need to ```make clean``` unless except in rare instances. Make should handle the dependencies for you.



Online CAGRA example for the E450 experiments:
--
``` ./bin/grutinizer -Hmg lib/libRCNPhistos.so config/online.cal config/rcnp_constants.cal "data/cagra/E450/run_1044.gtd03_000_0109" data/gr/E450/run1044.bld -o /dev/null ```

* Note: the path for the CAGRA data has to be given in quotation marks in contrast to the GR data!
* The file name ```run_1044.gtd03_000_0109``` is composed of:

  **- run number: 1044**

  **- IOC number: gtd03 -> IOC 3**

  **- file number: 000 (first file recorded for IOC 3)**

  **- boardID: 0109**

  **- after 2GB another file is opened, e.g., run_1044.gtd03_001_0109**

* To unpack all files for a given IOC replace ```run_1044.gtd03_000_0109``` by ```run_1044.gtd03_*```
* To unpack all files replace ```run_1044.gtd03_000_0109``` by ```run_1044.gtd0*```


Offline CAGRA example for the E450 experiments:
--

* For an offline analysis I'd suggest to unpack the data into a root file without creating histograms and opening the GUI. The option ```-qm``` takes care of this:

```./bin/grutinizer -qm lib/libRCNPhistos.so config/online.cal config/rcnp_constants.cal "data/cagra/E450/run_1044.gtd0*" data/gr/E450/run1044.bld -o run1044.root```

* The data from CAGRA and GR are now sorted into combined events. Create a root file ```hist1044.root``` including the histograms defined in ```histos/RCNPhistos.cxx```:

```./bin/grutinizer -qH lib/libRCNPhistos.so config/online.cal config/rcnp_constants.cal run1044.root```

* The ```hist1044.root``` can be loaded into the the GUI by:

```./bin/grutinizer -g hist1044.root```

* Use ```-Hg``` to simultaneously open the GUI while creating the histograms.

General remarks
--

* If the grutinizer crashes you can open a new terminal and log into miho-1, miho-2, aino-1 or aino-2:

```> ssh -Y cagragr@miho-1```

* Go to ```/home/cagragr/ana/jisaak/template_campaign``` and run the command ```grutsh``` before you can use grutinizer.
* Energy calibration parameters for the individual detectors are defined in ```config/online.cal``` or other similar named files. Please check or ask current shifts.

ToDo's during experiment shifts
--

* After a new run is started please make an entry into the logbook. Please make sure that you DO NOT use a run number in the CAGRA DAQ which already exists otherwise the DAQ crashes and has to be restarted hardware-wise.

* Start grutinizer for the new run and check at least if we have true coincidences between GR and CAGRA. Therefore, open in the GUI of grutinizer the directory ```COIN``` and check the time-difference spectra named ```Diff_CAGRA_GR``` and ```Diff_CAGRA_GR_Cond```. The first spectrum shows the time-difference without any conditions/cuts on anything, while the second spectrum applied a cut on the X_TH in the GR to clean-up the spectrum. In both spectra prompt peaks around bin=250 should be visible.

* Check from time to time if the rates shown in the VME windows of the CAGRA DAQ are updating and or not frozen or something like this.

* Check the Buffer of each VME. It happened sometimes that they ran full and the DAQ crashed and it had to be restarted hardware-wise in the experimental area. The Buffers usually oscillate between 390 and 400.

* Check the Rate Meter of GR for the rates and live times. Every 4-6 hours the GR DAQ suddenly has shows live times < 2%. This will be also visible in the CAGRA rates, because less coincidences are written, obviously. The reason for the sudden dead time rise is unclear at the moment.

* If one of the checks is negative stop the current run and start a new one. Please leave a comment in the logbook. If the problems still occur please contact one of the Core Member: Tamii, Kobayashi or Isaak.

* After you started a new run please unpack the previous one. The CAGRA and GR raw data will be sorted and events are created and saved into a root file. There exists a script in the working directory called ```sort_raw```. It's better to open a new terminal, e.g. by typing in the current console ```xterm &``` and run the script with the respective run number XXXX: ```./sort_raw XXXX```. The unpacked data are written to ```data/rootfiles/E450/RunXXXX```.

* You can check if there are runs which where not unpacked by now and make up for it. Open a new xterm and do as described above. I suggest to run at maximum 3 unpacking procedures simultaneously.



Key program options (flags) for running GRUTinizer:

Options:

   arg                   Input file(s)

  **-H [ --histos ]        attempt to run events through MakeHisto lib.**

  **-m [ --sort-multiple ] If passed multiple raw data files, treat them as one file.**

  **-s [ --sort ]          Attempt to loop through root files.**

  **--build-window arg     Build window, timestamp units**

  **-g [ --start-gui ]     Start the GUI**

  **-o [ --output ] arg    Root output file**

  **-S [ --gr-singles ]    Ignore GR timestamps and take singles.**

  **--hist-output arg      Output file for histograms**

  **-q [ --quit ]          Run in batch mode**

  **-h -? [ --help ]       Show this help message**




**INTERACTIVE ANALYSIS**
--


In addition to the many libraries in TGRUTAnalysis to make analysis easier, GRUTinizer also takes control of some of the behind the scene functions of ROOT to make analysis a bit easier. 

To take advantage of these features all one has to do is start GRUTinizer!  They are implemented when examining and natural root classes, whether the where draw from a tree, made fresh or load from a file.
 
  * **GlobalFunctions** <br>
    
  * **GCanvas** <br>  The most notable difference, is the replacement of the TCanvas with GCanvas.  This replacement is done naturally - no changes from either existing ROOT scripts or GRUTinizer scripts are needed to take advantage of the GCanvas default behavior.  

## Universal Commands

| Key | Action |
|:-----|:------|
| **F2** | Show/Hide Editor Tab |

## TH1 Commands
|  Key  | Action |
|:------|:------|
|        | All normal root commands/interactions.|
| **m**  | Toggle on/off marker mode; when on, the histogram will remember and display the last four clicks as marks on the histogram.|
| **p**  | If the 1d hist was made using the global ProjectionX/ProjectionY; gating the original 2D matrix this histogram came from is possible by placing markers around the gate and pressing p.  The gates spectra is immediately drawn. |
| **B** | Cycle through types of automatic background subtraction used when projecting with **p**.  Current types include: No subtraction, Fraction of the total, subtract gate from the 3rd marker (gate size set to the distance between marker 1 and 2). |
| **b** | Set the background, how it is set depends on **B**.|
| **n**  | Remove all markers / functions drawn on the histogram (not gates!).|
| **e**  | Expand the x-axis range between the last two markers.|
| **E**  | Bring up dialogue box used to set desired x-axis range.|
| **o**  | Unzoom the entire histogram.|
| **Arrow Left/Right** |  When zoomed in, mover the display region to the left/right by one half of the region currently displayed.|
| **Arrow Up/Down** |  Quickly display the next histogram stored in memory, especially useful when gating to go back and forth between gates and the total projection. (currently only available in GH1D)  |
| **f**  | Ryan D's TPeak Fit (proper skewd gaus for gamma-rays with automatic bg) with minimum output. |  
| **g**  | Gaus fit with linear background, displays results of the fit **RESULTS STILL NEED TO BE VERIFIED** |
| **i**  | Raw integral of counts between the two markers |
| **s**  | Show peak values. |
| **S**  | Remove peak values. |
| **l**  | Toggle y-axis linear. |


## GH2I Commands
|  Key  | Action |
|:------|:------|
|        | All normal root commands/interactions.|
| middle-click | Select the current pad, current pad is outlined by a red border. |
| **e**  | Expand the x-axis range between the two markers.|
| **g**  | Create a TCuG on the canvas, name scheme is _cut# where # is tracked from the start of the program.|
| **o**  | Unzoom the entire histogram, x and y.|
| **x**  | Make and display a total projection of the x-axis.|
| **X**  | Make and display a one bin projection of the x-axis, arrow up/down will cycle through all bins.|
| **y**  | Make and display a total projection of the y-axis.|
| **Y**  | Make and display a one bin projection of the y-axis, arrow up/down will cycle through all bins.|


## TGraph/TGraphErrors Commands
|  Key  | Action |:
|:------|:------|
|        | All normal root commands/interactions.|
| **p**  | Print the graph to the terminal.


``` make -j15 ```

Note, you do not need to ```make clean``` unless except in rare instances. Make should handle the dependencies for you.



Key program options (flags) for running GRUTinizer:

Options:

   arg                   Input file(s)

  **-H [ --histos ]        attempt to run events through MakeHisto lib.**

  **-m [ --sort-multiple ] If passed multiple raw data files, treat them as one file.**

  **-s [ --sort ]          Attempt to loop through root files.**

  **--build-window arg     Build window, timestamp units**

  **-g [ --start-gui ]     Start the GUI**

  **-o [ --output ] arg    Root output file**

  **-S [ --gr-singles ]    Ignore GR timestamps and take singles.**

  **--hist-output arg      Output file for histograms**

  **-q [ --quit ]          Run in batch mode**

  **-h -? [ --help ]       Show this help message**


Online CAGRA example:
--
``` ./bin/grutinizer -Hmg config/rcnpchannels.cal libraries/libRCNPhistos.so cagra_data/data1/run_1018.gtd01_000_010* -o /dev/null ```
* Above, -H specifies that a histogram library will be used. For each built event, the function MakeHistograms in ```./histos/RCNPhistos.cxx ``` will be called. Any new *.cxx files put in ./histos will be automatically compiled into a corresponding library and put in ./libraries.

* Above -m indicates that multiple raw files will be sorted simultaneously and attempted to event correlated. This flag should almost always be present when multiple raw files to be sorted.

* -g launches the graphical user interface. Some helper functions and key commands are described in the Interactive Analysis section below.

* -o /dev/null indicates that the output root tree will not be created, and only histograms will be made for online monitoring

Online Grand Raiden example:
--
``` ./bin/grutinizer -Hg config/rcnpchannels.cal libraries/libRCNPhistos.so online.bld -o /dev/null ```
* All options here are similar to those found above, but now we are using a special online.bld file to indicate that GRUTinizer should utilize the online functionality of the GRAnalyzer to retrieve GR event data. Note that online mode only works when you are running GRUTinizer from *aino-1*, since this is the computer on which the Grand Raiden DAQ is run.


Offline CAGRA example:
--
``` ./bin/grutinizer -Hmq config/rcnpchannels.cal libraries/libRCNPhistos.so cagra_data/data1/run_1018.gtd01_000_010* -o run1018.root ```

* The only difference here is that a -o specifies the name of the rootfile that will contain the built root tree from the raw data (ie, we are not skipping this step as in the previous example). This can be nice for calibrations as this root tree can be read in just like any other raw file and analyzed, without needing to resort all the data.

* -q specifies that we want to close GRUTinizer when analysis and sorting is done. Note that in this mode, the GUI is not opened.


Offline Grand Raiden example:
--
``` ./bin/grutinizer -SHg config/rcnpchannels.cal libraries/libRCNPhistos.so ~/data/run1016.bld  ```

* In this case, we are running the graphical interface in offline mode, while sorting data from a file that was produced at some point in the past. Since we did not specify an output filename, the file should detect the run number and name it appropriately (note that this doesnt yet work CAGRA, I think).

* -S here is a special flag to indicate GR singles mode. This is necessarry if for some reason, the MYRIAD timestamp module is not functioning, but looking at the raw singles data is still desired. This flag can be used in offline or online mode.


Coincident CAGRA + GR example:

``` ./bin/grutinizer -Hmg config/rcnpchannels.cal libraries/libRCNPhistos.so online.bld cagra_data/data1/run_1018.gtd01_000_010* -o /dev/null  ```

* Here we attach online to both GR and to raw CAGRA data
 
* In this setting -m is critical as it indicates that the different raw sources should be built together.
 



**INTERACTIVE ANALYSIS**
--


In addition to the many libraries in TGRUTAnalysis to make analysis easier, GRUTinizer also takes control of some of the behind the scene functions of ROOT to make analysis a bit easier. 

To take advantage of these features all one has to do is start GRUTinizer!  They are implemented when examining and natural root classes, whether the where draw from a tree, made fresh or load from a file.
 
  * **GlobalFunctions** <br>
    
  * **GCanvas** <br>  The most notable difference, is the replacement of the TCanvas with GCanvas.  This replacement is done naturally - no changes from either existing ROOT scripts or GRUTinizer scripts are needed to take advantage of the GCanvas default behavior.  

## Universal Commands

| Key | Action |
|:-----|:------|
| **F2** | Show/Hide Editor Tab |

## TH1 Commands
|  Key  | Action |
|:------|:------|
|        | All normal root commands/interactions.|
| **m**  | Toggle on/off marker mode; when on, the histogram will remember and display the last four clicks as marks on the histogram.|
| **p**  | If the 1d hist was made using the global ProjectionX/ProjectionY; gating the original 2D matrix this histogram came from is possible by placing markers around the gate and pressing p.  The gates spectra is immediately drawn. |
| **B** | Cycle through types of automatic background subtraction used when projecting with **p**.  Current types include: No subtraction, Fraction of the total, subtract gate from the 3rd marker (gate size set to the distance between marker 1 and 2). |
| **b** | Set the background, how it is set depends on **B**.|
| **n**  | Remove all markers / functions drawn on the histogram (not gates!).|
| **e**  | Expand the x-axis range between the last two markers.|
| **E**  | Bring up dialogue box used to set desired x-axis range.|
| **o**  | Unzoom the entire histogram.|
| **Arrow Left/Right** |  When zoomed in, mover the display region to the left/right by one half of the region currently displayed.|
| **Arrow Up/Down** |  Quickly display the next histogram stored in memory, especially useful when gating to go back and forth between gates and the total projection. (currently only available in GH1D)  |
| **f**  | Ryan D's TPeak Fit (proper skewd gaus for gamma-rays with automatic bg) with minimum output. |  
| **g**  | Gaus fit with linear background, displays results of the fit **RESULTS STILL NEED TO BE VERIFIED** |
| **i**  | Raw integral of counts between the two markers |
| **s**  | Show peak values. |
| **S**  | Remove peak values. |
| **l**  | Toggle y-axis linear. |


## GH2I Commands
|  Key  | Action |
|:------|:------|
|        | All normal root commands/interactions.|
| middle-click | Select the current pad, current pad is outlined by a red border. |
| **e**  | Expand the x-axis range between the two markers.|
| **g**  | Create a TCuG on the canvas, name scheme is _cut# where # is tracked from the start of the program.|
| **o**  | Unzoom the entire histogram, x and y.|
| **x**  | Make and display a total projection of the x-axis.|
| **X**  | Make and display a one bin projection of the x-axis, arrow up/down will cycle through all bins.|
| **y**  | Make and display a total projection of the y-axis.|
| **Y**  | Make and display a one bin projection of the y-axis, arrow up/down will cycle through all bins.|


## TGraph/TGraphErrors Commands
|  Key  | Action |:
|:------|:------|
|        | All normal root commands/interactions.|
| **p**  | Print the graph to the terminal.





## Add-on for the CAGRA+GR Campaign @ RCNP October 2016:
--

Most of the things explained above are still valid, such as the different options (flags) and universal commands. However, the syntax for the unpacking is a bit different. Furthermore, the working directory is:

It was supposed to be: /home/cagragr/ana/campaignOct2016/, BUT it turned out not to be :-D
So the working directory is actually is: /home/cagragr/ana/jisaak/template_campaign

Online CAGRA example for the E450 experiments:
--
``` ./bin/grutinizer -Hmg lib/libRCNPhistos.so config/online.cal config/rcnp_constants.cal "data/cagra/E450/run_1044.gtd03_000_0109" data/gr/E450/run1044.bld -o /dev/null ```

* Note: the path for the CAGRA data has to be given in quotation marks in contrast to the GR data!
* The file name ```run_1044.gtd03_000_0109``` is composed of:

  **- run number: 1044**

  **- IOC number: gtd03 -> IOC 3**

  **- file number: 000 (first file recorded for IOC 3)**

  **- boardID: 0109**

  **- after 2GB another file is opened, e.g., run_1044.gtd03_001_0109**

* To unpack all files for a given IOC replace ```run_1044.gtd03_000_0109``` by ```run_1044.gtd03_*```
* To unpack all files replace ```run_1044.gtd03_000_0109``` by ```run_1044.gtd0*```


Offline CAGRA example for the E450 experiments:
--

* For an offline analysis I'd suggest to unpack the data into a root file without creating histograms and opening the GUI. The option ```-qm``` takes care of this:

```./bin/grutinizer -qm lib/libRCNPhistos.so config/online.cal config/rcnp_constants.cal "data/cagra/E450/run_1044.gtd0*" data/gr/E450/run1044.bld -o run1044.root```

* The data from CAGRA and GR are now sorted into combined events. Create a root file ```hist1044.root``` including the histograms defined in ```histos/RCNPhistos.cxx```:

```./bin/grutinizer -qH lib/libRCNPhistos.so config/online.cal config/rcnp_constants.cal run1044.root```

* The ```hist1044.root``` can be loaded into the the GUI by:

```./bin/grutinizer -g hist1044.root```

* Use ```-Hg``` to simultaneously open the GUI while creating the histograms.

General remarks
--

* If the grutinizer crashes you can open a new terminal and log into miho-1, miho-2, aino-1 or aino-2:

```> ssh -Y cagragr@miho-1```

* Go to ```/home/cagragr/ana/jisaak/template_campaign``` and run the command ```grutsh``` before you can use grutinizer.
* Energy calibration parameters for the individual detectors are defined in ```config/online.cal``` or other similar named files. Please check or ask current shifts.

ToDo's during experiment shifts
--

* After a new run is started please make an entry into the logbook. Please make sure that you DO NOT use a run number in the CAGRA DAQ which already exists otherwise the DAQ crashes and has to be restarted hardware-wise.

* Start grutinizer for the new run and check at least if we have true coincidences between GR and CAGRA. Therefore, open in the GUI of grutinizer the directory ```COIN``` and check the time-difference spectra named ```Diff_CAGRA_GR``` and ```Diff_CAGRA_GR_Cond```. The first spectrum shows the time-difference without any conditions/cuts on anything, while the second spectrum applied a cut on the X_TH in the GR to clean-up the spectrum. In both spectra prompt peaks around bin=250 should be visible.

* Check from time to time if the rates shown in the VME windows of the CAGRA DAQ are updating and or not frozen or something like this.

* Check the Buffer of each VME. It happened sometimes that they ran full and the DAQ crashed and it had to be restarted hardware-wise in the experimental area. The Buffers usually oscillate between 390 and 400.

* Check the Rate Meter of GR for the rates and live times. Every 4-6 hours the GR DAQ suddenly has shows live times < 2%. This will be also visible in the CAGRA rates, because less coincidences are written, obviously. The reason for the sudden dead time rise is unclear at the moment.

* If one of the checks is negative stop the current run and start a new one. Please leave a comment in the logbook. If the problems still occur please contact one of the Core Member: Tamii, Kobayashi or Isaak.

* After you started a new run please unpack the previous one. The CAGRA and GR raw data will be sorted and events are created and saved into a root file. There exists a script in the working directory called ```sort_raw```. It's better to open a new terminal, e.g. by typing in the current console ```xterm &``` and run the script with the respective run number XXXX: ```./sort_raw XXXX```. The unpacked data are written to ```data/rootfiles/E450/RunXXXX```.

* You can check if there are runs which where not unpacked by now and make up for it. Open a new xterm and do as described above. I suggest to run at maximum 3 unpacking procedures simultaneously.
