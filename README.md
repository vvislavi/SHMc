# SHMc
SHMc + FastReso + corona spectra (+ratio + RAA) calculation, published in JHEP 07 (2021) 035.
***************************************
# 1. Quick guide
****************************************
## 1.1. Prerequisites
- To calculate BW+FastReso spectra, we need FastReso kernels. These can be calculated using FastReso package:
https://github.com/amazeliauskas/FastReso
There are sample kernels provided with the package for D0, D+, D*+, Ds, and Lc, but they might be not the most up-to-date of writing.
- For corona contribution, we need the measured pp data. This data _does not_ come with the package, but I add a bash script that fetches the latest (as of writing) relevant data from HEPData. To get it, run:
source DownloadHEPData.sh
- For drawing plots, we have an option to add AA data for comparisons (will be discussed later). These also need to be fetched from HEPData and are downloaded in the previous step.
- I would usually run the code in aliroot. At time of writing the README, I noticed that when running code in ROOT, I get a lot of warnings about interpolation going out of range. These are not present when running in aliroot though. The results, whether running with root or aliroot, do not differ despite the warnings. I'll address the source of this error little bit later.

## 1.2. Running the code
- To calculate different species spectra, RAA, and ratios for a given setup, run:
`aliroot -l SpectraCombination.cxx` Alternatively, one could also run it with root, but then there will be quite few warning messages regarding extrapolation going out of range. I'm currently looking into this.
- This calculates different variations of spectra, raa, and ratios to D0 based on given event parameters described in *Scripts/EventParameters.h* and set in *Scripts/Config.C*, and stores them in the output directory (`outputDir`), specified in *Scripts/EventParameters.h*. Then, the mean values and errors are calculated for each observable and stored in the same file. In addition, the mean values and errors for each observable are printed out into text files, found in outputDir/DataTables/.

- To draw "standard" plots (do not store them), run (both root and aliroot works):
`root -l DrawPlots.cxx\(0\)`
- To draw and save the plots, run:
`root -l DrawPlots.cxx\(1\)`
- The `DrawPlots(int SavePlots=kTRUE)` function defined in *DrawPlots.cxx* reads files from outputDir (defined in Scripts/EventParameters.h) and draws all plots in a similar style that was used for the paper. The legend and label positions are pre-set and so they might overlap with data or model points, and should be explicitly tuned in respective functions. If the `SavePlots=1`, these plots are stored in directory {outputDir}_Plots, where {outputDir} is defined in Scripts/EventParameters.h. One does not need to create these directories manually, as they get created during run-time. Also note that the DrawPlots.cxx macro *does not* recalculate the spectra, but instead uses spectra that have been precalculated in SpectraCombination.cxx macro.

## 1.3. Changing particle and/or event parameters

The main configuration is given by three files:
- *Scripts/ParticleParameters.h* contains the definitions of particles that we are/might be interested in. If we want to add another particle into the mix, we have to first define it in this file. The different members of the `struct Particle` are defined within the file.
- *Scripts/EventParameters.h* contains different event configuration settings, such as temperature, beta_max, freeze-out surface, etc. These members are described in the EventParameters.h file. One could use just a single configuration and change it's parameters, but I would strongly suggest against this to make book keeping easier
- *Scripts/Config.C* is the main setting file, defining:
   - `KernelFile`, `ThermalKernelFile`: paths to full and thermal kernel files
   - `gEv`: selected event configuration to be used in all the scripts. Set this to Central, Peripheral, CentralMusic, or whichever other configuration you define in EventParameters.h
   - `allParticles`: a vector of all the particles to consider. All of them have to be defined in ParticleParameters.h. For all these particles, spectra, RAA, and ratios to D0 (except for D0, naturally), are calculated.
   - `allSurf`: all surface types to be considered. The implementation of these surfaces are in Scripts/GetBWSpectra.C and allSurf are used for drawing different models
   - `surfLabels`: labels for each surface when drawn on a legend
   - `addHiddenStates`: a flag to add "with hidden states" to the "SHMc+FastReso+corona" label in legends. NOTE!!!! This flag doesn _not_ add additional states to anything; it only modifies the drawn legend. To modify additional states, one needs to set the scaling factor in ParticleParameters.h for each particle.

***************************************
# 2. Known issues
****************************************
- When running SpectraCombination.cxx in ROOT, one gets many warning messages about interpolating TH3D out-of-range. This is because for fast interpolation (when histogram bins are of the same size), the interpolation is done using the two adjacent bins and this creates a problem for boundary values. In particular, if the first bin is e.g. 0.1-0.2 and we are interpolating at 0.11, there is not bin below 0.1 to use for interpolation. This I am working to address.
- When quitting after:
`root -l DrawPlots.cxx`
there are warning messages about clearing lists. This is because of the way how root handles canvases and TF1, and in particular is because the TF1s drawn in `BWOnly_DrawRatiosToThermal` are destroyed when removing the canvas, but the destructor of TList of root objects inside the canvas is invoked and the same functions are "attempted" to be deleted again.
