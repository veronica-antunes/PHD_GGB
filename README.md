# PHD GGB paper

Here you can find some of the source codes described or used to get the results presented in the below publication. The codes are described in detail in the documentation provided. More information can be find in my PhD thesis or in the paper.

If you use any of these tools, please cite as:

Antunes V., Planès, T., Zahradník, J., Obermann, A., Alvizuri, C., Carrier, A.,& Lupi, M. (2020). 
Seismotectonics and 1-D velocity model of the Greater Geneva Basin, France–Switzerland. 
Geophysical Journal International, 221(3), 2026-2047. 
DOI: https://doi.org/10.1093/gji/ggaa129


Each script runs independently from the others and contains one or several functions that can be used to obtain the results. Together with each script, I provide a PDF document that can be consulted. This PDF has all information needed to use the tools.

  - Magnitude.py: Calculates the Magnitude of the events present in a catalogue. The input file needs to be in NORDIC format and follow a SEISAN database structure (Check SEISAN manual for more details - http://seisan.info/)

  - Picks_Quality.pt: Gives a weight number to each station pick. Like in Magnitude.py script, the input file needs to be in NORDIC format and follow a SEISAN database structure. The weight is based on the signal to noise ratio around the pick and the uncertainty interval. The Appendix section in Antunes et al., 2020 explains the process in detail.
  
  - Velest_tools.py: A set of tools that work together with VELEST output files and might be used to automatize VELEST. For more details on VELEST code check https://seg.ethz.ch/software/velest.html.

  - get_isolapz.py: A tool to help to create the input files necessary to use the ISOLA toolbox. The tool creates: i) instrument response file (pz file, with sensitivity, poles and zeros information) from a /several dataless or from an online ArcLink server; ii) converts the waveforms into SAC format. For more information on ISOLA toolbox package, check http://geo.mff.cuni.cz/~jz/for_ISOLAnews/.
