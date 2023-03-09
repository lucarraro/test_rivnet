# test_rivnet

Scripts supporting Carraro, L., "Seamless extraction and analysis of river networks in R:  the rivnet package"

 - `test_rivnet.R` is the main script file, and which produces Figs. 2, 3 and 4 of the manuscript.
 - `rivers4.csv` contains inputs to function `extract_river()` for the four basins used in Fig. 3.
 - `results_4rivers.rda` stores river data for Fig. 3. It is produced by `test_rivnet.R`.
 - `ARAcoord.csv`: contains data on wastewater treatment plants for the Thur catchment, and is used in Fig. 4.
 
 Note that the landcover map of Switzerland must be downloaded in order for the script to produce Fig. 4 (see lines 187-188 in `test_rivnet.R`).
