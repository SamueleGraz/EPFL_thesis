**EPFL_thesis**

[**Asaf Marco 's paper**]([URL]https://www.nature.com/articles/s41593-020-00717-0#data-availability) data analysis Linux part:
- Download the data -> 00_download.R
- STAR index and gene count -> 01_star.R
- Data wrangling in order to have a final table with the total counts per gene and condition -> 02_final_table.R
- deseq2 library usage -> 03_deseq2.R
- PCA and euclidean distances across the samples -> 04_sample_correlation.R
-  DEG according to specific contrasts. Here BasalVsEarly, BasalVsLate and BasalVsReactivated were done -> 05_DEGs_analysis.R
-  Random vulcano plot -> 06_vulcano_plots_RAW.R
- Dockerfile where the environments were created in order to work with conda in a reproducible way -> Dockerfile.txt
- folder with the environments created with conda -> envs
