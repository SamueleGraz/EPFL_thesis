**EPFL_thesis**

*RNA_seq_from_Asaf_Marco_paper folder*

[**Asaf Marco 's paper**]([URL]https://www.nature.com/articles/s41593-020-00717-0#data-availability) data analysis Linux part:
- 00_download.R = Download the data
- 01_star.R = STAR index and gene count
- 02_final_table.R = Data wrangling in order to have a final table with the total counts per gene and condition
- 03_deseq2.R = deseq2 library usage
- 04_sample_correlation.R = PCA and euclidean distances across the samples
- 05_DEGs_analysis.R = DEG according to specific contrasts. Here BasalVsEarly, BasalVsLate and BasalVsReactivated were done
- 06_vulcano_plots_RAW.R = Random vulcano plot
- Dockerfile.txt = Dockerfile where the environments were created in order to work with conda in a reproducible way
- envs = folder with the environments created with conda

*marco_asaf's_data_comparison folder*

- script_marco_data_comparison.R = I used the filter log2fold > 2 and padj < 0.01 like what they claimed in the paper. I also made the right contrasts accordingly to what they did:
  + BasalVsEarly, EarlyVsLate and LateVsReactivated
