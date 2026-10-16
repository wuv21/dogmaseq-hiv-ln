# Non-follicular T cells predominate in the lymph node HIV reservoir during antiretroviral therapy

This is the repo for the downstream analysis of the scDOGMAseq profiling of lymph node mononuclear cells from people without HIV and people with HIV (both viremic and ART-treated).

The raw and processed data can be downloaded in the NCBI GEO repository under accension number `GSE313072`.

To recreate the figures, the `main_preprocessing.R` and `main_downstream.R` files need to be run first within this R project (with the associated renv environment). In the `figures/` directory, all figures and all supplemental figures (with the exception of Supplemental Figure 1 and subplots which are generated with BioRender, FlowJo, or Illustrator). Note that there was an addition of a new supplemental figure during revisions so the prefix naming is changed starting with Supplemental Figure 3. Please see the below list for the script:figure pairing.

  - `figure1.R` - Figure 1, Supplemental Figure 2, and Supplemental Table 2
  - `figure2.R` - Figure 2 and Supplemental Figures 3,4 (note that these are saved with the prefix of supfig_3-1 and supfig_3-2, respectively)
  - `figure3.R` - Figure 3 and Supplemental Figures 5 (note that this is saved with the prefix of supfig_4)
  - `figure4.R` - Figure 4 and Supplemental Figures 6 (note that this is saved with the prefix of supfig_5)
  - `figure5.R` - Figure 5 and Supplemental Figures 7,8 (note that these are saved with the prefix of supfig_6 and supfig_7, respectively)
  
## Publication information
TBA

