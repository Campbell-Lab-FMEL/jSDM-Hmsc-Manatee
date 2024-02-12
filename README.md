[![DOI](https://zenodo.org/badge/639979997.svg)](https://zenodo.org/badge/latestdoi/639979997)

# jSDM-Hmsc-Manatee
This repository accompanies the manuscript **"Land use predicts proportion of West Nile virus vector-competent mosquitoes"** by Amely M. Bauer, Robert P. Guralnick, Shelley A. Whitehead, Narayani Barve, Julie M. Allen, and Lindsay P. Campbell, published in  Ecosphere. To read the open-access article, follow this link: [https://doi.org/10.1038/s41598-023-30751-4](https://doi.org/10.1002/ecs2.4771). 

We used a “Hierarchical model of species communities” (Hmsc) joint species distribution modeling (jSDM) approach that included mosquito species vector competency for West Nile virus (WNV) as a trait to understand community-level responses to land cover and predict joint species distributions in Manatee County, Florida. We assembled species presence/absence data of mosquito assamblages sampled from 2016 to 2020, and included percent land cover of five common land cover types as predictor variables to investigate patterns in species richness and community-weighted proportions of WNV-competent vector species. The discussion of our results highlights the value of community-level analyses to predict joint vector distributions that can inform where greatest transmission hazard may occur

## Code Files

The scripts document data preparation and analyses conducted for this study (using the `Hmsc` package [^1]) and contain example code for plotting model outputs. 

* **Wallace_Hmsc_Manatee_2016_2020_model_setup_Rscript.R**
  * Define and fit jSDM (run on HPC)
  
* **Wallace_Hmsc_Manatee_2016_2020_model_fit_Rscript.R** 
  * Evaluate model fit in terms of explanatory and predictive power (run on HPC)
 
* **Hmsc_Manatee_2016_2020_model_postprocessing_Rscript.R** 
  * Explore model convergence and fit
  * Post-processing steps: inspect model fit, compute variance partitioning, and make predictions 


<br/>
<br/>
<br/>

*This repository contains the version of the code and data files at the time of manuscript submission, which may undergo changes in future*

[^1]: Tikhonov G, Ovaskainen O, Oksanen J, de Jonge M, Opedal O, Dallas T (2022). Hmsc: Hierarchical Model of Species Communities. R package version 3.0-13, <https://CRAN.R-project.org/package=Hmsc>, <https://github.com/hmsc-r/HMSC>
