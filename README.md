# Shiny Immune Map Results

This repository contains pilot data and resources for building a **Shiny app** to communicate the results of our immune map paper dynamically. The aim is to create an interactive platform that visualises and explores cell proportion changes in the immune system across the body.

## About the Data

The provided data includes insights into immune cell proportions, and gene expression changes throughout the body. These data are foundational for developing a personalised medicine approach to immune system visualisation.

The project builds on findings from our preprint: [Preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3). We have revamped the manuscript and are preparing it for resubmission.

## Objectives and Tasks

1. **Interactive Tables**  
   - Develop tables to explore hypothesis test statistics.  
   - Example: The file `estimates_age_bins_effect_tibble_only.rds` provides an initial dataset for this feature.  

2. **Interactive Immune Composition Visualisation**  
   - Implement an interactive [anatogram](https://github.com/jespermaag/gganatogram).  
   - Users will be able to adjust (activate/deactivate) 3 sliders (age bins, sex, ethnicity) to create a personalised medicine representation of immune cell composition as a body heatmap.
   - The input data is [here](https://github.com/stemangiola/shiny_immune_map_results/blob/main/estimates_age_bins_effect_tibble_only.rds) 

   - **Additional Features:**  
     - Include a dynamic table displaying mean values and errors for cell proportions and RNA abundance.  
     - Use the [sccomp](https://github.com/MangiolaLaboratory/sccomp) backend to generate cell proportion predictions based on user inputs.  

   ![Body Heatmap Preview](https://github.com/user-attachments/assets/0bcd17e8-7665-423a-8993-9c0074068705)

3. **Gene Expression Interface**  
   - Design a parallel interface to visualise and interact with gene expression profiles and changes.  


### Prerequisites

- Familiarity with R and Shiny.
- Knowledge of `gganatogram` for anatomical visualisations.
- Understanding of the `sccomp` package for backend processing, and differential-expression models.

### How to Contribute

Please get in contact with us. 
