# shiny_immune_map_results
This is a repository that includes the pilot data necessary to start building a shiny app to communicate the results of our immune map paper dynamically

The data is for cell proportion changes in the immune system throughout the body. 

This is the [preprint](https://www.biorxiv.org/content/10.1101/2023.06.08.542671v3) that we revamped and soon ready to resubmit.

The tasks are

- lists the interactive tables (one provided here, estimates_age_bins_effect_tibble_only.rds) of hypothesis test statistics. 
- has an interactive [anatogram](https://github.com/jespermaag/gganatogram) with slides to create a personalised medicine representation of immune cell composition as a body heatmap. This can be accompained with dynamic table which visualises the mean and error values of the cell proportion and gene RNA abundance.

 <img width="326" alt="image" src="https://github.com/user-attachments/assets/0bcd17e8-7665-423a-8993-9c0074068705">

 For this task [sccomp](https://github.com/MangiolaLaboratory/sccomp) backend is needed to make the proportion predictions, inputted by the sliders.

We will need to build an analogue interface for the gene expression changes.
