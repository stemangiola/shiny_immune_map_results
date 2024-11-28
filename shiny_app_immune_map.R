# Prepare data for shiny APP

library(tidyverse)
library(targets)
library(tarchetypes)
library(sccomp)


tar_script({
  
  library(tidyverse)
  library(targets)
  library(tarchetypes)
  library(glue)
  library(crew)
  library(crew.cluster)
  
  age_range_sex_specific <- function(age_bins, sex) {
    # Initialise an empty vector to store the age ranges in years
    age_ranges <- integer(length(age_bins))
    
    # Define age thresholds for male, female, and unknown
    male_thresholds <- c(3, 13, 21, 40, 55)
    female_thresholds <- c(3, 13, 19, 36, 50)
    unknown_thresholds <- c(3, 13, 20, 38, 52)
    
    # Map age groups to approximate age ranges for each sex category
    for (i in seq_along(age_bins)) {
      if (sex[i] == "male") {
        age_ranges[i] <- dplyr::case_when(
          age_bins[i] == "Infancy" ~ male_thresholds[1] / 2,
          age_bins[i] == "Childhood" ~ (male_thresholds[1] + male_thresholds[2]) / 2,
          age_bins[i] == "Adolescence" ~ (male_thresholds[2] + male_thresholds[3]) / 2,
          age_bins[i] == "Young Adulthood" ~ (male_thresholds[3] + male_thresholds[4]) / 2,
          age_bins[i] == "Middle Age" ~ (male_thresholds[4] + 55) / 2,
          age_bins[i] == "Senior" ~ 55,
          TRUE ~ NA_integer_
        )
      } else if (sex[i] == "female") {
        age_ranges[i] <- dplyr::case_when(
          age_bins[i] == "Infancy" ~ female_thresholds[1] / 2,
          age_bins[i] == "Childhood" ~ (female_thresholds[1] + female_thresholds[2]) / 2,
          age_bins[i] == "Adolescence" ~ (female_thresholds[2] + female_thresholds[3]) / 2,
          age_bins[i] == "Young Adulthood" ~ (female_thresholds[3] + female_thresholds[4]) / 2,
          age_bins[i] == "Middle Age" ~ (female_thresholds[4] + 50) / 2,
          age_bins[i] == "Senior" ~ 50,
          TRUE ~ NA_integer_
        )
      } else if (sex[i] == "unknown") {
        age_ranges[i] <- dplyr::case_when(
          age_bins[i] == "Infancy" ~ unknown_thresholds[1] / 2,
          age_bins[i] == "Childhood" ~ (unknown_thresholds[1] + unknown_thresholds[2]) / 2,
          age_bins[i] == "Adolescence" ~ (unknown_thresholds[2] + unknown_thresholds[3]) / 2,
          age_bins[i] == "Young Adulthood" ~ (unknown_thresholds[3] + unknown_thresholds[4]) / 2,
          age_bins[i] == "Middle Age" ~ (unknown_thresholds[4] + 52) / 2,
          age_bins[i] == "Senior" ~ 52,
          TRUE ~ NA_integer_
        )
      } else {
        stop("Each element of 'sex' must be either 'male', 'female', or 'unknown'.")
      }
    }
    
    return(age_ranges)
  }
  
  tar_option_set(
    
    memory = "transient", 
    garbage_collection = 100, 
    storage = "worker", 
    retrieval = "worker", 
    error = "continue", 
    # cue = tar_cue(mode = "never"), 
    workspace_on_error = TRUE,
    format = "qs",
    
    #-----------------------#
    # SLURM
    #-----------------------#
    controller = crew_controller_group(
      
      crew_controller_slurm(
        name = "slurm_160",
        workers = 20,
        tasks_max = 20,
        seconds_idle = 30,
        crashes_error = 10,
        options_cluster = crew_options_slurm(
          memory_gigabytes_required = c(160), 
          cpus_per_task = c(2), 
          time_minutes = c(60*8),
          verbose = T
        )
      )
    ),
    # debug = "metadata_grouped_pseudobulk_processed_split_1",
    
  )
  
  list(
    tarchetypes::tar_file(estimates_file, "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/estimates_continuous_age_plus_age_bins.rds"),
    tarchetypes::tar_file(input_data_file, "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/cell_metadata_1_0_1_sccomp_input_counts.rds"),
    tar_target(input_data, readRDS(input_data_file)),
    tarchetypes::tar_group_by(
      formulae_and_new_data,
      { 
        estimate = readRDS(estimates_file)
        
        tribble(
          ~ name, ~ formula_composition, ~ new_data, ~ contrasts,
          
          # Tissue level, age
          "age",
          "~ 1 + age_days_scaled + age_bin_sex_specific + (1 + age_days_scaled + age_bin_sex_specific | tissue_groups)",  
          input_data |> 
            distinct(tissue_groups) |> 
            expand_grid(age_bin_sex_specific = input_data |> pull(age_bin_sex_specific) |> unique()) |> 
            mutate(sex = "unknown") |> 
            mutate(age_days = age_range_sex_specific(age_bin_sex_specific, sex)) |> 
            mutate(age_days_scaled = 
                     ( age_days -  (20 * 365)  ) / 
                     (input_data |> pull(age_days) |> sd())
            ) |>           
            mutate(n = 5000) |> 
            rowid_to_column("sample_id") |> 
            mutate(sample_id = sample_id |> as.character()) |> 
            
            # Drop combinations not present in the model
            inner_join(input_data |> distinct(tissue_groups, age_bin_sex_specific)) |> 
            mutate(ethnicity = "european", sex = "female"),
          
          # Contrasts
          inner_join(
            estimate |>
              filter(parameter |> str_detect("age_days_scaled___")) |>
              distinct(parameter) |>
              tidyr::extract(parameter, "tissue_groups", ".+___(.+)", remove = FALSE) |>
              extract(parameter, "fixed_continuous", "(.*)___.*", remove = FALSE) |> 
              dplyr::rename(parameter_continuous = parameter),
            
            estimate |>
              filter(parameter |> str_detect("age_bin_sex_specific[a-zA-Z ]+___")) |>
              distinct(parameter) |>
              tidyr::extract(parameter, "tissue_groups", ".+___(.+)", remove = FALSE) |>
              tidyr::extract(parameter, c("age_bin_groups"), "age_bin_sex_specific([a-zA-Z ]+)___.+", remove = FALSE) |>
              
              extract(parameter, "fixed_bins", "(.*)___.*", remove = FALSE) |> 
              dplyr::rename(parameter_bins = parameter) 
          ) |> 
            
            mutate(contrast = glue("`{fixed_continuous}` + `{parameter_continuous}` + `{fixed_bins}` + `{parameter_bins}`") |> as.character()) |>
            unite("age_bin_tissue_groups", age_bin_groups, tissue_groups) |> 
            select(age_bin_tissue_groups, contrast) |> 
            deframe( ),
          
          # Tissue level, sex
          "sex",
          "~ 1 + sex + (1 + sex | tissue_groups)",  
          input_data |> 
            distinct(tissue_groups) |> 
            expand_grid(sex = c("male", "female")) |> 
            
            mutate(n = 5000) |> 
            rowid_to_column("sample_id") |> 
            mutate(sample_id = sample_id |> as.character()) |> 
            
            # Drop combinations not present in the model
            inner_join(input_data |> distinct(tissue_groups, sex)) |> 
            mutate(ethnicity = "European", age_bin_sex_specific = "Adolescence"),
          
          # Contrasts
          estimate |>
            filter(parameter |> str_detect("^sexmale___")) |>
            distinct(parameter) |>
            tidyr::extract(parameter, "tissue_groups", ".+___(.+)", remove = FALSE) |>
            extract(parameter, "fixed", "(.*)___.*", remove = FALSE) |> 
            dplyr::rename(parameter = parameter) |> 
            
            mutate(contrast = glue("`{fixed}` + `{parameter}` ") |> as.character()) |>
            select(tissue_groups, contrast) |> 
            deframe( ),
          
          # Tissue level, ethnicity
          "ethnicity",
          "~ 1 + ethnicity_groups + (1 + ethnicity_groups | tissue_groups)",  
          input_data |> 
            distinct(tissue_groups) |> 
            expand_grid(ethnicity_groups = input_data |> pull(ethnicity_groups) |> unique()) |> 
            
            mutate(n = 5000) |> 
            rowid_to_column("sample_id") |> 
            mutate(sample_id = sample_id |> as.character()) |> 
            
            # Drop combinations not present in the model
            inner_join(input_data |> distinct(tissue_groups, ethnicity_groups)) ,
          
          # Contrasts
          estimate |>
            filter(parameter |> str_detect("^ethnicity_groups[a-zA-Z ]+___")) |>
            distinct(parameter) |>
            tidyr::extract(parameter, c("ethnicity_groups_tissue_groups"), "ethnicity_groups([a-zA-Z ]+___.+)", remove = FALSE) |>
            extract(parameter, "fixed", "(.*)___.*", remove = FALSE) |> 
            dplyr::rename(parameter = parameter) |> 
            
            mutate(contrast = glue("`{fixed}` + `{parameter}` ") |> as.character()) |>
            select(ethnicity_groups_tissue_groups, contrast) |> 
            deframe( ),
          
          # Tissue level, age + sex
          "age + sex",
          "~ 1 + age_days_scaled + age_bin_sex_specific*sex + (1 + age_days_scaled + age_bin_sex_specific*sex | tissue_groups)",  
          input_data |> 
            distinct(tissue_groups) |> 
            expand_grid(age_bin_sex_specific = input_data |> pull(age_bin_sex_specific) |> unique()) |> 
            expand_grid(sex = c("male", "female")) |> 
            mutate(age_days = age_range_sex_specific(age_bin_sex_specific, sex)) |> 
            #left_join(map_age_to_scaled_age)
            mutate(age_days_scaled = 
                     ( age_days -  (20 * 365)  ) / 
                     (input_data |> pull(age_days) |> sd())
            ) |> 
            mutate(n = 5000) |> 
            rowid_to_column("sample_id") |> 
            mutate(sample_id = sample_id |> as.character()) |> 
            
            # Drop combinations not present in the model
            inner_join(input_data |> distinct(tissue_groups, age_bin_sex_specific, sex)) ,
          
          # Contrasts
          list(estimate |>
              filter(parameter |> str_detect("age_days_scaled___")) |>
              distinct(parameter) |>
              tidyr::extract(parameter, "tissue_groups", ".+___(.+)", remove = FALSE) |>
              extract(parameter, "fixed_continuous", "(.*)___.*", remove = FALSE) |> 
              dplyr::rename(parameter_continuous = parameter),
            
            estimate |>
              filter(parameter |> str_detect("age_bin_sex_specific[a-zA-Z ]+___")) |>
              distinct(parameter) |>
              tidyr::extract(parameter, "tissue_groups", ".+___(.+)", remove = FALSE) |>
              tidyr::extract(parameter, c("age_bin_groups"), "age_bin_sex_specific([a-zA-Z ]+)___.+", remove = FALSE) |>
              
              extract(parameter, "fixed_bins", "(.*)___.*", remove = FALSE) |> 
              dplyr::rename(parameter_bins = parameter) ,
            
            estimate |>
              filter(parameter |> str_detect("age_bin_sex_specific[a-zA-Z ]+:sexmale___")) |>
              distinct(parameter) |>
              tidyr::extract(parameter, "tissue_groups", ".+___(.+)", remove = FALSE) |>
              tidyr::extract(parameter, c("age_bin_groups"), "age_bin_sex_specific([a-zA-Z ]+):sexmale___.+", remove = FALSE) |>
              
              extract(parameter, "fixed_bins_interaction", "(.*)___.*", remove = FALSE) |> 
              dplyr::rename(parameter_bins_interaction = parameter) 
          ) |>
            reduce(inner_join) |> 
            
            mutate(contrast = glue("`{fixed_continuous}` + `{parameter_continuous}` + `{fixed_bins}` + `{parameter_bins}` + `{fixed_bins_interaction}` + `{parameter_bins_interaction}`") |> as.character()) |>
            mutate(age_bin_tissue_groups = glue("{age_bin_groups} {tissue_groups} sexmale")) |> 
            select(age_bin_tissue_groups, contrast) |> 
            deframe( ),
          
          # Tissue level, age + ethnicity
          "age + ethnicity",
          "~ 1 + age_days_scaled + age_bin_sex_specific + ethnicity_groups + (1 + age_days_scaled + age_bin_sex_specific  + ethnicity_groups | tissue_groups)",  
          input_data |> 
            distinct(tissue_groups) |> 
            expand_grid(age_bin_sex_specific = input_data |> pull(age_bin_sex_specific) |> unique()) |> 
            expand_grid(ethnicity_groups = input_data |> pull(ethnicity_groups) |> unique()) |> 
            mutate(sex = "unknown") |> 
            mutate(age_days = age_range_sex_specific(age_bin_sex_specific, sex)) |> 
            mutate(age_days_scaled = 
                     ( age_days -  (20 * 365)  ) / 
                     (input_data |> pull(age_days) |> sd())
            ) |> 
            mutate(n = 5000) |> 
            rowid_to_column("sample_id") |> 
            mutate(sample_id = sample_id |> as.character()) |> 
            
            # Drop combinations not present in the model
            inner_join(input_data |> distinct(tissue_groups, age_bin_sex_specific, ethnicity_groups)) ,
          
          # Contrasts
          NA_character_,
          
          # Tissue level, sex + ethnicity
          "sex + ethnicity",
          "~ 1 + sex + ethnicity_groups + (1 + sex + ethnicity_groups | tissue_groups)",  
          input_data |> 
            distinct(tissue_groups) |> 
            expand_grid(age_bin_sex_specific = input_data |> pull(age_bin_sex_specific) |> unique()) |> 
            expand_grid(ethnicity_groups = input_data |> pull(ethnicity_groups) |> unique()) |> 
            mutate(sex = "unknown") |> 
            mutate(age_days = age_range_sex_specific(age_bin_sex_specific, sex)) |> 
            mutate(age_days_scaled = 
                     ( age_days -  (20 * 365)  ) / 
                     (input_data |> pull(age_days) |> sd())
            ) |> 
            mutate(n = 5000) |> 
            rowid_to_column("sample_id") |> 
            mutate(sample_id = sample_id |> as.character()) |> 
            
            # Drop combinations not present in the model
            inner_join(input_data |> distinct(tissue_groups, sex, age_bin_sex_specific, ethnicity_groups)) ,
          
          # Contrasts
          NA_character_
          
        ) },
      formula_composition,
      packages = c("dplyr", "tibble", "tidyr", "glue"), 
      deployment = "main"
    ),
    tar_target(
      predictions,
      formulae_and_new_data |> mutate(prediction = map2(
        formula_composition, new_data,
        ~ estimates_file |> 
          readRDS() |> 
          sccomp_predict(
            formula_composition = .x |> as.formula(),
            new_data = .y,
            mcmc_seed = 42,
            summary_instead_of_draws = TRUE
          )                               
        
      )) |> 
        
        # Make the object lighter
        select(-new_data),
      pattern = map(formulae_and_new_data),
      resources = tar_resources(crew = tar_resources_crew("slurm_160")), 
      packages = "sccomp"
    ),
    
    tar_target(
      hypothesis_testing,
      formulae_and_new_data |> mutate(prediction = map(
        contrasts, 
        ~ estimates_file |> 
          readRDS() |> 
          sccomp_test(contrasts = .x)                               
        
      )) |> 
        
        # Make the object lighter
        select(-new_data),
      pattern = map(formulae_and_new_data),
      resources = tar_resources(crew = tar_resources_crew("slurm_160")), 
      packages = "sccomp"
    )
  )
},
script = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/_targets_shiny.R",
ask=FALSE
)


job::job({
  
  tar_make(
    script = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/_targets_shiny.R", 
    store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/_targets_shiny", 
    reporter = "verbose"
  )
  
})

tar_read(
  predictions,
  store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/_targets_shiny", 
) |> 
  mutate(name = c("age" ,    "age + ethnicity",     "age + sex" ,           "ethnicity"  ,   "sex"      ,      "sex + ethnicity") ) |> 
  select(name, everything()) |> 
  saveRDS("~/PostDoc/immuneHealthyBodyMap/cell_proportion_for_shiny_app.rds", compress = "xz")

system("~/bin/rclone copy /vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1/cell_metadata_1_0_1_sccomp_input.rds box_adelaide:/Mangiola_ImmuneAtlas/taskforce_shared_folder/")


tar_workspace(
  predictions_f5799da0fa28762a,
  script = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/_targets_shiny.R", 
  store = "/vast/projects/mangiola_immune_map/PostDoc/immuneHealthyBodyMap/sccomp_on_cellNexus_1_0_1_regularised/_targets_shiny", 
)
