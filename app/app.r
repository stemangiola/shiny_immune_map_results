library(shiny)
library(gganatogram)
library(shinyBS)
library(shinyWidgets)
library(DT)

# Factors for easier table filtering
prop_data <- readRDS("cell_proportion_for_shiny_app.rds")
prop_data$cell_type_unified_ensemble <- factor(prop_data$cell_type_unified_ensemble)
prop_data$ethnicity_groups <- factor(prop_data$ethnicity_groups)
prop_data$sex <- factor(prop_data$sex)
prop_data$tissue_groups <- factor(prop_data$tissue_groups)

# Need to connect tissue_groups with appropriate anatomical organs
male_organ_list <- list(
  "blood" = c("leukocyte"),
  #"epithelium and mucosal tissues" = c("nose", "pleura", "endometrium", "throat", "nasal_pharynx", "nasal_septum", "esophagus"),
  "epithelium and mucosal tissues" = c("pleura"),
  "large intestine" = c("caecum", "rectum", "colon"),
  "respiratory system" = c("lung", "bronchus", "pulmonary_valve", "diaphragm"),
  "spleen" = c("spleen"),
  "gastrointestinal accessory organs" = c("liver", "pancreas", "gall_bladder"),
  "lymphatic system" = c("lymph_node", "tonsil"),
  "bone marrow" = c("bone_marrow"),
  "small intestine" = c("small_intestine", "ileum", "duodenum"),
  "renal system" = c("kidney", "renal_cortex", "urinary_bladder"),
  "nasal, oral, and pharyngeal regions" = c("nose", "nasal_pharynx", "throat", "tongue", "salivary_gland", "parotid_gland", "submandibular_gland"),
  "integumentary system (skin)" = c("skin"),
  "cerebral lobes and cortical areas" = c("frontal_cortex", "temporal_lobe", "prefrontal_cortex"),
  "endocrine system" = c("thyroid_gland", "pituitary_gland", "adrenal_gland"),
  "prostate" = c("prostate"),
  "trachea" = c("trachea"),
  "sensory-related structures" = c("nerve", "amygdala"),
  "stomach" = c("stomach"),
  "oesophagus" = c("esophagus", "gastroesophageal_junction"),
  "cardiovascular system" = c("heart", "aorta", "left_ventricle", "mitral_valve", "tricuspid_valve", "atrial_appendage", "coronary_artery"),
  "adipose tissue" = c("adipose_tissue"),
  "brainstem and cerebellar structures" = c("brain", "spinal_cord", "cerebellum", "cerebellar_hemisphere")
)

female_organ_list <- list(
  "blood" = c("leukocyte"),
  #"epithelium and mucosal tissues" = c("nose", "pleura", "endometrium", "throat", "nasal_pharynx", "nasal_septum", "esophagus"),
  "epithelium and mucosal tissues" = c("pleura"),
  "large intestine" = c("caecum", "rectum", "colon"),
  "respiratory system" = c("lung", "bronchus", "pulmonary_valve", "diaphragm"),
  "spleen" = c("spleen"),
  "gastrointestinal accessory organs" = c("liver", "pancreas", "gall_bladder"),
  "lymphatic system" = c("lymph_node", "tonsil"),
  "bone marrow" = c("bone_marrow"),
  "small intestine" = c("small_intestine", "ileum", "duodenum"),
  "renal system" = c("kidney", "renal_cortex", "urinary_bladder"),
  "nasal, oral, and pharyngeal regions" = c("nose", "nasal_pharynx", "nasal_septum", "throat", "salivary_gland", "parotid_gland", "submandibular_gland"),
  "integumentary system (skin)" = c("skin"),
  "cerebral lobes and cortical areas" = c("frontal_cortex", "temporal_lobe", "prefrontal_cortex"),
  "endocrine system" = c("thyroid_gland", "pituitary_gland", "adrenal_gland"),
  "breast" = c("breast"),
  "trachea" = c("trachea"),
  "sensory-related structures" = c("nerve"),
  "stomach" = c("stomach"),
  "oesophagus" = c("esophagus", "gastroesophageal_junction"),
  "female reproductive system" = c("ectocervix", "uterine_cervix", "vagina", "ovary", "uterus", "fallopian_tube", "endometrium"),
  "cardiovascular system" = c("heart", "aorta", "left_ventricle", "mitral_valve", "tricuspid_valve", "atrial_appendage", "coronary_artery"),
  "adipose tissue" = c("adipose_tissue"),
  "brainstem and cerebellar structures" = c("brain", "spinal_cord", "cerebellum", "cerebellar_hemisphere")
)

make_plot_data <- function(data, celltype, male_organ_map, female_organ_map, age = NULL, ethnicity = NULL, sex = "male") {
    if (!is.null(age)) {
        out <- data[data$age_bin_sex_specific == age, ]
    }

    if (!is.null(ethnicity)) {
        out <- data[data$ethnicity_groups == ethnicity, ]
    }

    out <- data[data$cell_type_unified_ensemble == celltype & data$sex == sex, ]

    if (sex == "male") {
        # Get value for each organ based on the matching name of the organ group it is a part of
        out_df <- data.frame(organ = unlist(male_organ_map))
        out_df$tissue_groups <- gsub("\\d", "", rownames(out_df))
    } else {
        out_df <- data.frame(organ = unlist(female_organ_map))
        out_df$tissue_groups <- gsub("\\d", "", rownames(out_df))
    }

    out_df$value <- out$proportion_mean[match(out_df$tissue_groups, out$tissue_groups)]

    out_df
}

ui <- navbarPage(
    "Immune Proportion Map",
    tabPanel(
        "Immune Proportions",
        sidebarLayout(
            sidebarPanel(width = 2,
                h3("Data Selection"),
                selectInput("cell_type", "Cell Type", choices = unique(prop_data$cell_type_unified_ensemble)),
                selectInput("ethnicity", "Ethnicity", choices = unique(prop_data$ethnicity_groups)),
                selectInput("age", "Age", choices = unique(prop_data$age_bin_sex_specific))
            ),
            mainPanel(
                fluidRow(
                    column(6, 
                        h3("Male Proportions"),
                        plotOutput("male_anatogram", width = "400px", height = "600px")
                    ),
                    column(6,
                        h3("Female Proportions"),
                        plotOutput("female_anatogram", width = "400px", height = "600px")
                    )
                )
            )
        )
    ),
    tabPanel(
        "Immune Proportions Table",
        br(),
        div(DTOutput("full_data"), style = "font-size:80%;")
    ),
    tabPanel(
        "Gene Expression",
        sidebarLayout(
            sidebarPanel(
                h3("Settings"),
                p("placeholder")
            ),
            mainPanel(
                h3("placeholder"),
                p("content")
            )
        )
    )
)

server <- function(input, output) {
    output$full_data <- renderDT({
        datatable(prop_data,
            rownames = FALSE,
            filter = "top",
            extensions = c("Buttons"),
            caption = "Full Immune Proportion Data",
            options = list(
                search = list(regex = TRUE),
                pageLength = 10,
                dom = "Blfrtip",
                buttons = c("copy", "csv", "excel", "pdf", "print")
            )
        ) %>% formatStyle(0, target = "row", lineHeight = "50%")
    })

    male_data <- reactive({make_plot_data(prop_data,
            celltype = input$cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            age = input$age,
            ethnicity = input$ethnicity,
            sex = "male"
        )
    })

    female_data <- reactive({make_plot_data(prop_data,
            celltype = input$cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            age = input$age,
            ethnicity = input$ethnicity,
            sex = "female"
        )
    })

    output$male_anatogram <- renderPlot({
        gganatogram(data = male_data(), sex = "male", fill = "value", organism = "human") + theme_void()
    })

    output$female_anatogram <- renderPlot({
        gganatogram(data = female_data(), sex = "female", fill = "value", organism = "human") + theme_void()
    })
}

shinyApp(ui = ui, server = server)
