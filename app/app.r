library(shiny)
library(gganatogram)
library(shinyBS)
library(shinyWidgets)
library(DT)
library(viridis)
library(colourpicker)
library(svglite)

# Factors for easier table filtering
prop_data <- readRDS("cell_proportion_for_shiny_app.rds")
prop_data$cell_type_unified_ensemble <- factor(prop_data$cell_type_unified_ensemble)
prop_data$ethnicity_groups <- factor(prop_data$ethnicity_groups)
prop_data$sex <- factor(prop_data$sex)
prop_data$tissue_groups <- factor(prop_data$tissue_groups)

# Need to connect tissue_groups with appropriate anatomical organs
male_organ_list <- list(
    "blood" = c("leukocyte"),
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
    "cerebral lobes and cortical areas" = c("brain", "frontal_cortex", "temporal_lobe", "prefrontal_cortex"),
    "endocrine system" = c("thyroid_gland", "pituitary_gland", "adrenal_gland"),
    "prostate" = c("prostate"),
    "trachea" = c("trachea"),
    "sensory-related structures" = c("nerve", "amygdala"),
    "stomach" = c("stomach"),
    "oesophagus" = c("esophagus", "gastroesophageal_junction"),
    "cardiovascular system" = c("heart", "aorta", "left_ventricle", "mitral_valve", "tricuspid_valve", "atrial_appendage", "coronary_artery"),
    "adipose tissue" = c("adipose_tissue"),
    "brainstem and cerebellar structures" = c("spinal_cord", "cerebellum", "cerebellar_hemisphere")
)

female_organ_list <- list(
    "blood" = c("leukocyte"),
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
    "cerebral lobes and cortical areas" = c("brain", "frontal_cortex", "temporal_lobe", "prefrontal_cortex"),
    "endocrine system" = c("thyroid_gland", "pituitary_gland", "adrenal_gland"),
    "breast" = c("breast"),
    "trachea" = c("trachea"),
    "sensory-related structures" = c("nerve"),
    "stomach" = c("stomach"),
    "oesophagus" = c("esophagus", "gastroesophageal_junction"),
    "female reproductive system" = c("ectocervix", "uterine_cervix", "vagina", "ovary", "uterus", "fallopian_tube", "endometrium"),
    "cardiovascular system" = c("heart", "aorta", "left_ventricle", "mitral_valve", "tricuspid_valve", "atrial_appendage", "coronary_artery"),
    "adipose tissue" = c("adipose_tissue"),
    "brainstem and cerebellar structures" = c("spinal_cord", "cerebellum", "cerebellar_hemisphere")
)

make_plot_data <- function(data, celltype, male_organ_map, female_organ_map, tissue_groups = NULL, age = NULL, ethnicity = NULL, sex = "male") {
    out <- data

    if (!is.null(age)) {
        out <- out[out$age_bin_sex_specific == age, ]
    }

    if (!is.null(ethnicity)) {
        out <- out[out$ethnicity_groups == ethnicity, ]
    }

    out <- out[out$cell_type_unified_ensemble == celltype & out$sex == sex, ]

    if (sex == "male") {
        # Get value for each organ based on the matching name of the organ group it is a part of
        out_df <- data.frame(organ = unlist(male_organ_map))

        # Order organs so plot order puts smaller stuff on top, generally
        out_df$organ <- factor(out_df$organ, levels = hgMale_key$organ)

        out_df$tissue_groups <- gsub("\\d", "", rownames(out_df))
    } else {
        out_df <- data.frame(organ = unlist(female_organ_map))
        out_df$organ <- factor(out_df$organ, levels = hgFemale_key$organ)

        out_df$tissue_groups <- gsub("\\d", "", rownames(out_df))
    }

    out_df$value <- out$proportion_mean[match(out_df$tissue_groups, out$tissue_groups)]
    out_df$CI_lower <- out$proportion_lower[match(out_df$tissue_groups, out$tissue_groups)]
    out_df$CI_upper <- out$proportion_upper[match(out_df$tissue_groups, out$tissue_groups)]

    if (!is.null(tissue_groups)) {
        out_df <- out_df[out_df$tissue_groups %in% tissue_groups, ]
    }

    out_df <- with(out_df, out_df[order(organ), ])

    out_df
}

ui <- navbarPage(
    "Immune Proportion Map",
    header = tags$head(
        # Note the wrapping of the string in HTML()
        tags$style(HTML("
        hr {
            margin-top: 10px;
            margin-bottom: 10px;
        }
        h3 {
            margin-top: 5px;
            margin-bottom: 5px;
        }"))
    ),
    tabPanel(
        "Proportions by Sex",
        sidebarLayout(
            sidebarPanel(
                width = 3,
                h4("Data Selection"),
                selectInput("t1_cell_type", "Cell Type",
                    choices = unique(prop_data$cell_type_unified_ensemble)
                ),
                selectInput("t1_ethnicity", "Ethnicity",
                    choices = unique(prop_data$ethnicity_groups),
                    selected = "European"
                ),
                sliderTextInput("t1_age", "Age",
                    choices = c("Infancy", "Childhood", "Adolescence", "Young Adulthood", "Middle Age", "Senior"),
                    grid = TRUE,
                    selected = "Adolescence"
                ),
                pickerInput("t1_tissue_groups", "Tissue Groups",
                    choices = unique(prop_data$tissue_groups),
                    multiple = TRUE,
                    selected = unique(prop_data$tissue_groups),
                    options = pickerOptions(
                        actionsBox = TRUE,
                        size = 10,
                        selectedTextFormat = "count > 3"
                    )
                ),
                hr(),
                h4("Plot Aesthetics"),
                fluidRow(
                    column(
                        6,
                        colourInput("t1_outline_colour", "Outline Colour", value = "lightgray"),
                        prettyCheckbox("t1_outline", "Plot Outline", value = TRUE),
                        selectInput("t1_dl_type", "Download Format",
                            choices = c(".png", ".pdf", ".svg")
                        )
                    ),
                    column(
                        6,
                        selectInput("t1_palette", "Palette",
                            choices = c("viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo")
                        ),
                        numericInput("t1_opacity", "Opacity",
                            value = 0.51, min = 0.01, max = 1, step = 0.1
                        ),
                        prettyCheckbox("t1_reverse", "Reverse Palette", value = FALSE)
                    )
                ),
                div(style = "text-align: center;", actionButton("t1_update", "Update Plots"))
            ),
            mainPanel(
                width = 9,
                fluidRow(
                    column(
                        6,
                        h3("Male Proportions"),
                        plotOutput("t1_male_anatogram", width = "320px", height = "450px"),
                        downloadButton("t1_male_dl", "Download Plot"),
                        hr(),
                        div(DTOutput("t1_male_props"), style = "font-size:70%;")
                    ),
                    column(
                        6,
                        h3("Female Proportions"),
                        plotOutput("t1_female_anatogram", width = "320px", height = "450px"),
                        downloadButton("t1_female_dl", "Download Plot"),
                        hr(),
                        div(DTOutput("t1_female_props"), style = "font-size:70%;")
                    )
                )
            )
        )
    ),
    tabPanel(
        "Immune Proportions Table",
        br(),
        div(DTOutput("full_data"), style = "font-size:70%;")
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
        ) %>%
            formatStyle(0, target = "row", lineHeight = "50%") %>%
            formatRound(c("value", "CI_lower", "CI_upper"), 4)
    })

    t1_male_data <- reactive({
        input$t1_update

        make_plot_data(prop_data,
            celltype = isolate(input$t1_cell_type),
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = isolate(input$t1_tissue_groups),
            age = isolate(input$t1_age),
            ethnicity = isolate(input$t1_ethnicity),
            sex = "male"
        )
    })

    t1_female_data <- reactive({
        input$t1_update

        make_plot_data(prop_data,
            celltype = isolate(input$t1_cell_type),
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = isolate(input$t1_tissue_groups),
            age = isolate(input$t1_age),
            ethnicity = isolate(input$t1_ethnicity),
            sex = "female"
        )
    })

    t1_max_val <- reactive({
        max(t1_female_data()$value, t1_male_data()$value, na.rm = TRUE)
    })

    t1_male_plot <- reactive({
        direc <- ifelse(isolate(input$t1_reverse), -1, 1)

        p <- gganatogram(
            data = t1_male_data(), sex = "male", fill = "value",
            organism = "human", outline = isolate(input$t1_outline),
            fillOutline = isolate(input$t1_outline_colour),
        ) + theme_void()

        p + scale_fill_viridis(
            option = isolate(input$t1_palette),
            alpha = isolate(input$t1_opacity),
            direction = direc,
            limits = c(0, t1_max_val())
        )
    })

    output$t1_male_anatogram <- renderPlot({
        t1_male_plot()
    })

    output$t1_male_dl <- downloadHandler(
        filename = function() { paste0("male_anatogram", isolate(input$t1_dl_type)) },
        content = function(file) {
            ggsave(file, t1_male_plot(), width = 6)
        }
    )

    output$t1_male_props <- renderDT(
        {
            input$update

            datatable(isolate(t1_male_data()),
                rownames = FALSE,
                filter = "top",
                extensions = c("Buttons"),
                options = list(
                    search = list(regex = TRUE),
                    pageLength = 10,
                    dom = "Blfrtip",
                    buttons = c("copy", "csv", "excel", "pdf", "print"),
                    autoWidth = FALSE,
                    columnDefs = list(list(width = "30%", targets = 1))
                )
            ) %>%
                formatStyle(0, target = "row", lineHeight = "50%") %>%
                formatRound(c("value", "CI_lower", "CI_upper"), 4)
        }
    )

    t1_female_plot <- reactive({
        direc <- ifelse(isolate(input$t1_reverse), -1, 1)

        p <- gganatogram(
            data = t1_female_data(), sex = "female", fill = "value",
            organism = "human", outline = isolate(input$t1_outline),
            fillOutline = isolate(input$t1_outline_colour),
        ) + theme_void()

        p + scale_fill_viridis(
            option = isolate(input$t1_palette),
            alpha = isolate(input$t1_opacity),
            direction = direc,
            limits = c(0, t1_max_val())
        )
    })

    output$t1_female_anatogram <- renderPlot({
        t1_female_plot()
    })

    output$t1_female_dl <- downloadHandler(
        filename = function() { paste0("female_anatogram", isolate(input$t1_dl_type)) },
        content = function(file) {
            ggsave(file, t1_female_plot(), width = 6)
        }
    )

    output$t1_female_props <- renderDT(
        {
            input$update

            datatable(isolate(t1_female_data()),
                rownames = FALSE,
                filter = "top",
                extensions = c("Buttons"),
                options = list(
                    search = list(regex = TRUE),
                    pageLength = 10,
                    dom = "Blfrtip",
                    buttons = c("copy", "csv", "excel", "pdf", "print"),
                    autoWidth = FALSE,
                    columnDefs = list(list(width = "30%", targets = 1))
                )
            ) %>%
                formatStyle(0, target = "row", lineHeight = "50%") %>%
                formatRound(c("value", "CI_lower", "CI_upper"), 4)
        }
    )
}

shinyApp(ui = ui, server = server)
