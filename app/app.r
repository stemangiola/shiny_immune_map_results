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
        }
        .shiny-split-layout {
            width: 98%;
        }
        "))
    ),
    tabPanel(
        "Comparisons by Sex",
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
                )
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
        "Comparisons by Age",
        sidebarLayout(
            sidebarPanel(
                width = 2,
                h4("Data Selection"),
                selectInput("t2_cell_type", "Cell Type",
                    choices = unique(prop_data$cell_type_unified_ensemble)
                ),
                selectInput("t2_ethnicity", "Ethnicity",
                    choices = unique(prop_data$ethnicity_groups),
                    selected = "European"
                ),
                pickerInput("t2_tissue_groups", "Tissue Groups",
                    choices = unique(prop_data$tissue_groups),
                    multiple = TRUE,
                    selected = unique(prop_data$tissue_groups),
                    options = pickerOptions(
                        actionsBox = TRUE,
                        size = 10,
                        selectedTextFormat = "count > 3"
                    )
                ),
                radioGroupButtons(
                    inputId = "t2_sex",
                    label = "Sex",
                    choices = c("male", "female")
                ),
                hr(),
                h4("Plot Aesthetics"),
                fluidRow(
                    column(
                        6,
                        colourInput("t2_outline_colour", "Outline Colour", value = "lightgray"),
                        prettyCheckbox("t2_outline", "Plot Outline", value = TRUE),
                        selectInput("t2_dl_type", "Download Format",
                            choices = c(".png", ".pdf", ".svg")
                        )
                    ),
                    column(
                        6,
                        selectInput("t2_palette", "Palette",
                            choices = c("viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo")
                        ),
                        numericInput("t2_opacity", "Opacity",
                            value = 0.51, min = 0.01, max = 1, step = 0.1
                        ),
                        prettyCheckbox("t2_reverse", "Reverse Palette", value = FALSE)
                    )
                )
            ),
            mainPanel(
                width = 10,
                fluidRow(
                    column(
                        2,
                        h3("Infancy"),
                        plotOutput("t2_infancy_anatogram", width = "250px", height = "360px"),
                        splitLayout(downloadButton("t2_infancy_dl", "Download"), actionButton("t2_infancy_show", "Show Data"))
                    ),
                    column(
                        2,
                        h3("Childhood"),
                        plotOutput("t2_childhood_anatogram", width = "250px", height = "360px"),
                        splitLayout(downloadButton("t2_childhood_dl", "Download"), actionButton("t2_childhood_show", "Show Data"))
                    ),
                    column(
                        2,
                        h3("Adolescence"),
                        plotOutput("t2_adolescence_anatogram", width = "250px", height = "360px"),
                        splitLayout(downloadButton("t2_adolescence_dl", "Download"), actionButton("t2_adolescence_show", "Show Data"))
                    ),
                    column(
                        2,
                        h3("Young Adulthood"),
                        plotOutput("t2_ya_anatogram", width = "250px", height = "360px"),
                        splitLayout(downloadButton("t2_ya_dl", "Download"), actionButton("t2_ya_show", "Show Data"))
                    ),
                    column(
                        2,
                        h3("Middle Age"),
                        plotOutput("t2_middleage_anatogram", width = "250px", height = "360px"),
                        splitLayout(downloadButton("t2_middleage_dl", "Download"), actionButton("t2_middleage_show", "Show Data"))
                    ),
                    column(
                        2,
                        h3("Senior"),
                        plotOutput("t2_senior_anatogram", width = "250px", height = "360px"),
                        splitLayout(downloadButton("t2_senior_dl", "Download"), actionButton("t2_senior_show", "Show Data"))
                    )
                ),
                fluidRow(
                    hr(),
                    div(DTOutput("t2_props"), style = "font-size:70%;")
                )
            )
        )
    ),
    tabPanel(
        "Comparisons by Ethnicity",
        sidebarLayout(
            sidebarPanel(
                width = 2,
                h4("Data Selection"),
                selectInput("t3_cell_type", "Cell Type",
                    choices = unique(prop_data$cell_type_unified_ensemble)
                ),
                selectInput("t3_age", "Age",
                    choices = unique(prop_data$age_bin_sex_specific),
                    selected = "European"
                ),
                pickerInput("t3_tissue_groups", "Tissue Groups",
                    choices = unique(prop_data$tissue_groups),
                    multiple = TRUE,
                    selected = unique(prop_data$tissue_groups),
                    options = pickerOptions(
                        actionsBox = TRUE,
                        size = 10,
                        selectedTextFormat = "count > 3"
                    )
                ),
                radioGroupButtons(
                    inputId = "t3_sex",
                    label = "Sex",
                    choices = c("male", "female")
                ),
                hr(),
                h4("Plot Aesthetics"),
                fluidRow(
                    column(
                        6,
                        colourInput("t3_outline_colour", "Outline Colour", value = "lightgray"),
                        prettyCheckbox("t3_outline", "Plot Outline", value = TRUE),
                        selectInput("t3_dl_type", "Download Format",
                            choices = c(".png", ".pdf", ".svg")
                        )
                    ),
                    column(
                        6,
                        selectInput("t3_palette", "Palette",
                            choices = c("viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo")
                        ),
                        numericInput("t3_opacity", "Opacity",
                            value = 0.51, min = 0.01, max = 1, step = 0.1
                        ),
                        prettyCheckbox("t3_reverse", "Reverse Palette", value = FALSE)
                    )
                )
            ),
            mainPanel(
                width = 10,
                fluidRow(
                    splitLayout(
                        tagList(
                            h3("European"),
                            plotOutput("t3_euro_anatogram", width = "210px", height = "300px"),
                            splitLayout(downloadButton("t3_euro_dl"), actionButton("t3_euro_show", "Show Data"))
                        ),
                        tagList(
                            h3("African"),
                            plotOutput("t3_afri_anatogram", width = "210px", height = "300px"),
                            splitLayout(downloadButton("t3_afri_dl"), actionButton("t3_afri_show", "Show Data"))
                        ),
                        tagList(
                            h3("East Asian"),
                            plotOutput("t3_easian_anatogram", width = "210px", height = "300px"),
                            splitLayout(downloadButton("t3_easian_dl"), actionButton("t3_easian_show", "Show Data"))
                        ),
                        tagList(
                            h3("Hisp./Latin Amer."),
                            plotOutput("t3_hisp_anatogram", width = "210px", height = "300px"),
                            splitLayout(downloadButton("t3_hisp_dl"), actionButton("t3_hisp_show", "Show Data"))
                        ),
                        tagList(
                            h3("Nat. Amer./Pac. Isl."),
                            plotOutput("t3_napi_anatogram", width = "210px", height = "300px"),
                            splitLayout(downloadButton("t3_napi_dl"), actionButton("t3_napi_show", "Show Data"))
                        ),
                        tagList(
                            h3("South Asian"),
                            plotOutput("t3_sasian_anatogram", width = "210px", height = "300px"),
                            splitLayout(downloadButton("t3_sasian_dl"), actionButton("t3_sasian_show", "Show Data"))
                        ),
                        tagList(
                            h3("Other/Unknown"),
                            plotOutput("t3_other_anatogram", width = "210px", height = "300px"),
                            splitLayout(downloadButton("t3_other_dl"), actionButton("t3_other_show", "Show Data"))
                        )
                    )
                ),
                fluidRow(
                    hr(),
                    div(DTOutput("t3_props"), style = "font-size:70%;")
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
            formatRound(c("proportion_mean", "proportion_lower", "proportion_upper"), 4)
    })

    ### Comparisons by Sex Tab

    t1_male_data <- reactive({
        make_plot_data(prop_data,
            celltype = input$t1_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t1_tissue_groups,
            age = input$t1_age,
            ethnicity = input$t1_ethnicity,
            sex = "male"
        )
    })

    t1_female_data <- reactive({
        make_plot_data(prop_data,
            celltype = input$t1_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t1_tissue_groups,
            age = input$t1_age,
            ethnicity = input$t1_ethnicity,
            sex = "female"
        )
    })

    t1_max_val <- reactive({
        max(t1_female_data()$value, t1_male_data()$value, na.rm = TRUE)
    })

    t1_male_plot <- reactive({
        direc <- ifelse(input$t1_reverse, -1, 1)

        p <- gganatogram(
            data = t1_male_data(), sex = "male", fill = "value",
            organism = "human", outline = input$t1_outline,
            fillOutline = input$t1_outline_colour,
        ) + theme_void()

        p + scale_fill_viridis(
            option = input$t1_palette,
            alpha = input$t1_opacity,
            direction = direc,
            limits = c(0, t1_max_val())
        )
    })

    output$t1_male_anatogram <- renderPlot({
        t1_male_plot()
    })

    output$t1_male_dl <- downloadHandler(
        filename = function() {
            paste0("male_anatogram", input$t1_dl_type)
        },
        content = function(file) {
            ggsave(file, t1_male_plot(), width = 6)
        }
    )

    output$t1_male_props <- renderDT({
        datatable(t1_male_data(),
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
    })

    t1_female_plot <- reactive({
        direc <- ifelse(input$t1_reverse, -1, 1)

        p <- gganatogram(
            data = t1_female_data(), sex = "female", fill = "value",
            organism = "human", outline = input$t1_outline,
            fillOutline = input$t1_outline_colour,
        ) + theme_void()

        p + scale_fill_viridis(
            option = input$t1_palette,
            alpha = input$t1_opacity,
            direction = direc,
            limits = c(0, t1_max_val())
        )
    })

    output$t1_female_anatogram <- renderPlot({
        t1_female_plot()
    })

    output$t1_female_dl <- downloadHandler(
        filename = function() {
            paste0("female_anatogram", input$t1_dl_type)
        },
        content = function(file) {
            ggsave(file, t1_female_plot(), width = 6)
        }
    )

    output$t1_female_props <- renderDT({
        datatable(t1_female_data(),
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
    })

    ### Comparisons by Age Tab

    t2_data <- reactive({
        infancy <- make_plot_data(prop_data,
            celltype = input$t2_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t2_tissue_groups,
            age = "Infancy",
            ethnicity = input$t2_ethnicity,
            sex = input$t2_sex
        )

        childhood <- make_plot_data(prop_data,
            celltype = input$t2_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t2_tissue_groups,
            age = "Childhood",
            ethnicity = input$t2_ethnicity,
            sex = input$t2_sex
        )

        adolescence <- make_plot_data(prop_data,
            celltype = input$t2_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t2_tissue_groups,
            age = "Adolescence",
            ethnicity = input$t2_ethnicity,
            sex = input$t2_sex
        )

        ya <- make_plot_data(prop_data,
            celltype = input$t2_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t2_tissue_groups,
            age = "Young Adulthood",
            ethnicity = input$t2_ethnicity,
            sex = input$t2_sex
        )

        middleage <- make_plot_data(prop_data,
            celltype = input$t2_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t2_tissue_groups,
            age = "Middle Age",
            ethnicity = input$t2_ethnicity,
            sex = input$t2_sex
        )

        senior <- make_plot_data(prop_data,
            celltype = input$t2_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t2_tissue_groups,
            age = "Senior",
            ethnicity = input$t2_ethnicity,
            sex = input$t2_sex
        )

        list(
            infancy = infancy,
            childhood = childhood,
            adolescence = adolescence,
            ya = ya,
            middleage = middleage,
            senior = senior
        )
    })

    t2_max_val <- reactive({
        max(t2_data()$infancy$value,
            t2_data()$childhood$value,
            t2_data()$adolescence$value,
            t2_data()$ya$value,
            t2_data()$middleage$value,
            t2_data()$senior$value,
            na.rm = TRUE
        )
    })

    t2_prop_data <- reactiveVal()

    output$t2_props <- renderDT({
        req(t2_prop_data)

        t2_prop_data()
    })

    # Observers to update the data table when the "Show Data" button is clicked
    observeEvent(input$t2_infancy_show, {
        t2_prop_data(t2_infancy_props())
    })

    observeEvent(input$t2_childhood_show, {
        t2_prop_data(t2_childhood_props())
    })

    observeEvent(input$t2_adolescence_show, {
        t2_prop_data(t2_adolescence_props())
    })

    observeEvent(input$t2_ya_show, {
        t2_prop_data(t2_ya_props())
    })

    observeEvent(input$t2_middleage_show, {
        t2_prop_data(t2_middleage_props())
    })

    observeEvent(input$t2_senior_show, {
        t2_prop_data(t2_senior_props())
    })

    observeEvent(t2_data(), {
        t2_prop_data(NULL)
    })

    t2_infancy_plot <- reactive({
        direc <- ifelse(input$t2_reverse, -1, 1)

        pdata <- t2_data()$infancy

        if (is.null(pdata) || nrow(pdata) == 0) {
            p <- ggplot() +
                theme_void() +
                theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
                geom_text(aes(x = 0.5, y = 0.5, label = "No results."),
                    inherit.aes = FALSE, check_overlap = TRUE
                )
        } else {
            p <- gganatogram(
                data = pdata, sex = input$t2_sex, fill = "value",
                organism = "human", outline = input$t2_outline,
                fillOutline = input$t2_outline_colour,
            ) + theme_void()

            p <- p + scale_fill_viridis(
                option = input$t2_palette,
                alpha = input$t2_opacity,
                direction = direc,
                limits = c(0, t2_max_val())
            )
        }

        p
    })

    output$t2_infancy_anatogram <- renderPlot({
        t2_infancy_plot()
    })

    output$t2_infancy_dl <- downloadHandler(
        filename = function() {
            paste0("infancy_anatogram", input$t2_dl_type)
        },
        content = function(file) {
            ggsave(file, t2_infancy_plot(), width = 6)
        }
    )

    t2_infancy_props <- reactive({
        dat <- t2_data()$infancy
        dat <- dat[!is.na(dat$value), ]

        datatable(dat,
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
    })

    t2_childhood_plot <- reactive({
        direc <- ifelse(input$t2_reverse, -1, 1)

        pdata <- t2_data()$childhood

        if (is.null(pdata) || nrow(pdata) == 0) {
            p <- ggplot() +
                theme_void() +
                theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
                geom_text(aes(x = 0.5, y = 0.5, label = "No results."),
                    inherit.aes = FALSE, check_overlap = TRUE
                )
        } else {
            p <- gganatogram(
                data = pdata, sex = input$t2_sex, fill = "value",
                organism = "human", outline = input$t2_outline,
                fillOutline = input$t2_outline_colour,
            ) + theme_void()

            p <- p + scale_fill_viridis(
                option = input$t2_palette,
                alpha = input$t2_opacity,
                direction = direc,
                limits = c(0, t2_max_val())
            )
        }

        p
    })

    output$t2_childhood_anatogram <- renderPlot({
        t2_childhood_plot()
    })

    output$t2_childhood_dl <- downloadHandler(
        filename = function() {
            paste0("childhood_anatogram", input$t2_dl_type)
        },
        content = function(file) {
            ggsave(file, t2_childhood_plot(), width = 6)
        }
    )

    t2_childhood_props <- reactive({
        dat <- t2_data()$childhood
        dat <- dat[!is.na(dat$value), ]

        datatable(dat,
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
    })

    t2_adolescence_plot <- reactive({
        direc <- ifelse(input$t2_reverse, -1, 1)

        pdata <- t2_data()$adolescence

        if (is.null(pdata) || nrow(pdata) == 0) {
            p <- ggplot() +
                theme_void() +
                theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
                geom_text(aes(x = 0.5, y = 0.5, label = "No results."),
                    inherit.aes = FALSE, check_overlap = TRUE
                )
        } else {
            p <- gganatogram(
                data = pdata, sex = input$t2_sex, fill = "value",
                organism = "human", outline = input$t2_outline,
                fillOutline = input$t2_outline_colour,
            ) + theme_void()

            p <- p + scale_fill_viridis(
                option = input$t2_palette,
                alpha = input$t2_opacity,
                direction = direc,
                limits = c(0, t2_max_val())
            )
        }

        p
    })

    output$t2_adolescence_anatogram <- renderPlot({
        t2_adolescence_plot()
    })

    output$t2_adolescence_dl <- downloadHandler(
        filename = function() {
            paste0("adolescence_anatogram", input$t2_dl_type)
        },
        content = function(file) {
            ggsave(file, t2_adolescence_plot(), width = 6)
        }
    )

    t2_adolescence_props <- reactive({
        dat <- t2_data()$adolescence
        dat <- dat[!is.na(dat$value), ]

        datatable(dat,
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
    })

    t2_ya_plot <- reactive({
        direc <- ifelse(input$t2_reverse, -1, 1)

        pdata <- t2_data()$ya

        if (is.null(pdata) || nrow(pdata) == 0) {
            p <- ggplot() +
                theme_void() +
                theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
                geom_text(aes(x = 0.5, y = 0.5, label = "No results."),
                    inherit.aes = FALSE, check_overlap = TRUE
                )
        } else {
            p <- gganatogram(
                data = pdata, sex = input$t2_sex, fill = "value",
                organism = "human", outline = input$t2_outline,
                fillOutline = input$t2_outline_colour,
            ) + theme_void()

            p <- p + scale_fill_viridis(
                option = input$t2_palette,
                alpha = input$t2_opacity,
                direction = direc,
                limits = c(0, t2_max_val())
            )
        }

        p
    })

    output$t2_ya_anatogram <- renderPlot({
        t2_ya_plot()
    })

    output$t2_ya_dl <- downloadHandler(
        filename = function() {
            paste0("ya_anatogram", input$t2_dl_type)
        },
        content = function(file) {
            ggsave(file, t2_ya_plot(), width = 6)
        }
    )

    t2_ya_props <- reactive({
        dat <- t2_data()$ya
        dat <- dat[!is.na(dat$value), ]

        datatable(dat,
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
    })

    t2_middleage_plot <- reactive({
        direc <- ifelse(input$t2_reverse, -1, 1)

        pdata <- t2_data()$middleage

        if (is.null(pdata) || nrow(pdata) == 0) {
            p <- ggplot() +
                theme_void() +
                theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
                geom_text(aes(x = 0.5, y = 0.5, label = "No results."),
                    inherit.aes = FALSE, check_overlap = TRUE
                )
        } else {
            p <- gganatogram(
                data = pdata, sex = input$t2_sex, fill = "value",
                organism = "human", outline = input$t2_outline,
                fillOutline = input$t2_outline_colour,
            ) + theme_void()

            p <- p + scale_fill_viridis(
                option = input$t2_palette,
                alpha = input$t2_opacity,
                direction = direc,
                limits = c(0, t2_max_val())
            )
        }

        p
    })

    output$t2_middleage_anatogram <- renderPlot({
        t2_middleage_plot()
    })

    output$t2_middleage_dl <- downloadHandler(
        filename = function() {
            paste0("middleage_anatogram", input$t2_dl_type)
        },
        content = function(file) {
            ggsave(file, t2_middleage_plot(), width = 6)
        }
    )

    t2_middleage_props <- reactive({
        dat <- t2_data()$middleage
        dat <- dat[!is.na(dat$value), ]

        datatable(dat,
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
    })

    t2_senior_plot <- reactive({
        direc <- ifelse(input$t2_reverse, -1, 1)

        pdata <- t2_data()$senior

        if (is.null(pdata) || nrow(pdata) == 0) {
            p <- ggplot() +
                theme_void() +
                theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
                geom_text(aes(x = 0.5, y = 0.5, label = "No results."),
                    inherit.aes = FALSE, check_overlap = TRUE
                )
        } else {
            p <- gganatogram(
                data = pdata, sex = input$t2_sex, fill = "value",
                organism = "human", outline = input$t2_outline,
                fillOutline = input$t2_outline_colour,
            ) + theme_void()

            p <- p + scale_fill_viridis(
                option = input$t2_palette,
                alpha = input$t2_opacity,
                direction = direc,
                limits = c(0, t2_max_val())
            )
        }

        p
    })

    output$t2_senior_anatogram <- renderPlot({
        t2_senior_plot()
    })

    output$t2_senior_dl <- downloadHandler(
        filename = function() {
            paste0("senior_anatogram", input$t2_dl_type)
        },
        content = function(file) {
            ggsave(file, t2_senior_plot(), width = 6)
        }
    )

    t2_senior_props <- reactive({
        dat <- t2_data()$senior
        dat <- dat[!is.na(dat$value), ]

        datatable(dat,
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
    })

    ### Comparisons by Ethnicities Tab
    t3_data <- reactive({
        euro <- make_plot_data(prop_data,
            celltype = input$t3_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t3_tissue_groups,
            age = input$t3_age,
            ethnicity = "European",
            sex = input$t3_sex
        )

        afri <- make_plot_data(prop_data,
            celltype = input$t3_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t3_tissue_groups,
            age = input$t3_age,
            ethnicity = "African",
            sex = input$t3_sex
        )

        easian <- make_plot_data(prop_data,
            celltype = input$t3_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t3_tissue_groups,
            age = input$t3_age,
            ethnicity = "East Asian",
            sex = input$t3_sex
        )

        hisp <- make_plot_data(prop_data,
            celltype = input$t3_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t3_tissue_groups,
            age = input$t3_age,
            ethnicity = "Hispanic/Latin American",
            sex = input$t3_sex
        )

        napi <- make_plot_data(prop_data,
            celltype = input$t3_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t3_tissue_groups,
            age = input$t3_age,
            ethnicity = "Native American & Pacific Islander",
            sex = input$t3_sex
        )

        other <- make_plot_data(prop_data,
            celltype = input$t3_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t3_tissue_groups,
            age = input$t3_age,
            ethnicity = "Other/Unknown",
            sex = input$t3_sex
        )

        sasian <- make_plot_data(prop_data,
            celltype = input$t3_cell_type,
            male_organ_map = male_organ_list,
            female_organ_map = female_organ_list,
            tissue_groups = input$t3_tissue_groups,
            age = input$t3_age,
            ethnicity = "South Asian",
            sex = input$t3_sex
        )

        list(
            euro = euro,
            afri = afri,
            easian = easian,
            hisp = hisp,
            napi = napi,
            other = other,
            sasian = sasian
        )
    })

    t3_max_val <- reactive({
        max(t3_data()$euro$value,
            t3_data()$afri$value,
            t3_data()$easian$value,
            t3_data()$hisp$value,
            t3_data()$napi$value,
            t3_data()$other$value,
            t3_data()$sasian$value,
            na.rm = TRUE
        )
    })

    t3_prop_data <- reactiveVal()

    output$t3_props <- renderDT({
        req(t3_prop_data)

        t3_prop_data()
    })

    t3_euro_plot <- reactive({
        direc <- ifelse(input$t3_reverse, -1, 1)

        pdata <- t3_data()$euro

        if (is.null(pdata) || nrow(pdata) == 0) {
            p <- ggplot() +
                theme_void() +
                theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
                geom_text(aes(x = 0.5, y = 0.5, label = "No results."),
                    inherit.aes = FALSE, check_overlap = TRUE
                )
        } else {
            p <- gganatogram(
                data = pdata, sex = input$t3_sex, fill = "value",
                organism = "human", outline = input$t3_outline,
                fillOutline = input$t3_outline_colour,
            ) + theme_void()

            p <- p + scale_fill_viridis(
                option = input$t3_palette,
                alpha = input$t3_opacity,
                direction = direc,
                limits = c(0, t3_max_val())
            )
        }

        p
    })

    # Observers to update the data table when the "Show Data" button is clicked
    observeEvent(input$t3_euro_show, {
        t3_prop_data(t3_euro_props())
    })

    observeEvent(input$t3_afri_show, {
        t3_prop_data(t3_afri_props())
    })

    observeEvent(input$t3_easian_show, {
        t3_prop_data(t3_easian_props())
    })

    observeEvent(input$t3_hisp_show, {
        t3_prop_data(t3_hisp_props())
    })

    observeEvent(input$t3_napi_show, {
        t3_prop_data(t3_napi_props())
    })

    observeEvent(input$t3_other_show, {
        t3_prop_data(t3_other_props())
    })

    observeEvent(input$t3_sasian_show, {
        t3_prop_data(t3_sasian_props())
    })

    observeEvent(t3_data(), {
        t3_prop_data(NULL)
    })

    output$t3_euro_anatogram <- renderPlot({
        t3_euro_plot()
    })

    output$t3_euro_dl <- downloadHandler(
        filename = function() {
            paste0("euro_anatogram", input$t3_dl_type)
        },
        content = function(file) {
            ggsave(file, t3_euro_plot(), width = 6)
        }
    )

    t3_euro_props <- reactive({
        dat <- t3_data()$euro
        dat <- dat[!is.na(dat$value), ]

        datatable(dat,
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
            formatStyle(0, target = "row", lineHeight = "70%") %>%
            formatRound(c("value", "CI_lower", "CI_upper"), 4)
    })

    t3_afri_plot <- reactive({
        direc <- ifelse(input$t3_reverse, -1, 1)

        pdata <- t3_data()$afri

        if (is.null(pdata) || nrow(pdata) == 0) {
            p <- ggplot() +
                theme_void() +
                theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
                geom_text(aes(x = 0.5, y = 0.5, label = "No results."),
                    inherit.aes = FALSE, check_overlap = TRUE
                )
        } else {
            p <- gganatogram(
                data = pdata, sex = input$t3_sex, fill = "value",
                organism = "human", outline = input$t3_outline,
                fillOutline = input$t3_outline_colour,
            ) + theme_void()

            p <- p + scale_fill_viridis(
                option = input$t3_palette,
                alpha = input$t3_opacity,
                direction = direc,
                limits = c(0, t3_max_val())
            )
        }

        p
    })

    output$t3_afri_anatogram <- renderPlot({
        t3_afri_plot()
    })

    output$t3_afri_dl <- downloadHandler(
        filename = function() {
            paste0("afri_anatogram", input$t3_dl_type)
        },
        content = function(file) {
            ggsave(file, t3_afri_plot(), width = 6)
        }
    )

    t3_afri_props <- reactive({
        dat <- t3_data()$afri
        dat <- dat[!is.na(dat$value), ]

        datatable(dat,
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
            formatStyle(0, target = "row", lineHeight = "70%") %>%
            formatRound(c("value", "CI_lower", "CI_upper"), 4)
    })

    t3_easian_plot <- reactive({
        direc <- ifelse(input$t3_reverse, -1, 1)

        pdata <- t3_data()$easian

        if (is.null(pdata) || nrow(pdata) == 0) {
            p <- ggplot() +
                theme_void() +
                theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
                geom_text(aes(x = 0.5, y = 0.5, label = "No results."),
                    inherit.aes = FALSE, check_overlap = TRUE
                )
        } else {
            p <- gganatogram(
                data = pdata, sex = input$t3_sex, fill = "value",
                organism = "human", outline = input$t3_outline,
                fillOutline = input$t3_outline_colour,
            ) + theme_void()

            p <- p + scale_fill_viridis(
                option = input$t3_palette,
                alpha = input$t3_opacity,
                direction = direc,
                limits = c(0, t3_max_val())
            )
        }

        p
    })

    output$t3_easian_anatogram <- renderPlot({
        t3_easian_plot()
    })

    output$t3_easian_dl <- downloadHandler(
        filename = function() {
            paste0("east_asian_anatogram", input$t3_dl_type)
        },
        content = function(file) {
            ggsave(file, t3_easian_plot(), width = 6)
        }
    )

    t3_easian_props <- reactive({
        dat <- t3_data()$easian
        dat <- dat[!is.na(dat$value), ]

        datatable(dat,
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
            formatStyle(0, target = "row", lineHeight = "70%") %>%
            formatRound(c("value", "CI_lower", "CI_upper"), 4)
    })

    t3_hisp_plot <- reactive({
        direc <- ifelse(input$t3_reverse, -1, 1)

        pdata <- t3_data()$hisp

        if (is.null(pdata) || nrow(pdata) == 0) {
            p <- ggplot() +
                theme_void() +
                theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
                geom_text(aes(x = 0.5, y = 0.5, label = "No results."),
                    inherit.aes = FALSE, check_overlap = TRUE
                )
        } else {
            p <- gganatogram(
                data = pdata, sex = input$t3_sex, fill = "value",
                organism = "human", outline = input$t3_outline,
                fillOutline = input$t3_outline_colour,
            ) + theme_void()

            p <- p + scale_fill_viridis(
                option = input$t3_palette,
                alpha = input$t3_opacity,
                direction = direc,
                limits = c(0, t3_max_val())
            )
        }

        p
    })

    output$t3_hisp_anatogram <- renderPlot({
        t3_hisp_plot()
    })

    output$t3_hisp_dl <- downloadHandler(
        filename = function() {
            paste0("hispanic_la_anatogram", input$t3_dl_type)
        },
        content = function(file) {
            ggsave(file, t3_hisp_plot(), width = 6)
        }
    )

    t3_hisp_props <- reactive({
        dat <- t3_data()$hisp
        dat <- dat[!is.na(dat$value), ]

        datatable(dat,
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
            formatStyle(0, target = "row", lineHeight = "70%") %>%
            formatRound(c("value", "CI_lower", "CI_upper"), 4)
    })

    t3_napi_plot <- reactive({
        direc <- ifelse(input$t3_reverse, -1, 1)

        pdata <- t3_data()$napi

        if (is.null(pdata) || nrow(pdata) == 0) {
            p <- ggplot() +
                theme_void() +
                theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
                geom_text(aes(x = 0.5, y = 0.5, label = "No results."),
                    inherit.aes = FALSE, check_overlap = TRUE
                )
        } else {
            p <- gganatogram(
                data = pdata, sex = input$t3_sex, fill = "value",
                organism = "human", outline = input$t3_outline,
                fillOutline = input$t3_outline_colour,
            ) + theme_void()

            p <- p + scale_fill_viridis(
                option = input$t3_palette,
                alpha = input$t3_opacity,
                direction = direc,
                limits = c(0, t3_max_val())
            )
        }

        p
    })

    output$t3_napi_anatogram <- renderPlot({
        t3_napi_plot()
    })

    output$t3_napi_dl <- downloadHandler(
        filename = function() {
            paste0("native_american_pacific_islander_anatogram", input$t3_dl_type)
        },
        content = function(file) {
            ggsave(file, t3_napi_plot(), width = 6)
        }
    )

    t3_napi_props <- reactive({
        dat <- t3_data()$napi
        dat <- dat[!is.na(dat$value), ]

        datatable(dat,
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
            formatStyle(0, target = "row", lineHeight = "70%") %>%
            formatRound(c("value", "CI_lower", "CI_upper"), 4)
    })

    t3_other_plot <- reactive({
        direc <- ifelse(input$t3_reverse, -1, 1)

        pdata <- t3_data()$other

        if (is.null(pdata) || nrow(pdata) == 0) {
            p <- ggplot() +
                theme_void() +
                theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
                geom_text(aes(x = 0.5, y = 0.5, label = "No results."),
                    inherit.aes = FALSE, check_overlap = TRUE
                )
        } else {
            p <- gganatogram(
                data = pdata, sex = input$t3_sex, fill = "value",
                organism = "human", outline = input$t3_outline,
                fillOutline = input$t3_outline_colour,
            ) + theme_void()

            p <- p + scale_fill_viridis(
                option = input$t3_palette,
                alpha = input$t3_opacity,
                direction = direc,
                limits = c(0, t3_max_val())
            )
        }

        p
    })

    output$t3_other_anatogram <- renderPlot({
        t3_other_plot()
    })

    output$t3_other_dl <- downloadHandler(
        filename = function() {
            paste0("other_unknown_anatogram", input$t3_dl_type)
        },
        content = function(file) {
            ggsave(file, t3_other_plot(), width = 6)
        }
    )

    t3_other_props <- reactive({
        dat <- t3_data()$other
        dat <- dat[!is.na(dat$value), ]

        datatable(dat,
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
            formatStyle(0, target = "row", lineHeight = "70%") %>%
            formatRound(c("value", "CI_lower", "CI_upper"), 4)
    })

    t3_sasian_plot <- reactive({
        direc <- ifelse(input$t3_reverse, -1, 1)

        pdata <- t3_data()$sasian

        if (is.null(pdata) || nrow(pdata) == 0) {
            p <- ggplot() +
                theme_void() +
                theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
                geom_text(aes(x = 0.5, y = 0.5, label = "No results."),
                    inherit.aes = FALSE, check_overlap = TRUE
                )
        } else {
            p <- gganatogram(
                data = pdata, sex = input$t3_sex, fill = "value",
                organism = "human", outline = input$t3_outline,
                fillOutline = input$t3_outline_colour,
            ) + theme_void()

            p <- p + scale_fill_viridis(
                option = input$t3_palette,
                alpha = input$t3_opacity,
                direction = direc,
                limits = c(0, t3_max_val())
            )
        }

        p
    })

    output$t3_sasian_anatogram <- renderPlot({
        t3_sasian_plot()
    })

    output$t3_sasian_dl <- downloadHandler(
        filename = function() {
            paste0("south_asian_anatogram", input$t3_dl_type)
        },
        content = function(file) {
            ggsave(file, t3_sasian_plot(), width = 6)
        }
    )

    t3_sasian_props <- reactive({
        dat <- t3_data()$sasian
        dat <- dat[!is.na(dat$value), ]

        datatable(dat,
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
            formatStyle(0, target = "row", lineHeight = "70%") %>%
            formatRound(c("value", "CI_lower", "CI_upper"), 4)
    })
}

shinyApp(ui = ui, server = server)
