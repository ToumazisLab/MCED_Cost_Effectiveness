# ui.R
source("Global.R")



cancer_types_overall <- c(
  "anus",
  "bladder",
  "breast",
  "cervix",
  "colon and rectum",
  "esophagus",
  "head and neck",
  "kidney and renal pelvis",
  "liver and intrahepatic bile duct",
  "lung and bronchus",
  "lymphoma",
  "ovarian",
  "pancreas",
  "prostate",
  "stomach",
  "uterus"
)

display_names <- str_to_title(cancer_types_overall)
# Fix small words (so “Colon and Rectum” instead of “Colon And Rectum”):
display_names <- str_replace_all(display_names,
                                 c("\\bAnd\\b" = "and", "\\bOf\\b" = "of", "\\bIn\\b" = "in")
)

# 3) Build a named vector: names = “Colon and Rectum” → value = “colon and rectum”
cancer_choices <- setNames(cancer_types_overall, display_names)

ss_choices <- c(
  "Slow"        = "1",
  "Medium"      = "3",
  "Medium fast" = "5",
  "Fast"        = "7"
)

survival_labels <- c(
  floor_max   = "Lead-time corrected",
  add_shift   = "Full stage-shift benefit",
  max_new_old = "Crossover-guarded",
  naive       = "No correction"
)

#### Sidebar ####
sidebar <- dashboardSidebar(
  # Always show the "landing" page:
  # Then, add a *second* sidebarMenu that only appears if consent is checked:
  conditionalPanel(
    condition = "output.userConsented",
    sidebarMenu(
      id = "main_menu",
      menuItem("Customize Run", tabName = "model_params", icon = icon("dashboard")),
      menuItem("Results", tabName = "results", icon = icon("th"))
    )
  )
)

#### Body ####
body <- dashboardBody(
  tags$head(
    # Existing tooltips/popovers
    tags$script(HTML("
    $(function () {
      $('[data-toggle=\"tooltip\"]').tooltip();
      $('[data-toggle=\"popover\"]').popover({
        trigger: 'focus',
        placement: 'right',
        html: true
      });
    });
  ")),
    
    # New handler for tab activation
    tags$script(HTML("
    Shiny.addCustomMessageHandler('activateTab', function(message) {
      var el = document.getElementById(message.id);
      if (el) el.click();
    });
  ")),
    tags$style(HTML("
  .box { margin-bottom: 10px !important; }
  .slider-input { margin-bottom: 10px; }
")),
    tags$style(HTML("
    
    .centered-block {
      display: flex;
      justify-content: center;
      align-items: center;
      gap: 40px; /* Optional: increases spacing */
    }
    
    .arrow-block {
      margin-top: 10px;
      margin-bottom: 10px;
    }
    
    .diagnostic-grid {
      flex-wrap: nowrap; /* Prevent line wrap */
    }

    .diagnostic-group {
      flex: 1;
    }
    
    .takehome-box-ce {
    border: 2px solid #4CAF50;
    border-radius: 8px;
    padding: 15px;
    background-color: #f6fff6;
    justify-content: center;
    color: #2e7d32;
    font-size: 20px;
    margin-bottom: 10px;
  }
  
  .takehome-box-notce {
    border: 2px solid #f44336;
    border-radius: 8px;
    padding: 15px;
    background-color: #fff6f6;
    justify-content: center;
    color: #b71c1c;
    font-size: 20px;
    margin-bottom: 10px;
  }
  
")),
    tags$style(HTML("
    #run_tabbox_wrapper .nav-tabs {
      position: sticky;
      top: 0;
      z-index: 999;
      background-color: white;
      border-bottom: 1px solid #ccc;
      padding-top: 8px;
    }
  ")),
    tags$style(HTML("
.badge { display:inline-block; padding:2px 8px; border-radius:10px; font-size:12px; margin-left:6px; }
.badge-ok { background:#e8f5e9; color:#2e7d32; }
.badge-default { background:#e3f2fd; color:#0d47a1; }
.badge-missing { background:#fff3e0; color:#e65100; }
.badge-invalid { background:#ffebee; color:#b71c1c; }
.helper { color:#666; font-size:12px; margin-top:6px; }
")),
    tags$style(HTML("
  .run-panel {
    background: #f9fafb;
    border: 1px solid #e5e7eb;
    border-radius: 8px;
    padding: 12px 12px 8px 12px;
  }
  /* Optional: make right panel stick as you scroll inside long box */
  @media (min-width: 992px) {
    .run-panel { position: sticky; top: 12px; }
  }
  /* Tidy up the checklist spacing */
  .helper { color:#6b7280; font-size: 90%; margin-top: 4px; }
  .badge { 
    display:inline-block; padding:2px 8px; border-radius:10px; font-size: 12px; 
    margin-left:6px; vertical-align: middle;
  }
  .badge-ok { background:#e6f4ea; color:#1e7e34; }
  .badge-default { background:#eef2ff; color:#3730a3; }
  .badge-missing { background:#fff1f2; color:#b91c1c; }
  .badge-invalid { background:#fef3c7; color:#92400e; }
"))
  ),
  tabItems(
    #### Model Parameters Tab ####
    tabItem(
      tabName = "model_params",
      
      ##### Row 1: Cancer selection ----
      fluidRow(
        box(
          title = tags$span(style = "font-size:20px;", strong("1. Select cancer types")),
          width = 12,
          div(style = "column-count:4; column-gap:20px; margin-top:8px;",
              awesomeCheckboxGroup(
                inputId = "selected_cancers",
                label   = NULL,
                choices = cancer_choices
              )
          ),
          div(style = "margin-top:10px;",
              actionButton("reset_all", "Reset All",  icon = icon("redo")),
              actionButton("select_all","Select All", icon = icon("check-square"))
          )
        )
      ),
      
      ##### Row 2: Two columns ----
      fluidRow(
        height = "500px",
        ## Left Column
        column(
          width = 4,
          
          # 2a. Cost Box
          box(
            width = NULL,
            height = "200px",
            title = strong("2a. Cost"),
            solidHeader = TRUE,
            numericInput("mced_test_cost", "MCD Test Cost ($)", value = 500, min = 0, max = 1e5, step = 50)
          ),
          
          # 2b. More Model Parameters Box
          box(
            width = NULL,
            height = "300px",
            title = strong("2b. More model parameters"),
            collapsible = TRUE,
            collapsed = FALSE,
            status = "primary",
            solidHeader = TRUE,
            p("These parameters let you fine-tune compliance, follow-up, and invasive rates."),
              div(
              style = "max-height: 300px; overflow-y: auto; padding-right: 10px;",
              sliderInput("mced_age_range", "MCD Age Range (yrs)", 40, 100, value = c(45, 80), step = 1),
              uiOutput("fp_track_slider_ui"),
              sliderInput("mced_compliance_percent", "MCD Compliance (%)", 0, 100, 75, step = 5),
              sliderInput("u_mced_fp_compliance", "Follow-up after positive MCD (%)", 0, 100, 100, step = 10),
              sliderInput("u_invasive_multi", "Increased invasive rate after dual FP (%)", 0, 100, 10, step = 10),

              radioButtons(
                inputId = "survival_method",
                label   = "Survival handling for stage-shifted cases",
                choices = c(
                  "Lead-time corrected (recommended)" = "floor_max",
                  "Full stage-shift benefit"          = "add_shift",
                  "Crossover-guarded only"            = "max_new_old"
                ),
                selected = "floor_max"
              ),

              radioButtons(
                inputId = "ss_scenario",
                label   = "Dwell-time scenario",
                choices = ss_choices,
                selected = "3"   # set to your current input_params$ss_scenario default
              )
            )
          )
        ),
        
        ## Right Column
        column(
          width = 8,
          box(
            width = NULL,
            height = "500px",
            title = strong("2c. Accuracy assumptions"),
            solidHeader = TRUE,
            tabsetPanel(
              id = "acc_tabs", type = "pills",
              
              tabPanel("Sensitivity",
                       div(style = "height: 10px;"),  # adds vertical space
                       p(strong("Adjust Sensitivity by Stage")),
                       p("Each row represents a selected cancer type. Enter the test sensitivity for each stage."),
                       div(
                         style = "max-height: 420px; overflow-y: auto; padding-right: 10px;",
                         uiOutput("sensitivity_or_placeholder")
                       )
              ),
              
              tabPanel("Specificity",
                       div(style = "height: 10px;"),
                       selectInput("spec_mode", "How do you want to set specificity?",
                                   choices = c("- Select One -" = "", "Set one value for all cancers" = "all", "Set separately by cancer" = "by_cancer"),
                                   selected = "all"),
                       div(style = "max-height: 420px; overflow-y: auto; padding-right: 10px;",
                           uiOutput("specificity_or_placeholder")
                       )
              )
            )
          )
        )
      ),
      
      ##### Row 3: Run model ----
      ##### Row 3: Review & Run ----
      fluidRow(
        column(
          width = 8, offset = 2,   # 2 + 8 + 2 = 12 → centered
          box(
            title = tags$span(style = "font-size:20px;", strong("3. Review & Run")),
            width = 12,            # fill the column
            status = "primary",
            solidHeader = TRUE,
            fluidRow(
              column(
                width = 7,  # use full 12 inside the box (e.g., 7/5)
                div(class = "run-panel", uiOutput("checklist_ui"))
              ),
              column(
                width = 5,
                div(class = "run-panel", uiOutput("run_model_ui"))
              )
            )
          )
        )
      )
      
    ),
    #### Second tab content ####
    tabItem(tabName = "results",
            h2("Results"),
            box(
              title = "Key Terms",
              width = 12,
              collapsible = TRUE,
              collapsed = TRUE,
              
              p(strong("Incremental Cost-Effectiveness Ratio (ICER):"),
                "The ICER is calculated as the ratio between the change in cost divided by the change in QALYs between strategies. It representes the additional cost incurred per QALY gained from adding the MCD test to UC.
                  It’s compared to a willingness-to-pay threshold to determine if an intervention is cost-effective. A strategy is cost-effective is the ICER is lower than the willingness-to-pay threshold."),
              
              p(strong("Possible Screening Outcomes:"),br(),
                "• True Positive (TP): Cancer is present and the test detects cancer",br(),
                "• False Positive (FP): Cancer is NOT present and the test mistakenly detects cancer",br(),
                "• True Negative (TN): Cancer is NOT present and the test does not detect cancer",br(),
                "• False Negative (FN): Cancer is present and the test does not detect cancer"),
              
              p(strong("Using Plotly:"),
                "You can hover to see tooltips, zoom in, or click items in the legend to hide/show data series."),
              
              p(strong("How the Results Are Organized:"),
                "The Cost-Effectiveness Plane and TP–FP Bubble present summary results across all model runs.
     To the right, the 'Results' section may contain up to 10 individual tabs—one for each model run—with detailed, run-specific outputs.
     This allows you to compare and review different scenarios or parameter sets.")
            ),
            fluidRow(
              # -------- LEFT COLUMN: Two stacked boxes -------------
              column(
                width = 6,
                
                # -- Top box: Overall results --
                box(
                  title = "Overall Results",
                  status = "primary",
                  width = 12,
                  selectInput(
                    inputId = "wtp_threshold",
                    label   = "Willingness-to-Pay Threshold ($/QALY)",
                    choices = c(
                      "$50,000"  = 50000,
                      "$75,000"  = 75000,
                      "$100,000" = 100000,
                      "$150,000" = 150000,
                      "$200,000" = 200000,
                      "$300,000" = 300000
                    ),
                    selected = 100000
                  ),
                  tabBox(
                    id = "results_plots",
                    width = 12,
                    tabPanel(
                      title = "Cost‑Effectiveness",
                      p(strong("Cost‑Effectiveness (CE) Plane:"),
                        "Each point corresponds to a specific run of the model. The x-axis shows incremental QALYs (i.e., benefit) gained by supplementing Usual Care (UC) with the specificed MCD test (MCD+UC), and the y-axis shows incremental costs incurred relative to UC alone.
                        The slope between points is the ICER (Incremental Cost-Effectiveness Ratio). The origin corresponds to the Usual Care (UC). Users have the option to select their willingness-to-pay threshold from the dropdown menu above."),
                      plotlyOutput("ce_plotly", height = "400px")
                    ),
                    tabPanel(
                      title = "TP‑FP Bubble",
                      p(strong("TP–FP Bubble Plot:"),
                        "Displays tradeoffs between true positives and false positives across all model runs.
                        Each bubble corresponds to a specific run of the model and it size represents the total cost of the MCD+UC arm."),
                      plotlyOutput("bubble_plotly", height = "400px")
                    )
                  )
                ),
                
                # -- Bottom box: Additional content --
                uiOutput("ce_summary_box")
              ),
              
              # -------- RIGHT COLUMN: Dynamic Run Tabs ------------
              column(
                width = 6,
                uiOutput("dynamic_runs_tabbox_ui")
              )
            ),
            
            absolutePanel(
              id = "add_run_button_panel",
              bottom = 5, right = 30, width = "auto", draggable = FALSE,
              actionButton("add_new_run", "+ Add New Run", 
                           style = "color: white; background-color: #007bff; font-weight: bold; 
                        border-radius: 20px; padding: 10px 20px; font-size: 16px;")
            )
    )
  )
)
#### Dashboard Page ####
dashboardPage(
  dashboardHeader(title = "Multi-Cancer Detection"),
  sidebar,
  body
)