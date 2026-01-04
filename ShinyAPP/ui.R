# Load the required libraries
library("shiny")
library("shinythemes")
library("shinyWidgets")
library("shinymaterial")
library("markdown")
library("shinyFiles")
library("ff")
library("bit")


options(shiny.maxRequestSize = 1024^3)

# Define the UI
ui <- fluidPage(
  tags$head(
    tags$script(
      HTML("
        document.title = 'OpEnCAST';
      ")
    )
  ),
  theme = shinytheme("superhero"),
  tags$head(
    tags$style(HTML(
      "
      body {
        background:
    linear-gradient(
        to right,
        #1A0E06 0%,   /* dark brown */
        #4A2A18 50%,  /* warm brown center glow */
        #1A0E06 100%  /* dark brown */
    ),
    radial-gradient(circle at center, rgba(255, 255, 255, 0.08), transparent 70%);
background-blend-mode: overlay;

      
     .navbar {
    background-color: #001a12;  /* Very dark black-green */
    color: #fff;
    font-size: 19px;
}

.navbar-brand {
    color: #A8FFBF !important;  /* soft mint green accent */
    font-size: 36px;
    font-weight: bold;
}

.nav-link {
    color: #C7FFD8 !important;  /* pale green for readability */
    font-size: 24px;
}

.nav-link:hover {
    color: #7CFFAF !important;  /* brighter mint on hover */
}



      .jumbotron {
        background-color: #f8f9fa;
        padding: 40px;
        margin-bottom: 30px;
      }
      
      .jumbotron h2 {
        font-weight: bold;
        color: #333;
      }
      
      .jumbotron p {
        color: #555;
      }
      
      .panel-heading {
        background-color: #333;
        color: #fff;
        padding: 10px;
        border-radius: 4px;
      }
      
      .panel-body {
        padding: 20px;
      }
      
      .btn-primary {
        background-color: #007bff;
        border-color: #007bff;
      }
      
      .btn-primary:hover {
        background-color: #0069d9;
        border-color: #0062cc;
      }
      
      .custom-sidebar {
       background-color: rgba(0,0,0, 0.5);
        width: 400px;
      }
      /* Custom Yellow Button */
.btn-yellow {
  background-color: #f2c200 !important;   /* golden yellow */
  border-color: #b68410 !important;
  color:  #3b2200 !important;              /* dark green text */
  font-weight: 600;
  font-size: 16px;
}
.btn-yellow:hover {
  background-color: #ffd629 !important;
  border-color: #c89615 !important;
  color: #000 !important;
}

      
      .custom-align {
        display: flex;
        align-items: center;
        justify-content: center;
      }
      
      .predict-panel {
        background-color: rgba(0,0,0, 0.5);
        padding: 20px;
        height: 310px;
        border-radius: 10px;
        box-shadow: 0 0 15px rgba(255, 215, 0, 0.8);
        overflow-y: auto;
      }
      /* Center only the table output inside the prediction panel */
.predict-panel table {
  margin-left: auto !important;
  margin-right: auto !important;
  text-align: center;
}

/* Optional: center the table cells too */
.predict-panel td, 
.predict-panel th {
  text-align: center !important;
}

body { padding-bottom: 10px; }

      
      .predict-description {
        color: #fff;
        text-align: center;
        margin-bottom: 30px;
        font-size: 28px;
      }
      
      .predict-button {
        display: flex;
        justify-content: center;
        margin-top: 30px;
      }
      /* Increase font size of output tables */
  .predict-panel table {
  font-size: 14px !important;   /* Change size here: 20px, 24px, 28px */
  font-weight: 300;              /* Optional: make it semi-bold */
  text-align: center !important;
  }

      /* Align radio input and image nicely */
  .choice-with-img {
  display: inline-flex;
  align-items: center;
  gap: 8px;
  font-size: 18px;
  }

/* make the whole label clickable and a bit larger */
  .radio label {
  cursor: pointer;
  padding: 6px 8px;
  display: inline-flex;
  align-items: center;
}
  .italic-plant {
    font-style: italic;
}

      "
    ))
  ),
  
  navbarPage(
    title = NULL,
    collapsible = TRUE,
    theme = "bootstrap",
    tabPanel("Home", includeHTML("www/aboutPMeth.html")),
    
    
    tabPanel(
      "Predict",
      sidebarPanel(
        tags$label(h2("Input Details", align = "center")),
        #textInput("dir_name", "Enter Directory Name"),
        #actionButton("save_button", "Save Files"),
        fileInput("fasta_file", "Upload Fasta file",accept = c(".fasta", ".fa", ".FASTA")),
       
        radioButtons(
          inputId = "Reference_type",
          label = "Choose Reference Model:",
          choiceNames = list(
            # Monocot choice (image + caption)
            tagList(
              tags$div(class = "choice-with-img",
                       tags$img(src = "monocot.png", alt = "Monocot (Oryza)", style = "height:48px; margin-right:20px; vertical-align:middle;"),
                       HTML("<span>Monocot (<i>Oryza sp.</i>)</span>")
              )
            ),
            # Dicot choice
            tagList(
              tags$div(class = "choice-with-img",
                       tags$img(src = "dicot.png", alt = "Dicot (Arabidopsis)", style = "height:48px; margin-right:20px; vertical-align:middle;"),
                       HTML("<span>Dicot (<i>Arabidopsis sp.</i>)</span>")
              )
            ),
            # Dual choice (you can use a small combined icon or simple text)
            tagList(
              tags$div(class = "choice-with-img",
                       tags$img(src = "both.png", alt = "Both Models", style = "height:48px; margin-right:20px; vertical-align:middle;"),
                       HTML("<span>Both References (<i>Oryza sp.</i>+<i>Arabidopsis sp.</i>)</span>")
              )
            )
          ),
          choiceValues = c("Rice", "Arabidopsis", "both"),
          selected = "both"
        )
        ,
        
        actionButton("submitbutton", "Run Analysis", class = "btn btn-yellow"),
      
        class = "custom-sidebar"
      ),
      mainPanel(
        div(
          class = "predict-panel",
          style = "margin-bottom: 30px;",
          p(class = "predict-description", "DNA Methylation Status for Monocot will appear here"),
          h3("DNA Methylation Status"),
          tableOutput("msmono"),
          div(class = "predict-button",
          downloadButton("downloadmonocot", "Download Predictions", class = "btn btn-primary")
          ))),
        mainPanel(
          div(
            class = "predict-panel",
            style = "margin-bottom: 30px;",
            p(class = "predict-description", "DNA Methylation Status for Dicot  will appear here"),
            h3("DNA Methylation Status"),
          tableOutput("msdi"),
          div(class = "predict-button",
              
              downloadButton("downloaddicot", "Download Predictions", class = "btn btn-primary")
          )
        )
      ), includeMarkdown("predict.md")
    ),
    tabPanel(
      "Resource Availability",
      titlePanel("Resources"),
      div(
        tabsetPanel(
          tabPanel(
            "Example Data",
            includeMarkdown("example_data.md"),
            br(),
            downloadButton("downloadButton1", "Download Zipped File", class = "btn btn-yellow")
          ),
          tabPanel(
            "Sequence Manipulation Tools",
            # short description + external link button
            tags$div(style = "padding: 18px;",
                     tags$h4("Sequence Manipulation Suite — Split FASTA"),
                     tags$p("Use the online Split FASTA tool to split large FASTA files into smaller chunks of 41bp as required. Click the button below to open the tool in a new tab.Enter sequence length as 41bp in the Split Fasta tool."),
                     tags$a(
                       "Open Split FASTA (Sequence Manipulation Suite)",
                       href = "https://www.bioinformatics.org/sms2/split_fasta.html",
                       target = "_blank",
                       rel = "noopener noreferrer",
                       class = "btn btn-yellow"
                     ),
                     
                     tags$p(style = "margin-top:8px; color:#ddd;"),
                     tags$small("Note: This opens an external website. Your data will be handled by that external site.")
            )
          )
        ),
        align = "justify"
      ),
    includeMarkdown("predict.md")
  ),
    
    tabPanel("Developers", includeHTML("www/developersSweep.html")),
    tabPanel("Contact Details", includeHTML("www/contact-details.html")),
 
  
  ),
  
  tags$div(
    id = "globalFooter",
    style="
    position: fixed;
    bottom: 0;
    left: 0;
    width: 100%;
    background-color: #001a12;
    padding: 8px 12px;
    z-index: 9999;
  ",
    tags$div(
      style="display:flex; justify-content:center; align-items:center; gap:20px;
           color:white; font-weight:bold;",
      
      HTML("Copyright © ICAR–Indian Agricultural Statistics Research Institute,
          New Delhi–110012. All rights reserved."),
      
      tags$span(id="liveSlot"),
      tags$span(id="totalSlot"),
      tags$span(HTML("<span id='clock'></span>"))
    )
  ),
  
  # ---- JS that injects counters + clock ----
  tags$script(HTML("
  setInterval(function(){

    var c = document.getElementById('clock');
    if(c){
      var now = new Date();
      c.textContent = now.toLocaleString('en-IN',{
        day:'2-digit',month:'short',year:'numeric',
        hour:'2-digit',minute:'2-digit',second:'2-digit'
      });
    }

    if(window.Shiny){
  Shiny.onInputChange('__ping__', Math.random());
  Shiny.setInputValue('__req__', Math.random());
}

  },1000);
")),
  
  tags$script(HTML("
  Shiny.addCustomMessageHandler('counts', function(m){
    var l = document.getElementById('liveSlot');
    var t = document.getElementById('totalSlot');
    if(l) l.textContent = m.live;
    if(t) t.textContent = m.total;
  });
"))
  
)
