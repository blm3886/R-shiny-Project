library(shiny)
#127.0.0.1:7344
#options(shiny.host = '127.0.0.1')
#options(shiny.port = 7344)
#runApp('shinyapp')
#library(shinythemes)
StudyClass = list(`Ectoderm` = list("TCGA Head and Neck Cancer","TCGA Glioblastoma","TCGA Lower Grade Glioma","TCGA Pheochromocytoma & Paraganglioma", #Brain/CNS
                                    "TCGA Melanoma", #Skin
                                    "TCGA Ocular melanomas" #Eye
),
`Mesoderm` = list("TCGA Adrenocortical Cancer", #Adernal Gland
                  "TCGA Kidney Chromophobe","TCGA Kidney Clear Cell Carcinoma","TCGA Kidney Papillary Cell Carcinoma", #Kidney
                  "TCGA Large B-cell Lymphoma", #Lymph
                  "TCGA Acute Myeloid Leukemia", #Blood
                  "TCGA Uterine Carcinosarcoma", "TCGA Endometrioid Cancer", #Uterus
                  "TCGA Sarcoma",
                  "TCGA Ovarian Cancer", "TCGA Testicular Cancaer", #Reproductive System
                  "TCGA Breast Cancer" 
),
`Endoderm` = list("TCGA Bladder Cancer", #Bladder/Urinary Tract
                  "TCGA thymoma", #Thymus
                  "TCGA Lung Adenocarcinoma", "TCGA Lung Squamous Cell Carcinoma", "TCGA Mesothelioma", #Lung
                  "TCGA Pancreatic Cancer", #Pancreas
                  "TCGA Prostate Cancer", #Prostate
                  "TCGA Liver Cancer", #Liver
                  "TCGA Thyroid Cancer", #Thyroid
                  "TCGA Esophageal Cancer","TCGA Stomach Cancer","TCGA Colon Cancer"
)
)
load("UI DATA.RData")
Gene_List = GENES[order(GENES)]
#Cancer_studies = dir("../PANCAN_TOOL/")#c("TCGA Glioblastoma","TCGA Ovarian Cancer","TCGA Lung Adenocarcinoma","TCGA Lung Squamous Cell Carcinoma","TCGA Prostate Cancer","TCGA Endometrioid Cancer","TCGA Bladder Cancer","TCGA Testicular Cancaer","TCGA Esophageal Caner","TCGA Pancreatic Cancer","TCGA Kidney Papillary Cell Carcinoma","TCGA Liver Cancer","TCGA Cervical Cancer","TCGA Sarcoma","TCGA Breast Cancer","TCGA thymoma","TCGA Mesothelioma","TCGA Colon Cancer","TCGA Stomach Cancer","TCGA Melanoma","TCGA Bile Duct Cancer","TCGA Kidney Clear Cell Carcinoma","TCGA Formalin Fixed Paraffin-Embedded Pilot Phase II","TCGA Thyroid Cancer","TCGA Head and Neck Cancer","TCGA Acute Myeloid Leukemia","TCGA Rectal Cancer","TCGA Lower Grade Glioma","TCGA Large B-cell Lymphoma","TCGA Kidney Chromophobe","TCGA Uterine Carcinosarcoma","TCGA Adrenocortical Cancer","TCGA Pheochromocytoma & Paraganglioma","TCGA Ocular melanomas")
#Cancer_studies = sort(Cancer_studies);
#T = read.table("GeneList.txt",sep="\t");
#Gene_List = T$V1#GENES[order(GENES)]
#GraphHeight = "600px";
#GraphWidth = "100%";
#BoxplotOptions = c("Select one: " = "","Age","AJCC Cancer Stage","DFS Status",
#                   "Metastasis Stage","Neoadjuvent","Nodal Stage","Normal vs Tumor",
#                   #"Pathological Tumor Stage",
#                   "Pharmaceutical Therapy","Race","Radiation Treatment","Sex")
#BoxplotOptions = BoxplotOptions[order(BoxplotOptions)];
#P_BoxplotOptions = c("Select one: " = "","Age","DFS Status",
#                   "Metastasis Stage","Neoadjuvent","Nodal Stage",
#                   #"Pathological Tumor Stage",
#                   "Pharmaceutical Therapy","Race","Radiation Treatment","Sex")
#P_BoxplotOptions = P_BoxplotOptions[order(P_BoxplotOptions)];
#SurvOptions = c("Select one: " = "","Age",
#                "Expression","Pharmaceutical Therapy",
#                "Race","Radiation Treatment","Sex")
shinyUI(fluidPage(
  #theme = shinytheme("paper"),
                  #titlePanel("Gupta Lab"),
                  navbarPage(div("Gupta Lab",tags$img(src="Logo.png", style="float:left; margin:10px;margin-top:-12px",height = 50, width = 90)),
                             tabPanel("Home"),
                             tabPanel("DNA"),
                             tabPanel("RNA",
                                      sidebarLayout(
                                        sidebarPanel(
                                          width = "3",
                                          selectInput(inputId = "StudyType",label = "Cancer Study:",choices = c("Select one:"="","One Study","Multiple Studies")
                                          ),
                                          conditionalPanel(condition = "input.StudyType == 'One Study'",
                                                           selectInput(inputId = "Cancerstudy",label = "Select Study:",choices = c("Select one: " = "",StudyClass))),
                                          conditionalPanel(condition = "input.StudyType == 'Multiple Studies'",
                                                           selectInput(inputId = "M_Cancerstudy",label = "Select Study:",choices = c("Select any: " = "",StudyClass),multiple = TRUE)),
                                          selectizeInput(inputId = "Genes",label = "Gene of Interest:",choices = as.character(Gene_List),multiple = TRUE, options = list(maxItems = 6)),
                                          conditionalPanel(condition = "input.StudyType != ''",
                                                           selectInput(inputId = "Output_opt",
                                                                       label = "Analysis:",
                                                                       choices = c("Select one: " = "","Heatmap","Boxplot","Survival")),
                                                           #### HEAT ####
                                                           conditionalPanel(
                                                             condition = "input.Output_opt == 'Heatmap' && input.Genes.length > 1 && input.StudyType == 'One Study'",
                                                             selectInput(inputId = "HeatType",
                                                                         label = "Heatmap for:",
                                                                         choices = c("Select one: " = "",
                                                                                     "Normal Vs Tumor",
                                                                                     "Pharmaceutical Therapy",
                                                                                     "Radiation Treatment")
                                                             ),
                                                             conditionalPanel(condition = "input.HeatType == 'Radiation Treatment' || input.HeatType == 'Pharmaceutical Therapy'",
                                                                              selectInput(inputId = "ZScore", label = "Z-Score",choices = c(0.5,1,1.5,2),selected = 1))
                                                           ),
                                                           
                                                           #### BOX #####
                                                           conditionalPanel(
                                                             condition = "input.Output_opt == 'Boxplot' && input.StudyType == 'One Study'",
                                                             selectInput(inputId = "BoxType",
                                                                         label = "Boxplot for:",
                                                                         choices = BoxplotOptions
                                                             )  
                                                           ),
                                                           conditionalPanel(
                                                             condition = "input.Output_opt == 'Boxplot' && input.StudyType == 'Multiple Studies'",
                                                             selectInput(inputId = "M_BoxType",
                                                                         label = "Boxplot for:",
                                                                         choices = c("Select one: " = "","Normal Samples","Tumor Samples")
                                                             )  
                                                           ),
                                                           #### SURVIVAL ####
                                                           conditionalPanel(
                                                             condition = "input.Output_opt == 'Survival' && input.StudyType == 'One Study' ",
                                                             selectInput(inputId = "SurvType",
                                                                         label = "Survival wrt:",
                                                                         choices = SurvOptions
                                                             ),
                                                             conditionalPanel(condition = "input.SurvType == 'Expression'",
                                                                              selectInput(inputId = "Z_Score", label = "Z-Score",choices = c(0.5,1,1.5,2),selected = 1)),
                                                             selectInput(inputId = "SurvivalTime", 
                                                                         label = "Survival Time",
                                                                         choices = c("Median Survival","Time of your choice"),
                                                                         selected = "Median Survival"),
                                                             conditionalPanel(condition="input.SurvivalTime == 'Time of your choice'",
                                                                              sliderInput(inputId = "Time", label = "Time (in Days): ", min = 1,max = 100,value = 50,step = 2))
                                                           )
                                                           #####
                                          ),
                                          actionButton("AnalyseRNA","Analyze!")
                                        ),
                                        mainPanel(
                                          width = "9",
                                          #### HEAT OUT ####
                                          conditionalPanel(condition = "input.Output_opt=='Heatmap' && (input.Cancerstudy != '' || input.M_Cancerstudy != '' ) && input.Genes.length > 1",
                                                           conditionalPanel(condition = "input.HeatType == 'Pharmaceutical Therapy' || input.HeatType == 'Radiation Treatment'",plotOutput("R_Heat_T", width = GraphWidth, height = GraphHeight)),
                                                           #conditionalPanel(condition = "input.HeatType == 'Radiation Treatment'",plotOutput("HeatRT", width = GraphWidth, height = GraphHeight)),
                                                           conditionalPanel(condition = "input.HeatType == 'Normal Vs Tumor'",plotOutput("HeatNT" , width = GraphWidth, height = GraphHeight))
                                          ),
                                          conditionalPanel(condition = "input.Output_opt=='Heatmap' && (input.Cancerstudy != '' || input.M_Cancerstudy != '' ) && input.Genes.length <= 1",
                                                           textOutput("HeatText")
                                          ),
                                          #### BOX & SURV OUT ####
                                          conditionalPanel(condition = "(input.Output_opt=='Boxplot' || input.Output_opt=='Survival')  && input.Genes.length != 0 && (input.Cancerstudy != '' || input.M_Cancerstudy != '' )",
                                                           uiOutput("OutPlot")
                                          ),
                                          #####
                                          conditionalPanel(condition = "input.Output_opt=='Boxplot' && input.Genes.length != 0 && input.M_Cancerstudy != ''",
                                                           conditionalPanel(condition = "input.M_BoxType == 'Normal Samples' || input.M_BoxType == 'Tumor Samples'",plotOutput("BoxplotALL", width = GraphWidth, height = GraphHeight))
                                          )
                                        )
                                      )
                             ),
                             tabPanel("Protein",
                                      sidebarLayout(
                                        sidebarPanel(
                                          width = "3",
                                          selectInput(inputId = "P_StudyType",label = "Cancer Study:",choices = "One Study"),
                                          conditionalPanel(condition = "input.P_StudyType == 'One Study'",
                                                           selectInput(inputId = "P_Cancerstudy",label = "Select Study:",choices = c("Select one: " = "",StudyClass))),
                                          #conditionalPanel(condition = "input.P_StudyType == 'Multiple Studies'",
                                          #                selectInput(inputId = "P_M_Cancerstudy",label = "Select Study:",choices = c("Select any: " = "",Cancer_studies),multiple = TRUE)),
                                          selectInput(inputId = "P_Genes",label = "Gene of Interest:",choices = as.character(Gene_List),multiple = TRUE
                                                      #list(`E` = list("V","B","N"),`H` = list("A","S","D")),
                                          ),
                                          conditionalPanel(condition = "input.P_StudyType != ''",
                                                           selectInput(inputId = "P_Output_opt",
                                                                       label = "Analysis:",
                                                                       choices = c("Select one: " = "","Heatmap","Survival")),
                                                           #### HEAT ####
                                                           conditionalPanel(
                                                             condition = "input.P_Output_opt == 'Heatmap' && input.P_Genes.length > 1",
                                                             selectInput(inputId = "P_HeatType",
                                                                         label = "Heatmap for:",
                                                                         choices = c("Select one: " = "", "All",
                                                                                     "Pharmaceutical Therapy",
                                                                                     "Radiation Treatment")
                                                             ),
                                                             selectInput(inputId = "P_ZScore", label = "Z-Score",choices = c(0.5,1,1.5,2),selected = 1)
                                                           ),
                                                           
                                                           #### BOX #####
                                                           #conditionalPanel(
                                                          #   condition = "input.P_Output_opt == 'Boxplot' && input.P_StudyType == 'One Study'",
                                                          #   selectInput(inputId = "P_BoxType",
                                                          #               label = "Boxplot for:",
                                                          #               choices = P_BoxplotOptions
                                                          #   )  
                                                          # ),
                                                           #conditionalPanel(
                                                           #  condition = "input.P_Output_opt == 'Boxplot' && input.P_StudyType == 'Multiple Studies'",
                                                           #  selectInput(inputId = "P_M_BoxType",
                                                           #              label = "Boxplot for:",
                                                           #              choices = c("Select one: " = "","Normal Samples","Tumor Samples")
                                                           #  )  
                                                           #),
                                                           #### SURVIVAL ####
                                                           conditionalPanel(
                                                             condition = "input.P_Output_opt == 'Survival' && input.P_StudyType == 'One Study'",
                                                             selectInput(inputId = "P_SurvType",
                                                                         label = "Survival wrt:",
                                                                         choices = SurvOptions
                                                             ),
                                                             conditionalPanel(condition = "input.P_SurvType == 'Expression'",
                                                                              selectInput(inputId = "P_Z_Score", label = "Z-Score",choices = c(0.5,1,1.5,2),selected = 1)),
                                                             selectInput(inputId = "P_SurvivalTime", 
                                                                         label = "Survival Time",
                                                                         choices = c("Median Survival","Time of your choice"),
                                                                         selected = "Median Survival"),
                                                             conditionalPanel(condition="input.P_SurvivalTime == 'Time of your choice'",
                                                                              sliderInput(inputId = "Time", label = "Time (in Days): ", min = 1,max = 100,value = 50,step = 2))
                                                           )
                                                           #####
                                          )
                                        ),
                                        
                                        mainPanel(
                                          width = "9",
                                          #### HEAT OUT ####
                                          conditionalPanel(condition = "input.P_Output_opt=='Heatmap' && input.P_Cancerstudy != '' && input.P_Genes.length > 1",
                                                           conditionalPanel(condition = "input.P_HeatType == 'Pharmaceutical Therapy' || input.P_HeatType == 'Radiation Treatment' || input.P_HeatType == 'All'",plotOutput("P_Heat_T", width = GraphWidth, height = GraphHeight))
                                          ),
                                          conditionalPanel(condition = "input.P_Output_opt=='Heatmap' && input.P_Cancerstudy != '' && input.P_Genes.length <= 1",
                                                           textOutput("P_HeatText")
                                          ),
                                          #### BOX & SURV OUT ####
                                          conditionalPanel(condition = "(input.P_Output_opt=='Boxplot' || input.P_Output_opt=='Survival')  && input.P_Genes.length != 0 && input.P_Cancerstudy != ''",
                                                           uiOutput("P_OutPlot")
                                          )
                                          #####
                                        )
                                      )
                             )
                  )
)
)