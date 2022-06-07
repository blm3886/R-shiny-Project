library(shiny)
library(gplots)
library(KMsurv)
library(survival)
library(survminer)
library(gridExtra)
library(cgdsr)

Gene_List = read.table("H2A Isoform.txt",sep="\t");
Gene_List = Gene_List$V1;
#Gene_List = GENES[order(GENES)]
load("ServerData.RData")
GetNormal = function(GOI, Type){
  Normal = data.frame()
  for(i in 1:length(GOI)){
    ROWIND = which(Gene_Info == GOI[i],arr.ind = TRUE)[1]
    ROWIND = GIndexes[ROWIND]
    FilePath = gsub(pattern = "TYPE",replacement = Type,x = "../PANCAN_TOOL/TYPE/Expression Data/Expression_TYPE_CT_GENEINT.csv");
    FilePath = gsub(pattern = "GENEINT",replacement = ROWIND,x = FilePath);
    if(!is.na(ROWIND)){
      Nor = read.csv(file = sub("CT","Normal",FilePath),header = TRUE,row.names = 1)
      if(i == 1){
        Normal = Nor[GOI[i]]
      }else{
        Normal = cbind(Normal,Nor[GOI[i]])
      }
    }
  }
  Normal = `row.names<-`(Normal,row.names(Nor[1]))
  return(Normal)
  
}
GetTumor = function(GOI, Type){
  Tumor = data.frame()
  for(i in 1:length(GOI)){
    ROWIND = which(Gene_Info == GOI[i],arr.ind = TRUE)[1]
    ROWIND = GIndexes[ROWIND]
    FilePath = gsub(pattern = "TYPE",replacement = Type,x = "../PANCAN_TOOL/TYPE/Expression Data/Expression_TYPE_CT_GENEINT.csv");
    FilePath = gsub(pattern = "GENEINT",replacement = ROWIND,x = FilePath);
    if(!is.na(ROWIND)){
      Tum = read.csv(file = sub("CT","Tumor",FilePath),header = TRUE,row.names = 1)
      if(i == 1){
        Tumor = Tum[GOI[i]]
      }else{
        Tumor = cbind(Tumor,Tum[GOI[i]])
      }
    }
  }
  Tumor = `row.names<-`(Tumor,row.names(Tum[1]))
  return(Tumor)
}
GetGeneticProfile = function(GOI, Type){
  GP = data.frame()
  Selected_Study = as.character(ProCT$cBioStudy[ProCT$Cancer_Type == Type])
  for(i in 1:length(GOI)){
    ROWIND = which(Gene_Info == GOI[i],arr.ind = TRUE)[1]
    ROWIND = GIndexes[ROWIND]
    if(!is.na(ROWIND)){
      FilePath = gsub(pattern = "FILE",replacement = Selected_Study,x = "../PANCAN_TOOL/TYPE/Expression Data/cBio-Expression_GENE_FILE_ZScore.csv");
      FilePath = gsub(pattern = "TYPE",replacement = Type,x = FilePath);
      FilePath = gsub(pattern = "GENE",replacement = ROWIND,x = FilePath);
      gp = read.csv(file = FilePath,header = TRUE,row.names = 1)
      if(i == 1){
        GP = gp[GOI[i]]
      }else{
        GP = cbind(GP,gp[GOI[i]])
      }
    }
  }
  GP = `row.names<-`(GP,row.names(gp[1]))
  return(GP)
}
P_GetGeneticProfile = function(GOI, Type){
  GP = data.frame()
  Selected_Study = as.character(ProCT$cBioStudy[ProCT$Cancer_Type == Type])
  if(length(grep("Protein",dir(gsub(pattern = "TYPE",Type,"../PANCAN_TOOL/TYPE/")))) != 0){
    for(i in 1:length(GOI)){
      ROWIND = which(Gene_Info == GOI[i],arr.ind = TRUE)[1]
      ROWIND = GIndexes[ROWIND]
      if(!is.na(ROWIND)){
        FilePath = gsub(pattern = "TYPE",replacement = Type,x = "../PANCAN_TOOL/TYPE/Protein/");
        Pat = gsub(pattern = "SELSTUDY",replacement = Selected_Study,x = "^SELSTUDY.*_ROW\\.csv")
        Pat = gsub(pattern = "ROW",replacement = ROWIND,x = Pat)
        FileTOR = dir(FilePath)[grep(Pat,dir(FilePath))]
        gp = read.csv(file = paste0(FilePath,FileTOR),header = TRUE,row.names = 1)
        if(match(GOI[i],colnames(gp),nomatch = 0) != 0){
          if(i == 1){
            GP = gp[GOI[i]]
          }else{
            GP = cbind(GP,gp[GOI[i]])
          }
        }
      }
    }
    GP = `row.names<-`(GP,row.names(gp[1]))
  }else GP = NULL
  return(GP)
}
GetExpression = function(GOI, Type){
  Expression = data.frame()
  for(i in 1:length(GOI)){
    ROWIND = which(Gene_Info == GOI[i],arr.ind = TRUE)[1]
    ROWIND = GIndexes[ROWIND]
    EFilePath = gsub(pattern = "TYPE",replacement = Type,x = "../PANCAN_TOOL/TYPE/Expression Data/Expression_TYPE_GENEINT.csv");
    EFilePath = gsub(pattern = "GENEINT",replacement = ROWIND,x = EFilePath);
    if(!is.na(ROWIND)){
      Exp = read.csv(file = EFilePath,header = TRUE,row.names = 1);
      if(i == 1){
        Expression = Exp[GOI[i]]
      }else{
        Expression = cbind(Expression, Exp[GOI[i]])
      }
    }
  }
  Expression = `row.names<-`(Expression,row.names(Exp[1]))
  return(Expression)
}
GetDataCBP = function(DataCode,Selected_Study,Gene_List){
  mycgds = CGDS("http://www.cbioportal.org/");
  AllCaselist = getCaseLists(mycgds,Selected_Study);
  AllGeneticProfile = getGeneticProfiles(mycgds,Selected_Study)
  if(DataCode == 1 || DataCode == 2){
    Case_row = grep("RNA Seq V2",ignore.case = TRUE, AllCaselist$case_list_name);
    if(length(Case_row)>0){
      mycaselist = AllCaselist[Case_row,1];
      #print(paste('Getting Clinical Data for the case list: ', mycaselist));
      myclinicaldata = getClinicalData(mycgds,mycaselist);
      GenePro_row = grep("rna_seq_v2_mrna_median_Zscores",ignore.case = TRUE, AllGeneticProfile$genetic_profile_id);
      if(length(GenePro_row)>0 && DataCode == 2){
        mygeneticprofile = AllGeneticProfile[GenePro_row,1];
        #print(paste("Getting Expression data: " , mygeneticprofile));
        GeneticProfile = getProfileData(mycgds,Gene_List,mygeneticprofile,mycaselist);
      }#else GeneticProfile = NULL
    }#else myclinicaldata = NULL
  }
  if(DataCode == 3 || DataCode ==4){
    Prot_Case_row = grep("RPPA Data",ignore.case = TRUE, AllCaselist$case_list_name);
    if(length(Prot_Case_row) != 0){
      Prot_mycaselist = AllCaselist[Prot_Case_row,1];
      P_myclinicaldata_Prot = getClinicalData(mycgds,Prot_mycaselist);
      Pro_row = grep("protein_quantification_zscores",ignore.case = TRUE, AllGeneticProfile$genetic_profile_id);
      if(length(Pro_row) == 0)Pro_row = grep("rppa_zscores",ignore.case = TRUE, AllGeneticProfile$genetic_profile_id)
      if(length(Pro_row) != 0){
        P_mygeneticprofile = AllGeneticProfile[Pro_row,1];
        P_GeneticProfile = getProfileData(mycgds,Gene_List,P_mygeneticprofile,Prot_mycaselist);
      }#else P_GeneticProfile = NULL
    }#else P_myclinicaldata_Prot = NULL
  }
  if(DataCode == 1)Data = myclinicaldata
  if(DataCode == 2)Data = GeneticProfile
  if(DataCode == 3)Data = P_myclinicaldata_Prot
  if(DataCode == 4)Data = P_GeneticProfile
  return(Data)
}
Generate_surv = function(KM,survival_data,GOI,Legend_Title,Legend_label,Time){
  X = as.data.frame(x = KM$strata);
  Plot_Question = 1
  if(Time == 0){
    if(length(KM$surv[(KM$surv < 0.50) == TRUE])>0){
      Plot = ggsurvplot(KM, data = survival_data, 
                        title = paste('Kaplan-Meier Estimate for',GOI),
                        surv.median.line = 'hv',
                        font.main = 22, font.x = c(18,"bold"), font.y = c(18,"bold"),
                        linetype = "solid", censor=TRUE, censor.shape=134,censor.size = 6.5, 
                        palette= "Set1",
                        pval = TRUE, pval.coord = c(ceiling(max(survival_data$OS_MONTHS,na.rm = T))/2,0.97), pval.size= 6,
                        legend.title= Legend_Title, legend.labs=Legend_label,
                        axes.offset = FALSE, xlim = c(0,(ceiling(max(survival_data$OS_MONTHS,na.rm = T))+5)),
                        ggtheme = theme_survminer(font.tickslab = c(16,"bold","black"),font.legend = 16) + theme(
                          plot.title = element_text(hjust = 0.5),
                          plot.margin = margin(25,25,25,25),
                          plot.subtitle = element_text(hjust = 0.5,size = 14),
                          panel.border = element_rect(linetype = "solid",fill = NA),
                          #legend.box.background = element_rect(fill=NA,linetype="blank",size=0),
                          legend.margin = margin(-5,0,0,0),
                          legend.justification = c("right", "top"),
                          panel.grid.major = element_line(linetype = "dashed",colour = "grey85"),
                          panel.grid.minor = element_line(linetype = "dashed",colour = "grey90"),
                          axis.ticks.length = unit(.25, "cm")
                        )
      );
    }else{
      Plot = ggsurvplot(KM, data = survival_data, 
                        title = paste('Kaplan-Meier Estimate for',GOI),
                        #surv.median.line = 'hv',
                        font.main = 22, font.x = c(18,"bold"), font.y = c(18,"bold"),
                        linetype = "solid", censor=TRUE, censor.shape=134,censor.size = 6.5, 
                        palette= "Set1",
                        pval = TRUE, pval.coord = c(ceiling(max(survival_data$OS_MONTHS,na.rm = T))/2,0.97), pval.size= 6,
                        legend.title= Legend_Title, legend.labs=Legend_label,
                        axes.offset = FALSE, xlim = c(0,(ceiling(max(survival_data$OS_MONTHS,na.rm = T))+5)),
                        ggtheme = theme_survminer(font.tickslab = c(16,"bold","black"),font.legend = 16) + theme(
                          plot.title = element_text(hjust = 0.5),
                          plot.margin = margin(25,25,25,25),
                          plot.subtitle = element_text(hjust = 0.5,size = 14),
                          panel.border = element_rect(linetype = "solid",fill = NA),
                          #legend.box.background = element_rect(fill=NA,linetype="blank",size=0),
                          legend.margin = margin(-5,0,0,0),
                          legend.justification = c("right", "top"),
                          panel.grid.major = element_line(linetype = "dashed",colour = "grey85"),
                          panel.grid.minor = element_line(linetype = "dashed",colour = "grey90"),
                          axis.ticks.length = unit(.25, "cm")
                        )
      );
      
    }
  }else{
    if(Time > 0){
      surv_prob = as.vector(summary(KM,time = Time)$surv);
      prob_df = data.frame(x1 = rep(Time,length(surv_prob)) , x2 = rep(Time,length(surv_prob)) , y1 = 0 , y2 = surv_prob);
      Label = "(X,Y)";
      for(i in 1:length(surv_prob)){
        Label[i] = sub("X",Time,"(X,Y)");
        Label[i] = sub("Y",signif(surv_prob[i],2),Label[i]);
      }
      Plot = ggsurvplot(KM, data = survival_data, 
                        title = paste('Kaplan-Meier Estimate for',GOI),
                        font.main = 22, font.x = c(18,"bold"), font.y = c(18,"bold"),
                        linetype = "solid", censor=TRUE, censor.shape=134,censor.size = 6.5, 
                        palette= "Set1",
                        pval = TRUE, pval.coord = c(ceiling(max(survival_data$OS_MONTHS,na.rm = T))/2,0.97), pval.size= 6,
                        legend.title= Legend_Title, legend.labs=Legend_label,
                        axes.offset = FALSE, xlim = c(0,(ceiling(max(survival_data$OS_MONTHS,na.rm = T))+5)),
                        ggtheme = theme_survminer(font.tickslab = c(16,"bold","black"),font.legend = 16) + theme(
                          plot.title = element_text(hjust = 0.5),
                          plot.margin = margin(25,25,25,25),
                          plot.subtitle = element_text(hjust = 0.5,size = 14),
                          panel.border = element_rect(linetype = "solid",fill = NA),
                          #legend.box.background = element_rect(fill=NA,linetype="blank",size=0),
                          legend.margin = margin(-5,0,0,0),
                          legend.justification = c("right", "top"),
                          panel.grid.major = element_line(linetype = "dashed",colour = "grey85"),
                          panel.grid.minor = element_line(linetype = "dashed",colour = "grey90"),
                          axis.ticks.length = unit(.25, "cm")
                        )
      );
      Plot$plot = Plot$plot + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
                                           linetype = "dashed", size = 0.6, color="black",data = prob_df)+ #Vertical
        geom_segment(aes(x = x1, y = y2, xend = 0,  yend = y2), linetype = "dashed", size = 0.6,color="black",data = prob_df)+ #Horizontal
        geom_text(mapping = aes(x = x1,y = y2),data=prob_df,size = 4.5,color = "black",label = Label,
                  nudge_x = rep(4.5,length(surv_prob)),nudge_y = rep(0.02,length(surv_prob)))+
        geom_point(mapping = aes(x=x1,y=y2),data = prob_df,color="green");
    }
  }
  print(Plot$plot)
}

#Gene_List = GENES
shinyServer(function(input, output, session) {
######## RNA ######
  GOI <- eventReactive(input$AnalyseRNA,{
    validate(
      need(expr = input$Genes != "",'Please select a gene of interest')
    )
    input$Genes
  })
  Type <- eventReactive(input$AnalyseRNA,{
    validate(need(expr = input$Cancerstudy != "",message = "Please Select a Study"))
    input$Cancerstudy
  })
  M_Type <- eventReactive(input$AnalyseRNA,{
    validate(need(expr = input$M_Cancerstudy != "",message = "Please Select a Study"))
    input$M_Cancerstudy
  })
  myclinicaldata <- eventReactive(input$AnalyseRNA,{
    Type = input$Cancerstudy
    Selected_Study = as.character(ProCT$cBioStudy[ProCT$Cancer_Type == Type])
    FilePath = gsub(pattern = "TYPE",replacement = Type,x = "../PANCAN_TOOL/TYPE/Clinical Data/cBio-Clinical_FILE.csv");
    myclinicaldata = read.csv(file = sub("FILE",Selected_Study,FilePath),header = TRUE,row.names = 1)
    return(myclinicaldata);
  })
  Time <- reactive({
    if(input$SurvivalTime == "Time of your choice"){
      myclinicaldata = myclinicaldata()
      switch (input$SurvType,
              "Age" = {
                myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = "",Code_type = 3,Cutoff = 0);
                UN = unique(myclinicaldata$AGE_GROUP)
                if(match(1,UN,0) != 0)M1 = max(myclinicaldata$OS_MONTHS[myclinicaldata$AGE_GROUP == 1],na.rm = TRUE)else M1 = NA
                if(is.na(M1))M2 = max(myclinicaldata$OS_MONTHS[myclinicaldata$AGE_GROUP == 2],na.rm = TRUE)else M2 = NA
                if(is.na(M2))M3 = max(myclinicaldata$OS_MONTHS[myclinicaldata$AGE_GROUP == 3],na.rm = TRUE)else M3 = NA
                if(is.na(M3))M4 = max(myclinicaldata$OS_MONTHS[myclinicaldata$AGE_GROUP == 4],na.rm = TRUE)else M4 = NA
                if(is.na(M4))M5 = max(myclinicaldata$OS_MONTHS[myclinicaldata$AGE_GROUP == 5],na.rm = TRUE)else M5 = NA
                MaxVal = as.integer(max(c(M1,M2,M3,M4,M5),na.rm = TRUE)) 
              },
              "Expression" = {
                #myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = 1,Cutoff = as.numeric(input$Z_Score),GOI = GOI);
                #High = grep("High",ignore.case = TRUE, myclinicaldata$Expression)
                #Low = grep("Low",ignore.case = TRUE, myclinicaldata$Expression)
                
              },
              "Pharmaceutical Therapy" = {
                myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = "",Cutoff = 0);
                Yes = grep("Yes",ignore.case = TRUE, myclinicaldata$PHARMACEUTICAL_TX_ADJUVANT);
                No = grep("No",ignore.case = TRUE, myclinicaldata$PHARMACEUTICAL_TX_ADJUVANT);
                if(length(Yes)>0&&length(No)>0){
                  MaxVal = min(c(as.integer(max(myclinicaldata$OS_MONTHS[Yes],na.rm = TRUE)),as.integer(max(myclinicaldata$OS_MONTHS[No],na.rm = TRUE))))
                }
                else MaxVal = NA
              },
              "Race" = {
                myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,"",0);
                RACES = as.character(unique(myclinicaldata$RACE));
                Index_RACE = c();
                for(Pat in RACES){
                  if(Pat != ""){
                    Index_RACE = c(Index_RACE,grep(pattern = Pat,x = myclinicaldata$RACE));
                  }
                }
                MaxVal = as.integer(min(c(myclinicaldata$OS_MONTHS[Index_RACE],1),na.rm = TRUE))
              },
              "Radiation Treatment" = ,
              "Sex" = {
                myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = 2,Cutoff = 0);
                UN = unique(myclinicaldata$SEX)
                if(match(1,UN,0) != 0)M1 = max(myclinicaldata$OS_MONTHS[myclinicaldata$SEX == 1],na.rm = TRUE)else M1 = NA
                if(match(2,UN,0) != 0)M2 = max(myclinicaldata$OS_MONTHS[myclinicaldata$SEX == 2],na.rm = TRUE)else M2 = NA
                MaxVal = as.integer(min(c(M1,M2),na.rm = TRUE)) 
              },
              {MaxVal = as.integer(max(myclinicaldata$OS_MONTHS,na.rm = TRUE))}
      )
      updateSliderInput(session, "Time", value = input$Time,
                        min = 1, max = MaxVal, step = 1)
      Time = input$Time;
    }else Time = 0
  })
  Normal <- eventReactive(input$AnalyseRNA,{
    Type = Type()
    GOI = GOI()
    Normal = GetNormal(GOI,Type);
    return(Normal);
  })
  Tumor <- eventReactive(input$AnalyseRNA,{
    Type = Type()
    GOI = GOI()
    Tumor = GetTumor(GOI,Type);
    return(Tumor)
  })
  Expression <- eventReactive(input$AnalyseRNA,{
    Type = Type()
    GOI = GOI()
    Expression = GetExpression(GOI,Type);
    return(Expression)
  })
  GeneticProfile <- eventReactive(input$AnalyseRNA,{
    Type = Type()
    GOI = GOI()
    GeneticProfile = GetGeneticProfile(GOI,Type);
    return(GeneticProfile)
  })
  #### TABS ####
  pObjects <- eventReactive(input$AnalyseRNA,{
    J <- length(GOI())
    Tabnames <- paste0(GOI()) 
    Plot_type <- input$Output_opt
    Study <- input$StudyType
    Analysis = "";
    Type = Type();
    if(Plot_type == "Boxplot")Analysis <- input$BoxType
    if(Plot_type == "Survival")Analysis <- input$SurvType
    list(Gene_Length=J, Gene_List = Tabnames, Plot = Plot_type, Analysis = Analysis,Type = Type, Study = Study)
  })
  outputNodes <- reactive({ # output node names
    pobjects <- pObjects()
    if (is.null(pobjects)) return(NULL)  
    J <- pobjects$Gene_Length
    list(pnodes=paste0("pnode", LETTERS[1:J])) # plot outputs
  })
  observe({ 
    pobjects <- pObjects()
    if (!is.null(pobjects)) {
      outnodes <- outputNodes()
      pnodes <- outnodes$pnodes
      tests <- pobjects$Gene_List
      Plot <- pobjects$Plot
      Study <- pobjects$Study
      J <- pobjects$Gene_Length
      Analysis <- pobjects$Analysis
      Type <- pobjects$Type;
      ## tab 1, 2, ..., J
      I <- input$tab1
      GOI <- tests[as.numeric(I)]
      if(!is.null(Plot) && !is.null(GOI) && length(GOI) != 0 && !is.null(Type)){
        #### BOXPLOT ####
        if(Plot == "Boxplot"){
          if(Study == "One Study"){
            #FilePath = gsub(pattern = "TYPE",replacement = Type,x = "../PANCAN_TOOL/TYPE/Expression Data/Expression_TYPE_CT.csv");
            Normal = Normal()#read.csv(file = sub("CT","Normal",FilePath),header = TRUE,row.names = 1)
            Tumor = Tumor()#read.csv(file = sub("CT","Tumor",FilePath),header = TRUE,row.names = 1)
            Expression = Expression()#read.csv(file = EFilePath,header = TRUE,row.names = 1);
            Selected_Study = as.character(ProCT$cBioStudy[ProCT$Cancer_Type == Type])
            myclinicaldata = myclinicaldata()#GetDataCBP(1,Selected_Study,Gene_List = GOI);
            if(!is.null(myclinicaldata)){
              #### AGE ####
              if(Analysis == "Age"){
                myclinical = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = "", Code_type = 4);
                AG = sort(unique(myclinical$AGE_GROUP));
                AG1 = row.names(myclinical[grep(1,myclinical$AGE_GROUP),]);
                if(length(AG1)!= 0)AG1 = Search_Clin_Exp(AG1,Expression)else rm(AG1)
                AG2 = row.names(myclinical[grep(2,myclinical$AGE_GROUP),]); 
                if(length(AG2)!= 0)AG2 = Search_Clin_Exp(AG2,Expression)else rm(AG2)
                if(J>0){
                  if(match(GOI,colnames(Tumor),nomatch = 0) != 0){
                    if(length(AG1)!= 0 & length(AG2)!= 0){
                      for(i in 1:J){
                        if(i == as.numeric(I)){
                          output[[pnodes[i]]] <- renderPlot({ # plot in each tab
                            AG11 = AG1[GOI]
                            AG22 = AG2[GOI]
                            BoxData = c(AG11,AG22);
                            BoxNames = c(paste("0-45\nSample-size: ",length(AG11[,1]),sep = ""),paste("46-100\nSample-size: ",length(AG22[,1]),sep =""));
                            BoxPlotName = paste(GOI,Type,"Boxplot AGE GROUP.tiff");
                            BoxTitle = paste(Type,"\nAge Group Boxplot")
                            BoxLegend = Calc_P_val_TwoSample(Data1 = AG11,Data2 = AG22,"0-45 Vs 46-100");
                            BoxColor = c("green","red");
                            Generate_Box(GOI,BoxPlotName,BoxData,BoxNames,BoxColor,BoxTitle,BoxLegend);
                          },width = PlotWidth,height = PlotHeight)
                        }
                      }
                    }else{
                      for(i in 1:J){ 
                        output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                      }
                    }
                  }else{
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              #### CANCER STAGE ####
              if(Analysis == "AJCC Cancer Stage"){
                myclinical = myclinicaldata
                if(match("AJCC_PATHOLOGIC_TUMOR_STAGE",colnames(myclinical),nomatch = 0) !=0 ){myclinical$Cancer_stage = as.character(myclinical$AJCC_PATHOLOGIC_TUMOR_STAGE)
                }else{
                  if(match("CLINICAL_STAGE",colnames(myclinical),nomatch = 0) !=0 ) myclinical$Cancer_stage = as.character(myclinical$CLINICAL_STAGE)
                }
                myclinical = tryCatch({
                  Cancer_stage_handle(Dataframe = myclinical,
                                      Pattern = c("Stage I[A|B|C]","Stage II[A|B|C]","Stage III[A|B|C]","Stage IV[A|B|C]"),
                                      Replacement = c("Stage I","Stage II","Stage III","Stage IV"))
                },
                error = function(cond){
                  textOutput("Error in Stage Data")
                },
                warning = function(cond){
                  textOutput("Warning: Stage Data")
                }
                )
                Cancer_Stage = c("Stage I", "Stage II", "Stage III", "Stage IV");
                Indexes = c();
                Complete_indexes = c();
                Index_ID = c();
                Exp_Indexes = c();
                for(Q in 1:length(Cancer_Stage)){
                  if(Cancer_Stage[Q]!=""){
                    for(ZX in 1:length(myclinical[,1])){
                      if(as.character(myclinical$Cancer_stage[ZX]) == as.character(Cancer_Stage[Q])){
                        Indexes[length(Indexes)+1] = ZX;
                        Match_exp = grep(as.character(row.names(myclinical[ZX,])),row.names(Expression));
                        if(length(Match_exp)>0){
                          Exp_Indexes[length(Exp_Indexes)+1] = Match_exp;
                          Complete_indexes[length(Complete_indexes)+1] = Match_exp;
                        }
                      }
                    }
                    if(length(Indexes)>0){
                      assign(paste(Cancer_Stage[Q],sep = ""),Expression[Exp_Indexes,])
                    }
                    Index_ID[length(Index_ID)+1] = length(Exp_Indexes);
                    Indexes = c();
                    Exp_Indexes = c();
                  }
                }
                ##OP##
                if(J>0){
                  if(match(GOI,colnames(Tumor),nomatch = 0) != 0){
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot({
                        BoxPlotName = paste(GOI,Type,"Boxplot Stagewise.tiff");
                        BoxData = c(Normal[GOI],`Stage I`[GOI],`Stage II`[GOI],`Stage III`[GOI],`Stage IV`[GOI]);
                        BoxNames = c(paste("Normal\n(",length(Normal[,GOI]),")",sep = ""),
                                     paste("Stage I\n(",length(`Stage I`[,GOI]),")",sep = ""),
                                     paste("Stage II\n(",length(`Stage II`[,GOI]),")",sep = ""),
                                     paste("Stage III\n(",length(`Stage III`[,GOI]),")",sep = ""),
                                     paste("Stage IV\n(",length(`Stage IV`[,GOI]),")",sep = "")
                        );
                        BoxColor = rainbow(7);
                        BoxTitle = paste(Type,"\nCancer Stagewise Boxplot");
                        BoxLegend = c(Calc_P_val_TwoSample(Data1 = Normal[GOI],Data2 = `Stage I`[GOI],"Normal Vs Stage I"),
                                      Calc_P_val_TwoSample(Data1 = Normal[GOI],Data2 = `Stage II`[GOI],"Normal Vs Stage II"),
                                      Calc_P_val_TwoSample(Data1 = Normal[GOI],Data2 = `Stage III`[GOI],"Normal Vs Stage III"),
                                      Calc_P_val_TwoSample(Data1 = Normal[GOI],Data2 = `Stage IV`[GOI],"Normal Vs Stage IV"))
                        Generate_Box(GOI,BoxPlotName,BoxData,BoxNames,BoxColor,BoxTitle,BoxLegend);
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              #### DFS ####
              if(Analysis == "DFS Status"){
                myclinical = myclinicaldata
                if(length(grep("DFS_STATUS",x = colnames(myclinical),ignore.case = TRUE))>0){
                  Yes = row.names(myclinical[grep("DiseaseFree",ignore.case = TRUE, myclinical$DFS_STATUS),]);
                  No = row.names(myclinical[grep("Recurred/Progressed",ignore.case = TRUE, myclinical$DFS_STATUS),]);
                  if(length(Yes)!= 0)Yes = Search_Clin_Exp(Yes,Expression)
                  if(length(No)!= 0)No = Search_Clin_Exp(No,Expression)
                  if(match(GOI,colnames(Tumor),nomatch = 0) != 0){
                    if(length(Yes[GOI]) > 0 & length(No[GOI]) > 0){
                      for(i in 1:J){ 
                        output[[pnodes[i]]] <- renderPlot({ # plot in each tab
                          BoxData = c(No[GOI],Yes[GOI]);
                          BoxLegend = Calc_P_val_TwoSample(Data1 = No[GOI],Data2 = Yes[GOI],"Recurred/Progressed Vs Disease Free");
                          BoxTitle = paste(Type,"\nDFS Status Boxplot")
                          BoxNames = c(paste("Recurred/Progressed\nSample-size: ",length(No[,GOI]),sep = ""),paste("Disease Free\nSample-size: ",length(Yes[,GOI]),sep = ""))
                          BoxColor = c("red","green");
                          BoxPlotName = paste(GOI,Type,"Boxplot Radiation treatment.tiff");
                          Generate_Box(GOI,BoxPlotName,BoxData,BoxNames,BoxColor,BoxTitle,BoxLegend);
                        },width = PlotWidth,height = PlotHeight)
                      }  
                    }else{
                      for(i in 1:J){ 
                        output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                      }
                    }
                  }else{
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              #### METASTASIS STAGE ####
              if(Analysis == "Metastasis Stage"){
                myclinical = myclinicaldata
                cM0 = grep(pattern = "cM0",x = myclinical$AJCC_METASTASIS_PATHOLOGIC_PM);
                M0 = grep(pattern = "M0",x = myclinical$AJCC_METASTASIS_PATHOLOGIC_PM);
                M0 = M0[!M0 %in% cM0];
                M1 = grep("M1",myclinical$AJCC_METASTASIS_PATHOLOGIC_PM);
                M0 = row.names( myclinical[M0,]);
                if(length(M0)!= 0)M0 = Search_Clin_Exp(M0,Expression)
                M1 = row.names( myclinical[M1,]);
                if(length(M1)!= 0)M1 = Search_Clin_Exp(M1,Expression)
                if(J>0){
                  if(match(GOI,colnames(Tumor),nomatch = 0) != 0){
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot({
                        BoxPlotName = paste(GOI,Type,"Boxplot Metastasis.tiff");
                        BoxData = c(M0[GOI],M1[GOI]);
                        BoxTitle = paste(Type,"\nMetastasis Boxplot")
                        BoxLegend = Calc_P_val_TwoSample(Data1 = M0[GOI],Data2 = M1[GOI],"No distant metastasis Vs Distant metastasis");
                        BoxNames = c(paste("No distant metastasis\nSample-size: ",length(M0[,GOI]),sep = ""),paste("Distant metastasis\nSample-size: ",length(M1[,GOI]),sep = ""))
                        BoxColor = c("green","red");
                        Generate_Box(GOI,BoxPlotName,BoxData,BoxNames,BoxColor,BoxTitle,BoxLegend);
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              #### NEOADJUVENT ####
              if(Analysis == "Neoadjuvent"){
                myclinical = myclinicaldata
                if(length(grep("HISTORY_NEOADJUVANT",x = colnames(myclinical),ignore.case = TRUE))>0){
                  Yes = row.names(myclinical[grep("Yes",ignore.case = TRUE, myclinical$HISTORY_NEOADJUVANT_TRTYN),]);
                  No = row.names(myclinical[grep("No",ignore.case = TRUE, myclinical$HISTORY_NEOADJUVANT_TRTYN),]);
                  if(length(Yes)!= 0)Yes = Search_Clin_Exp(Yes,Expression)
                  if(length(No)!= 0)No = Search_Clin_Exp(No,Expression)
                  if(match(GOI,colnames(Tumor),nomatch = 0) != 0){
                    if(length(Yes[GOI]) >0 & length(No[GOI])>0){
                      for(i in 1:J){ 
                        output[[pnodes[i]]] <- renderPlot({
                          BoxData = c(No[GOI],Yes[GOI]);
                          BoxLegend = Calc_P_val_TwoSample(Data1 = No[GOI],Data2 = Yes[GOI],"No Vs Yes");
                          BoxTitle = paste(Type,"\nHistory of Neo Adjuvent Treatment Boxplot")
                          BoxNames = c(paste("No\nSample-size: ",length(No[,GOI]),sep = ""),paste("Yes\nSample-size: ",length(Yes[,GOI]),sep = ""))
                          BoxColor = c("green","red");
                          BoxPlotName = paste(GOI,Type,"Boxplot Radiation treatment.tiff");
                          Generate_Box(GOI,BoxPlotName,BoxData,BoxNames,BoxColor,BoxTitle,BoxLegend);
                        },width = PlotWidth,height = PlotHeight)
                      }
                    }else{
                      for(i in 1:J){ 
                        output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                      }
                    }
                  }else{
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              #### NODAL STAGE ####
              if(Analysis == "Nodal Stage"){
                myclinical = myclinicaldata
                M0 = grep(pattern = "N0",x = myclinical$AJCC_NODES_PATHOLOGIC_PN);
                M1 = grep("[N1|N2|N3]",myclinical$AJCC_METASTASIS_PATHOLOGIC_PM);
                M0 = row.names( myclinical[M0,]);
                if(length(M0)!= 0)M0 = Search_Clin_Exp(M0,Expression)
                M1 = row.names( myclinical[M1,]);
                if(length(M1)!= 0)M1 = Search_Clin_Exp(M1,Expression)
                if(J > 0){
                  if(match(GOI,colnames(Tumor),nomatch = 0) != 0){
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot({
                        BoxPlotName = paste(GOI,Type,"Boxplot Nodal Stage.tiff");
                        BoxData = c(M0[GOI],M1[GOI]);
                        BoxTitle = paste(Type,"\nNodal Stage Boxplot")
                        BoxLegend = Calc_P_val_TwoSample(Data1 = M0[GOI],Data2 = M1[GOI],"No cancer found in the lymph nodes Vs Involvement of regional lymph nodes");
                        BoxNames = c(paste("no cancer found in the lymph nodes\nSample-size: ",length(M0[,GOI]),sep = ""),paste("Involvement of regional lymph nodes\nSample-size: ",length(M1[,GOI]),sep = ""))
                        BoxColor = c("green","red");
                        Generate_Box(GOI,BoxPlotName,BoxData,BoxNames,BoxColor,BoxTitle,BoxLegend);
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              #### NORMAL VS TUMOR ####
              if(Analysis == "Normal vs Tumor"){
                if(J>0){
                  if(!is.null(Normal) || length(Normal) != 0){
                    #if(match(GOI,colnames(Tumor),nomatch = 0) != 0){
                      for(i in 1:J){ 
                        output[[pnodes[i]]] <- renderPlot({
                          BoxPlotName = paste(GOI,Type,"Boxplot NormaVSTumor.tiff");
                          BoxLegend = Calc_P_val_TwoSample(Data1 = Normal[GOI],Data2 = Tumor[GOI],String = "Normal Vs Tumor");
                          BoxData = c(Normal[GOI],Tumor[GOI]);
                          BoxNames = c(paste("Normal\nSample size: ",length(Normal[,GOI]),sep = ""),
                                       paste("Tumor\nSample size:",length(Tumor[,GOI]),sep = ""));
                          BoxColor = c("green","red");
                          BoxTitle = paste(Type,"\nNormal VS Tumor Boxplot");
                          Generate_Box(GOI,BoxPlotName,BoxData,BoxNames,BoxColor,BoxTitle,BoxLegend);
                        },width = PlotWidth,height = PlotHeight)
                      }
                    #}else{
                    #  for(i in 1:J){ 
                    #    output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    #  }
                    #}
                  }else{
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(Type))
                    }
                  }
                }
              }
              #### PHARMACEUTICAL THERAPY ####
              if(Analysis == "Pharmaceutical Therapy"){
                myclinical = myclinicaldata
                if(length(grep("Pharmaceutical_tx",x = colnames(myclinical),ignore.case = TRUE))>0){
                  Yes = row.names(myclinical[grep("Yes",ignore.case = TRUE, myclinical$PHARMACEUTICAL_TX_ADJUVANT),]);
                  No = row.names(myclinical[grep("No",ignore.case = TRUE, myclinical$PHARMACEUTICAL_TX_ADJUVANT),]);
                  if(length(Yes)!= 0)Yes = Search_Clin_Exp(Yes,Expression)
                  if(length(No)!= 0)No = Search_Clin_Exp(No,Expression)
                  if(J>0){
                    if(match(GOI,colnames(Tumor),nomatch = 0) != 0){
                      if(length(Yes[GOI]) > 0 & length(No[GOI])>0){
                        for(i in 1:J){ 
                          output[[pnodes[i]]] <- renderPlot({
                            BoxData = c(No[GOI],Yes[GOI]);
                            BoxLegend = Calc_P_val_TwoSample(Data1 = No[GOI],Data2 = Yes[GOI],"No Vs Yes")
                            BoxTitle = paste(Type,"\nPharmaceutical Therapy Boxplot")
                            BoxNames = c(paste("No\nSample-size: ",length(No[,GOI]),sep = ""),paste("Yes\nSample-size: ",length(Yes[,GOI]),sep = ""))
                            BoxColor = c("green","red");
                            BoxPlotName = paste(GOI,Type,"Boxplot Pharmaceutical.tiff");
                            Generate_Box(GOI,BoxPlotName,BoxData,BoxNames,BoxColor,BoxTitle,BoxLegend);
                          },width = PlotWidth,height = PlotHeight)
                        }
                      }else{
                        for(i in 1:J){ 
                          output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                        }
                      }
                    }else{
                      for(i in 1:J){ 
                        output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                      }
                    }
                  }
                }else{
                  for(i in 1:J){ 
                    output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                  }
                }
              }
              #### RACE ####
              if(Analysis == "Race"){
                myclinical = myclinicaldata
                if(length(grep("RACE",colnames(myclinical),ignore.case = T))>0){
                  RACES = as.character(sort(unique(myclinical$RACE)));
                  for(Pat in RACES){
                    if(Pat != ""){
                      RInd = row.names(myclinical[grep(pattern = Pat,x = myclinical$RACE,ignore.case = TRUE),]);
                      if(length(RInd)!= 0)assign(paste(Pat,sep = ""),Search_Clin_Exp(RInd,Expression))
                      rm(RInd);
                    }
                  }
                  if(J>0){
                    if(match(GOI,colnames(Tumor),nomatch = 0) != 0){
                      for(i in 1:J){ 
                        output[[pnodes[i]]] <- renderPlot({
                          BoxData = c(Normal[GOI]);
                          BoxNames = c(paste("Normal\n(",length(Normal[,GOI]),")",sep = ""));
                          BoxLegend = c();
                          for(Data in RACES){
                            if(Data != ""){
                              if(length(grep(pattern = paste("",Data,sep = ""),x = ls()))>0){
                                print(Data)
                                BoxData = c(BoxData,get(x = Data)[GOI]);
                                BoxNames = c(BoxNames,paste(sub(" OR ","\n",Data,15),"\n(",dim(get(x = Data))[1],")",sep = "")) 
                                BoxLegend = c(BoxLegend,as.character(Calc_P_val_TwoSample(Data1 = Normal[GOI],Data2 = get(x = Data)[GOI],
                                                                                          paste("Normal Vs ",Data,sep = "") ) ) )
                              }
                            }
                          }
                          BoxPlotName = paste(GOI,Type,"Boxplot Race.tiff");
                          BoxColor = rainbow(7);
                          BoxTitle = paste(Type,"\nRace Boxplot");
                          Generate_Box(GOI = GOI,BoxPlotName = BoxPlotName,BoxData = BoxData,BoxNames = BoxNames,
                                       BoxColor =   BoxColor,BoxTitle = BoxTitle,BoxLegend = BoxLegend);
                        },width = PlotWidth,height = PlotHeight)
                      }
                    }else{
                      for(i in 1:J){ 
                        output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                      }
                    }
                  }
                }else{
                  for(i in 1:J){ 
                    output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                  }
                }
              }
              #### RADIATION TREATMENT ####
              if(Analysis == "Radiation Treatment"){
                myclinical = myclinicaldata
                if(length(grep("Radiation_treatment",x = colnames(myclinical),ignore.case = TRUE))>0){
                  Yes = row.names(myclinical[grep("Yes",ignore.case = TRUE, myclinical$RADIATION_TREATMENT_ADJUVANT),]);
                  No = row.names(myclinical[grep("No",ignore.case = TRUE, myclinical$RADIATION_TREATMENT_ADJUVANT),]);
                  if(length(Yes)!= 0)Yes = Search_Clin_Exp(Yes,Expression)
                  if(length(No)!= 0)No = Search_Clin_Exp(No,Expression)
                  if(match(GOI,colnames(Tumor),nomatch = 0) != 0){
                    if(length(Yes[GOI]) > 0 & length(No[GOI])> 0){
                      if(J>0){
                        for(i in 1:J){ 
                          output[[pnodes[i]]] <- renderPlot({
                            BoxData = c(No[GOI],Yes[GOI]);
                            BoxLegend = Calc_P_val_TwoSample(Data1 = No[GOI],Data2 = Yes[GOI],"No Vs Yes");
                            BoxTitle = paste(Type,"\nRadiation Treatment Boxplot")
                            BoxNames = c(paste("No\nSample-size: ",length(No[,GOI]),sep = ""),paste("Yes\nSample-size: ",length(Yes[,GOI]),sep = ""))
                            BoxColor = c("green","red");
                            BoxPlotName = paste(GOI,Type,"Boxplot Radiation treatment.tiff");
                            Generate_Box(GOI,BoxPlotName,BoxData,BoxNames,BoxColor,BoxTitle,BoxLegend);
                          },width = PlotWidth,height = PlotHeight )
                        }
                      }
                    }else{
                      for(i in 1:J){ 
                        output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                      }
                    }
                  }else{
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              #### SEX/GENDER ####
              if(Analysis == "Sex"){
                myclinical = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = "", Code_type = 2);
                Male = row.names( myclinical[grep(1,myclinical$SEX),]);
                if(length(Male)!= 0)Male = Search_Clin_Exp(Male,Expression)
                Female = row.names( myclinical[grep(2,myclinical$SEX),]);
                if(length(Female)!= 0)Female = Search_Clin_Exp(Female,Expression)
                if(length(Male)!= 0 & length(Female)!= 0){
                  if(J>0){
                    if(match(GOI,colnames(Tumor),nomatch = 0) != 0){
                      for(i in 1:J){ 
                        output[[pnodes[i]]] <- renderPlot({ # plot in each tab
                          BoxPlotName = paste(GOI,Type,"Boxplot GENDER.tiff");
                          BoxData = c(Male[GOI],Female[GOI]);
                          BoxTitle = paste(Type,"\nGender Boxplot")
                          BoxLegend = Calc_P_val_TwoSample(Data1 = Male[GOI],Data2 = Female[GOI],"Male Vs Female");
                          BoxNames = c(paste("Male\nSample-size: ",length(Male[,GOI]),sep = ""),paste("Female\nSample-size: ",length(Female[,GOI]),sep = ""))
                          BoxColor = c("green","red");
                          Generate_Box(GOI,BoxPlotName,BoxData,BoxNames,BoxColor,BoxTitle,BoxLegend);
                        },width = PlotWidth,height = PlotHeight)
                      }
                    }else{
                      for(i in 1:J){ 
                        output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                      }
                    }
                  }
                }
              }
              
            }else{
              for(i in 1:J){ 
                output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
              }
            }
          }else{
            if(Study == "Multiple Studies"){
              if(Analysis == "Normal Samples" || Analysis == "Tumor Samples"){
                if(J>0){
                  for(i in 1:J){
                    if(i == as.numeric(I)){
                      output[[pnodes[i]]] <- renderPlot({
                        BOXGOI = GOI;
                        #BOXGOI = BOXGOI[1]  #Single Gene Input
                        Types = input$M_Cancerstudy#M_Type();
                        BoxData = c();
                        BoxNames = c();
                        BoxType = Analysis;
                        BoxType = substr(x = BoxType,start = 1,stop = as.numeric(regexpr("Samples",BoxType))-2);
                        for(Type in Types){
                          BoxFiles = dir(paste("../PANCAN_TOOL/",Type,"/Expression Data/",sep=""));
                          if(length(grep("Normal",BoxFiles))!= 0){
                            INTIALS = unique(substr(BOXGOI,1,1))
                            FilePath = gsub(pattern = "TYPE",replacement = Type,x = paste0("../PANCAN_TOOL/TYPE/Expression Data/Expression_TYPE_",BoxType,"_GENEINT.csv"));
                            FilePath = gsub(pattern = "GENEINT",replacement = INTIALS,x = FilePath);
                            BoxN = tryCatch({
                              read.csv(file = FilePath, header = TRUE, row.names = 1);
                            },
                            error = function(cond){
                              message("Error in Reading file",cond)
                            },
                            warning = function(cond){
                              message("Warning: File reading",cond)
                            }
                            )
                            BoxN = BoxN[BOXGOI]
                            Names = gsub(pattern = " ",replacement = "\n",x = paste(substring(text = Type, regexpr(pattern = "TCGA",Type)+5)));
                            #else {BoxNormal = NA}
                            BoxData = c(BoxData,BoxN);
                            BoxNames = c(BoxNames,Names);
                          }
                          BoxPlotName = paste(BOXGOI,"Boxplot Samples.tiff");
                          BoxColor = rainbow(length(BoxData));
                          BoxTitle = paste0(BoxType," Samples Boxplot");
                          BoxLegend = "NO";
                        }
                        Generate_Box(BOXGOI,BoxPlotName,BoxData,BoxNames,BoxColor,BoxTitle,BoxLegend);
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
            }
          }
        }
        #### SURVIVAL ####
        if(Plot == "Survival"){
          myclinicaldata = myclinicaldata()
          Type = input$Cancerstudy# Type();
          #if(Type != ""){
          Selected_Study = as.character(ProCT$cBioStudy[ProCT$Cancer_Type == Type])
          GeneticProfile = GeneticProfile() 
          GG = c();
          if(!is.null(GeneticProfile)){
            for(GL in colnames(GeneticProfile)){
              if(match(x = GL,table = Gene_List,nomatch = 0) != 0){
                GG = c(GG,GL)
              }
            }
            GeneticProfile = GeneticProfile[GG];
            if(!is.null(myclinicaldata)){
              if(Analysis == "Age"){
                if(J>0){
                  if(match(x = GOI,table = colnames(GeneticProfile),nomatch = 0) != 0){
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot({ # plot in each tab
                        myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = 3,Cutoff = 0);
                        KM <- survfit(Surv(OS_MONTHS, VITAL_STATUS) ~ AGE_GROUP, type="kaplan-meier", conf.type="log", data = myclinicaldata)
                        survdiff <- survdiff(Surv(OS_MONTHS, VITAL_STATUS) ~ AGE_GROUP, data = myclinicaldata) #Tests if there is a difference between two or more survival curves
                        Legend_label = c("0-20","20-40","40-60","60-80","80-120")[sort(unique(myclinicaldata$AGE_GROUP))];
                        survival_data = myclinicaldata;
                        Legend_Title = "Age Group";
                        Generate_surv(KM,survival_data,GOI,Legend_Title,Legend_label,Time())
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              if(Analysis == "Expression"){
                if(J>0){
                  if(match(x = GOI,table = colnames(GeneticProfile),nomatch = 0) != 0){
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot({
                        myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = 1,Cutoff = as.numeric(input$Z_Score),GOI = GOI);
                        High = grep("High",ignore.case = TRUE, myclinicaldata$Expression)
                        Low = grep("Low",ignore.case = TRUE, myclinicaldata$Expression)
                        if(length(High)>0 & length(Low)>0){
                          survival_data = rbind(myclinicaldata[High,],myclinicaldata[Low,]);
                          KM <- survfit(Surv(OS_MONTHS, VITAL_STATUS) ~ Expression, type="kaplan-meier", conf.type="log", data = survival_data)
                          survdiff <- survdiff(Surv(OS_MONTHS, VITAL_STATUS) ~ Expression, data = survival_data) #Tests if there is a difference between two or more survival curves
                          Legend_label = c(paste('High',GOI),paste('Low',GOI));
                          Legend_Title = "Expression";
                          Generate_surv(KM,survival_data,GOI,Legend_Title,Legend_label,Time());
                        }else{DataNA(GOI)}
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              if(Analysis == "Pharmaceutical Therapy"){
                if(J>0){
                  if(match(x = GOI,table = colnames(GeneticProfile),nomatch = 0) != 0){
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot({ # plot in each tab
                        if(length(grep("Pharmaceutical_tx",x = colnames(myclinicaldata),ignore.case = TRUE))>0){
                          myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = "",Cutoff = 0);
                          Yes = grep("Yes",ignore.case = TRUE, myclinicaldata$PHARMACEUTICAL_TX_ADJUVANT);
                          No = grep("No",ignore.case = TRUE, myclinicaldata$PHARMACEUTICAL_TX_ADJUVANT);
                          survival_data = rbind(myclinicaldata[Yes,],myclinicaldata[No,]);
                          KM <- survfit(Surv(OS_MONTHS, VITAL_STATUS) ~ PHARMACEUTICAL_TX_ADJUVANT, type="kaplan-meier", conf.type="log", data = survival_data)
                          survdiff <- survdiff(Surv(OS_MONTHS, VITAL_STATUS) ~ PHARMACEUTICAL_TX_ADJUVANT, data = survival_data) #Tests if there is a difference between two or more survival curves
                          Legend_label = sub("PHARMACEUTICAL_TX_ADJUVANT=","",names(KM$strata));
                          Legend_Title = "Pharmaceutical Therapy";
                          if(length(Yes)>0 & length(No)>0){
                            Generate_surv(KM,survival_data,GOI,Legend_Title,Legend_label,Time());
                          }
                        }
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              if(Analysis == "Race"){
                if(J>0){
                  if(match(x = GOI,table = colnames(GeneticProfile),nomatch = 0) != 0){
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot({ # plot in each tab
                        myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,"",0);
                        RACES = as.character(unique(myclinicaldata$RACE));
                        Index_RACE = c();
                        for(Pat in RACES){
                          if(Pat != ""){
                            Index_RACE = c(Index_RACE,grep(pattern = Pat,x = myclinicaldata$RACE));
                          }
                        }
                        survival_data = myclinicaldata[Index_RACE,];
                        KM <- survfit(Surv(OS_MONTHS, VITAL_STATUS) ~ RACE, type="kaplan-meier", conf.type="log", data = survival_data)
                        survdiff <- survdiff(Surv(OS_MONTHS, VITAL_STATUS) ~ RACE, data = survival_data) 
                        Legend_label = sub("RACE=","",names(KM$strata));
                        Legend_Title = "Race";
                        Generate_surv(KM,survival_data,GOI,Legend_Title,Legend_label,Time())
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              if(Analysis == "Radiation Treatment"){
                if(J>0){
                  if(match(x = GOI,table = colnames(GeneticProfile),nomatch = 0) != 0){
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot({ # plot in each tab
                        if(length(grep("Radiation_treatment",x = colnames(myclinicaldata),ignore.case = TRUE))>0){
                          myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = "",Cutoff = 0);
                          Yes = grep("Yes",ignore.case = TRUE, myclinicaldata$RADIATION_TREATMENT_ADJUVANT)
                          No = grep("No",ignore.case = TRUE, myclinicaldata$RADIATION_TREATMENT_ADJUVANT)
                          survival_data = rbind(myclinicaldata[Yes,],myclinicaldata[No,]);
                          KM <- survfit(Surv(OS_MONTHS, VITAL_STATUS) ~ RADIATION_TREATMENT_ADJUVANT, type="kaplan-meier", conf.type="log", data = survival_data)
                          survdiff <- survdiff(Surv(OS_MONTHS, VITAL_STATUS) ~ RADIATION_TREATMENT_ADJUVANT, data = survival_data) #Tests if there is a difference between two or more survival curves
                          Legend_label = sub("RADIATION_TREATMENT_ADJUVANT=","",names(KM$strata));
                          Legend_Title = "Radiation Treatment";
                          if(length(Yes)>0 & length(No)>0){
                            Generate_surv(KM,survival_data,GOI,Legend_Title,Legend_label,Time());
                          }
                        }
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              if(Analysis == "Sex"){
                if(J>0){
                  if(match(x = GOI,table = colnames(GeneticProfile),nomatch = 0) != 0){
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot({ # plot in each tab
                        myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = 2,Cutoff = 0);
                        KM <- survfit(Surv(OS_MONTHS, VITAL_STATUS) ~ SEX, type="kaplan-meier", conf.type="log", data = myclinicaldata)
                        survdiff <- survdiff(Surv(OS_MONTHS, VITAL_STATUS) ~ SEX, data = myclinicaldata) #Tests if there is a difference between two or more survival curves
                        Legend_label = c("Male","Female");
                        Legend_Title = "Gender";
                        survival_data = myclinicaldata;
                        Generate_surv(KM,survival_data,GOI,Legend_Title,Legend_label,Time())
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:J){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
            }else{
              for(i in 1:J){ 
                output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
              }
            }
          }else{
            for(i in 1:J){ 
              output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
            }
          }
        }
      }
    }
  })
  output$OutPlot = renderUI({
    tabs <- list(NULL)
    ## temporary firsttab (disappears after data selection) :
    tabs[[1]] <- tabPanel("Data",value="0")
    ## permanent tabs : 1, 2, ..., J
    pobjects <- pObjects()
    if (!is.null(pobjects)) { 
      outnodes <- outputNodes()
      pnodes <- outnodes$pnodes
      tabnames <- pobjects$Gene_List
      J <- length(tabnames)
      if(J>0){
        if(!is.null(tabnames)){
          for(i in 1:J){
            tabs[[i]] <- tabPanel(tabnames[i],plotOutput(pnodes[i]),value=i)
          }
        } 
        tabs$id <- "tab1"
        do.call(tabsetPanel, tabs)
      }
    }
  })
  #### BOXPLOT ####
  output$BoxplotALL <- renderPlot({
    BOXGOI = GOI();
    BOXGOI = BOXGOI[1]  #Single Gene Input
    Types = M_Type();
    BoxData = c();
    BoxNames = c();
    BoxType = input$M_BoxType;
    BoxType = substr(x = BoxType,start = 1,stop = as.numeric(regexpr("Samples",BoxType))-2);
    if(length(BOXGOI == 1)){
      for(Type in Types){
        BoxFiles = dir(paste("../PANCAN_TOOL/",Type,"/Expression Data",sep=""));
        if(length(grep(BoxType,BoxFiles))!= 0){
          INTIALS = unique(substr(BOXGOI,1,1))
          if(BoxType == "Normal")BoxN = GetNormal(GOI = BOXGOI,Type = Type)
          if(BoxType == "Tumor")BoxN = GetTumor(GOI = BOXGOI,Type = Type)
          BoxN = BoxN[BOXGOI]
          Names = gsub(pattern = " ",replacement = "\n",x = paste(substring(text = Type, regexpr(pattern = "TCGA",Type)+5)));
          #else {BoxNormal = NA}
          BoxData = c(BoxData,BoxN);
          BoxNames = c(BoxNames,Names);
        }
        BoxPlotName = paste(BOXGOI,"Boxplot Samples.tiff");
        BoxColor = rainbow(length(BoxData));
        BoxTitle = paste0(BoxType," Samples Boxplot");
        BoxLegend = "NO";
      }
      Generate_Box(BOXGOI,BoxPlotName,BoxData,BoxNames,BoxColor,BoxTitle,BoxLegend);
    }
  })
  #### HEATMAPS ####
  output$HeatText <- renderText("Select at least 2 gene of interest.")
  output$HeatNT <- renderPlot({
    Type = Type();
    GOI = GOI()
    if(length(GOI)>1){
      #INTIALS = unique(substr(GOI,1,1))
      #Normal = data.frame()
      #Tumor = data.frame()
      #for(i in 1:length(INTIALS)){
      #  FilePath = gsub(pattern = "TYPE",replacement = Type,x = "../PANCAN_TOOL/TYPE/Expression Data/Expression_TYPE_CT_GENEINT.csv");
      #  FilePath = gsub(pattern = "GENEINT",replacement = INTIALS[i],x = FilePath);
        #Nor = read.csv(file = sub("CT","Normal",FilePath),header = TRUE,row.names = 1)
        #Tum = read.csv(file = sub("CT","Tumor",FilePath),header = TRUE,row.names = 1)
        #GOGG = GOI[grep(paste0("^",INTIALS[i]),GOI)]
        #if(i == 1){
          #Normal = Nor[GOGG]
        #  Tumor = Tum[GOGG]
        #}else{
          #Normal = cbind(Normal,Nor[GOGG])
        #  Tumor = cbind(Tumor,Tum[GOGG])
        #}
      #}
#      Normal = `row.names<-`(Normal,row.names(Nor[1]))
      Tumor = Tumor()#`row.names<-`(Tumor,row.names(Tum[1]))
      Normal = Normal()
      Heat = data.frame(t(Normal[,GOI]),t(Tumor[,GOI]));
      Column_Color= c(rep("gray",length(Normal[,1])),rep("black",length(Tumor[,1])));
      HeatmapName = sub("HEAT",Type,"HEAT-heatmap.tiff");
      Title = paste(substring(text = Type, regexpr(pattern = "TCGA",HeatmapName)+5), "Correlation");
      heatmap.2(as.matrix(Heat),
                main = Title,         # heat map title
                notecol="black",      # change font color of cell labels to black
                density.info="none",  # turns off density plot inside color legend
                trace="none",         # turns off trace lines inside the heat map
                margins =c(8,10),     # widens margins around plot
                col=my_palette,       # use on color palette defined earlier
                dendrogram = "none",
                Colv = "NA",
                Rowv = "NA",
                keysize = 1,
                #breaks = col_breaks,
                ColSideColors = Column_Color,
                cexRow = 1.5,
                lwid = c(0.1,0.75,4),
                lhei = c(0.25,0.5,0.1,1),
                lmat=rbind(c(4,4,4),c(3,5,0),c(0,1,1),c(0,2,2))          
      )# +
      legend("topright",                                 # location of the legend on the heatmap plot
             legend = c("Normal", "Tumor"), # category labels
             col = c("gray", "black"),                  # color key
             lty= 1,                                            # line style
             lwd = 8)                                          # line width
    }
  })
  output$R_Heat_T <- renderPlot({
    Type = Type();
    GOI = GOI()
    HEAT_TYPE = input$HeatType
    Selected_Study = as.character(ProCT$cBioStudy[ProCT$Cancer_Type == Type])
    myclinicaldata = myclinicaldata()#GetDataCBP(1,Selected_Study,Gene_List = GOI)
    Selected_Study = as.character(ProCT$cBioStudy[ProCT$Cancer_Type == Type])
    FilePath = gsub(pattern = "TYPE",replacement = Type,x = "../PANCAN_TOOL/TYPE/Expression Data/cBio-Expression_FILE_ZScore.csv");
    GeneticProfile = GeneticProfile()#tryCatch({GetDataCBP(DataCode = 2,Selected_Study,Gene_List = GOI)})
    if(HEAT_TYPE == "Pharmaceutical Therapy"){
      if(length(grep("Pharmaceutical",x = colnames(myclinicaldata),ignore.case = TRUE))>0){
        myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = "",Cutoff = 0);
        Yes = grep("Yes",ignore.case = TRUE, myclinicaldata$PHARMACEUTICAL_TX_ADJUVANT)
        No = grep("No",ignore.case = TRUE, myclinicaldata$PHARMACEUTICAL_TX_ADJUVANT)
        Yes = match(row.names(myclinicaldata[Yes,]),row.names(GeneticProfile))
        No = match(row.names(myclinicaldata[No,]),row.names(GeneticProfile))
        Title = paste(Type,"-Pharmaceutical Therapy",sep = "");
        Max_val = abs(as.numeric(input$ZScore));
        Min_val = (-1*Max_val);
        my_palette <- colorRampPalette(c("green", "black", "red"))(n = 299);
        ColFontSize = 0.7;        #Font Size of Row labels
        col_breaks = c(seq(-5,Min_val,length=100),seq((Min_val+0.04),(Max_val-0.04),length=100),seq(Max_val,5,length=100));
        Column_Color= c(rep("gray",length(No)),rep("black",length(Yes)))
        HEATDATA = as.matrix(t(rbind(GeneticProfile[No,GOI],GeneticProfile[Yes,GOI])))
      }
    }
    if(HEAT_TYPE == "Radiation Treatment"){
      if(length(grep("Radiation_treatment",x = colnames(myclinicaldata),ignore.case = TRUE))>0){
        myclinicaldata4 = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = "",Cutoff = 0);
        Yes = grep("Yes",ignore.case = TRUE, myclinicaldata4$RADIATION_TREATMENT_ADJUVANT)
        No = grep("No",ignore.case = TRUE, myclinicaldata4$RADIATION_TREATMENT_ADJUVANT)
        Yes = match(row.names(myclinicaldata[Yes,]),row.names(GeneticProfile))
        No = match(row.names(myclinicaldata[No,]),row.names(GeneticProfile))
        Title = paste(Type,"-Radiation Treatment",sep = "");
        Max_val = abs(as.numeric(input$ZScore));
        Min_val = (-1*Max_val);
        my_palette <- colorRampPalette(c("green", "black", "red"))(n = 299);
        ColFontSize = 0.7;        #Font Size of Row labels
        col_breaks = c(seq(-5,Min_val,length=100),seq((Min_val+0.04),(Max_val-0.04),length=100),seq(Max_val,5,length=100));
        Column_Color= c(rep("gray",length(No)),rep("black",length(Yes)))
        HEATDATA = as.matrix(t(rbind(GeneticProfile[No,GOI],GeneticProfile[Yes,GOI])))
      }
    }
    if(HEAT_TYPE == "Radiation Treatment" || HEAT_TYPE == "Pharmaceutical Therapy"){
      heatmap.2(x = HEATDATA,na.color = "white",
                main = Title,         # heat map title
                notecol = "black",      # change font color of cell labels to black
                density.info = "none",  # turns off density plot inside color legend
                trace = "none",         # turns off trace lines inside the heat map
                margins = c(8,10),     # widens margins around plot
                col = my_palette,       # use on color palette defined earlier
                dendrogram = "none",
                Colv = "NA",
                Rowv = "NA",
                keysize = 1,
                breaks = col_breaks,
                ColSideColors = Column_Color,
                key.xlab = "Z-score",
                cexCol = ColFontSize,
                cexRow = 1.5,
                na.rm = 0,
                cex.main = 4,
                cex = 1,
                lwid = c(0.1,0.75,4),
                lhei = c(0.25,0.5,0.1,1),
                lmat=rbind(c(4,4,4),c(3,5,0),c(0,1,1),c(0,2,2))
      );
      legend("topright",                                 # location of the legend on the heatmap plot
             legend = c("No", "Yes"), # category labels
             col = c("gray", "black"),inset = c(0,0.1),                  # color key
             lty= 1,  #cex.09,                                          # line style
             lwd = 8);
    }
  })
####### PROTEIN   #####
  P_Type <- reactive({
    validate(need(expr = input$P_Cancerstudy != "",message = "Please Select a Study"))
    input$P_Cancerstudy
  })
  P_myclinicaldata <- reactive({
    Type = P_Type();
    Selected_Study = as.character(ProCT$cBioStudy[ProCT$Cancer_Type == Type])
    if(length(grep("Protein",dir(gsub(pattern = "TYPE",Type,"../PANCAN_TOOL/TYPE/")))) != 0){
      FilePath = gsub(pattern = "TYPE",replacement = Type,x = "../PANCAN_TOOL/TYPE/Protein/");
      P_myclinicaldata = read.csv(file = paste0(FilePath,dir(FilePath)[grep("Clinical",dir(FilePath))]),header = TRUE,row.names = 1)
    }
    return(P_myclinicaldata);
  })
  P_Time <- reactive({
    if(input$SurvivalTime == "Time of your choice"){
      myclinicaldata = myclinicaldata()
      switch (input$SurvType,
              "Age" = {
                myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = "",Code_type = 3,Cutoff = 0);
                UN = unique(myclinicaldata$AGE_GROUP)
                if(match(1,UN,0) != 0)M1 = max(myclinicaldata$OS_MONTHS[myclinicaldata$AGE_GROUP == 1],na.rm = TRUE)else M1 = NA
                if(is.na(M1))M2 = max(myclinicaldata$OS_MONTHS[myclinicaldata$AGE_GROUP == 2],na.rm = TRUE)else M2 = NA
                if(is.na(M2))M3 = max(myclinicaldata$OS_MONTHS[myclinicaldata$AGE_GROUP == 3],na.rm = TRUE)else M3 = NA
                if(is.na(M3))M4 = max(myclinicaldata$OS_MONTHS[myclinicaldata$AGE_GROUP == 4],na.rm = TRUE)else M4 = NA
                if(is.na(M4))M5 = max(myclinicaldata$OS_MONTHS[myclinicaldata$AGE_GROUP == 5],na.rm = TRUE)else M5 = NA
                MaxVal = as.integer(max(c(M1,M2,M3,M4,M5),na.rm = TRUE)) 
              },
              "Expression" = {
                #myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = 1,Cutoff = as.numeric(input$Z_Score),GOI = GOI);
                #High = grep("High",ignore.case = TRUE, myclinicaldata$Expression)
                #Low = grep("Low",ignore.case = TRUE, myclinicaldata$Expression)
                
              },
              "Pharmaceutical Therapy" = {
                myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = "",Cutoff = 0);
                Yes = grep("Yes",ignore.case = TRUE, myclinicaldata$PHARMACEUTICAL_TX_ADJUVANT);
                No = grep("No",ignore.case = TRUE, myclinicaldata$PHARMACEUTICAL_TX_ADJUVANT);
                if(length(Yes)>0&&length(No)>0){
                  MaxVal = min(c(as.integer(max(myclinicaldata$OS_MONTHS[Yes],na.rm = TRUE)),as.integer(max(myclinicaldata$OS_MONTHS[No],na.rm = TRUE))))
                }
                else MaxVal = NA
              },
              "Race" = {
                myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,"",0);
                RACES = as.character(unique(myclinicaldata$RACE));
                Index_RACE = c();
                for(Pat in RACES){
                  if(Pat != ""){
                    Index_RACE = c(Index_RACE,grep(pattern = Pat,x = myclinicaldata$RACE));
                  }
                }
                MaxVal = as.integer(min(c(myclinicaldata$OS_MONTHS[Index_RACE],1),na.rm = TRUE))
              },
              "Radiation Treatment" = ,
              "Sex" = {
                myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = 2,Cutoff = 0);
                UN = unique(myclinicaldata$SEX)
                if(match(1,UN,0) != 0)M1 = max(myclinicaldata$OS_MONTHS[myclinicaldata$SEX == 1],na.rm = TRUE)else M1 = NA
                if(match(2,UN,0) != 0)M2 = max(myclinicaldata$OS_MONTHS[myclinicaldata$SEX == 2],na.rm = TRUE)else M2 = NA
                MaxVal = as.integer(min(c(M1,M2),na.rm = TRUE)) 
              },
              {MaxVal = as.integer(max(myclinicaldata$OS_MONTHS,na.rm = TRUE))}
      )
      updateSliderInput(session, "Time", value = input$Time,
                        min = 1, max = MaxVal, step = 1)
      Time = input$Time;
    }else Time = 0
  })
  P_GOI <- reactive({
    validate(
      need(expr = input$P_Genes != "",'Please select a gene of interest')
    )
    input$P_Genes
  })
  P_GeneticProfile <- reactive({
    Type = P_Type()
    GOI = P_GOI()
    GeneticProfile = P_GetGeneticProfile(GOI,Type);
    print(head(GeneticProfile))
    return(GeneticProfile)
  })
  output$P_Heat_T <- renderPlot({
    Type = P_Type();
    HEAT_TYPE = input$P_HeatType
    GOI = P_GOI()
    Selected_Study = as.character(ProCT$cBioStudy[ProCT$Cancer_Type == Type])
    myclinicaldata = P_myclinicaldata()#GetDataCBP(3,Selected_Study,Gene_List = GOI)
    Selected_Study = as.character(ProCT$cBioStudy[ProCT$Cancer_Type == Type])
    #print(colnames(myclinicaldata))
    #print(grep("Pharmaceutical",x = colnames(myclinicaldata),ignore.case = TRUE))
    #print(grep("Radiation",x = colnames(myclinicaldata),ignore.case = TRUE))
    GeneticProfile = P_GeneticProfile()
    if(!is.null(myclinicaldata) && !is.null(GeneticProfile) && match(x = GOI,table = colnames(GeneticProfile),nomatch = 0) != 0){
      HEATDATA = 0;
      if(HEAT_TYPE == "Pharmaceutical Therapy"){
        if(length(grep("Pharmaceutical",x = colnames(myclinicaldata),ignore.case = TRUE))>0){
          myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = "",Cutoff = 0);
          Yes = grep("Yes",ignore.case = TRUE, myclinicaldata$PHARMACEUTICAL_TX_ADJUVANT)
          No = grep("No",ignore.case = TRUE, myclinicaldata$PHARMACEUTICAL_TX_ADJUVANT)
          Yes = match(row.names(myclinicaldata[Yes,]),row.names(GeneticProfile))
          No = match(row.names(myclinicaldata[No,]),row.names(GeneticProfile))
          Title = paste(Type,"-Pharmaceutical Therapy",sep = "");
          Column_Color= c(rep("gray",length(No)),rep("black",length(Yes)))
          HEATDATA = as.matrix(t(rbind(GeneticProfile[No,],GeneticProfile[Yes,])))
        }
      }
      if(HEAT_TYPE == "Radiation Treatment"){
        if(length(grep("Radiation_treatment",x = colnames(myclinicaldata),ignore.case = TRUE))>0){
          myclinicaldata4 = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = "",Cutoff = 0);
          Yes = grep("Yes",ignore.case = TRUE, myclinicaldata4$RADIATION_TREATMENT_ADJUVANT)
          No = grep("No",ignore.case = TRUE, myclinicaldata4$RADIATION_TREATMENT_ADJUVANT)
          Yes = match(row.names(myclinicaldata[Yes,]),row.names(GeneticProfile))
          No = match(row.names(myclinicaldata[No,]),row.names(GeneticProfile))
          Title = paste(Type,"-Radiation Treatment",sep = "");
          Column_Color= c(rep("gray",length(No)),rep("black",length(Yes)))
          HEATDATA = as.matrix(t(rbind(GeneticProfile[No,],GeneticProfile[Yes,])))
        }
      }
      if(HEAT_TYPE == "All"){
        Title = Type
        Column_Color= rep("white",length(GeneticProfile[,1]))
        HEATDATA = as.matrix(t(GeneticProfile))
        #print(length(HEATDATA))
      }
      Max_val = abs(as.numeric(input$P_ZScore));
      Min_val = (-1*Max_val);
      my_palette <- colorRampPalette(c("green", "black", "red"))(n = 299);
      ColFontSize = 0.7;        #Font Size of Row labels
      col_breaks = c(seq(-5,Min_val,length=100),seq((Min_val+0.04),(Max_val-0.04),length=100),seq(Max_val,5,length=100));
      if(is.matrix(HEATDATA)){
        heatmap.2(x = HEATDATA,na.color = "white",
                  main = Title,         # heat map title
                  notecol = "black",      # change font color of cell labels to black
                  density.info = "none",  # turns off density plot inside color legend
                  trace = "none",         # turns off trace lines inside the heat map
                  margins = c(8,10),     # widens margins around plot
                  col = my_palette,       # use on color palette defined earlier
                  dendrogram = "none",
                  Colv = "NA",
                  Rowv = "NA",
                  keysize = 1,
                  breaks = col_breaks,
                  ColSideColors = Column_Color,
                  key.xlab = "Z-score",
                  cexCol = ColFontSize,
                  cexRow = 1.5,
                  na.rm = FALSE,
                  cex.main = 4,
                  cex = 1,
                  lwid = c(0.1,0.75,4),
                  lhei = c(0.25,0.5,0.1,1),
                  lmat=rbind(c(4,4,4),c(3,5,0),c(0,1,1),c(0,2,2))
        );
        legend("topright",                                 # location of the legend on the heatmap plot
               legend = c("No", "Yes"), # category labels
               col = c("gray", "black"),inset = c(0,0.1),                  # color key
               lty= 1,  #cex.09,                                          # line style
               lwd = 8);
      }else{DataNA(Type)}
    }
  })
  Prot_Objects <- reactive({
    J <- length(P_GOI())
    Tabnames <- paste0(P_GOI()) 
    Plot_type <- input$P_Output_opt
    Study <- input$P_StudyType
    Analysis = "";
    Type = input$P_Cancerstudy
    #if(Plot_type == "Boxplot")Analysis <- input$BoxType
    if(Plot_type == "Survival")Analysis <- input$P_SurvType
    list(Gene_Length=J, Gene_List = Tabnames, Plot = Plot_type, Analysis = Analysis,Type = Type, Study = Study)
  })
  Prot_outputNodes <- reactive({ # output node names
    pobjects <- Prot_Objects()
    if (is.null(pobjects)) return(NULL)  
    J <- pobjects$Gene_Length
    list(pnodes=paste0("pnode", LETTERS[1:J])) # plot outputs
  })
  observe({ 
    pobjects <- Prot_Objects()
    if (!is.null(pobjects)) {
      outnodes <- Prot_outputNodes()
      pnodes <- outnodes$pnodes
      Genes <- pobjects$Gene_List
      Plot <- pobjects$Plot
      Study <- pobjects$Study
      Gene_Length <- pobjects$Gene_Length
      Analysis <- pobjects$Analysis
      ## tab 1, 2, ..., J
      I <- input$tabP1
      GOI <- Genes[as.numeric(I)]
      if(!is.null(Plot) && !is.null(GOI)){
        myclinicaldata = P_myclinicaldata()#GetDataCBP(DataCode = 3,Selected_Study = Selected_Study,GOI)#myclinicaldata()
        GeneticProfile = P_GeneticProfile()
        #### SURVIVAL ####
        if(Plot == "Survival"){
          print(colnames(GeneticProfile))
          message(GOI);
          Type = input$Cancerstudy# Type();
          #if(Type != ""){
          Selected_Study = as.character(ProCT$cBioStudy[ProCT$Cancer_Type == Type])
          GG = c();
          if(!is.null(GeneticProfile)){
            if(!is.null(myclinicaldata)){
              if(Analysis == "Age"){
                if(Gene_Length>0){
                  if(match(x = GOI,table = colnames(GeneticProfile),nomatch = 0) != 0){
                    for(i in 1:Gene_Length){ 
                      output[[pnodes[i]]] <- renderPlot({ # plot in each tab
                        myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = 3,Cutoff = 0);
                        KM <- survfit(Surv(OS_MONTHS, VITAL_STATUS) ~ AGE_GROUP, type="kaplan-meier", conf.type="log", data = myclinicaldata)
                        survdiff <- survdiff(Surv(OS_MONTHS, VITAL_STATUS) ~ AGE_GROUP, data = myclinicaldata) #Tests if there is a difference between two or more survival curves
                        Legend_label = c("0-20","20-40","40-60","60-80","80-120")[sort(unique(myclinicaldata$AGE_GROUP))];
                        survival_data = myclinicaldata;
                        Legend_Title = "Age Group";
                        Generate_surv(KM,survival_data,GOI,Legend_Title,Legend_label,Time())
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:Gene_Length){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              if(Analysis == "Expression"){
                if(Gene_Length>0){
                  if(match(x = GOI,table = colnames(GeneticProfile),nomatch = 0) != 0){
                    for(i in 1:Gene_Length){ 
                      output[[pnodes[i]]] <- renderPlot({
                        myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = 1,Cutoff = as.numeric(input$P_Z_Score),GOI = GOI);
                        High = grep("High",ignore.case = TRUE, myclinicaldata$Expression)
                        Low = grep("Low",ignore.case = TRUE, myclinicaldata$Expression)
                        if(length(High)>0 & length(Low)>0){
                          survival_data = rbind(myclinicaldata[High,],myclinicaldata[Low,]);
                          KM <- survfit(Surv(OS_MONTHS, VITAL_STATUS) ~ Expression, type="kaplan-meier", conf.type="log", data = survival_data)
                          survdiff <- survdiff(Surv(OS_MONTHS, VITAL_STATUS) ~ Expression, data = survival_data) #Tests if there is a difference between two or more survival curves
                          Legend_label = c(paste('High',GOI),paste('Low',GOI));
                          Legend_Title = "Expression";
                          Generate_surv(KM,survival_data,GOI,Legend_Title,Legend_label,Time());
                        }else{
                          DataNA(GOI)
                        }
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:Gene_Length){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              if(Analysis == "Pharmaceutical Therapy"){
                if(Gene_Length>0){
                  if(match(x = GOI,table = colnames(GeneticProfile),nomatch = 0) != 0){
                    for(i in 1:Gene_Length){ 
                      output[[pnodes[i]]] <- renderPlot({ # plot in each tab
                        if(length(grep("Pharmaceutical_tx",x = colnames(myclinicaldata),ignore.case = TRUE))>0){
                          myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = "",Cutoff = 0);
                          Yes = grep("Yes",ignore.case = TRUE, myclinicaldata$PHARMACEUTICAL_TX_ADJUVANT);
                          No = grep("No",ignore.case = TRUE, myclinicaldata$PHARMACEUTICAL_TX_ADJUVANT);
                          if(length(Yes)>0 & length(No)>0){
                            survival_data = rbind(myclinicaldata[Yes,],myclinicaldata[No,]);
                            KM <- survfit(Surv(OS_MONTHS, VITAL_STATUS) ~ PHARMACEUTICAL_TX_ADJUVANT, type="kaplan-meier", conf.type="log", data = survival_data)
                            survdiff <- survdiff(Surv(OS_MONTHS, VITAL_STATUS) ~ PHARMACEUTICAL_TX_ADJUVANT, data = survival_data) #Tests if there is a difference between two or more survival curves
                            Legend_label = sub("PHARMACEUTICAL_TX_ADJUVANT=","",names(KM$strata));
                            Legend_Title = "Pharmaceutical Therapy";
                            Generate_surv(KM,survival_data,GOI,Legend_Title,Legend_label,Time());
                          }else{DataNA(GOI)}
                        }
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:Gene_Length){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              if(Analysis == "Race"){
                if(Gene_Length>0){
                  if(match(x = GOI,table = colnames(GeneticProfile),nomatch = 0) != 0){
                    for(i in 1:Gene_Length){ 
                      output[[pnodes[i]]] <- renderPlot({ # plot in each tab
                        myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,"",0);
                        RACES = as.character(unique(myclinicaldata$RACE));
                        Index_RACE = c();
                        for(Pat in RACES){
                          if(Pat != ""){
                            Index_RACE = c(Index_RACE,grep(pattern = Pat,x = myclinicaldata$RACE));
                          }
                        }
                        survival_data = myclinicaldata[Index_RACE,];
                        KM <- survfit(Surv(OS_MONTHS, VITAL_STATUS) ~ RACE, type="kaplan-meier", conf.type="log", data = survival_data)
                        survdiff <- survdiff(Surv(OS_MONTHS, VITAL_STATUS) ~ RACE, data = survival_data) 
                        Legend_label = sub("RACE=","",names(KM$strata));
                        Legend_Title = "Race";
                        Generate_surv(KM,survival_data,GOI,Legend_Title,Legend_label,Time())
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:Gene_Length){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              if(Analysis == "Radiation Treatment"){
                if(Gene_Length>0){
                  if(match(x = GOI,table = colnames(GeneticProfile),nomatch = 0) != 0){
                    for(i in 1:Gene_Length){ 
                      output[[pnodes[i]]] <- renderPlot({ # plot in each tab
                        if(length(grep("Radiation_treatment",x = colnames(myclinicaldata),ignore.case = TRUE))>0){
                          myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = "",Cutoff = 0);
                          Yes = grep("Yes",ignore.case = TRUE, myclinicaldata$RADIATION_TREATMENT_ADJUVANT)
                          No = grep("No",ignore.case = TRUE, myclinicaldata$RADIATION_TREATMENT_ADJUVANT)
                          survival_data = rbind(myclinicaldata[Yes,],myclinicaldata[No,]);
                          KM <- survfit(Surv(OS_MONTHS, VITAL_STATUS) ~ RADIATION_TREATMENT_ADJUVANT, type="kaplan-meier", conf.type="log", data = survival_data)
                          survdiff <- survdiff(Surv(OS_MONTHS, VITAL_STATUS) ~ RADIATION_TREATMENT_ADJUVANT, data = survival_data) #Tests if there is a difference between two or more survival curves
                          Legend_label = sub("RADIATION_TREATMENT_ADJUVANT=","",names(KM$strata));
                          Legend_Title = "Radiation Treatment";
                          if(length(Yes)>0 & length(No)>0){
                            Generate_surv(KM,survival_data,GOI,Legend_Title,Legend_label,Time());
                          }
                        }
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:Gene_Length){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
              if(Analysis == "Sex"){
                if(Gene_Length>0){
                  if(match(x = GOI,table = colnames(GeneticProfile),nomatch = 0) != 0){
                    for(i in 1:Gene_Length){ 
                      output[[pnodes[i]]] <- renderPlot({ # plot in each tab
                        myclinicaldata = Code_Clinical(Clinical_File = myclinicaldata, Exp_File = GeneticProfile,Code_type = 2,Cutoff = 0);
                        KM <- survfit(Surv(OS_MONTHS, VITAL_STATUS) ~ SEX, type="kaplan-meier", conf.type="log", data = myclinicaldata)
                        survdiff <- survdiff(Surv(OS_MONTHS, VITAL_STATUS) ~ SEX, data = myclinicaldata) #Tests if there is a difference between two or more survival curves
                        Legend_label = c("Male","Female");
                        Legend_Title = "Gender";
                        survival_data = myclinicaldata;
                        Generate_surv(KM,survival_data,GOI,Legend_Title,Legend_label,Time())
                      },width = PlotWidth,height = PlotHeight)
                    }
                  }else{
                    for(i in 1:Gene_Length){ 
                      output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
                    }
                  }
                }
              }
            }else{
              for(i in 1:Gene_Length){ 
                output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
              }
            }
          }else{
            for(i in 1:Gene_Length){ 
              output[[pnodes[i]]] <- renderPlot(DataNA(GOI),width = PlotWidth,height = PlotHeight)
            }
          }
        }
      }
    }
  })
  output$P_OutPlot <- renderUI({
    tabs <- list(NULL)
    ## temporary firsttab (disappears after data selection) :
    tabs[[1]] <- tabPanel("Data",value="0")
    ## permanent tabs : 1, 2, ..., J
    pobjects <- Prot_Objects()
    if (!is.null(pobjects)) { 
      outnodes <- Prot_outputNodes()
      pnodes <- outnodes$pnodes
      tabnames <- pobjects$Gene_List
      Gene_Length <- length(tabnames)
      if(Gene_Length>0){
        if(!is.null(tabnames)){
          for(i in 1:Gene_Length){
            tabs[[i]] <- tabPanel(tabnames[i],plotOutput(pnodes[i]),value=i)
          }
        } 
        tabs$id <- "tabP1"
        do.call(tabsetPanel, tabs)
      }
    }
    
  })
  })