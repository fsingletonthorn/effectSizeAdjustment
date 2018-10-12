library(plyr)
library(tidyverse)
library(reshape)


         studies = 16
         analysis.type = 1
         Nmin.raw  = 30
         Nmin.cond = 15
         subset    = c("all","WEIRD","NON-WEIRD")[1]
         rootdir   =normalizePath(paste0(getwd(),"/Data/EmbargoFolder"))
         indir     = list(RAW.DATA = "RAW.DATA.PRIVATE",MASTERKEY = "MASTERKEY", SOURCEINFO = "SOURCEINFO")
         outdir    = list(ROBJECTS = "ROBJECTS",RESULTS.RDS = "RESULTS.RDS")
         onlineTables = TRUE
  
  
  paths <- unique(c(paste0(rootdir,"/",indir), paste0(rootdir,"/",outdir)))
  for(d in paths[!dir.exists(paths)]){
    dir.create(d, showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  tp  <- analysis.type
  wop <- options(warn=-1, expressions=10000)
  
  # Load Key Table
  
  if(indir$MASTERKEY==""|onlineTables){
    ML2.key <- get.GoogleSheet(data='ML2masteRkey')$df
    disp(paste("Downloaded keytable Googlesheet: ML2_masteRkey [https://docs.google.com/spreadsheets/d/1fqK3WHwFPMIjNVVvmxpMEjzUETftq_DmP5LzEhXxUHA/]"), header = "get.analyses", footer = FALSE)
  } else {
    ML2.key <- read.xlsx(file.path(rootdir,indir$MASTERKEY,"ML2_masteRkey.xlsx"),"ML2masteRkey")
    disp(paste0("Loaded keytable from disk: ML2_masteRkey.xlsx [",file.path(rootdir,indir$MASTERKEY),"]"), header = "get.analyses", footer = FALSE)
  }
  
  ML2.key <- ML2.key[!is.na(ML2.key$unique.id),]
  
  # Load data
  if(indir$RAW.DATA==""|onlineTables){
    ML2.key <- get.GoogleSheet(data='ML2masteRkey')$df
    disp(paste("Downloaded data from OSF: 'ML2_RawData_S1.rds' and 'ML2_RawData_S2.rds'"), header = FALSE, footer = FALSE)
  } 
    load(file.path(rootdir,indir$RAW.DATA,"ML2_S1.rda", fsep = "\\"))
    ML2.S1 <- ML2_S1
    load(file.path(rootdir,indir$RAW.DATA,"ML2_S2.rda", fsep = "\\"))
    ML2.S2 <- ML2_S2
    disp(paste0("Loaded data from disk: 'ML2_RawData_S1.rds' and 'ML2_RawData_S1.rds'[",file.path(rootdir,indir$RAW.DATA,"ML2_S2.rda", fsep = "\\")), header = FALSE, footer = FALSE)

  # Load information about sources
  if(indir$SOURCEINFO==""|onlineTables){
    SourceInfoTable    <- get.GoogleSheet(url = "https://docs.google.com/spreadsheets/d/1Qn_kVkVGwffBAmhAbpgrTjdxKLP1bb2chHjBMVyGl1s/pub?gid=1435507167&single=true&output=csv")$df
    disp(paste("Downloaded information about the data sources from Googlesheet: 'ML2_SourceInfo.xlsx' [https://docs.google.com/spreadsheets/d/1Qn_kVkVGwffBAmhAbpgrTjdxKLP1bb2chHjBMVyGl1s/]"), header = FALSE, footer = FALSE)
  } else {
    SourceInfoTable    <- readRDS(file.path(rootdir,indir$SOURCEINFO,"ML2_SourceInfo.xlsx"))
    disp(paste0("Loaded information about the data sources from disk: 'MML2_SourceInfo.xlsx' [",file.path(rootdir,indir$SOURCEINFO),"]"), header = FALSE, footer = FALSE)
  }
  
  # Decide which analyses to run on which groups
  toRun  <- decide.analysis(ML2.key, studies, tp)
  

  # Prepare list objects
  ML2.data           <- vector("list", length = length(ML2.key$study.analysis))
  names(ML2.data)    <- ML2.key$study.analysis
  ML2.output         <- vector("list", length = length(ML2.key$study.analysis))
  names(ML2.output)  <- ML2.key$study.analysis
  ML2.rawdata        <- vector("list", length = length(ML2.key$study.analysis))
  names(ML2.rawdata) <- ML2.key$study.analysis
  #ML2.descriptivesDump <- list()
  
  studs <- toRun$studiess
  cnt   <- 0
  
  # START STUDIES ----------------------------------
  
  for(s in studs){
    
    # Get the correct slate according to info in ML2.key['study.slate']
    if(ML2.key[s,'study.slate'] == 1){ML2.df <- ML2.S1
    } else {
      ML2.df <- ML2.S2
    }
    
    # Add a unique ID
    ML2.df$uID = seq(1, nrow(ML2.df))
    
    # Get info to create a dataset for the current study
    # keytable <- ML2.key[s,]
    ML2.in <- get.info(ML2.key[s, ], colnames(ML2.df), subset)
    
    
    # Generate chain to select variables for the data frame and create a filter chain for the variables to use for analysis
    # Info based on KeyTable information in study.vars, cases.include, site.include, params.NA
    ML2.id <- get.chain(ML2.in)
    
    # Apply the df chain to select relevant subset of variables
    ML2.df <- eval(parse(text=paste("ML2.df", ML2.id$df)))
    
    
    if(NROW(ML2.df)>0){
      
      ML2.df$study.order <- NA
      stmp <- strsplit(ML2.df$StudyOrderN,"[|]")
      
      Stud <- ML2.key$study.name[[s]]
      if(Stud%in%"Tversky"){Stud <- "Tversky.Gati"}
      if(Stud%in%"Rottenstreich"){Stud <- "Rottenstrich"}
      if(Stud%in%"Ross"&(ML2.key[s,'study.slate'] == 1)){Stud <- "Ross.Slate1"}
      if(Stud%in%"Ross"&(ML2.key[s,'study.slate'] == 2)){Stud <- "Ross.Slate2"}
      if(Stud%in%"vanLange"){Stud <- "VanLange"}
      if(Stud%in%"Giessner"){Stud <- "Geissner"}
      
      ML2.df$study.order <- laply(seq_along(stmp), function(o){which(grepl(Stud,stmp[[o]]))%00%NA})
      
      # Loop over groups within study
      ugroup       <- sort(na.exclude(unique(eval(parse(text=toRun$ugroup)))))
      tp           <- toRun$tp
      ML2.sr       <- list()
      ML2.var      <- list()
      outputSource <- list()
      dataSource   <- list()
      raw.df       <- list()
      clean.df     <- list()
      testVarEqual <- ML2.in$stat.params$var.equal
      
      
      # runGroups <- sort(na.exclude(unique(ML2.df$study.order)))
      
      cnt          <- cnt + 1
      
      ifelse(tp[cnt]==1,
             runGroups <- "all",
             runGroups <- ugroup
      )
      
      disp(paste(s, ML2.key$study.analysis[[s]],"- START"), header = toupper(ML2.key$study.analysis[[s]]), footer = FALSE)
      cat("\n")
      
      # START GROUPS ----------------------------------------------
      
      
      for(g in seq_along(runGroups)){
        
        listIT     <- FALSE
        nMin1      <- FALSE
        nMin2      <- FALSE
        compN <- compN1 <- compN2 <- 0
        
        if(tp[cnt]<4){
          if(runGroups[g]=="all"){gID <- rep(TRUE, nrow(ML2.df))} else {gID <- ML2.df$source%in%runGroups[g]}
        } else {
          gID <-  ML2.df$study.order%in%runGroups[g]
        }
        
        # Check nMin
        if(sum(gID, na.rm=TRUE) >= Nmin.raw){
          nMin1 <- TRUE
          # Get a list containing the data frames to be used in the analysis
          ML2.sr[[g]] <- get.sourceData(ML2.id, ML2.df[gID, ], ML2.in)
        }
        
        # Double-check nMin
        if(nMin1){
          compN  <- ML2.sr[[g]]$N
          compN1 <- sum(ML2.sr[[g]]$RawDataFilter[[1]]$Included, na.rm = TRUE)
          compN2 <- sum(ML2.sr[[g]]$RawDataFilter[[2]]$Included, na.rm = TRUE)
          if(any(compN >= Nmin.raw)&(all(compN1>=Nmin.cond, compN2>=Nmin.cond))){nMin2 <- TRUE}
        }
        
        # START ANALYSIS ----------------------------------------
        
        if(all(nMin1,nMin2)){
          
          # Organize and Calculate variables for the analysis using function according to ML2.info: 'stat.vars'
          ML2.var[[g]] <- eval(parse(text=paste0(ML2.key[s,'stat.vars'],'(ML2.sr[[g]])',collapse="")))
          
          # Check equal variance assumption
          if(!is.na(testVarEqual)){
            if(testVarEqual){
              logtxt <- paste(s,ML2.key$study.analysis[[s]],'-', runGroups[g])
              ML2.in$stat.params$var.equal <- decide.EqualVar(ML2.var[[g]],ML2.in$study.vars.labels, ML2.key[s, ], group = logtxt)
            }}
          
          # Run the analysis according to ML2.key: 'stat.test'
          stat.params <<- ML2.in$stat.params
          stat.test   <- try.CATCH(with(ML2.var[[g]],eval(parse(text = ML2.key$stat.test[[s]]))))
          
          # if(grepl("approximation", stat.test$warning)){stat.test$warning <- NULL}
          
          if(all(is.null(stat.test$warning), grepl("simpleWarning",stat.test$warning),
                 !grepl("Error", stat.test$value[[1]]),
                 !grepl("message", names(unlist(stat.test))[1]))){
            stat.test  <- stat.test$value
            ConsoleOut <- paste(capture.output(print(stat.test)),collapse="\n")
            listIT     <- TRUE
          }
          
          # START RECORD DATA -------------------------------------
          
          if(listIT){
            
            describe <- get.descriptives(stat.test = stat.test,
                                         vars      = ML2.var[[g]],
                                         keytable  = ML2.key[s,])
            
            var.lor <- ifelse(grepl("OR",describe$test$estype),
                              sum(1/(table(ML2.var[[g]]$Condition,ML2.var[[g]]$Response)), na.rm = TRUE),
                              NA)
            
            ESCI  <-   generateOutput(describe        = describe,
                                      var.lor         = var.lor,
                                      runningGroup    = runGroups[g],
                                      runningAnalysis = paste(s,ML2.key$study.analysis[[s]]))
            
            # Raw and clean datasets
            
            if(length(ML2.sr[[g]]$RawDataFilter)>1){
              case.include <- ML2.sr[[g]]$RawDataFilter[[1]]$Included|ML2.sr[[g]]$RawDataFilter[[2]]$Included
              df1 <- ML2.sr[[g]]$RawDataFilter[[1]][ ,-which(colnames( ML2.sr[[g]]$RawDataFilter[[1]])=="Included")]
              raw.df[[g]] <-  cbind.data.frame(df1, analysis.type = c("Global","Primary","Secondary","Order")[tp[cnt]],subset=subset,case.include = case.include)
            } else {
              case.include <- ML2.sr[[g]]$RawDataFilter[[1]]$Included
              df1 <- ML2.sr[[g]]$RawDataFilter[[1]][ ,-which(colnames( ML2.sr[[g]]$RawDataFilter[[1]])=="Included")]
              raw.df[[g]] <-  cbind.data.frame(df1, analysis.type = c("Global","Primary","Secondary","Order")[tp[cnt]],subset=subset,cases.include = case.include)
            }
            
            
            if(tp<4){
              if(runGroups[g]!="all"){
                fID <- unique(ML2.df$.id[ML2.df$source==runGroups[g]])
                sID <- SourceInfoTable$Source%in%runGroups[g]&SourceInfoTable$Filename%in%fID
                if(sum(sID)==1){
                  SourceInfo1 <- SourceInfoTable[sID, ]
                  SourceInfo2 <- raw.df[[g]] %>% filter(case.include) %>% group_by(source) %>%
                    summarise(
                      N.sources.global    = length(unique(Source.Global)),
                      N.sources.primary   = length(unique(Source.Primary)),
                      N.sources.secondary = length(unique(Source.Secondary)),
                      N.countries         = length(unique(Country)),
                      N.locations         = length(unique(Location)),
                      N.languages         = length(unique(Language)),
                      Pct.WEIRD           = mean(Weird, na.rm=TRUE)*100,
                      Tbl.Execution       = paste0(capture.output(table(Execution)),collapse="\n"),
                      Tbl.subjectpool     = paste0(capture.output(table(SubjectPool)),collapse="\n"),
                      Tbl.setting       = paste0(capture.output(table(Setting)),collapse="\n"),
                      Tbl.Tablet        = paste0(capture.output(table(Tablet)),collapse="\n"),
                      Tbl.Pencil        = paste0(capture.output(table(Pencil)),collapse="\n"),
                      N.studyorders1    = length(unique(StudyOrderN)),
                      N.IDiffOrderN     = length(unique(IDiffOrderN)),
                      N.uIDs            = length(unique(uID)),
                      N.studyorders2    = length(unique(study.order)),
                      Tbl.analysistype  = paste0(capture.output(table(analysis.type)),collapse="\n"),
                      Tbl.subset        = paste0(capture.output(table(subset)),collapse="\n") #,
                     # N.cases.included  = sum(case.include, na.rm=TRUE),
                     # N.cases.excluded  = sum(raw.df[[g]]$case.include==FALSE,na.rm=TRUE)
                    )
                  SourceInfo<-cbind(SourceInfo1,SourceInfo2)
                  # colnames(SourceInfo) <- c("name","name.Global",colnames(SourceInfoTable)[3:NCOL(SourceInfoTable)])
                }
              } else {
                
                SourceInfo <- raw.df[[g]] %>% filter(case.include) %>%
                  summarise(
                    N.sources.global    = length(unique(Source.Global)),
                    N.sources.primary   = length(unique(Source.Primary)),
                    N.sources.secondary = length(unique(Source.Secondary)),
                    N.countries         = length(unique(Country)),
                    N.locations         = length(unique(Location)),
                    N.languages         = length(unique(Language)),
                    Pct.WEIRD           = mean(Weird, na.rm=TRUE)*100,
                    Tbl.Execution       = paste0(capture.output(table(Execution)),collapse="\n"),
                    Tbl.subjectpool     = paste0(capture.output(table(SubjectPool)),collapse="\n"),
                    Tbl.setting       = paste0(capture.output(table(Setting)),collapse="\n"),
                    Tbl.Tablet        = paste0(capture.output(table(Tablet)),collapse="\n"),
                    Tbl.Pencil        = paste0(capture.output(table(Pencil)),collapse="\n"),
                    N.studyorders1    = length(unique(StudyOrderN)),
                    N.IDiffOrderN     = length(unique(IDiffOrderN)),
                    N.uIDs            = length(unique(uID)),
                    N.studyorders2    = length(unique(study.order)),
                    Tbl.analysistype  = paste0(capture.output(table(analysis.type)),collapse="\n"),
                    Tbl.subset        = paste0(capture.output(table(subset)),collapse="\n")
                  )
                    #N.cases.included  = n(),
                    #N.cases.excluded  = sum(raw.df[[g]]$case.include==FALSE,na.rm=TRUE)
              }
            } else {
              SourceInfo <- raw.df[[g]] %>% filter(case.include) %>%
                summarise(
                  N.sources.global    = length(unique(Source.Global)),
                  N.sources.primary   = length(unique(Source.Primary)),
                  N.sources.secondary = length(unique(Source.Secondary)),
                  N.countries         = length(unique(Country)),
                  N.locations         = length(unique(Location)),
                  N.languages         = length(unique(Language)),
                  Pct.WEIRD           = mean(Weird, na.rm=TRUE)*100,
                  Tbl.Execution       = paste0(capture.output(table(Execution)),collapse="\n"),
                  Tbl.subjectpool     = paste0(capture.output(table(SubjectPool)),collapse="\n"),
                  Tbl.setting       = paste0(capture.output(table(Setting)),collapse="\n"),
                  Tbl.Tablet        = paste0(capture.output(table(Tablet)),collapse="\n"),
                  Tbl.Pencil        = paste0(capture.output(table(Pencil)),collapse="\n"),
                  N.studyorders1    = length(unique(StudyOrderN)),
                  N.IDiffOrderN     = length(unique(IDiffOrderN)),
                  N.uIDs            = length(unique(uID)),
                  N.studyorders2    = length(unique(study.order)),
                  Tbl.analysistype  = paste0(capture.output(table(analysis.type)),collapse="\n"),
                  Tbl.subset        = paste0(capture.output(table(subset)),collapse="\n")
                  # N.cases.included  = n(),
                  #  N.cases.excluded  = sum(raw.df[[g]]$case.include==FALSE,na.rm=TRUE)
                )
            }
            
            rownames(SourceInfo) <- NULL
            
            test  <- describe$test
            descr <- describe$descr.raw
            outputSource[[g]] <- get.output(key      = ML2.key[s,],
                                            vars     = ML2.var[[g]],
                                            descr    = descr,
                                            group    = runGroups[g],
                                            analysis = c("Global","Primary","Secondary","Order")[tp[cnt]],
                                            varEqual = stat.params$var.equal,
                                            test     = test,
                                            ESCI     = ESCI,
                                            test.ConsoleOutput = ConsoleOut,
                                            SourceInfo = SourceInfo,
                                            stat.test = stat.test)
            
            # Data list for output to spreadsheet
            dataSource[[g]] <- list(
              study.id      = ML2.key$study.id[[s]],
              study.slate   = ML2.key$study.slate[[s]],
              study.name    = ML2.key$study.name[[s]],
              study.source  = runGroups[g],
              analysis.type = c("Global","Primary","Secondary","Order")[tp[cnt]],
              analysis.name = ML2.key$study.analysis[[s]],
              subset        = subset,
              stat.info     = ML2.in,
              stat.data.cleanchain = ML2.id,
              stat.data.raw       = raw.df[[g]],
              stat.data.cleaned   = ML2.sr[[g]][1:length(ML2.sr[[g]])-1],
              stat.data.analysed  = ML2.var[[g]][1:length(ML2.var[[g]])-1],
              stat.test = stat.test)
            
            suppressMessages(clean.df[[g]] <- ldply(dataSource[[g]]$stat.data.analysed,melt))
            colnames(clean.df[[g]])[colnames(clean.df[[g]])==".id"] <- "Condition"
            
            rm(stat.params)
          } else {# LISTIT
            cat("\nListIT = FALSE\n")
            if(grepl("observations",as.character(stat.test$value))){
              disp(paste(s, ML2.key$study.analysis[[s]],'-',
                         runGroups[g],'>> Not enough observations'),
                   header = FALSE, footer = FALSE)
            } else {
              disp(paste(s,ML2.key$study.analysis[[s]],'-', runGroups[g],'>> stat.test failed:'),
                   header = FALSE, footer = FALSE)
              # disp(paste('value: ',stat.test$value),
              #      header = FALSE, footer = FALSE)
              disp(paste('warning:',stat.test$warning),
                   header = FALSE, footer = FALSE)
            }
            ConsoleOut <- paste(gsub("[[:punct:]]", "", stat.test$warning, perl = TRUE), collapse="\n")
            NN <- lengths(ML2.var[[g]])
            NN <- NN[!names(NN)=="N"]
            N  <- rep(ML2.var[[g]]$N,length.out = length(NN))
            ML2.rnd <- llply(seq_along(NN), function(nn) rnorm(N[nn]))#eval(parse(text=paste(names(NN[nn])," = rnorm(101)"))))
            names(ML2.rnd) <- names(NN)
            stat.test  <- try.CATCH(with(ML2.var[[g]],eval(parse(text = ML2.key$stat.test[[s]]))))
            stat.test  <- stat.test$value
            #RANDOMdata <- TRUE
            # Analysis error, there may be descriptives, but set ESCI to NA
            ESCI[1:length(ESCI)] <- rep(NA,length(ESCI))
          }#Listit = FALSE
        } # all nMin 1,2
        
        # Report on errors
        
        if(!nMin1){
          disp(paste0(s,' ',ML2.key$study.analysis[[s]],' - ',
                      runGroups[g],' not included in results >> Cases in source file (',
                      sum(gID, na.rm = TRUE),') < Nmin.raw (',Nmin.raw,')'),
               header = FALSE, footer = FALSE)
        } # Check nMin 1}
        if(!nMin2){
          disp(paste0(s,' ',ML2.key$study.analysis[[s]],' - ',
                      runGroups[g],' not included in results >> Valid cases after varfun (n',
                      c(1,2)[compN < Nmin.cond],"=", compN[compN < Nmin.cond],') < Nmin.cond (',Nmin.cond,')'),
               header = FALSE, footer = FALSE)
        } # Double-check nMin
        
      } # iterate groups
      
      disp(paste(s, ML2.key$study.analysis[[s]],"- COMPLETED"), header = FALSE)
      
      ML2.output[[s]]  <- ldply(outputSource)
      ML2.rawdata[[s]] <- ldply(raw.df)
      
      if(outdir$ROBJECTS!=""){
        save(dataSource, file = file.path(rootdir,outdir$ROBJECTS,paste0(ML2.key$study.analysis[[s]],"_",c("Global","Primary","Secondary","Order")[tp[cnt]],".RData")))
      }
      
      rm(ML2.in, ML2.var, ML2.id, ML2.df, ML2.sr, outputSource, dataSource, raw.df, clean.df, descr, SourceInfo, nMin1, nMin2, listIT)
      
    } else { # if nrow > 0
      
      disp(paste(s, ML2.key$study.analysis[[s]],"- SKIPPED"), header = FALSE)
      
      ML2.output[[s]]  <- NULL
      ML2.rawdata[[s]] <- NULL
      
      rm(ML2.in, ML2.var, ML2.id, ML2.df, ML2.sr)
    }
    
  } # for s i studies
  
  options(wop)
 # return(list(raw.case   = ML2.rawdata,
  #            aggregated = ML2.output)
 # ) # ), descriptivesDump = ML2.descriptivesDump))
  