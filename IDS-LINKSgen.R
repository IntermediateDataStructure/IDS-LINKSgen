
  #########################
  ####Short description####
  #########################
  
  # This script transforms LINKS-IDS from the IDS format into the pedigree format. It was written by Rick Mourits, Ingrid van Dijk, and Kees Mandemakers.
  # For further information please see: Mourits, R.J., Van Dijk, I.K. & Mandemakers (2020). From matched certificates to related persons: Building a dataset 
  # from Links-Zeeland 2017.02. Historical Life Course Studies.
  
  # Version:              1.0
  # Date:                 2 September 2020
  # Purpose:              Transform LINKS-IDS into LINKS-gen
  # Rationale:            The conversion script is developed to transform LINKS-IDS into one rectangular dataframe that is readily available for statistical analysis.
  # Background:           The tool was developed in the context of the Genes, germs, and Resources project to be able to analyse LINKS-Zeeland 
  
  # Compatability:        The script only compatible with LINKS. Minor recoding is required to apply the script to other IDS data bases.
  # Input:                INDIVIDUAL and INDIV_INDIV from LINKS-Zeeland v2017.02 (to run script)
  #                       HISCO - US SES.txt (to enrich occupational data with OCC1950, OCCSCORE, PRESGL, SEI, and NPBOSS)
  # Outputs:              Running this script will produce two rectangular tables: Pedigree and Occupation. For a description of the data bases, see:
  #                             Mourits, R.J., Van Dijk, I.K. & Mandemakers, K. (2020). From matched certificates to related persons: Building a dataset 
  #                             from Links-Zeeland. Historical Life Course Studies.
  
  # Requirements:         In order for the program to be run successfully, the user needs to install R and the R-packages: data.table, dplyr, and lubridate.
  #                       These R-packages can be installed by running line 41 INSTALL.PACKAGES(....) without the #
  
  # Contact information:  Rick Mourits, r.j.mourits@uu.nl, Utrecht University, Utrecht, the Netherlands
  
  
  #####################################
  ####Step 1. Import and check files###
  #####################################
  
  # This section opens all required data files. To proceed please set the library where files are located using the function: SETWD()
  
  # Before you run the script, make sure that the checks in step 1 guarantee that: 
  #   1. Records in Individual are not duplicated
  #   2. Relations in Indiv_Indiv are not duplicated
  #   3. Each Id_I_2 in Indiv_Indiv is related to no more than 1 father and 1 mother
  
  
  #load required packages
    #install.packages("data.table"); install.packages("dplyr") #drop # to install data.table, dplyr, and lubridate
    library("data.table"); library("dplyr"); library("lubridate")
  
  #set working directory
    setwd("C:/Surfdrive/Data/LINKS-IDS/") #set library where all files are located
    getwd()
    
  #main files, CONTEXT is not required to construct Pedigree file
    Individual <- fread("INDIVIDUAL.txt")
  #relation files, CONTEXT_CONTEXT & INDIV_CONTEXT are not required to construct Pedigree file
    Indiv_Indiv <- fread("INDIV_INDIV.txt")
    
    
  #check whether entries in Individual, Context, and Indiv_Indiv are unique
    #1. Records in Individual are not duplicated, based on: person identifier, value, type, and date stamp
      length(which(duplicated(Individual[,c("Id_I", "Value", "Type", "Day", "Month", "Year")]))) #if not 0, do NOT run the script
    #2. Relations in Indiv_Indiv are not duplicated, based on: person identifiers person 1 and person identifier person 2
      length(which(duplicated(Indiv_Indiv[,c("Id_I_1", "Id_I_2")]))) #if not 0, do NOT run the script, relations should be the same across certificates (child, father, mother, bride, groom, etc.)
    #3. Each Id_I_2 in Indiv_Indiv is related to no more than 1 father and 1 mother
      length(which(duplicated(Indiv_Indiv[which(Indiv_Indiv$Relation=="Father"),c("Id_I_2")]))) #if not 0, do NOT run the script
      length(which(duplicated(Indiv_Indiv[which(Indiv_Indiv$Relation=="Mother"),c("Id_I_2")]))) #if not 0, do NOT run the script
      
  #load functions
      safe.ifelse <- function(cond, yes, no) structure(ifelse(cond, yes, no), class = class(yes))
      
    
  
  #############################################################################################################
  ###Step 2. Construct basic pedigree file with identifier + sex for ego, and identifier for mother & father###
  #############################################################################################################
  
  #This section gathers the identifier + sex of all mentioned persons in one data frame named Pedigree, and add relations to parents through the parental identifiers: Id_I_1.
  # 2a. lists all persons in LINKS_IDS and their sex
  # 2b. adds the person's father from LINKS_IDS to Pedigree
  # 2c. adds the person's mother from LINKS_IDS to Pedigree
      
  #2a. use all persons ever mentioned
   #with known sex
    #from Individual table with known sex
      Pedigree <- Individual[which(Individual$Type=="SEX"), c("Id_I", "Value", "Year")]
      colnames(Pedigree) <- c("Id_person", "Sex", "Year")
    #from Indiv_Indiv table with implicitly known sex
      #women
      Pedigree2 <- Indiv_Indiv[which(Indiv_Indiv$Relation=="Bride" | Indiv_Indiv$Relation=="Wife" | Indiv_Indiv$Relation=="Mother" | Indiv_Indiv$Relation=="Mother-in-law"), c("Id_I_1", "Relation", "Year")]
      colnames(Pedigree2) <- c("Id_person", "Sex", "Year")
      Pedigree2$Sex <- "Female"
      #men
      Pedigree3 <- Indiv_Indiv[which(Indiv_Indiv$Relation=="Groom" | Indiv_Indiv$Relation=="Husband" | Indiv_Indiv$Relation=="Father" | Indiv_Indiv$Relation=="Father-in-law"), c("Id_I_1", "Relation", "Year")]
      colnames(Pedigree3) <- c("Id_person", "Sex", "Year")
      Pedigree3$Sex <- "Male"
    #bind and select first known sex, as it is most reliable
      Pedigree <- rbind(Pedigree, Pedigree2, Pedigree3)
      Pedigree <- Pedigree %>% arrange(Id_person, Year, Sex) #sort by date
      Pedigree <- rbind(Pedigree[which(Pedigree$Sex!="Unknown"),],Pedigree[which(Pedigree$Sex=="Unknown"),]) #prioritize female/male over unknown 
      Pedigree <- Pedigree[which(!duplicated(Pedigree$Id_person)),] #select first mentioned age
      Pedigree$Year <- NULL
      rm(Pedigree2, Pedigree3)
  #clean environment
      length(Pedigree$Id_person)
      
  #2b. add identifier father
    #select identifier father
      Father <- Indiv_Indiv[which(Indiv_Indiv$Relation=="Father"),c("Id_I_1", "Id_I_2")]
      colnames(Father) <- c("Id_father", "Id_person")
    #merge
      Pedigree <- merge(Pedigree, Father, by="Id_person", all.x=T)
      
  #2c. add identifier mother
      Mother <- Indiv_Indiv[which(Indiv_Indiv$Relation=="Mother"),c("Id_I_1", "Id_I_2")]
      colnames(Mother) <- c("Id_mother", "Id_person")
    #merge
      Pedigree <- merge(Pedigree, Mother, by="Id_person", all=T)
      
  #2d. clean environment
      rm(Father, Mother)
      
      
  
  ################################################################################################
  ###Step 3. Add identifier partner + location + date for all known marriages to Pedigree file ###
  ################################################################################################
  
  #This section adds marital information to the basic Pedigree file established in step 2
  # 3a. Selects all known marriages in the data frame Marriages
  # 3b. Checks whether information on the marriage death in Marriages is not missing
  # 3c. Adds the identifier of the marriage partner to Marriages
  # 3d. Adds the marriage partner, date, and location from Marriages chronologically to Pedigree
  # 3e. Counts the number of known marriages for each person in Pedigree and adds a flag on the availability of marital information
      
  #make backup
      Backup <- Pedigree
      
      
  #3a. Select Entries
      Marriages <- Individual[which(Individual$Type=="MARRIAGE_LOCATION"),] 
      
  #3b. Inspect Marriages
    #number of marriages
      length(which(Individual$Type=="MARRIAGE_LOCATION")); length(Marriages$Type); length(which(Marriages$Estimation=="Exact")) #if not identical, fix before continuing the script
    #check whether dates are missing
      length(which(Marriages$Day==0 | is.na(Marriages$Day))) #if not 0, fix before continuing the script
      length(which(Marriages$Month==0 | is.na(Marriages$Month))) #if not 0, fix before continuing the script
      length(which(Marriages$Year==0 | is.na(Marriages$Year))) #if not 0, fix before continuing the script
    #check for duplicates
      length(which(duplicated(Marriages[,c("Id_I", "Year", "Month", "Day")])))  #if not 0, fix before continuing the script
      
  #3c. Add identifier marriage partners
      Marriage_ID <- Indiv_Indiv[which(Indiv_Indiv$Relation=="Bride" | Indiv_Indiv$Relation=="Groom"), c("Id_I_1", "Id_I_2", "Day", "Month", "Year")]
      colnames(Marriage_ID) <- c("Id_I", "Id_partner", "Day", "Month", "Year")
      Marriages <- merge(Marriages, Marriage_ID, by=c("Id_I", "Day", "Month", "Year"), all.x=T)
      
  #3d. Construct variables
    #Compile marriage date
      Marriages$M_date <- paste(Marriages$Year, Marriages$Month, Marriages$Day, sep="-")
      Marriages$M_date <- ymd(Marriages$M_date)
    #Determine marriage order & 
      Marriages <- Marriages %>% arrange(Id_I, M_date) %>% group_by(Id_I) %>% mutate(Marriage_N=n(), M_number=row_number()) %>% ungroup()
    #Determine marriage location
      Marriages$M_location <- Marriages$Value_Id_C
    #Add Id_partner, M_date & M_location to Pedigree
      #use repeat function to link marriage 1 to marriage N (max = currently 5)
      x <- 1
      end <- max(Marriages$Marriage_N)
      repeat{
        #filter
        df <- Marriages[which(Marriages$M_number==x),c("Id_I", "Id_partner", "M_date", "M_location")]
        colnames(df)[2:4] <- paste(c("Id_partner", "M_date", "M_location"), x, sep="_")
        Pedigree <- merge(Pedigree, df, by.x="Id_person", by.y="Id_I", all=T)
        if (x==end){
          break}
        x <- x+1
      }
      rm(x,end,df)
      
  #3e. Flag available information
    #Add Marriage_N
      Marriages <- Marriages[which(!duplicated(Marriages$Id_I)),c("Id_I", "Marriage_N")]
      Pedigree <- merge(Pedigree, Marriages, by.x="Id_person", by.y="Id_I", all=T)
      Pedigree$Marriage_N[is.na(Pedigree$Marriage_N)] <- 0
    #Flag marriage certificate
      Pedigree$M <- safe.ifelse(Pedigree$Marriage_N>0, 1, 0)
    #Flag marriage certificate parents
      M_parents <- Indiv_Indiv[which(Indiv_Indiv$Relation=="Bride"), c("Id_I_1", "Id_I_2")]
      colnames(M_parents) <- c("Id_mother", "Id_father")
      M_parents$M_parents <- 1
      Pedigree <- merge(Pedigree, M_parents, by=c("Id_mother", "Id_father"), all.x=T)
      Pedigree$M_parents[is.na(Pedigree$M_parents)] <- 0
      
  #3f. clean environment
      rm(Marriage_ID, Marriages, M_parents)
      
  
  
  ###############################################################################
  ###Step 4. Add information on birth/death date and location to pedigree file###
  ###############################################################################
  
  #This section adds information on the date and location of birth and/or death to the Pedigree file.
  # 4a. Adds the date of death to Pedigree; flags the availability of death information in Pedigree
  # 4b. Adds the location of death to Pedigree
  # 4c. Adds exact birth dates to Pedigree and estimates birth dates from available time ranges; flags the availability of death information in Pedigree
  # 4d. Adds the location of birth to Pedigree
      
  #make backup
      Backup2 <- Pedigree
      
      
  #4a. Construct death date
      Death <- Individual[which(Individual$Type=="DEATH_DATE" | Individual$Type=="STILLBIRTH_DATE"),]
    #check whether dates are missing
      length(which(Death$Day==0 | is.na(Death$Day))) #if not 0, fix before continuing the script
      length(which(Death$Month==0 | is.na(Death$Month))) #if not 0, fix before continuing the script
      length(which(Death$Year==0 | is.na(Death$Year))) #if not 0, fix before continuing the script
    #check for duplicates
      length(which(duplicated(Death[,c("Id_I", "Year", "Month", "Day")])))  #if not 0, fix before continuing the script
    #Compile death date
      Death$D_date <- paste(Death$Year, Death$Month, Death$Day, sep="-")
      Death$D_date <- ymd(Death$D_date)
    #flag available dead certificates
      Death$D <- 1
    #flag stillbirths
      Death$D_deadonregistration <- safe.ifelse(Death$Type=="STILLBIRTH_DATE", 1, 0)
    #add to Pedigree
      Death <- Death[,c("Id_I", "D_date", "D", "D_deadonregistration")]
      Pedigree <- merge(Pedigree, Death, by.x="Id_person", by.y="Id_I", all=T)
      
  #4b. Add death location
      Death <- Individual[which(Individual$Type=="DEATH_LOCATION" | Individual$Type=="STILLBIRTH_LOCATION"),]
    #rename death location
      Death$D_location <- Death$Value_Id_C
    #add to pedigree
      Death <- Death[,c("Id_I", "D_location")]
      Pedigree <- merge(Pedigree, Death, by.x="Id_person", by.y="Id_I", all=T)
      Pedigree$D[is.na(Pedigree$D)] <- 0
      
  #4c. Construct birth date
      Birth <- Individual[which(Individual$Type=="BIRTH_DATE"),]
      length(which(duplicated(Birth$Id_I))) #if not 0, fix before continuing the script
    #select reported dates
      Birth_exact <- Birth[which(Birth$Estimation=="Exact"),]
        #check whether exact dates are missing
          length(which(Birth_exact$Day==0 | is.na(Birth_exact$Day))) #if not 0, fix before continuing the script
          length(which(Birth_exact$Month==0 | is.na(Birth_exact$Month))) #if not 0, fix before continuing the script
          length(which(Birth_exact$Year==0 | is.na(Birth_exact$Year))) #if not 0, fix before continuing the script
        #check for duplicates
          length(which(duplicated(Birth_exact[,c("Id_I", "Year", "Month", "Day")])))  #if not 0, fix before continuing the script
        #Compile Birth_exact date
          Birth_exact$B_date <- paste(Birth_exact$Year, Birth_exact$Month, Birth_exact$Day, sep="-")
          Birth_exact$B_date <- ymd(Birth_exact$B_date)
        #Assign birth range
          Birth_exact$B_min_date <- Birth_exact$B_date
          Birth_exact$B_max_date <- Birth_exact$B_date
        #Compile
          Birth_exact <- Birth_exact[,c("Id_I", "B_date", "B_min_date", "B_max_date")]
        #Flag available birth certificate
          Birth_exact$B <- 1
    #select age-based estimates
      Birth_age_based <- Birth[which(Birth$Estimation=="Age_based"),]
        #Check whether exact dates are missing
          length(which(Birth_age_based$Day==0 | is.na(Birth_age_based$Day))) #if not 0, fix before continuing the script
          length(which(Birth_age_based$Month==0 | is.na(Birth_age_based$Month))) #if not 0, fix before continuing the script
          length(which(Birth_age_based$Year==0 | is.na(Birth_age_based$Year))) #if not 0, fix before continuing the script
        #Check for duplicates
          length(which(duplicated(Birth_age_based[,c("Id_I", "Year", "Month", "Day")])))  #if not 0, fix before continuing the script
        #Compile max birth date
          Birth_age_based$B_max_date <- paste(Birth_age_based$Year, Birth_age_based$Month, Birth_age_based$Day, sep="-")
          Birth_age_based$B_max_date <- ymd(Birth_age_based$B_max_date)
        #Compute min birth date
          Birth_age_based$B_min_date <- Birth_age_based$B_max_date %m-% years(1) %m+% days(1)
        #Compute mean birth date
          Birth_age_based$B_date <- Birth_age_based$B_max_date %m-% days(182)
        #Compile
          Birth_age_based <- Birth_age_based[,c("Id_I", "B_date", "B_min_date", "B_max_date")]
        #Flag unavailable birth certificate
          Birth_age_based$B <- 0
    #select rule-based estimates
      Birth_rule_based <- Birth[which(Birth$Estimation=="Rule_based"),]
        #Check whether exact dates are missing
          length(which(Birth_rule_based$Start_day==0 | is.na(Birth_rule_based$Start_day))) #if not 0, fix before continuing the script
          length(which(Birth_rule_based$End_day==0 | is.na(Birth_rule_based$End_day))) #if not 0, fix before continuing the script
          length(which(Birth_rule_based$Start_month==0 | is.na(Birth_rule_based$Start_month))) #if not 0, fix before continuing the script
          length(which(Birth_rule_based$End_month==0 | is.na(Birth_rule_based$End_month))) #if not 0, fix before continuing the script
          length(which(Birth_rule_based$Start_year==0 | is.na(Birth_rule_based$Start_year))) #if not 0, fix before continuing the script
          length(which(Birth_rule_based$End_year==0 | is.na(Birth_rule_based$End_year))) #if not 0, fix before continuing the script
        #Compile min birth date
          Birth_rule_based$B_min_date <- paste(Birth_rule_based$Start_year, Birth_rule_based$Start_month, Birth_rule_based$Start_day, sep="-")
          Birth_rule_based$B_min_date <- ymd(Birth_rule_based$B_min_date)
        #Compile max birth date
          Birth_rule_based$B_max_date <- paste(Birth_rule_based$End_year, Birth_rule_based$End_month, Birth_rule_based$End_day, sep="-")
          Birth_rule_based$B_max_date <- ymd(Birth_rule_based$B_max_date)
        #Check whether range is logically possible
          length(which(Birth_rule_based$B_max_date<Birth_rule_based$B_min_date)) #if not 0, fix before continuing the script
        #Calculate mean birth date
          Birth_rule_based$B_date <- as.numeric(Birth_rule_based$B_max_date - Birth_rule_based$B_min_date)
          Birth_rule_based$B_date <- safe.ifelse(Birth_rule_based$B_date>0 & Birth_rule_based$B_date<=1095, #interval between 1 day and 3 years
                                                 Birth_rule_based$B_max_date %m-% days(round(Birth_rule_based$B_date/2)), NA)
        #Compile
          Birth_rule_based <- Birth_rule_based[,c("Id_I", "B_date", "B_min_date", "B_max_date")]
        #Flag unavailable birth certificate
          Birth_rule_based$B <- 0
    #Compile estimates from stillbirths 
      Stillbirths <- Individual[which(Individual$Type=="STILLBIRTH_DATE"),]
        #B_date
          Stillbirths$B_date <- paste(Stillbirths$Year, Stillbirths$Month, Stillbirths$Day, sep="-")
          Stillbirths$B_date <- ymd(Stillbirths$B_date)
        #B_min_date
          Stillbirths$B_min_date <- Stillbirths$B_date
        #B_max_date
          Stillbirths$B_max_date <- Stillbirths$B_date
        #Compile
          Stillbirths <- Stillbirths[,c("Id_I", "B_date", "B_min_date", "B_max_date")]
        #Flag unavailable birth certificate
          Stillbirths$B <- 'Deadonregistration'
    #compile birth dates
      Birth <- rbind(Birth_exact, Birth_age_based, Birth_rule_based, Stillbirths)
    #check for duplicates
      length(which(duplicated(Birth$Id_I))) #ok if not 0
      Birth <- Birth[which(!duplicated(Birth$Id_I)),]
    #merge
      Pedigree <- merge(Pedigree, Birth, by.x="Id_person", by.y="Id_I", all=T)
      Pedigree$B[is.na(Pedigree$B)] <- 0
      
  #4d. Birth locations
      Birth <- Individual[which(Individual$Type=="BIRTH_LOCATION" | Individual$Type=="STILLBIRTH_LOCATION"),]
      Birth <- Birth[,c("Id_I", "Value_Id_C")]
    #check for duplicates
      length(which(duplicated(Birth$Id_I))) #ok if not 0
      Birth <- Birth[which(!duplicated(Birth$Id_I)),]
    #merge
      colnames(Birth) <- c("Id_person", "B_location")
      Pedigree <- merge(Pedigree, Birth, by="Id_person", all=T)
      
  #4e. clean environment
      rm(Birth, Birth_exact, Birth_age_based, Birth_rule_based, Death, Stillbirths)
      
      
  
  ######################################
  ###Step 5. add age parents at birth###
  ######################################
  
  #This section compares a person's and his parents' dates of birth to compute the age of his parents at birth
  # 5a. calculates the mother's age at birth for the persons in Pedigree
  # 5b. calculates the father's age at birth for the persons in Pedigree
      
  #make backup
      Backup3 <- Pedigree
      
  #5a. Add age of the mother at birth
    #Add birth date mother
      Agemother <- Pedigree[, c("Id_person", "B_date")]
      colnames(Agemother)[2] <- "B_date_mother"
    #Calculate age mother at birth
      Pedigree <- merge(Pedigree, Agemother, by.x="Id_mother", by.y="Id_person", all.x=T)
      Pedigree$B_age_mother <- as.numeric(difftime(Pedigree$B_date, Pedigree$B_date_mother, units="days"))
      Pedigree$B_age_mother <- Pedigree$B_age_mother * 0.00273790926 #convert days to years
    #Clean environment
      Pedigree$B_date_mother <- NULL
      
  #5b. #Add age of the father at birth
    #Add birth date father
      Agefather <- Pedigree[, c("Id_person", "B_date")]
      colnames(Agefather)[2] <- "B_date_father"
    #Calculate age father at birth
      Pedigree <- merge(Pedigree, Agefather, by.x="Id_father", by.y="Id_person", all.x=T)
      Pedigree$B_age_father <- as.numeric(difftime(Pedigree$B_date, Pedigree$B_date_father, units="days"))
      Pedigree$B_age_father <- Pedigree$B_age_father * 0.00273790926 #convert days to years
    #Clean environment
      Pedigree$B_date_father <- NULL
      
  #5c.Clean environment
      rm(Agemother, Agefather)
      
      
      
  ##############################################################
  ###Step 6. clean fields that contain impossible information###
  ##############################################################
  
  #This section checks whether the established life course and relations for all individuals in the pedigree file are logically possible.
  # 6a. A mother has to be at least 14 years old when she got her child, if NOT all parental information is cleaned from the pedigree file
  # 6b. A father has to be at least 14 years old when he got his child, if NOT all parental information is cleaned from the pedigree file
  # 6c. A person has to be alive at marriage, if NOT all marital information is cleaned from the pedigree file
  # 6d. A person has to be at least 14 years old at marriage, if NOT marital information is cleaned from the pedigree file
  # 6e. A person cannot be born after his death, if NOT birth information is cleaned from the pedigree file
      
  #make backup
      Backup4 <- Pedigree
      
  #6a. Negative age mother
      as.data.frame(table(round(Pedigree[which(Pedigree$B_age_mother<14),"B_age_mother"])))
      #set parents info to NA
      Pedigree$M_parents[Pedigree$B_age_mother<0] <- 0
      Pedigree$Id_mother[Pedigree$B_age_mother<0] <- NA
      Pedigree$Id_father[Pedigree$B_age_mother<0] <- NA
      Pedigree$B_age_father[Pedigree$B_age_mother<0] <- NA
      Pedigree$B_age_mother[Pedigree$B_age_mother<0] <- NA
      
  #6b. Negative age father
      as.data.frame(table(round(Pedigree[which(Pedigree$B_age_father<14),"B_age_father"])))
      #set parents info to NA
      Pedigree$M_parents[Pedigree$B_age_father<0] <- 0
      Pedigree$Id_mother[Pedigree$B_age_father<0] <- NA
      Pedigree$Id_father[Pedigree$B_age_father<0] <- NA
      Pedigree$B_age_mother[Pedigree$B_age_father<0] <- NA
      Pedigree$B_age_father[Pedigree$B_age_father<0] <- NA
      
  #6c. Marriage date after deathdate
      length(which(Pedigree$M_date_1>Pedigree$D_date))
      length(which(Pedigree$M_date_2>Pedigree$D_date))
      #set deathdate to NA
      Pedigree$D[Pedigree$M_date_1>Pedigree$D_date] <- 0
      Pedigree$D_date[Pedigree$M_date_1>Pedigree$D_date] <- NA
      Pedigree$D_date[Pedigree$M_date_2>Pedigree$D_date] <- NA
      
  #6d. Marriage before ego was 14
      Pedigree$M_age <- as.numeric(difftime(Pedigree$M_date_1,Pedigree$B_date,units="days"))
      Pedigree$M_age <- Pedigree$M_age * 0.00273790926 #convert days to years
      head(as.data.frame(table(Pedigree$M_age)),20)
      #Set Marriage1 to NA
      Pedigree$M[Pedigree$M_age<14] <- 0
      Pedigree$M_date_1[Pedigree$M_age<14] <- NA
      Pedigree$M_location_1[Pedigree$M_age<14] <- NA
      Pedigree$Id_partner_1[Pedigree$M_age<14] <- NA
      Pedigree$Marriage_N <- safe.ifelse(Pedigree$M_age<14, Pedigree$Marriage_N-1, Pedigree$Marriage_N)
      #Set Marriage2 to Marriage1
      Pedigree$M <- safe.ifelse(is.na(Pedigree$Id_partner_1) & !is.na(Pedigree$Id_partner_2), 1, Pedigree$M)
      Pedigree$M_date_1 <- safe.ifelse(is.na(Pedigree$Id_partner_1) & !is.na(Pedigree$Id_partner_2), Pedigree$M_date_2, Pedigree$M_date_1)
      Pedigree$M_location_1 <- safe.ifelse(is.na(Pedigree$Id_partner_1) & !is.na(Pedigree$Id_partner_2), Pedigree$M_location_2, Pedigree$M_location_1)
      Pedigree$Id_partner_1 <- safe.ifelse(is.na(Pedigree$Id_partner_1) & !is.na(Pedigree$Id_partner_2), Pedigree$Id_partner_2, Pedigree$Id_partner_1)
      #rm M_age
      Pedigree$M_age <- NULL
      
  #6e. Marriage twice to the same person
      Pedigree$M_date_2[Pedigree$Id_partner_1==Pedigree$Id_partner_2] <- NA
      Pedigree$M_location_2[Pedigree$Id_partner_1==Pedigree$Id_partner_2] <- NA
      Pedigree$Id_partner_2[Pedigree$Id_partner_1==Pedigree$Id_partner_2] <- NA
      Pedigree$Marriage_N <- safe.ifelse(Pedigree$Id_partner_1==Pedigree$Id_partner_2, Pedigree$Marriage_N-1, Pedigree$Marriage_N)
      
  #6f. Deathdate before birthdate
      length(which(Pedigree$D_date<Pedigree$B_date))
      Pedigree$B_date <- safe.ifelse(is.na(Pedigree$D_date), Pedigree$B_date, 
                                     safe.ifelse(Pedigree$D_date<Pedigree$B_date & Pedigree$D_date>=Pedigree$B_min_date, Pedigree$B_min_date, Pedigree$B_date)
                                     )
      Pedigree$B[Pedigree$D_date<Pedigree$B_date] <- 0
      Pedigree$B_min_date[Pedigree$D_date<Pedigree$B_date] <- NA
      Pedigree$B_max_date[Pedigree$D_date<Pedigree$B_date] <- NA
      Pedigree$B_date[Pedigree$D_date<Pedigree$B_date] <- NA
      
  
  
  ###############################################
  ###Step 7. improve quality family structures###
  ###############################################
  
  #This section further improves the quality of life course reconstructions within reconstructed families. Within some families, birth and death certificates were succesfully linked
  #to the same parents, but not to each other. In these families, observations on the same individual were falsely recorded as observations on different individuals within the same family.
      
  #To solve this problem, we added unmatched death certificates to birth and marriage certificates of children without a death certificate, provided 
  #that the estimated birth date overlapped with the actual date of the birth certificate, and that an unmatched death certificates was relatable to only one child in the family.
      
  # 7a. Add all unmatched death certificates without birth information, chronologically to Pedigree
  # 7b. Filter possible matches if: (1) no death certificate was already linked, and 
  #                                 (2) the estimated birth range of a death certificate contains the established birth date on a birth certificate. 
  # 7c. List all possible matches in long format and select unique new matches
  # 7d. Add information from the row with available death information to the row with available birth information
  # 7e. Mark the rows with improved family structures
  # 7f. Delete the now redundant rows with available death information
  # 7g. In Indiv_Indiv & Individuals, replace the person identifiers from the removed rows with the new person identifiers
      
      
  #make backup
      Backup5 <- Pedigree
      
  #7a. Add unlinked death certificates to siblings
    #Select unlinked death records
      D <- Pedigree[which(B==0 & D==1 & !is.na(Id_mother) & !is.na(Id_father)), c("Sex", "Id_mother", "Id_father", "Id_person", "B_min_date", "B_max_date", "Id_partner_1", "M_date_1", "D_date")]
    #Drop cases with a birth range of over 3 years
      D <- D[which(as.numeric(D$B_max_date - D$B_min_date)<=1095),] #1095 = 3*365
    #Order chronologically witin family
      D <- D %>% arrange(Id_mother, Id_father, B_max_date)
    #Filter unmatched entries
      #use repeat function to link marriage 1 to marriage N
      x <- 1
      end <- max(as.data.frame(table(paste(D$Id_mother, D$Id_father, sep="")))$Freq) #maximum number of unmatched certificates in a family
      repeat{
        #filter
        df <- D %>% group_by(Id_mother, Id_father) %>% filter(row_number()==x)
        colnames(df)[4:9] <- paste(c("Id_person", "B_min_date", "B_max_date", "Id_partner_1", "M_date_1", "D_date"), x, sep="_")
        Pedigree <- merge(Pedigree, df, by=c("Id_mother", "Id_father", "Sex"), all=T)
        if (x==end){
          break}
        x <- x+1
      }
      rm(x,end,df)
      
  #7b. select possible matches
    #use repeat function to filter matches where B_min or B_max falls within estimated birth range 
      x <- 1
      end <- max(as.data.frame(table(paste(D$Id_mother, D$Id_father, sep="")))$Freq) #maximum number of unmatched certificates in a family
      repeat{
        #select cases if:
        # 1. Birth certificate is available
        # 2. Death certificates is not available
        # 3. Known birthdate falls within estimated birthrange of an unlinked Death certificate in the family
        Pedigree[[paste0("Id_person_", x)]] <- safe.ifelse(Pedigree$B==1 & Pedigree$D==0 & Pedigree$B_date>=Pedigree[[paste0("B_min_date_", x)]] & Pedigree$B_date<=Pedigree[[paste0("B_max_date_", x)]], #B_max_date is within estimated birth range
                                                        Pedigree[[paste0("Id_person_", x)]], NA)
        #drop cases with death before marriage
        Pedigree[[paste0("Id_person_", x)]][Pedigree$M_date_1>Pedigree[[paste0("D_date_", x)]]] <- NA #M_date_1 
        Pedigree[[paste0("Id_person_", x)]][Pedigree$M_date_2>Pedigree[[paste0("D_date_", x)]]] <- NA #M_date_2
        Pedigree[[paste0("Id_person_", x)]][Pedigree$M_date_3>Pedigree[[paste0("D_date_", x)]]] <- NA #M_date_3
        Pedigree[[paste0("Id_person_", x)]][Pedigree$M_date_4>Pedigree[[paste0("D_date_", x)]]] <- NA #M_date_4
        Pedigree[[paste0("Id_person_", x)]][Pedigree$M_date_5>Pedigree[[paste0("D_date_", x)]]] <- NA #M_date_5
        if (x==end){
          break}
        x <- x+1
      }
      rm(x,end)
    #clean other variables
      #use repeat function 
      x <- 1
      end <- max(as.data.frame(table(paste(D$Id_mother, D$Id_father, sep="")))$Freq) #maximum number of unmatched certificates in a family
      repeat{
        Pedigree[[paste0("B_min_date_", x)]][is.na(Pedigree[[paste0("Id_person_", x)]])] <- NA
        Pedigree[[paste0("B_max_date_", x)]][is.na(Pedigree[[paste0("Id_person_", x)]])] <- NA
        Pedigree[[paste0("Id_partner_1_", x)]][is.na(Pedigree[[paste0("Id_person_", x)]])] <- NA
        Pedigree[[paste0("M_date_1_", x)]][is.na(Pedigree[[paste0("Id_person_", x)]])] <- NA
        Pedigree[[paste0("D_date_", x)]][is.na(Pedigree[[paste0("Id_person_", x)]])] <- NA
        if (x==end){
          break}
        x <- x+1
      }
      rm(D,x,end)
  
  #7c. List possible matches
    #split new entries & pedigree file
      Pedigree2 <- cbind(Pedigree[,4],Pedigree[,34:length(colnames(Pedigree))])
      Pedigree <- Pedigree[,1:33]
    #list Pedigree2 in long format
      x <- 2
      y <- 7
      end <- length(colnames(Pedigree2))
      colnames(Pedigree2)[2:end] <- rep(c("Id_person_x", "B_min_date_x", "B_max_date_x", "Id_partner_1_x", "M_date_x", "D_date_x"),(end-1)/6)
      a <- cbind(Pedigree2[,1],Pedigree2[,x:y])
      repeat{
        x <- x+6
        y <- y+6
        a <- rbind(a,
                   cbind(Pedigree2[,1],Pedigree2[,x:y]))
        if (y==end){
          break}
      }
      Pedigree2 <- a
      rm(a,x,y,end)
    #filter valid entries
      Pedigree2 <- Pedigree2[which(!is.na(Pedigree2$Id_person_x)),]
    #drop multiple entries
      Pedigree2 <- Pedigree2 %>% group_by(Id_person_x) %>% filter(n()==1) %>% ungroup()
      Pedigree2 <- Pedigree2 %>% group_by(Id_person) %>% filter(n()==1) %>% ungroup()
      
  #7d. Add new information
      Pedigree <- merge(Pedigree, Pedigree2, by="Id_person", all.x=T)
      Pedigree$D_date <- safe.ifelse(is.na(Pedigree$D_date), Pedigree$D_date_x, Pedigree$D_date)
      Pedigree$M_date_2 <- safe.ifelse(!is.na(Pedigree$M_date_1) & !is.na(Pedigree$M_date_x), Pedigree$M_date_x, Pedigree$M_date_2)
      Pedigree$Id_partner_2 <- safe.ifelse(!is.na(Pedigree$Id_partner_2) & !is.na(Pedigree$Id_partner_1_x), Pedigree$Id_partner_1_x, Pedigree$Id_partner_2)
      Pedigree$M_date_1 <- safe.ifelse(is.na(Pedigree$M_date_2) & !is.na(Pedigree$M_date_x), Pedigree$M_date_x, Pedigree$M_date_1)
      Pedigree$Id_partner_1 <- safe.ifelse(is.na(Pedigree$Id_partner_2) & !is.na(Pedigree$Id_partner_1_x), Pedigree$Id_partner_1_x, Pedigree$Id_partner_1)
      
  #7e.Flag newly established relations
      Pedigree$Postlink_D <- safe.ifelse(!is.na(Pedigree$Id_person_x), 1, 0)
      
  #7f. Remove Id_person_x from Pedigree
      Pedigree <- Pedigree[which(!(Pedigree$Id_person %in% Pedigree2$Id_person_x)),]
      Pedigree <- Pedigree[,c(1:33,40)]
      
  #7g. Replace Id_I in Indiv_Indiv & Individual for deleted entries 
      Pedigree2 <- Pedigree2[,1:2]
      colnames(Pedigree2) <- c("Newid", "Oldid")
    #replace Individual$Id_I
      Individual2 <- merge(Individual, Pedigree2, by.x="Id_I", by.y="Oldid", all.x=T)
      Individual2$Id_I <- ifelse(!is.na(Individual2$Newid), Individual2$Newid, Individual2$Id_I)
      Individual2$Newid <- NULL
    #replace Indiv_Indiv$Id_I_1
      Indiv_Indiv2 <- merge(Indiv_Indiv, Pedigree2, by.x="Id_I_1", by.y="Oldid", all.x=T)
      Indiv_Indiv2$Id_I_1 <- ifelse(!is.na(Indiv_Indiv2$Newid), Indiv_Indiv2$Newid, Indiv_Indiv2$Id_I_1)
      Indiv_Indiv2$Newid <- NULL
    #replace Indiv_Indiv$Id_I_2
      Indiv_Indiv2 <- merge(Indiv_Indiv, Pedigree2, by.x="Id_I_2", by.y="Oldid", all.x=T)
      Indiv_Indiv2$Id_I_2 <- ifelse(!is.na(Indiv_Indiv2$Newid), Indiv_Indiv2$Newid, Indiv_Indiv2$Id_I_2)
      Indiv_Indiv2$Newid <- NULL
      
      
      
  #####################################
  ###Step 8. compute extra variables###
  #####################################
      
  #This section computes the last variables required for the Pedigree file:
  # 8a. Twins
  # 8b. Flag availability of parental birth and death certificate
  # 8c. Compute age at death
  # 8d. Compute age at marriage
  # 8e. Identify the last date + age when a person was alive, and flag certificate
      
  
  #make backup
      Backup6 <- Pedigree
      
  #8a. Flag twins
      Pedigree <- Pedigree %>% group_by(Id_mother, Id_father, B_date) %>% mutate(Twin=n()) %>% ungroup()
      Pedigree$Twin[is.na(Pedigree$Id_mother) | is.na(Pedigree$B_date)] <- NA
      Pedigree$Twin <- safe.ifelse(Pedigree$Twin>1, 1, 0) 
      
  #8b. Flag available certificates parents
    #Flag birth and death mother
      BD <- Pedigree[,c("Id_person", "B", "D")]
      colnames(BD) <- c("Id_mother", "B_mother", "D_mother")
      Pedigree <- merge(Pedigree, BD, by="Id_mother", all.x=T)
      Pedigree$B_mother <- safe.ifelse(is.na(Pedigree$B_mother), 0, Pedigree$B_mother)
      Pedigree$D_mother <- safe.ifelse(is.na(Pedigree$D_mother), 0, Pedigree$D_mother)
    #Flag birth and death father
      BD <- Pedigree[,c("Id_person", "B", "D")]
      colnames(BD) <- c("Id_father", "B_father", "D_father")
      Pedigree <- merge(Pedigree, BD, by="Id_father", all.x=T)
      Pedigree$B_father <- safe.ifelse(is.na(Pedigree$B_father), 0, Pedigree$B_father)
      Pedigree$D_father <- safe.ifelse(is.na(Pedigree$D_father), 0, Pedigree$D_father)
    #Clean environment
      rm(BD)
     
  #8c. Compute age at death
    #make function
      compute_age <- function(datevar1, datevar2){
        as.numeric(difftime(datevar1, datevar2 ,units="days")) * 0.00273790926 #convert days to years
      }
    #run function
      Pedigree$D_age <- compute_age(Pedigree$D_date, Pedigree$B_date)
      
  #8d. Compute age at marriage
      #use repeat function for marriage 1 to marriage N
      x <- 1
      end <- length(which(grepl("M_date", colnames(Pedigree))))
      repeat{
        #filter
        Pedigree[[paste0("M_age_", x)]] <- compute_age(Pedigree[[paste0("M_date_", x)]], Pedigree$B_date)
        if (x==end){
          break}
        x <- x+1
      }
      rm(x,end)
      
  #8e. Mark last observation when alive
    #make table with entries
      Pedigree <- as.data.frame(Pedigree)
      Entries <- cbind(Pedigree[,c("Id_person", "B_date", "D_date")], Pedigree[,grepl("M_date_", names(Pedigree))]) #birth, date, and marriage dates
    #retrieve last childbirth
     #mothers
      LastBirth_f <- Pedigree %>% filter(!is.na(Id_father) & !is.na(B_date)) %>% group_by(Id_mother) %>% summarise(LastBirth_f=max(B_date)) %>% ungroup()
      Entries <- merge(Entries, LastBirth_f, by.x="Id_person", by.y="Id_mother", all.x=T)
     #fathers
      LastBirth_m <- Pedigree %>% filter(!is.na(Id_father) & !is.na(B_date)) %>% group_by(Id_father) %>% summarise(LastBirth_m=max(B_date)) %>% ungroup()
      LastBirth_m$LastBirth_m <- LastBirth_m$LastBirth_m %m-% months(9) #father was alive at conception
      Entries <- merge(Entries, LastBirth_m, by.x="Id_person", by.y="Id_father", all.x=T)
     #cast into one variable
      Entries$LastBirth <- safe.ifelse(!is.na(Entries$LastBirth_f), Entries$LastBirth_f, Entries$LastBirth_m)
     #clean Entries
      Entries$LastBirth_f <- NULL
      Entries$LastBirth_m <- NULL
    #find last entry date
      colnames(Entries) #check entries
      LastEntryDate <- apply(Entries[,2:length(Entries)], 1, max, na.rm=T)
      LastEntryDate[LastEntryDate=="-Inf"] <- NA
      Entries <- cbind(Entries, LastEntryDate) #Link to Id_Person
      Entries$LastEntryDate <- ymd(Entries$LastEntryDate)
     #remove cases where mentioned parent is dead
      length(which(Entries$LastEntryDate>Entries$D_date))
      Entries$LastEntryDate <- safe.ifelse(is.na(Entries$D_date), Entries$LastEntryDate,
                                           safe.ifelse(Entries$LastEntryDate>Entries$D_date, Entries$D_date, Entries$LastEntryDate))
    #calculate last entry age
      Entries$LastEntryAge <- compute_age(Entries$LastEntryDate, Entries$B_date)
    #mark certificate of last entry
     #set function
      assign_cert <- function(varname, varcat){
        safe.ifelse(is.na(varname), Entries$LastEntryCert,
                    safe.ifelse(Entries$LastEntryDate==varname, varcat, Entries$LastEntryCert))
      }
     #run function
      Entries$LastEntryCert <- NA
      Entries$LastEntryCert <- assign_cert(Entries$B_date, "Birth")
      Entries$LastEntryCert <- assign_cert(Entries$M_date_1, "First marriage")
      Entries$LastEntryCert <- assign_cert(Entries$LastBirth, "Birth child")
      Entries$LastEntryCert <- assign_cert(Entries$D_date, "Death")
      Entries$LastEntryCert <- safe.ifelse(is.na(Entries$LastEntryCert) & !is.na(Entries$LastEntryDate), "Second marriage", Entries$LastEntryCert)
    #merge
      Entries <- Entries[,c("Id_person", "LastEntryCert", "LastEntryAge", "LastEntryDate")]
      Pedigree <- merge(Pedigree, Entries, by="Id_person", all.x=T)
    #clean envrionment
      rm(Entries, LastEntryDate, assign_cert)
      
      
      
  ############################
  ###Step 9. Make SES table###
  ############################
      
  #This section makes a new data frame containing all occupational information. Data is stored in long format by Id_person & Date
  # 9a. Make a basic data set containing all mentioned occupations, the corresponding person identifier, and the year, month, and day of observation
  # 9b. Compile a date variable
  # 9c. Construct an age variable for individuals with a known or estimated date of birth
  # 9d. Match the location and source certificate
  # 9e. Add occupational codes
  # 9f. Deal with occupational codes that indicate missings
  # 9g. Add American occupational codes
      
      
  #9a. make basic file
      Occupation <- Individual2[which(Individual2$Type=="OCCUPATION_STANDARD"),c("Id_I", "Value", "Day", "Month", "Year")]
      Occupation <- as.data.frame(Occupation)
      colnames(Occupation)[1:2] <- c("Id_person", "Occupation")
      #remove repeated entries
      Occupation <- Occupation[which(!duplicated(Occupation[,c("Id_person", "Day", "Month", "Year")])),] #combination person identifier and time stamp not unique
      
  #9b. construct date 
      Occupation$Date <- paste(Occupation$Year, Occupation$Month, Occupation$Day, sep="-")
      Occupation$Date <- ymd(Occupation$Date)
      Occupation$Year <- NULL; Occupation$Month <- NULL; Occupation$Day <- NULL
      
  #9c. add birthdate & construct age
      Occupation <- merge(Occupation, Pedigree[,c("Id_person","B_date")], by="Id_person", all.x=T)
      Occupation$Age <- as.numeric(difftime(Occupation$Date, Occupation$B_date, units="days"))
      Occupation$Age <- Occupation$Age * 0.00273790926 #convert days to years
      Occupation$B_date <- NULL
      
  #9d. add location + certificate from certificates of ego
    #find all available certificates entries of ego
      Locations <- cbind(Pedigree[,c("Id_person", "D_date", "D_location")], 
                         Pedigree[,grepl("M_date", names(Pedigree))], #marriage dates
                         Pedigree[,grepl("M_location", names(Pedigree))] ) #marriage locations
    #merge
      Occupation <- merge(Occupation, Locations, by="Id_person", all.x=T)
      rm(Locations)
    #add location + certificate
     #set function
      assign_cert_loc <- function(varname, varcat, varlocation){
        #Occupation$Cert
          Occupation$Cert <<- safe.ifelse(is.na(varname), Occupation$Cert,
                                               safe.ifelse(varname==Occupation$Date, varcat, Occupation$Cert))
         #Occupation$Location
          Occupation$Location <<- safe.ifelse(is.na(varname), Occupation$Location,
                                                   safe.ifelse(varname==Occupation$Date, varlocation, Occupation$Location))
      }
     #run function
      Occupation$Cert <- NA ; Occupation$Location <- NA
      assign_cert_loc(Occupation$D_date, "Death", Occupation$D_location)
      #within loop for marriage certificates
      x <- 1
      end <- length(which(grepl("M_date", colnames(Pedigree))))
      repeat{
        #filter
        assign_cert_loc(Occupation[[paste0("M_date_", x)]], paste0("Marriage",x), Occupation[[paste0("M_location_", x)]])
        if (x==end){
          break}
        x <- x+1
      }
      #remove repeated entries
      Occupation <- Occupation[which(!duplicated(Occupation[,c("Id_person", "Date")])),]
      rm(x,end)
  #9d. add location + certificate from certificates of children
    #split Occupations
      Occupation2 <- Occupation[which(is.na(Occupation$Cert)),c("Id_person", "Occupation", "Date", "Age")]
      Occupation <- Occupation[which(!is.na(Occupation$Cert)),c("Id_person", "Occupation", "Date", "Age", "Location", "Cert")]
    #find all available certificates of children
      Locations_father <- cbind(Pedigree[which(!is.na(Pedigree$Id_father)),c("Id_father", "B_date", "B_location", "D_date", "D_location")],
                                Pedigree[which(!is.na(Pedigree$Id_father)),grepl("M_date", names(Pedigree))], #marriage dates
                                Pedigree[which(!is.na(Pedigree$Id_father)),grepl("M_location", names(Pedigree))] ) #marriage locations
      Locations_mother <- cbind(Pedigree[which(!is.na(Pedigree$Id_mother)),c("Id_mother", "B_date", "B_location", "D_date", "D_location")],
                                Pedigree[which(!is.na(Pedigree$Id_mother)),grepl("M_date", names(Pedigree))], #marriage dates
                                Pedigree[which(!is.na(Pedigree$Id_mother)),grepl("M_location", names(Pedigree))] ) #marriage locations
      #bind observations
        colnames(Locations_father) <- paste0(colnames(Locations_father),"_child")
        colnames(Locations_father)[1] <- "Id_person"
        colnames(Locations_mother) <- paste0(colnames(Locations_mother),"_child")
        colnames(Locations_mother)[1] <- "Id_person"
        Locations <- rbind(Locations_father, Locations_mother)
        rm(Locations_father, Locations_mother)
      #order and number
        Locations <- Locations %>% arrange(Id_person, B_date_child) %>% group_by(Id_person) %>% mutate(child=row_number()) %>% ungroup()
    #merge
      Occupation3 <- merge(Occupation2, Locations, by="Id_person", all.x=T)
      rm(Locations)
    #add location + certificate
     #set function
      #set function
      assign_cert_loc <- function(varname, varcat, varlocation){#Occupation3$Cert
          Occupation3$Cert <<- safe.ifelse(is.na(varname), Occupation3$Cert,
                                          safe.ifelse(varname==Occupation3$Date, varcat, Occupation3$Cert))
          #Occupation3$Location
          Occupation3$Location <<- safe.ifelse(is.na(varname), Occupation3$Location,
                                              safe.ifelse(varname==Occupation3$Date, varlocation, Occupation3$Location))
      }
     #run function
      Occupation3$Cert <- NA ; Occupation3$Location <- NA
      assign_cert_loc(Occupation3$B_date_child, paste0("B", Occupation3$child), Occupation3$B_location_child)
      assign_cert_loc(Occupation3$D_date_child, paste0("D", Occupation3$child), Occupation3$D_location_child)
      #within loop for marriage certificates
      x <- 1
      end <- length(which(grepl("M_date", colnames(Pedigree))))
      repeat{
        #filter
        assign_cert_loc(Occupation3[[paste0("M_date_", x, "_child")]], paste0("M",Occupation3$child), Occupation3[[paste0("M_location_", x, "_child")]])
        if (x==end){
          break}
        x <- x+1
      }
      rm(x,end)
    #bind information from certificates of ego + children
      Occupation3 <- Occupation3[which(!is.na(Occupation3$Cert)),c("Id_person", "Occupation", "Date", "Age", "Location", "Cert")]
      Occupation3 <- Occupation3[which(!duplicated(Occupation3[,c("Id_person", "Date")])),]
      Occupation <- rbind(Occupation, Occupation3)
    #add certificates for which the source is unclear
      Occupation2$Location <- NA; Occupation2$Cert <- NA
      Occupation <- rbind(Occupation, Occupation2)
      Occupation <- Occupation[which(!duplicated(Occupation[,c("Id_person", "Date")])),]
      rm(Occupation2, Occupation3)
      
  #9e. add HISCO, HISCAM, HISCLASS, SOCPO
    #set function
      add_occ_scheme <- function(scheme){
        df <- as.data.frame(Individual2[which(Individual2$Type==paste0("OCCUPATION_", scheme)),c("Id_I", "Value", "Day", "Month", "Year")])
        colnames(df)[1:2] <- c("Id_person", scheme)
       #add date 
        df$Date <- paste(df$Year, df$Month, df$Day, sep="-")
        df$Date <- ymd(df$Date)
        df$Year <- NULL; df$Month <- NULL; df$Day <- NULL
       #merge
        merge(Occupation, df, by=c("Id_person", "Date"), all.x=T)
      }
    #run function
      Occupation <- add_occ_scheme("HISCO")
      Occupation <- add_occ_scheme("HISCAM_U1")
      Occupation <- add_occ_scheme("HISCLASS")
      Occupation <- add_occ_scheme("SOCPO")
      Occupation <- Occupation[which(!duplicated(Occupation[,c("Id_person", "Date")])),]
      
  #9f. recode occupational data
    #HISCO
     #Recode 0: not yet coded
      Occupation$HISCO[Occupation$HISCO==0] <- NA
     #Recode -1: occupational status
      Occupation$HISCO[Occupation$HISCO=="-1"] <- "Occupational status, classifiable"
     #Recode -2: zonder beroep
      Occupation$HISCO[Occupation$HISCO=="-2"] <- "Without occupation / unemployed"
     #Recode -3: no occupational information
      Occupation$HISCO[Occupation$HISCO=="-3"] <- "Entry error / Unintelligible / Unclassifiable"
     #Recode 99999: unintelligible or not classifiable
      Occupation$HISCO[Occupation$HISCO=="-3"] <- "Entry error / Unintelligible / Unclassifiable"
    #HISCLASS, HISCAM, SOCPO
      Occupation$HISCLASS[Occupation$HISCLASS<=0] <- NA
      Occupation$HISCAM[Occupation$HISCAM<=0] <- NA 
      Occupation$SOCPO[Occupation$SOCPO<=0] <- NA
      Occupation$SOCPO[Occupation$SOCPO==999] <- NA
      
  #9g. add OCC1950, OCCSCORE, NPBOSS, PRESGL, SEI
      US <- read.csv("HISCO - US SES.txt", header=T, quote="", sep="\t", stringsAsFactors=F)
      Occupation <- merge(Occupation, US, by="HISCO", all.x=T)
      rm(US)
      
      
      
  ###############################
  ###Step 10. Export data file###
  ###############################
  
  #This section orders and saves the Pedigree & Occupational data frames
  # 10a. Order Pedigree
  # 10b. Order Occupation
  # 10c. Save Pedigree
  # 10d. Save Occupation
      
  #Backup
      Backup7 <- Pedigree
      
  #10a. Order Pedigree
      Pedigree <- Pedigree[,c("Id_person", "Id_mother", "Id_father",
                              "Sex", "Twin",
                              "LastEntryDate", "LastEntryAge", "LastEntryCert",
                              "B", "M","D", "Postlink_D", "M_parents", "B_mother", "B_father", "D_mother", "D_father",
                              "B_date", "B_max_date", "B_min_date", "B_location", "B_age_mother", "B_age_father",
                              "D_age", "D_date", "D_location", "D_deadonregistration",
                              "Marriage_N",
                              "Id_partner_1", "M_age_1", "M_date_1", "M_location_1", 
                              "Id_partner_2", "M_age_2", "M_date_2", "M_location_2", 
                              "Id_partner_3", "M_age_3", "M_date_3", "M_location_3", 
                              "Id_partner_4", "M_age_4", "M_date_4", "M_location_4", 
                              "Id_partner_5", "M_age_5", "M_date_5", "M_location_5")]
  
  #10b. Order Occupation
      Occupation <- Occupation[,c("Id_person", "Date", "Age", "Location", "Cert",
                                  "Occupation", 
                                  "HISCO", "HISCLASS", "SOCPO", "HISCAM_U1",
                                  "OCC1950", "OCCSCORE", "PRESGL", "SEI", "NPBOSS")]
      
  #10c. Save Pedigree
      write.table(Pedigree, file="Pedigree.txt", quote=F, sep ="\t", col.names=T, row.names=F)
      
  #10d. Save Occupations
      write.table(Occupation, file="Occupation.txt", quote=F, sep ="\t", col.names=T, row.names=F)
      
  
  