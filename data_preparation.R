
######################################################################################
#                                                                                     
#   Filename        : data_preparation.R	   												  
#                                                                                     
#   Project         : Consultation Arp Schnittger and Joke De Jaeger-Breat, Meiosis / Heat stress analysis                                                                       
#   Author          : Linda Krause                                                                    
#   Purpose         : Data preparation (from the raw data Excel file to Stata file which can be used in survival analysis)                                                                       
#   Date            : 22.06.2021                                                                       
#   Program Version : V1											  
#   R Version       : 3.5.3                                                                       
#
#
#   Input data files  : RawDataImaging_28022020_SPLIT_220920_revised091220_LK.xlsx and RawData_mutants_completed_FL_080121_LK.xlsx                                                               
#   Output data files : sum_values.dta and sum_values_mutations.dta                                                               
#                                                                                     
#     
#   Required packages: xlsx
#
#   Associated Programs:                                                              
#    - Include  :                                                          
#    - Run the present program before  : Statistical analysis by Anika Buchholz                                                                 
#    - Run the present program after   :                                                     
#                                                                                     
#######################################################################################


# to avoid java heap space error
options(java.parameters = "-Xmx8000m")
options(stringsAsFactors = F)



# Packages
library(xlsx)



# set path to raw data in Excel formal
namedat1 = "../Data/RawDataImaging_28022020_SPLIT_220920_revised091220_LK.xlsx"
namedat2 = "../Data/RawData_mutants_completed_FL_080121_LK.xlsx"




# Function to calculate time in each stage for each cell

calculateTimeinEachStage = function(ord_values, treatment){
   # for each column in ord_values we create a row in sum_values with the time in this stage
   # each stage has 3 columns: one with the time spend in the stage and then 2 describing whether the time ended or started with NA
   # there are 7 stages of interest (+ all stages together + 3 new timeranges), so 34 columns needed (one ID column)
   
   # prepare resulting df
   sum_values = data.frame(matrix(NA, ncol = 25, nrow = dim(ord_values)[2]-1))
   # adjust colnames
   temp = c("P", "H", "F", "M", "I", "M2", "T", "PT", "HT", "FT", "HM2", "MM2")
   temp2 = c("time_l", "time_u")
   colnames(sum_values) = c("ID", paste(rep(temp, each = length(temp2)), temp2, sep = "_"))
   # save ID
   sum_values$ID = colnames(ord_values)[1:(dim(ord_values)[2]-1)]
   
   # for loop through all cells
   for (i in 1:(dim(ord_values)[2]-1)){
      
      # get measurements for that cell
      x = as.character(ord_values[,i])
      
      # replace M2 and Mi since it cannot be differentiated from M when counting in concatenated string
      # M2 becomes R
      idx = which(x == "M2")
      x[idx] = "R"
      # Mi becomes Q
      idx = which(x == "Mi")
      x[idx] = "Q"
      
      # paste measurement into one string
      y = paste(as.character(x), collapse = "")
      
      
      ### go through all stages P-H-F-M-I-M2-T, get lower and upper end of interval
      
      
      
      ## SEVERAL STAGES
      
      ### H -- M2 ( = R)
      
      # positions of first and last stage in character vector x
      temp1 = which(x == "H")
      temp2 = which(x == "R")
      
      # only if both stages are present: save otherwise leave as NA
      if (length(temp1) > 0 & length(temp2) > 0){
         
         # duration is max of last stage minus min of first stage plus 1
         count = max(temp2) - min(temp1) + 1
         
         # check whether beginning and end is known (not censored) >> time_l und time_u are the same
         if (length(grep("PH", y)) > 0 & length(grep("RT", y)) > 0){
            sum_values$HM2_time_l[i] = count
            sum_values$HM2_time_u[i] = count
         }
         
         # check whether stage before or next stage is not present in string >> right censored data >> time_l time and time_u is NA
         # positions of stage before and next stage in character vector x
         temp3 = which(x == "P")
         temp4 = which(x == "T")
         if (length(temp3) == 0 | length(temp4) == 0){
            sum_values$HM2_time_l[i] = count
            sum_values$HM2_time_u[i] = NA  
         }else{
            # interval censored
            # time_l is number of characters of this stage in the character string
            sum_values$HM2_time_l[i] = count
            # time_u is difference between index of next stage 1st character (min) the last character auf stage before (max) -1
            sum_values$HM2_time_u[i] = min(temp4) - max(temp3) - 1
         }
      }
      
      
      ### H -- T 
      
      # positions of first and last stage in character vector x
      temp1 = which(x == "H")
      temp2 = which(x == "T")
      
      # only if both stages are present: save otherwise leave as NA
      if (length(temp1) > 0 & length(temp2) > 0){
         
         # duration is max of last stage minus min of first stage plus 1
         count = max(temp2) - min(temp1) + 1
         
         # check whether beginning and end is known (not censored) >> time_l und time_u are the same
         if (length(grep("PH", y)) > 0 & length(grep("TQ", y)) > 0){
            sum_values$HT_time_l[i] = count
            sum_values$HT_time_u[i] = count
         }
         
         # check whether stage before or next stage is not present in string >> right censored data >> time_l time and time_u is NA
         # positions of stage before and next stage in character vector x
         temp3 = which(x == "P")
         temp4 = which(x == "Q")
         if (length(temp3) == 0 | length(temp4) == 0){
            sum_values$HT_time_l[i] = count
            sum_values$HT_time_u[i] = NA  
         }else{
            # interval censored
            # time_l is number of characters of this stage in the character string
            sum_values$HT_time_l[i] = count
            # time_u is difference between index of next stage 1st character (min) the last character auf stage before (max) -1
            sum_values$HT_time_u[i] = min(temp4) - max(temp3) - 1
         }
      }
      
      
      ### F -- T
      
      # positions of first and last stage in character vector x
      temp1 = which(x == "F")
      temp2 = which(x == "T")
      
      # only if both stages are present: save otherwise leave as NA
      if (length(temp1) > 0 & length(temp2) > 0){
         
         # duration is max of last stage minus min of first stage plus 1
         count = max(temp2) - min(temp1) + 1
         
         # check whether beginning and end is known (not censored) >> time_l und time_u are the same
         if (length(grep("HF", y)) > 0 & length(grep("TQ", y)) > 0){
            sum_values$FT_time_l[i] = count
            sum_values$FT_time_u[i] = count
         }
         
         # check whether stage before or next stage is not present in string >> right censored data >> time_l time and time_u is NA
         # positions of stage before and next stage in character vector x
         temp3 = which(x == "H")
         temp4 = which(x == "Q")
         if (length(temp3) == 0 | length(temp4) == 0){
            sum_values$FT_time_l[i] = count
            sum_values$FT_time_u[i] = NA  
         }else{
            # interval censored
            # time_l is number of characters of this stage in the character string
            sum_values$FT_time_l[i] = count
            # time_u is difference between index of next stage 1st character (min) the last character auf stage before (max) -1
            sum_values$FT_time_u[i] = min(temp4) - max(temp3) - 1
         }
      }
      
      
      
      ### M -- M2 (= R)
      
      # positions of first and last stage in character vector x
      temp1 = which(x == "M")
      temp2 = which(x == "R")
      
      # only if both stages are present: save otherwise leave as NA
      if (length(temp1) > 0 & length(temp2) > 0){
         
         # duration is max of last stage minus min of first stage plus 1
         count = max(temp2) - min(temp1) + 1
         
         # check whether beginning and end is known (not censored) >> time_l und time_u are the same
         if (length(grep("FM", y)) > 0 & length(grep("RT", y)) > 0){
            sum_values$MM2_time_l[i] = count
            sum_values$MM2_time_u[i] = count
         }
         
         # check whether stage before or next stage is not present in string >> right censored data >> time_l time and time_u is NA
         # positions of stage before and next stage in character vector x
         temp3 = which(x == "F")
         temp4 = which(x == "T")
         if (length(temp3) == 0 | length(temp4) == 0){
            sum_values$MM2_time_l[i] = count
            sum_values$MM2_time_u[i] = NA  
         }else{
            # interval censored
            # time_l is number of characters of this stage in the character string
            sum_values$MM2_time_l[i] = count
            # time_u is difference between index of next stage 1st character (min) the last character auf stage before (max) -1
            sum_values$MM2_time_u[i] = min(temp4) - max(temp3) - 1
         }
      }
      
      
      ### P -- T
      
      # positions of first and last stage in character vector x
      temp1 = which(x == "P")
      temp2 = which(x == "T")
      
      # only if both stages are present: save otherwise leave as NA
      if (length(temp1) > 0 & length(temp2) > 0){
         
         # duration is max of last stage minus min of first stage plus 1
         count = max(temp2) - min(temp1) + 1
         
         # check whether beginning and end is known (not censored) >> time_l und time_u are the same
         if (length(grep("XP", y)) > 0 & length(grep("TQ", y)) > 0){
            sum_values$PT_time_l[i] = count
            sum_values$PT_time_u[i] = count
         }
         
         # check whether stage before or next stage is not present in string >> right censored data >> time_l time and time_u is NA
         # positions of stage before and next stage in character vector x
         temp3 = which(x == "X")
         temp4 = which(x == "Q")
         if (length(temp3) == 0 | length(temp4) == 0){
            sum_values$PT_time_l[i] = count
            sum_values$PT_time_u[i] = NA  
         }else{
            # interval censored
            # time_l is number of characters of this stage in the character string
            sum_values$PT_time_l[i] = count
            # time_u is difference between index of next stage 1st character (min) the last character auf stage before (max) -1
            sum_values$PT_time_u[i] = min(temp4) - max(temp3) - 1
         }
      }
      
      
      
      #### ONE STAGE AT A TIME #####
      
      
      ### stage P
      
      # get number of characters of this stage in the string 
      count = stringr::str_count(y, "P")
      
      # only if stage is present: save the number otherwise leave as NA
      if (count > 0){
         
         # check whether beginning and end is known (not censored) >> time_l und time_u are the same
         if (length(grep("XP", y)) > 0 & length(grep("PH", y)) > 0){
            sum_values$P_time_l[i] = count
            sum_values$P_time_u[i] = count
         }
         
         # check whether stage before or next stage is not present in string >> right censored data >> time_l time and time_u is NA
         # positions of stage before and next stage in character vector x
         temp1 = which(x == "X")
         temp2 = which(x == "H")
         if (length(temp1) == 0 | length(temp2) == 0){
            sum_values$P_time_l[i] = count
            sum_values$P_time_u[i] = NA  
         }else{
            # interval censored
            # time_l is number of characters of this stage in the character string
            sum_values$P_time_l[i] = count
            # time_u is difference between index of next stage 1st character (min) the last character auf stage before (max) -1
            sum_values$P_time_u[i] = min(temp2) - max(temp1) - 1
         }
      }
      
      
      ### stage H
      
      # get number of characters of this stage in the string 
      count = stringr::str_count(y, "H")
      
      # only if stage is present: save the number otherwise leave as NA
      if (count > 0){
         
         # check whether beginning and end is known (not censored) >> time_l und time_u are the same
         if (length(grep("PH", y)) > 0 & length(grep("HF", y)) > 0){
            sum_values$H_time_l[i] = count
            sum_values$H_time_u[i] = count
         }
         
         # check whether stage before or next stage is not present in string >> right censored data >> time_l time and time_u is NA
         # positions of stage before and next stage in character vector x
         temp1 = which(x == "P")
         temp2 = which(x == "F")
         if (length(temp1) == 0 | length(temp2) == 0){
            sum_values$H_time_l[i] = count
            sum_values$H_time_u[i] = NA  
         }else{
            # interval censored
            # time_l is number of characters of this stage in the character string
            sum_values$H_time_l[i] = count
            # time_u is difference between index of next stage 1st character (min) the last character auf stage before (max) -1
            sum_values$H_time_u[i] = min(temp2) - max(temp1) - 1
         }
      }
      
      
      
      ### stage F
      
      # get number of characters of this stage in the string 
      count = stringr::str_count(y, "F")
      
      # only if stage is present: save the number otherwise leave as NA
      if (count > 0){
         
         # check whether beginning and end is known (not censored) >> time_l und time_u are the same
         if (length(grep("HF", y)) > 0 & length(grep("FM", y)) > 0){
            sum_values$F_time_l[i] = count
            sum_values$F_time_u[i] = count
         }
         
         # check whether stage before or next stage is not present in string >> right censored data >> time_l time and time_u is NA
         # positions of stage before and next stage in character vector x
         temp1 = which(x == "H")
         temp2 = which(x == "M")
         if (length(temp1) == 0 | length(temp2) == 0){
            sum_values$F_time_l[i] = count
            sum_values$F_time_u[i] = NA  
         }else{
            # interval censored
            # time_l is number of characters of this stage in the character string
            sum_values$F_time_l[i] = count
            # time_u is difference between index of next stage 1st character (min) the last character auf stage before (max) -1
            sum_values$F_time_u[i] = min(temp2) - max(temp1) - 1
         }
      }
      
      
      
      ### stage M
      
      # get number of characters of this stage in the string 
      count = stringr::str_count(y, "M")
      
      # only if stage is present: save the number otherwise leave as NA
      if (count > 0){
         
         # check whether beginning and end is known (not censored) >> time_l und time_u are the same
         if (length(grep("FM", y)) > 0 & length(grep("MI", y)) > 0){
            sum_values$M_time_l[i] = count
            sum_values$M_time_u[i] = count
         }
         
         # check whether stage before or next stage is not present in string >> right censored data >> time_l time and time_u is NA
         # positions of stage before and next stage in character vector x
         temp1 = which(x == "F")
         temp2 = which(x == "I")
         if (length(temp1) == 0 | length(temp2) == 0){
            sum_values$M_time_l[i] = count
            sum_values$M_time_u[i] = NA  
         }else{
            # interval censored
            # time_l is number of characters of this stage in the character string
            sum_values$M_time_l[i] = count
            # time_u is difference between index of next stage 1st character (min) the last character auf stage before (max) -1
            sum_values$M_time_u[i] = min(temp2) - max(temp1) - 1
         }
      }
      
      
      ### stage I
      
      # get number of characters of this stage in the string 
      count = stringr::str_count(y, "I")
      
      # only if stage is present: save the number otherwise leave as NA
      if (count > 0){
         
         # check whether beginning and end is known (not censored) >> time_l und time_u are the same
         if (length(grep("MI", y)) > 0 & length(grep("IR", y)) > 0){
            sum_values$I_time_l[i] = count
            sum_values$I_time_u[i] = count
         }
         
         # check whether stage before or next stage is not present in string >> right censored data >> time_l time and time_u is NA
         # positions of stage before and next stage in character vector x
         temp1 = which(x == "M")
         temp2 = which(x == "R")
         if (length(temp1) == 0 | length(temp2) == 0){
            sum_values$I_time_l[i] = count
            sum_values$I_time_u[i] = NA  
         }else{
            # interval censored
            # time_l is number of characters of this stage in the character string
            sum_values$I_time_l[i] = count
            # time_u is difference between index of next stage 1st character (min) the last character auf stage before (max) -1
            sum_values$I_time_u[i] = min(temp2) - max(temp1) - 1
         }
      }
      
      
      ### stage M2 == R
      
      # get number of characters of this stage in the string 
      count = stringr::str_count(y, "R")
      
      # only if stage is present: save the number otherwise leave as NA
      if (count > 0){
         
         # check whether beginning and end is known (not censored) >> time_l und time_u are the same
         if (length(grep("IR", y)) > 0 & length(grep("RT", y)) > 0){
            sum_values$M2_time_l[i] = count
            sum_values$M2_time_u[i] = count
         }
         
         # check whether stage before or next stage is not present in string >> right censored data >> time_l time and time_u is NA
         # positions of stage before and next stage in character vector x
         temp1 = which(x == "I")
         temp2 = which(x == "T")
         if (length(temp1) == 0 | length(temp2) == 0){
            sum_values$M2_time_l[i] = count
            sum_values$M2_time_u[i] = NA  
         }else{
            # interval censored
            # time_l is number of characters of this stage in the character string
            sum_values$M2_time_l[i] = count
            # time_u is difference between index of next stage 1st character (min) the last character auf stage before (max) -1
            sum_values$M2_time_u[i] = min(temp2) - max(temp1) - 1
         }
      }
      
      
      ### stage T
      
      # get number of characters of this stage in the string 
      count = stringr::str_count(y, "T")
      
      # only if stage is present: save the number otherwise leave as NA
      if (count > 0){
         
         # check whether beginning and end is known (not censored) >> time_l und time_u are the same
         if (length(grep("RT", y)) > 0 & length(grep("TQ", y)) > 0){
            sum_values$T_time_l[i] = count
            sum_values$T_time_u[i] = count
         }
         
         # check whether stage before or next stage is not present in string >> right censored data >> time_l time and time_u is NA
         # positions of stage before and next stage in character vector x
         temp1 = which(x == "R")
         temp2 = which(x == "Q")
         if (length(temp1) == 0 | length(temp2) == 0){
            sum_values$T_time_l[i] = count
            sum_values$T_time_u[i] = NA  
         }else{
            # interval censored
            # time_l is number of characters of this stage in the character string
            sum_values$T_time_l[i] = count
            # time_u is difference between index of next stage 1st character (min) the last character auf stage before (max) -1
            sum_values$T_time_u[i] = min(temp2) - max(temp1) - 1
         }
      }
      
   }
   # get anther ID using regular expression (stuff between both "_")
   pattern = "_\\s*(.*?)\\s*_"
   result = regmatches(sum_values$ID, regexec(pattern, sum_values$ID))
   
   # add anther
   sum_values$anther = unlist(lapply(result, function(x){x[2]}))
   
   # add treatment
   sum_values$treatment = treatment
   
   # rearrage order of columns
   sum_values = sum_values[, c(1, 26, 27, 2:25)]
   
   # return result
   return(sum_values)
}



# Wild-type data

## Control

### load data
dat = read.xlsx(namedat1, sheetIndex = 1)

# check and remove complete NA row
idx = which(apply(is.na(dat), 1, sum) == dim(dat)[2])
if (length(idx) > 0){
   dat = dat[-idx, ]
}

### add rowID
dat$uID = paste0("ID_", dat$ID, "_", dat$meiocyte.nr)
# reorder to put uID in the 7th column
dat = dat[, c(1:6, dim(dat)[2], 7:(dim(dat)[2]-1))]


# split measurements from information
info = dat[, 1:7]
values = dat[, 5:dim(dat)[2]]


### fill "n" if possible
# if state before and after "n" is the same, fill it with this value


# go through all rows
for (i in 1:dim(values)[1]){
   # go through all columns with measurements
   for (j in 4:dim(values)[2]){
      if (!is.na(values[i, j])){
         # check if it is an "n"
         if (values[i, j] == "n"){
            idx = match(values[i, j-1], values[i, j:dim(values)[2]])
            # if the index is not NA, then a match was found
            if (!is.na(idx)){
               values[i, j] = values[i, j-1]  
            }
         } 
      }
   }
}


### change C and S to P
idx = which(values == "C", arr.ind = T)
values[idx] = "P"
idx = which(values == "S", arr.ind = T)
values[idx] = "P"



### set to ordered factor
# All "n" are automatically set to NA if "n" is none of the factor levels.
temp = as.data.frame(t(values[, 4:dim(values)[2]]))
colnames(temp) = values$uID
ord_values = data.frame(lapply(temp, function(x){factor(x, ordered = T, levels = c("X", "P", "H", "F", "M", "I", "M2", "T", "Mi"))}))



### calculate time in each stage for each cell
sum_values_control = calculateTimeinEachStage(ord_values, "control")




## Treatment 1 week 30 degress

### load data
dat = read.xlsx(namedat1, sheetIndex = 4)

# check and remove complete NA row
idx = which(apply(is.na(dat), 1, sum) == dim(dat)[2])
if (length(idx) > 0){
   dat = dat[-idx, ]
}


### add rowID
dat$uID = paste0("ID_", dat$ID, "_", dat$meiocyte.nr)
# reorder to put uID in the 7th column
dat = dat[, c(1:6, dim(dat)[2], 7:(dim(dat)[2]-1))]


# split measurements from information
info = dat[, 1:7]
values = dat[, 5:dim(dat)[2]]


### fill "n" if possible
# if state before and after "n" is the same, fill it with this value
# go through all rows
for (i in 1:dim(values)[1]){
   # go through all columns with measurements
   for (j in 4:dim(values)[2]){
      if (!is.na(values[i, j])){
         # check if it is an "n"
         if (values[i, j] == "n"){
            idx = match(values[i, j-1], values[i, j:dim(values)[2]])
            # if the index is not NA, then a match was found
            if (!is.na(idx)){
               values[i, j] = values[i, j-1]  
            }
         } 
      }
   }
}


### change C and S to P
idx = which(values == "C", arr.ind = T)
values[idx] = "P"
idx = which(values == "S", arr.ind = T)
values[idx] = "P"



### set to ordered factor
# All "n" are automatically set to NA if "n" is none of the factor levels.
temp = as.data.frame(t(values[, 4:dim(values)[2]]))
colnames(temp) = values$uID
ord_values = data.frame(lapply(temp, function(x){factor(x, ordered = T, levels = c("X", "P", "H", "F", "M", "I", "M2", "T", "Mi"))}))



### calculate time in each stage for each cell
sum_values_oneweek30 = calculateTimeinEachStage(ord_values, "oneweek30")




## Treatment heatshock 30 degrees - stimulation in time A

### load data
dat = read.xlsx(namedat1, sheetIndex = 2)

# check and remove complete NA row
idx = which(apply(is.na(dat), 1, sum) == dim(dat)[2])
if (length(idx) > 0){
   dat = dat[-idx, ]
}


# get early (indicated as A)
idx = which(dat$early_late == "A")
idx.c = which(colnames(dat) == "early_late")
dat = dat[idx, -7]


### add rowID
dat$uID = paste0("ID_", dat$ID, "_", dat$meiocyte.nr)
# reorder to put uID in the 7th column
dat = dat[, c(1:6, dim(dat)[2], 7:(dim(dat)[2]-1))]


# split measurements from information
info = dat[, 1:7]
values = dat[, 5:dim(dat)[2]]


### fill "n" if possible
# if state before and after "n" is the same, fill it with this value
# go through all rows
for (i in 1:dim(values)[1]){
   # go through all columns with measurements
   for (j in 4:dim(values)[2]){
      if (!is.na(values[i, j])){
         # check if it is an "n"
         if (values[i, j] == "n"){
            idx = match(values[i, j-1], values[i, j:dim(values)[2]])
            # if the index is not NA, then a match was found
            if (!is.na(idx)){
               values[i, j] = values[i, j-1]  
            }
         } 
      }
   }
}


### change C and S to P
idx = which(values == "C", arr.ind = T)
values[idx] = "P"
idx = which(values == "S", arr.ind = T)
values[idx] = "P"



### set to ordered factor
# All "n" are automatically set to NA if "n" is none of the factor levels.
temp = as.data.frame(t(values[, 4:dim(values)[2]]))
colnames(temp) = values$uID
ord_values = data.frame(lapply(temp, function(x){factor(x, ordered = T, levels = c("X", "P", "H", "F", "M", "I", "M2", "T", "Mi"))}))



### calculate time in each stage for each cell
sum_values_HS30_A = calculateTimeinEachStage(ord_values, "HS30_A")





## Treatment heatshock 30 degrees - stimulation in time B

### load data
dat = read.xlsx(namedat1, sheetIndex = 2)

# check and remove complete NA row
idx = which(apply(is.na(dat), 1, sum) == dim(dat)[2])
if (length(idx) > 0){
   dat = dat[-idx, ]
}

# get late (indicated as B)
idx = which(dat$early_late == "B")
idx.c = which(colnames(dat) == "early_late")
dat = dat[idx, -7]

### add rowID
dat$uID = paste0("ID_", dat$ID, "_", dat$meiocyte.nr)
# reorder to put uID in the 7th column
dat = dat[, c(1:6, dim(dat)[2], 7:(dim(dat)[2]-1))]


# split measurements from information
info = dat[, 1:7]
values = dat[, 5:dim(dat)[2]]


### fill "n" if possible
# if state before and after "n" is the same, fill it with this value
# go through all rows
for (i in 1:dim(values)[1]){
   # go through all columns with measurements
   for (j in 4:dim(values)[2]){
      if (!is.na(values[i, j])){
         # check if it is an "n"
         if (values[i, j] == "n"){
            idx = match(values[i, j-1], values[i, j:dim(values)[2]])
            # if the index is not NA, then a match was found
            if (!is.na(idx)){
               values[i, j] = values[i, j-1]  
            }
         } 
      }
   }
}


### change C and S to P
idx = which(values == "C", arr.ind = T)
values[idx] = "P"
idx = which(values == "S", arr.ind = T)
values[idx] = "P"



### set to ordered factor
# All "n" are automatically set to NA if "n" is none of the factor levels.
temp = as.data.frame(t(values[, 4:dim(values)[2]]))
colnames(temp) = values$uID
ord_values = data.frame(lapply(temp, function(x){factor(x, ordered = T, levels = c("X", "P", "H", "F", "M", "I", "M2", "T", "Mi"))}))



### calculate time in each stage for each cell
sum_values_HS30_B = calculateTimeinEachStage(ord_values, "HS30_B")




## Treatment heatshock 34 degrees - stimulation in time A

### load data
dat = read.xlsx(namedat1, sheetIndex = 3)


# check and remove complete NA row
idx = which(apply(is.na(dat), 1, sum) == dim(dat)[2])
if (length(idx) > 0){
   dat = dat[-idx, ]
}


# get early (indicated as A)
idx = which(dat$early_late == "A")
idx.c = which(colnames(dat) == "early_late")
dat = dat[idx, -7]



### add rowID
dat$uID = paste0("ID_", dat$ID, "_", dat$meiocyte.nr)
# reorder to put uID in the 7th column
dat = dat[, c(1:6, dim(dat)[2], 7:(dim(dat)[2]-1))]


# split measurements from information
info = dat[, 1:7]
values = dat[, 5:dim(dat)[2]]


### fill "n" if possible
# if state before and after "n" is the same, fill it with this value

# go through all rows
for (i in 1:dim(values)[1]){
   # go through all columns with measurements
   for (j in 4:dim(values)[2]){
      if (!is.na(values[i, j])){
         # check if it is an "n"
         if (values[i, j] == "n"){
            idx = match(values[i, j-1], values[i, j:dim(values)[2]])
            # if the index is not NA, then a match was found
            if (!is.na(idx)){
               values[i, j] = values[i, j-1]  
            }
         } 
      }
   }
}


### change C and S to P
idx = which(values == "C", arr.ind = T)
values[idx] = "P"
idx = which(values == "S", arr.ind = T)
values[idx] = "P"



### set to ordered factor
# All "n" are automatically set to NA if "n" is none of the factor levels.
temp = as.data.frame(t(values[, 4:dim(values)[2]]))
colnames(temp) = values$uID
ord_values = data.frame(lapply(temp, function(x){factor(x, ordered = T, levels = c("X", "P", "H", "F", "M", "I", "M2", "T", "Mi"))}))



### calculate time in each stage for each cell
sum_values_HS34_A = calculateTimeinEachStage(ord_values, "HS34_A")




## Treatment heatshock 34 degrees - stimulation in time B

### load data
dat = read.xlsx(namedat1, sheetIndex = 3)


# check and remove complete NA row
idx = which(apply(is.na(dat), 1, sum) == dim(dat)[2])
if (length(idx) > 0){
   dat = dat[-idx, ]
}


# get late (indicated as B)
idx = which(dat$early_late == "B")
idx.c = which(colnames(dat) == "early_late")
dat = dat[idx, -7]



### add rowID
dat$uID = paste0("ID_", dat$ID, "_", dat$meiocyte.nr)
# reorder to put uID in the 7th column
dat = dat[, c(1:6, dim(dat)[2], 7:(dim(dat)[2]-1))]


# split measurements from information
info = dat[, 1:7]
values = dat[, 5:dim(dat)[2]]


### fill "n" if possible
# if state before and after "n" is the same, fill it with this value
# go through all rows
for (i in 1:dim(values)[1]){
   # go through all columns with measurements
   for (j in 4:dim(values)[2]){
      if (!is.na(values[i, j])){
         # check if it is an "n"
         if (values[i, j] == "n"){
            idx = match(values[i, j-1], values[i, j:dim(values)[2]])
            # if the index is not NA, then a match was found
            if (!is.na(idx)){
               values[i, j] = values[i, j-1]  
            }
         } 
      }
   }
}


### change C and S to P
idx = which(values == "C", arr.ind = T)
values[idx] = "P"
idx = which(values == "S", arr.ind = T)
values[idx] = "P"



### set to ordered factor
# All "n" are automatically set to NA if "n" is none of the factor levels.
temp = as.data.frame(t(values[, 4:dim(values)[2]]))
colnames(temp) = values$uID
ord_values = data.frame(lapply(temp, function(x){factor(x, ordered = T, levels = c("X", "P", "H", "F", "M", "I", "M2", "T", "Mi"))}))



### calculate time in each stage for each cell
sum_values_HS34_B = calculateTimeinEachStage(ord_values, "HS34_B")



## combine all sum values in one df and export for stata and R


df = rbind(sum_values_control, sum_values_HS30_A, sum_values_HS30_B, sum_values_HS34_A, sum_values_HS34_B, sum_values_oneweek30)
# export for stata
haven::write_dta(df, path = "sum_values.dta")
# export for R
saveRDS(df, file = "sum_values.Rds")








# Mutation data

## DMC1 control

### load data
dat = read.xlsx(namedat2, sheetIndex = 1)

# check and remove complete NA row
idx = which(apply(is.na(dat), 1, sum) == dim(dat)[2])
if (length(idx) > 0){
   dat = dat[-idx, ]
}



### add rowID
dat$uID = paste0("ID_", dat$ID, "_", dat$cells)
# reorder to put uID in the 6th column
dat = dat[, c(1:5, dim(dat)[2], 6:(dim(dat)[2]-1))]


# split measurements from information
info = dat[, 1:6]
values = dat[, 6:dim(dat)[2]]


### fill "n" if possible
#if state before and after "n" is the same, fill it with this value

# go through all rows
for (i in 1:dim(values)[1]){
   # go through all columns with measurements
   for (j in 2:dim(values)[2]){
      if (!is.na(values[i, j])){
         # check if it is an "n"
         if (values[i, j] == "n"){
            idx = match(values[i, j-1], values[i, j:dim(values)[2]])
            # if the index is not NA, then a match was found
            if (!is.na(idx)){
               values[i, j] = values[i, j-1]  
            }
         } 
      }
   }
}


### change C and S to P
idx = which(values == "C", arr.ind = T)
values[idx] = "P"
idx = which(values == "S", arr.ind = T)
values[idx] = "P"



### set to ordered factor
# All "n" are automatically set to NA if "n" is none of the factor levels.
temp = as.data.frame(t(values[, 2:dim(values)[2]]))
colnames(temp) = values$uID
ord_values = data.frame(lapply(temp, function(x){factor(x, ordered = T, levels = c("X", "P", "H", "F", "M", "I", "M2", "T", "Mi"))}))



### calculate time in each stage for each cell
sum_values_dmc1_control = calculateTimeinEachStage(ord_values, "DMC_control")




## DMC1 HS34

### load data
dat = read.xlsx(namedat2, sheetIndex = 2)


# check and remove complete NA row
idx = which(apply(is.na(dat), 1, sum) == dim(dat)[2])
if (length(idx) > 0){
   dat = dat[-idx, ]
}



### add rowID
dat$uID = paste0("ID_", dat$ID, "_", dat$cells)
# reorder to put uID in the 6th column
dat = dat[, c(1:5, dim(dat)[2], 6:(dim(dat)[2]-1))]


# split measurements from information
info = dat[, 1:6]
values = dat[, 6:dim(dat)[2]]


### fill "n" if possible
# if state before and after "n" is the same, fill it with this value

# go through all rows
for (i in 1:dim(values)[1]){
   # go through all columns with measurements
   for (j in 2:dim(values)[2]){
      if (!is.na(values[i, j])){
         # check if it is an "n"
         if (values[i, j] == "n"){
            idx = match(values[i, j-1], values[i, j:dim(values)[2]])
            # if the index is not NA, then a match was found
            if (!is.na(idx)){
               values[i, j] = values[i, j-1]  
            }
         } 
      }
   }
}


### change C and S to P
idx = which(values == "C", arr.ind = T)
values[idx] = "P"
idx = which(values == "S", arr.ind = T)
values[idx] = "P"



### set to ordered factor
# All "n" are automatically set to NA if "n" is none of the factor levels.
temp = as.data.frame(t(values[, 2:dim(values)[2]]))
colnames(temp) = values$uID
ord_values = data.frame(lapply(temp, function(x){factor(x, ordered = T, levels = c("X", "P", "H", "F", "M", "I", "M2", "T", "Mi"))}))




### calculate time in each stage for each cell
sum_values_dmc1_HS34 = calculateTimeinEachStage(ord_values, "DMC_HS34")



## SPO11 control

### load data
dat = read.xlsx(namedat2, sheetIndex = 3)


# check and remove complete NA row
idx = which(apply(is.na(dat), 1, sum) == dim(dat)[2])
if (length(idx) > 0){
   dat = dat[-idx, ]
}



### add rowID
dat$uID = paste0("ID_", dat$ID, "_", dat$cells)
# reorder to put uID in the 6th column
dat = dat[, c(1:5, dim(dat)[2], 6:(dim(dat)[2]-1))]


# split measurements from information
info = dat[, 1:6]
values = dat[, 6:dim(dat)[2]]


### fill "n" if possible
# if state before and after "n" is the same, fill it with this value

# go through all rows
for (i in 1:dim(values)[1]){
   # go through all columns with measurements
   for (j in 2:dim(values)[2]){
      if (!is.na(values[i, j])){
         # check if it is an "n"
         if (values[i, j] == "n"){
            idx = match(values[i, j-1], values[i, j:dim(values)[2]])
            # if the index is not NA, then a match was found
            if (!is.na(idx)){
               values[i, j] = values[i, j-1]  
            }
         } 
      }
   }
}


### change C and S to P
idx = which(values == "C", arr.ind = T)
values[idx] = "P"
idx = which(values == "S", arr.ind = T)
values[idx] = "P"



### set to ordered factor
# All "n" are automatically set to NA if "n" is none of the factor levels.
temp = as.data.frame(t(values[, 2:dim(values)[2]]))
colnames(temp) = values$uID
ord_values = data.frame(lapply(temp, function(x){factor(x, ordered = T, levels = c("X", "P", "H", "F", "M", "I", "M2", "T", "Mi"))}))





### calculate time in each stage for each cell
sum_values_spo11_control = calculateTimeinEachStage(ord_values, "SPO11_control")




## SPO11 HS34

### load data
dat = read.xlsx(namedat2, sheetIndex = 4)


# check and remove complete NA row
idx = which(apply(is.na(dat), 1, sum) == dim(dat)[2])
if (length(idx) > 0){
   dat = dat[-idx, ]
}



### add rowID
dat$uID = paste0("ID_", dat$ID, "_", dat$cells)
# reorder to put uID in the 6th column
dat = dat[, c(1:5, dim(dat)[2], 6:(dim(dat)[2]-1))]


# split measurements from information
info = dat[, 1:6]
values = dat[, 6:dim(dat)[2]]


### fill "n" if possible
# if state before and after "n" is the same, fill it with this value


# go through all rows
for (i in 1:dim(values)[1]){
   # go through all columns with measurements
   for (j in 2:dim(values)[2]){
      if (!is.na(values[i, j])){
         # check if it is an "n"
         if (values[i, j] == "n"){
            idx = match(values[i, j-1], values[i, j:dim(values)[2]])
            # if the index is not NA, then a match was found
            if (!is.na(idx)){
               values[i, j] = values[i, j-1]  
            }
         } 
      }
   }
}


### change C and S to P
idx = which(values == "C", arr.ind = T)
values[idx] = "P"
idx = which(values == "S", arr.ind = T)
values[idx] = "P"



### set to ordered factor
# All "n" are automatically set to NA if "n" is none of the factor levels.
temp = as.data.frame(t(values[, 2:dim(values)[2]]))
colnames(temp) = values$uID
ord_values = data.frame(lapply(temp, function(x){factor(x, ordered = T, levels = c("X", "P", "H", "F", "M", "I", "M2", "T", "Mi"))}))




### calculate time in each stage for each cell
sum_values_spo11_HS34 = calculateTimeinEachStage(ord_values, "SPO11_HS34")




## MSH control

### load data
dat = read.xlsx(namedat2, sheetIndex = 5)


# check and remove complete NA row
idx = which(apply(is.na(dat), 1, sum) == dim(dat)[2])
if (length(idx) > 0){
   dat = dat[-idx, ]
}



### add rowID
dat$uID = paste0("ID_", dat$ID, "_", dat$cells)
# reorder to put uID in the 6th column
dat = dat[, c(1:5, dim(dat)[2], 6:(dim(dat)[2]-1))]


# split measurements from information
info = dat[, 1:6]
values = dat[, 6:dim(dat)[2]]


### fill "n" if possible
# if state before and after "n" is the same, fill it with this value

# go through all rows
for (i in 1:dim(values)[1]){
   # go through all columns with measurements
   for (j in 2:dim(values)[2]){
      if (!is.na(values[i, j])){
         # check if it is an "n"
         if (values[i, j] == "n"){
            idx = match(values[i, j-1], values[i, j:dim(values)[2]])
            # if the index is not NA, then a match was found
            if (!is.na(idx)){
               values[i, j] = values[i, j-1]  
            }
         } 
      }
   }
}


### change C and S to P
idx = which(values == "C", arr.ind = T)
values[idx] = "P"
idx = which(values == "S", arr.ind = T)
values[idx] = "P"



### set to ordered factor
# All "n" are automatically set to NA if "n" is none of the factor levels.
temp = as.data.frame(t(values[, 2:dim(values)[2]]))
colnames(temp) = values$uID
ord_values = data.frame(lapply(temp, function(x){factor(x, ordered = T, levels = c("X", "P", "H", "F", "M", "I", "M2", "T", "Mi"))}))



### calculate time in each stage for each cell
sum_values_msh_control = calculateTimeinEachStage(ord_values, "MSH_control")



## MSH HS34

### load data
dat = read.xlsx(namedat2, sheetIndex = 6)


# check and remove complete NA row
idx = which(apply(is.na(dat), 1, sum) == dim(dat)[2])
if (length(idx) > 0){
   dat = dat[-idx, ]
}



### add rowID
dat$uID = paste0("ID_", dat$ID, "_", dat$cells)
# reorder to put uID in the 6th column
dat = dat[, c(1:5, dim(dat)[2], 6:(dim(dat)[2]-1))]


# split measurements from information
info = dat[, 1:6]
values = dat[, 6:dim(dat)[2]]


### fill "n" if possible
# if state before and after "n" is the same, fill it with this value

# go through all rows
for (i in 1:dim(values)[1]){
   # go through all columns with measurements
   for (j in 2:dim(values)[2]){
      if (!is.na(values[i, j])){
         # check if it is an "n"
         if (values[i, j] == "n"){
            idx = match(values[i, j-1], values[i, j:dim(values)[2]])
            # if the index is not NA, then a match was found
            if (!is.na(idx)){
               values[i, j] = values[i, j-1]  
            }
         } 
      }
   }
}


### change C and S to P
idx = which(values == "C", arr.ind = T)
values[idx] = "P"
idx = which(values == "S", arr.ind = T)
values[idx] = "P"



### set to ordered factor
# All "n" are automatically set to NA if "n" is none of the factor levels.
temp = as.data.frame(t(values[, 2:dim(values)[2]]))
colnames(temp) = values$uID
ord_values = data.frame(lapply(temp, function(x){factor(x, ordered = T, levels = c("X", "P", "H", "F", "M", "I", "M2", "T", "Mi"))}))



### calculate time in each stage for each cell
sum_values_msh_HS34 = calculateTimeinEachStage(ord_values, "MSH_HS34")



## ATM control

### load data
dat = read.xlsx(namedat2, sheetIndex = 7)


# check and remove complete NA row
idx = which(apply(is.na(dat), 1, sum) == dim(dat)[2])
if (length(idx) > 0){
   dat = dat[-idx, ]
}



### add rowID
dat$uID = paste0("ID_", dat$ID, "_", dat$cells)
# reorder to put uID in the 6th column
dat = dat[, c(1:5, dim(dat)[2], 6:(dim(dat)[2]-1))]


# split measurements from information
info = dat[, 1:6]
values = dat[, 6:dim(dat)[2]]


### fill "n" if possible
# if state before and after "n" is the same, fill it with this value

# go through all rows
for (i in 1:dim(values)[1]){
   # go through all columns with measurements
   for (j in 2:dim(values)[2]){
      if (!is.na(values[i, j])){
         # check if it is an "n"
         if (values[i, j] == "n"){
            idx = match(values[i, j-1], values[i, j:dim(values)[2]])
            # if the index is not NA, then a match was found
            if (!is.na(idx)){
               values[i, j] = values[i, j-1]  
            }
         } 
      }
   }
}


### change C and S to P
idx = which(values == "C", arr.ind = T)
values[idx] = "P"
idx = which(values == "S", arr.ind = T)
values[idx] = "P"



### set to ordered factor
# All "n" are automatically set to NA if "n" is none of the factor levels.
temp = as.data.frame(t(values[, 2:dim(values)[2]]))
colnames(temp) = values$uID
ord_values = data.frame(lapply(temp, function(x){factor(x, ordered = T, levels = c("X", "P", "H", "F", "M", "I", "M2", "T", "Mi"))}))



### calculate time in each stage for each cell
sum_values_atm_control = calculateTimeinEachStage(ord_values, "ATM_control")



## ATM HS34

### load data
dat = read.xlsx(namedat2, sheetIndex = 8)


# check and remove complete NA row
idx = which(apply(is.na(dat), 1, sum) == dim(dat)[2])
if (length(idx) > 0){
   dat = dat[-idx, ]
}



### add rowID
dat$uID = paste0("ID_", dat$ID, "_", dat$cells)
# reorder to put uID in the 6th column
dat = dat[, c(1:5, dim(dat)[2], 6:(dim(dat)[2]-1))]


# split measurements from information
info = dat[, 1:6]
values = dat[, 6:dim(dat)[2]]


### fill "n" if possible
# if state before and after "n" is the same, fill it with this value

# go through all rows
for (i in 1:dim(values)[1]){
   # go through all columns with measurements
   for (j in 2:dim(values)[2]){
      if (!is.na(values[i, j])){
         # check if it is an "n"
         if (values[i, j] == "n"){
            idx = match(values[i, j-1], values[i, j:dim(values)[2]])
            # if the index is not NA, then a match was found
            if (!is.na(idx)){
               values[i, j] = values[i, j-1]  
            }
         } 
      }
   }
}


### change C and S to P
idx = which(values == "C", arr.ind = T)
values[idx] = "P"
idx = which(values == "S", arr.ind = T)
values[idx] = "P"



### set to ordered factor
# All "n" are automatically set to NA if "n" is none of the factor levels.
temp = as.data.frame(t(values[, 2:dim(values)[2]]))
colnames(temp) = values$uID
ord_values = data.frame(lapply(temp, function(x){factor(x, ordered = T, levels = c("X", "P", "H", "F", "M", "I", "M2", "T", "Mi"))}))



### calculate time in each stage for each cell
sum_values_atm_HS34 = calculateTimeinEachStage(ord_values, "ATM_HS34")



## combine all sum values in one df and export for stata
df = rbind(sum_values_spo11_control, sum_values_spo11_HS34, sum_values_dmc1_HS34, sum_values_dmc1_control, 
           sum_values_msh_control, sum_values_msh_HS34, sum_values_atm_control, sum_values_atm_HS34)
# export for stata
haven::write_dta(df, path = "sum_values_mutations.dta")
# export for R
saveRDS(df, file = "sum_values_mutations.Rds")


