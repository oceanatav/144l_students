#####EEMB 144L Intro to R CalFire #####
#Oceana Tavasieff
#7 Oct 2020

library(tidyverse)
library(readxl)
# to read in excel packages
#### Load Data ####

#Remember excel files can have multiple sheets!! ( unlike csv)
excel_sheets("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx")
# Returns the names of the sheets, of which there are two

cal_fire_metadata = read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet=1)
cal_fire_data= read_excel("Input_Data/week1/2013_2019_CALFIRE_Redbook.xlsx", sheet=2)

#argument sheet in read_excel can accept integer input (number of sheet) or character ("Name of sheet")

