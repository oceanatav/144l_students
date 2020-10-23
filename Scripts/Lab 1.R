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

#"alt -" is shortcut for <- 

##### Lab 2 #####

####  Initial Data Exploration ####

names(cal_fire_data)
dim(cal_fire_data)
class(cal_fire_data)
head(cal_fire_data)
str(cal_fire_data)

# ?function give R documentation
# ??function suggest all function containing "function"

county <- cal_fire_data$County_Unit

max_acres <- max(cal_fire_data$Total_Acres_Burned, na.rm=T) #max() returns maximum value
max_acres

max(cal_fire_data$Structures_Destroyed, na.rm = T) #max() will break if NA's exist, use argument na.rm=TRUE or = T

#### Basic Data Wrangling (dplyr) Functions ####

df1= select(cal_fire_data, County_Unit:Controlled_Date, Total_Acres_Burned:Civil_Fatalities) #to subset columns, remove 2 variables 
# select() pulls our variables

df2 = filter(df1, County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "ORANGE", "SAN DIEGO") & Total_Acres_Burned >= 500 | Fire_Name == "THOMAS") #to look at SoCal and acres burned greater than or equal to 500 acres, OR the thomas fire

#long-form of above: df2=filet(df1, County Unity="SANTA BARBARA), County_Unit="VENTURA", ...)

# filter() pulls out rows
# %in% = includes

df3 = arrange(df2, desc(Start_Date), Total_Acres_Burned)

# arrange() sorts dataframe. desc=descending order of (variable), then by Total Acress Burned

df4 = mutate_at(df3, vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0)
  
#replace NA's with 0's 

df5 = mutate(df4, Fatalities=Fire_Fatalities+Civil_Fatalities)
#total fatalities. adds new column named "Fatalities"

#### Data Wrangling; messing with time ####

#install.packages("lubridate")
library(lubridate)

df6 = mutate(df5,
             interv=interval(Start_Date, Controlled_Date), 
             dur=as.duration(interv), 
             days=as.numeric(dur, "days"))

#as.duration finds distance between dates prints to dur column
# as.numeric says we want to convert interval data to numeric

##used 15 lines to make 5 dataframes
#there's a better way, ka "piping"

#### Intro to Piping ####

##Goal: Restrict data to SoCal Coast, exlcudef ries burning <500 acres, add column showing sum of fatalities, change NA's to 0, arrange data and dates

# PIPE OPERATOR %>% = "and then" (ctrl + shift + m)
# %in% ?

socal_fires <- cal_fire_data %>% 
  filter(County_Unit %in% c("SANTA BARBARA","VENTURA","LOS ANGELES", "ORANGE", "SAN DIEGO") & Total_Acres_Burned >= 500 | Fire_Name == "THOMAS") %>%
  mutate_at(vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0) %>%
  mutate(Fatalities=Fire_Fatalities+Civil_Fatalities,
         interv=interval(Start_Date, Controlled_Date), dur=as.duration(interv), days=as.numeric(dur, "days"), County_Unit = ifelse(County_Unit == "VENTURA/SANTA BARBARA", "VENTURA", County_Unit)) %>%
      arrange(desc(Start_Date), Total_Acres_Burned)


#Thomas fire is saved as ventura/santa barbara. change a specific set of characters (ventura/ santa barbara) into Venture with ifelse()

#### GGPlot ####

socal.plot = socal_fires %>% rename(start=Start_Date, acres= Total_Acres_Burned, county = County_Unit) %>%  ggplot(aes(x=start, y=acres)) + geom_point(aes(color = county)) + ggtitle ("California SoCal Major Fires 2013-2018") + xlab("Date") + ylab("Acres Burned")

socal.plot + facet_grid(~county)




