library(tidyverse)

library(data.table)
library(rmarkdown)
library(kableExtra)

Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc")


setwd('C:\\Users/gqc954/PycharmProjects/Igenomix_development/')


FILENAME='HQ_families_40families_mzscoreMin075_18_1_2021.csv'

#READ DATA
df=fread(FILENAME,sep=',') 


chrOrder <-c(1:22,"X","Y")

df=df%>%
  filter(Chr!=0 & Chr!='XY') %>% 
  mutate(Chr=factor(Chr,levels=chrOrder, ordered=TRUE)) %>%
  arrange(Chr)




for (family in unique(df$familyid))
{  
  print(family)
  render("RSCRIPTS/Tables.Rmd",
           output_file=paste0(family,"_mzscore_NNCN_075.html"),
           params=list(new_title=paste("Output - thr 075 for all", family)))
  
}

