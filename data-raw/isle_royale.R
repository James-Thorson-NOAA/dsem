## code to prepare `DATASET` dataset goes here

data_dir =  R'(C:\Users\James.Thorson\Desktop\Git\dsem\data-raw)'

# Load and format
Data = read.csv( file.path(data_dir,"Data_wolves_moose_Isle_Royale_June2019--Tab1.csv"), skip=1 )
isle_royale = na.omit( Data[,c('year','wolves','moose')] )

# Export
setwd( R'(C:\Users\James.Thorson\Desktop\Git\dsem)' )
usethis::use_data( isle_royale, overwrite=TRUE )
