###################
# isle_royale
###################

data_dir =  R'(C:\Users\James.Thorson\Desktop\Git\dsem\data-raw)'

# Load and format
Data = read.csv( file.path(data_dir,"Data_wolves_moose_Isle_Royale_June2019--Tab1.csv"), skip=1 )
isle_royale = na.omit( Data[,c('year','wolves','moose')] )

# Export
setwd( R'(C:\Users\James.Thorson\Desktop\Git\dsem)' )
usethis::use_data( isle_royale, overwrite=TRUE )

##################
# lake_washington
##################

root_dir = R'(C:\Users\James.Thorson\Desktop\Git\DSEM-varying-paths)'
lakeWAplankton = read.csv( file.path(root_dir, "Data", "lakeWAplankton.csv"))
lake_washington =lakeWAplankton[,c("Year","Month","Temp","Daphnia","Leptodora","Cryptomonas")]

# Export
setwd( R'(C:\Users\James.Thorson\Desktop\Git\dsem)' )
usethis::use_data( lake_washington, overwrite=TRUE )

##################
# paramesium_didinium
##################

root_dir = R'(C:\Users\James.Thorson\Desktop\Git\DSEM-varying-paths)'
data_dir = file.path( root_dir, "Data" )
paramesium_didinium = read.csv( file.path(data_dir, "VeilleuxMS_fig11a.csv"), skip = 4, nrows = 71 )
paramesium_didinium = paramesium_didinium[,c('time','paramecium','didinium')]

# Export
setwd( R'(C:\Users\James.Thorson\Desktop\Git\dsem)' )
usethis::use_data( paramesium_didinium, overwrite=TRUE )

##################
# hare_lynx
##################

root_dir = R'(C:\Users\James.Thorson\Desktop\Git\DSEM-varying-paths)'
data_dir = file.path( root_dir, "Data" )
hare_lynx = read.csv( file.path(data_dir,"hare_lynx.csv") )

# Export
setwd( R'(C:\Users\James.Thorson\Desktop\Git\dsem)' )
usethis::use_data( hare_lynx, overwrite=TRUE )

##################
# pdo_departure_bay
##################

root_dir = R'(C:\Users\James.Thorson\Desktop\Git\DSEM-varying-paths)'
data_dir = file.path( root_dir, "Data" )

CSV = read.csv( file.path(data_dir,"Departure_Bay_PBS_-_Average_Monthly_Sea_Surface_Temperatures_1914-2025.csv"), skip = 1 )
Temp = as.matrix(CSV)
Temp = ifelse( Temp == 999.99, NA, Temp )

#pdo_dir = R'(C:\Users\james\OneDrive\Desktop\Work files (backup)\Collab-2018\2018 -- Pacific biological oscillation\PDO data)'
pdo_dir = data_dir
PDO = read.csv( file.path(pdo_dir,"PDO.csv"), col.names = colnames(Temp) )

#
pdo_departure_bay = merge( Temp[,1:2], PDO[,1:2], by = "YEAR" )
pdo_departure_bay = setNames( pdo_departure_bay, c("year", "departure_bay", "PDO") )

# Export
setwd( R'(C:\Users\James.Thorson\Desktop\Git\dsem)' )
usethis::use_data( pdo_departure_bay, overwrite=TRUE )

#################
# Compile help files
#################

setwd( R'(C:\Users\James.Thorson\Desktop\Git\dsem\src)' )
TMB::compile('dsem.cpp', framework = "TMBad" )
setwd( R'(C:\Users\James.Thorson\Desktop\Git\dsem)' )
devtools::document()
