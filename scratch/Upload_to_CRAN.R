

library(devtools)
setwd( R'(C:\Users\James.Thorson\Desktop\Git\dsem)' )

# Compile
if( FALSE ){
  document()
}

# Test install
install_local(force=TRUE, dep=TRUE, build_vignettes=TRUE, upgrade=FALSE)
#install_local(force=TRUE, dep=TRUE, build_vignettes=FALSE, upgrade=FALSE)
browseVignettes("dsem")

#
if( FALSE ){
  library(TMB)
  setwd( R'(C:\Users\James.Thorson\Desktop\Git\dsem\src)' )
  compile("dsem.cpp")
}

# Try building vignetttes
if( FALSE ){
  library(rmarkdown)
  setwd( R'(C:\Users\James.Thorson\Desktop\Git\dsem\vignettes)' )
  devtools::build_rmd("vignettes/vignette.Rmd")
  render( file.path(getwd(),"vignette.Rmd"), pdf_document())
}

# Try mapping dependencies
if( FALSE ){
  # Simple way
  library(renv)
  x = dependencies()

  # Listed dependencies
  tools::package_dependencies("dsem")
  tools::package_dependencies("dynlm")

  # All
  pack <- available.packages()
  pack["dynlm","Depends"]
  packrat:::recursivePackageDependencies("dynlm", ignore = "", lib.loc = .libPaths()[1], fields="Imports")
}

# Run checks ... doesn't seem to work
file.remove( file.path("vignettes","vignette.pdf") )
#check( remote = TRUE, manual=TRUE )
check( manual=TRUE )

# Check manual
if( FALSE ){
  tools::texi2pdf
}

# Check online but document first!
document()
check_win_devel()

# Submit to CRAN via devtools .. not preferred!
if( FALSE ){
  file.remove( file.path("vignettes","vignette.pdf") )
  release()
}

# Build for uploading via web interface
root_dir = R'(C:\Users\James.Thorson\Desktop\Git\dsem\scratch)'
# https://cran.r-project.org/submit.html
build( path=root_dir, manual=TRUE )


