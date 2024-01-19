## Startup code

# Check dependencies
# See: https://mail.google.com/mail/u/0/#inbox/QgrcJHsNjpqZNGJhNMVHcfGFDLLMfrvqqHl
.onLoad <- function(libname, pkgname) {
  #checkDepPackageVersion(dep_pkg="TMB", this_pkg="dsem")
  checkDepPackageVersion(dep_pkg="Matrix", this_pkg="TMB")
}

.onUnload <- function(libpath) {
  library.dynam.unload("dsem", libpath)
}