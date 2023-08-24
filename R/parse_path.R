
#' Parse path
#'
#' \code{parse_path} is copied from \code{sem::parse.path}
#'
#' Copied with permission from John Fox under licence GPL (>= 2)
#'
#' @return Tagged-list defining variables and direction for a specified path coefficient
#'
#' @param path text to parse
parse_path <-
function( path ){
  path.1 <- gsub("-", "", gsub(" ", "", path))
  direction <- if(regexpr("<>", path.1) > 0){
    2
  }else if(regexpr("<", path.1) > 0){
    -1
  }else if(regexpr(">", path.1) > 0){
    1
  }else{
    stop(paste("ill-formed path:", path))
  }
  path.1 <- strsplit(path.1, "[<>]")[[1]]
  out = list(first = path.1[1], second = path.1[length(path.1)], direction = direction)
  return(out)
}
