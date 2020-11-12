#' shortcut for paste0()
#' @export
p <- paste0

#' Update ORFik version
#' @param branch default: "master"
#' @param user default: Roleren
#' @export
updateORFik <- function(branch = "master", user = "Roleren")  {
  system(p("rm -r ", .libPaths()[1], "/00LOCK-ORFik"))
  devtools::install_github(paste0(user, "/ORFik"), ref = branch, )
}
