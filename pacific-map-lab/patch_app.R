args <- commandArgs(trailingOnly = TRUE)
app_path <- if (length(args)) args[[1]] else "/app/app.R"
if (!file.exists(app_path) || file.info(app_path)$size < 1) {
  stop("Pacific Map Lab app.R was not reconstructed correctly.")
}
message("Pacific Map Lab web-lite archive is ready at ", app_path)
