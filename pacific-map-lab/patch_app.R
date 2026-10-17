args <- commandArgs(trailingOnly = TRUE)
app_path <- if (length(args)) args[[1]] else "/app/app.R"
lines <- readLines(app_path, warn = FALSE)

replace_fixed <- function(old, new, expected = 1L) {
  hits <- grepl(old, lines, fixed = TRUE)
  n <- sum(hits)
  if (!identical(n, as.integer(expected))) {
    stop(sprintf("Patch target count mismatch for %s: expected %d, found %d", old, expected, n))
  }
  lines[hits] <<- sub(old, new, lines[hits], fixed = TRUE)
}

replace_fixed('APP_VERSION <- "2.3.2-web"', 'APP_VERSION <- "2.3.3-web-lowmem"')
replace_fixed('MAX_GRID_CELLS <- 9e5', 'MAX_GRID_CELLS <- 1.5e5')
replace_fixed('PREVIEW_RASTER_LIMIT <- 4.5e5', 'PREVIEW_RASTER_LIMIT <- 6e4')
replace_fixed('CONTOUR_RASTER_LIMIT <- 3.5e5', 'CONTOUR_RASTER_LIMIT <- 4e4')
replace_fixed('"Coastal California" = list(lon = -121.45, lat = 37.25, width = 940, height = 1175, res = 1)', '"Coastal California" = list(lon = -121.45, lat = 37.25, width = 940, height = 1175, res = 4)')
replace_fixed('"Custom" = list(lon = -121.45, lat = 37.25, width = 940, height = 1175, res = 1)', '"Custom" = list(lon = -121.45, lat = 37.25, width = 940, height = 1175, res = 4)')
replace_fixed('selected = "1"', 'selected = "4"')
replace_fixed('selected = "medium"', 'selected = "small"')
replace_fixed('sliderInput("densify_m", "Coastline vertex spacing after crop (m)", min = 50, max = 2000, value = 500, step = 50)', 'sliderInput("densify_m", "Coastline vertex spacing after crop (m)", min = 50, max = 2000, value = 1200, step = 50)')
replace_fixed('safe_num(input$bathy_res, 1)', 'safe_num(input$bathy_res, 4)', expected = 3L)
replace_fixed('safe_num(isolate(input$bathy_res), 1)', 'safe_num(isolate(input$bathy_res), 4)')
replace_fixed('input$ne_scale %||% "medium"', 'input$ne_scale %||% "small"', expected = 2L)
replace_fixed('safe_num(input$densify_m, 500)', 'safe_num(input$densify_m, 1200)')
replace_fixed('safe_num(isolate(input$densify_m), 500)', 'safe_num(isolate(input$densify_m), 1200)')

writeLines(lines, app_path, useBytes = TRUE)
message("Applied Pacific Map Lab low-memory hosting patch to ", app_path)
