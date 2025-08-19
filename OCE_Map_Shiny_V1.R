# Shiny app to create, edit, and export high-resolution 
# bathymetry maps of the Santa Barbara Channel using queried NOAA ETOPO (2022).

# Michael R. Cornish
# cornish@ucsb.edu
# Code and updates found at: https://github.com/mrcornish/OCE_Maps/tree/main

#--------------------------------------------------------------------------
req_pkgs <- c(
  "shiny","ggplot2","sf","marmap","rnaturalearth","rnaturalearthdata",
  "ggspatial","ggrepel","dplyr","units","tibble","viridisLite",
  "colourpicker","ragg","metR","scales","grid"
)
new_pkgs <- setdiff(req_pkgs, rownames(installed.packages()))
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org")
invisible(lapply(req_pkgs, library, character.only = TRUE))

options(timeout = max(300, getOption("timeout", 60)))

# 1) Helpers -----------------------------------------------------------

deg_per_km_lat <- function() 1 / 110.574     # degrees latitude per km
deg_per_km_lon <- function(lat) 1 / (111.320 * cospi(lat / 180))

utm_epsg_from_lonlat <- function(lon, lat) {
  zone <- floor((lon + 180) / 6) + 1
  zone <- max(1, min(60, zone))
  if (lat >= 0) 32600 + zone else 32700 + zone
}

make_bbox_from_center_km <- function(clon, clat, width_km, height_km) {
  dx <- (width_km  * deg_per_km_lon(clat)) / 2
  dy <- (height_km * deg_per_km_lat())    / 2
  xmin <- max(-180, clon - dx); xmax <- min(180, clon + dx)
  ymin <- max(-90,  clat - dy); ymax <- min(90,  clat + dy)
  list(
    bbox = sf::st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), crs = sf::st_crs(4326)),
    xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax
  )
}

make_depth_palette <- function(kind, n, custom_cols, reverse = FALSE) {
  cols <- switch(kind,
                 "Blues"    = colorRampPalette(c("#ffffff","#dff3ff","#b9e0ff","#8fc7ff",
                                                 "#5ea8ff","#2b86ff","#0b5bd6","#063b8c","#03274d"))(n),
                 "Viridis"  = viridisLite::viridis(n),
                 "Magma"    = viridisLite::magma(n),
                 "CyanBlue" = colorRampPalette(c("white","cyan","blue"))(n),
                 "Custom"   = if (length(custom_cols) >= 2) colorRampPalette(custom_cols)(n) else viridisLite::viridis(n),
                 viridisLite::viridis(n)
  )
  if (reverse) rev(cols) else cols
}

parse_manual_isobaths <- function(txt) {
  if (is.null(txt) || !nzchar(txt)) return(numeric(0))
  out <- suppressWarnings(as.numeric(strsplit(gsub("\\s", "", txt), ",")[[1]]))
  out[is.finite(out)]
}

# Robust NOAA fetcher with retries + resolution fallback
fetch_bathy_resilient <- function(lon_min, lon_max, lat_min, lat_max, res_pref,
                                  max_retries = 4, sleep_base = 1) {
  res_chain <- unique(c(as.numeric(res_pref), 0.25, 0.5, 1.0))
  last_err <- NULL
  for (res in res_chain) {
    for (i in seq_len(max_retries)) {
      ok <- try({
        b <- marmap::getNOAA.bathy(lon1 = lon_min, lon2 = lon_max,
                                   lat1 = lat_min, lat2 = lat_max,
                                   resolution = res)
        attr(b, "resolution_tried") <- res
        b
      }, silent = TRUE)
      if (!inherits(ok, "try-error")) return(ok)
      last_err <- ok
      Sys.sleep(sleep_base * 2^(i - 1))  # exponential backoff
    }
  }
  stop("NOAA bathymetry could not be fetched after retries/fallbacks. Last error: ",
       conditionMessage(attr(last_err, "condition")))
}

# Degree–minute tick labels (ASCII apostrophe to avoid “boxes”)
to_dm_lon <- function(x){
  deg  <- floor(abs(x)); mins <- round((abs(x) - deg) * 60)
  bump <- mins == 60
  if (any(bump)) { deg[bump] <- deg[bump] + 1; mins[bump] <- 0 }
  paste0(deg, "\u00B0", sprintf("%02d", mins), "'", ifelse(x < 0, "W", "E"))
}
to_dm_lat <- function(x){
  deg  <- floor(abs(x)); mins <- round((abs(x) - deg) * 60)
  bump <- mins == 60
  if (any(bump)) { deg[bump] <- deg[bump] + 1; mins[bump] <- 0 }
  paste0(deg, "\u00B0", sprintf("%02d", mins), "'", ifelse(x < 0, "S", "N"))
}

# Paper-grain underlay for vintage theme
make_paper_grain <- function(lon_min, lon_max, lat_min, lat_max, nx = 600, ny = 600) {
  g <- expand.grid(lon = seq(lon_min, lon_max, length.out = nx),
                   lat = seq(lat_min, lat_max, length.out = ny))
  n <- matrix(runif(nx * ny), nrow = ny, ncol = nx)
  k <- matrix(1, nrow = 5, ncol = 5) / 25
  pad <- matrix(0, nrow = ny + 4, ncol = nx + 4)
  pad[3:(ny+2), 3:(nx+2)] <- n
  sm <- matrix(0, nrow = ny, ncol = nx)
  for (i in 1:5) for (j in 1:5) sm <- sm + k[i, j] * pad[i:(i+ny-1), j:(j+nx-1)]
  g$grain <- (as.vector(sm) - min(sm)) / (max(sm) - min(sm))
  g
}

# Nautical-miles scale bar (vintage)
make_nmi_scalebar <- function(lon0, lat0, nmi = 5, segments = 5,
                              height_deg = 0.025, label = "5 nmi") {
  lon_per_nmi <- (1/60) / cos(lat0 * pi/180)
  rects <- tibble::tibble(
    xmin = lon0 + (0:(segments-1)) * lon_per_nmi,
    xmax = lon0 + (1:segments)     * lon_per_nmi,
    ymin = lat0, ymax = lat0 + height_deg,
    fill = rep(c("#1c4a64","#eae2cf"), length.out = segments)
  )
  ticks <- tibble::tibble(
    x = lon0 + (0:segments) * lon_per_nmi,
    y0 = lat0, y1 = lat0 + height_deg * c(1, rep(0.6, segments-1), 1)
  )
  lbl <- tibble::tibble(
    x = lon0 + nmi * lon_per_nmi / 2,
    y = lat0 + height_deg * 1.25,
    lab = label
  )
  list(rects = rects, ticks = ticks, lbl = lbl)
}

# Build ggplot figure from precomputed pieces
build_map <- function(sea_df, bathy_df, land_clip, coast_clip, isobaths, bbox, opts) {
  # Common pieces
  pal <- make_depth_palette(opts$palette, n = opts$pal_n,
                            custom_cols = opts$pal_custom, reverse = opts$pal_reverse)
  
  # Branch by theme ----------------------------------------------------
  if (identical(opts$map_theme, "Vintage nautical")) {
    # Bounds & ticks
    lon_min <- bbox$xmin; lon_max <- bbox$xmax
    lat_min <- bbox$ymin; lat_max <- bbox$ymax
    lon_breaks <- seq(lon_min, lon_max, length.out = 6)
    lat_breaks <- seq(lat_min, lat_max, length.out = 6)
    
    # Depth in fathoms for fill
    sea_df <- sea_df |>
      dplyr::mutate(depth_fm = depth_m / 1.8288)
    
    # Fathom breaks (labels & lines) in meters (negative)
    fm_to_m <- function(fm) -(fm * 1.8288)
    iso_fm_shallow <- c(3, 6, 10, 20, 30, 40)
    # deep labels use finer lines; go until ~max
    max_fm <- floor(max(sea_df$depth_fm, na.rm = TRUE) / 50) * 50
    iso_fm_deep <- seq(50, max(50, max_fm), by = 50)
    iso_shallow_m <- fm_to_m(iso_fm_shallow)
    iso_deep_m    <- fm_to_m(iso_fm_deep)
    
    # Water color steps (no legend)
    water_cols <- c("#bfe3f1","#d4edf7","#ecf5fa","#ffffff")
    
    # Paper grain
    grain_df <- make_paper_grain(lon_min, lon_max, lat_min, lat_max, nx = 600, ny = 600)
    
    # Scale bar position (BL)
    sb_lat <- lat_min + 0.08
    sb_lon <- lon_min + 0.12
    sb <- make_nmi_scalebar(sb_lon, sb_lat, nmi = 5, segments = 5, height_deg = 0.03,
                            label = "5 nautical miles")
    
    # Neatlines
    frame <- tibble::tibble(xmin = lon_min, xmax = lon_max, ymin = lat_min, ymax = lat_max)
    frame_inner <- tibble::tibble(xmin = lon_min + 0.01, xmax = lon_max - 0.01,
                                  ymin = lat_min + 0.01, ymax = lat_max - 0.01)
    
    p <- ggplot() +
      # Paper grain
      geom_raster(data = grain_df, aes(lon, lat, alpha = grain),
                  inherit.aes = FALSE, fill = "#3b2f20") +
      scale_alpha(range = c(0.01, 0.05), guide = "none") +
      
      # Water raster (fathoms)
      geom_raster(data = sea_df, aes(lon, lat, fill = depth_fm), interpolate = TRUE) +
      
      # Isobaths (two weights)
      geom_contour(data = bathy_df, aes(lon, lat, z = z),
                   breaks = iso_shallow_m, linewidth = 0.35, color = "#2c6486", alpha = 0.9) +
      geom_contour(data = bathy_df, aes(lon, lat, z = z),
                   breaks = iso_deep_m,    linewidth = 0.25, color = "#6f8ea3", alpha = 0.8) +
      
      # Fathom labels on contours
      metR::geom_text_contour(
        data = bathy_df,
        aes(lon, lat, z = z, label = after_stat(round(abs(level)/1.8288))),
        breaks = iso_shallow_m, size = 3, colour = "#2c6486", check_overlap = TRUE
      ) +
      metR::geom_text_contour(
        data = bathy_df,
        aes(lon, lat, z = z, label = after_stat(round(abs(level)/1.8288))),
        breaks = iso_deep_m, size = 2.7, colour = "#6f8ea3", check_overlap = TRUE
      ) +
      
      # Land & coastline
      geom_sf(data = land_clip,  fill = "#eae2cf", color = NA) +
      geom_sf(data = coast_clip, color = "#2f2c24", linewidth = 0.55,
              lineend = "round", linejoin = "round") +
      
      # Optional labels (cities/islands) if provided
      { if (isTRUE(opts$show_cities) && !is.null(opts$cities_df))
        geom_point(data = opts$cities_df, aes(lon, lat), size = opts$city_pt_size, color = "#333333")
        else list() } +
      { if (isTRUE(opts$show_cities) && !is.null(opts$cities_df))
        ggrepel::geom_text_repel(
          data = opts$cities_df, aes(lon, lat, label = name),
          size = opts$city_lab_size, family = "serif", seed = 42,
          box.padding = 0.25, point.padding = 0.1, min.segment.length = 0,
          segment.size = 0.25, color = "#2b2b2b"
        )
        else list() } +
      { if (isTRUE(opts$show_islands) && !is.null(opts$islands_df))
        ggrepel::geom_text_repel(
          data = opts$islands_df, aes(lon, lat, label = name),
          size = opts$island_lab_size, family = "serif", seed = 42,
          box.padding = 0.25, point.padding = 0.1, min.segment.length = 0,
          segment.size = 0.25, color = "#2b2b2b"
        )
        else list() } +
      
      # Scale bar
      annotate("rect", xmin = sb$rects$xmin, xmax = sb$rects$xmax,
               ymin = sb$rects$ymin, ymax = sb$rects$ymax,
               fill = sb$rects$fill, color = "#1b2a34", linewidth = 0.25) +
      annotate("segment", x = sb$ticks$x, xend = sb$ticks$x,
               y = sb$ticks$y0, yend = sb$ticks$y1,
               color = "#1b2a34", linewidth = 0.25) +
      annotate("text", x = sb$lbl$x, y = sb$lbl$y, label = sb$lbl$lab,
               family = "serif", size = 3.2, color = "#1b2a34") +
      
      # Compass rose
      ggspatial::annotation_north_arrow(
        location = opts$north_loc, which_north = "true",
        style = ggspatial::north_arrow_fancy_orienteering
      ) +
      
      # Neatlines
      annotate("rect", xmin = frame$xmin, xmax = frame$xmax, ymin = frame$ymin, ymax = frame$ymax,
               fill = NA, color = "#2e2c24", linewidth = 0.7) +
      annotate("rect", xmin = frame_inner$xmin, xmax = frame_inner$xmax,
               ymin = frame_inner$ymin, ymax = frame_inner$ymax,
               fill = NA, color = "#2e2c24", linewidth = 0.4) +
      
      # Axes/palette/theme
      coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max),
               expand = FALSE, crs = sf::st_crs(4326)) +
      scale_fill_gradientn(colours = water_cols,
                           values  = scales::rescale(c(0, 5, 10, 30, max_fm)),
                           guide   = "none") +
      scale_x_continuous(breaks = lon_breaks, labels = to_dm_lon) +
      scale_y_continuous(breaks = lat_breaks, labels = to_dm_lat) +
      labs(
        title    = ifelse(nzchar(opts$title), opts$title, "Santa Barbara Channel"),
        subtitle = ifelse(nzchar(opts$subtitle), opts$subtitle,
                          "Vintage-style nautical chart • Depths in fathoms"),
        caption  = opts$caption
      ) +
      theme_minimal(base_size = opts$base_size) +
      theme(
        text             = element_text(family = "serif", color = "#2b2b2b"),
        plot.background  = element_rect(fill = "#f6f2e8", color = NA),
        panel.background = element_rect(fill = "#ffffff", color = NA),
        panel.grid       = element_blank(),
        axis.title       = element_blank(),
        axis.text        = element_text(color = "#5a564b"),
        plot.title       = element_text(size = 16, face = "bold"),
        plot.subtitle    = element_text(size = 11, margin = margin(b = 6)),
        plot.caption     = element_text(size = 9, color = "#514f47", margin = margin(t = 6)),
        legend.position  = "none",
        plot.margin      = margin(8, 12, 10, 12)
      )
    
    return(p)
  }
  
  # --- Modern theme (original) --------------------------------------
  p <- ggplot() +
    geom_raster(
      data = sea_df,
      aes(x = lon, y = lat, fill = depth_m),
      interpolate = isTRUE(opts$interpolate)
    )
  
  if (isTRUE(opts$show_contours) && length(isobaths) > 0) {
    p <- p + geom_contour(
      data = bathy_df,
      aes(x = lon, y = lat, z = z),
      breaks = isobaths,
      color  = opts$contour_color,
      size   = opts$contour_width,
      alpha  = opts$contour_alpha
    )
  }
  
  p <- p +
    geom_sf(data = land_clip,  fill = opts$land_fill, color = NA) +
    geom_sf(data = coast_clip, color = opts$coast_color, linewidth = opts$coast_width)
  
  if (isTRUE(opts$show_cities) && !is.null(opts$cities_df)) {
    p <- p + geom_point(data = opts$cities_df, aes(lon, lat), size = opts$city_pt_size) +
      ggrepel::geom_text_repel(
        data = opts$cities_df, aes(lon, lat, label = name),
        size = opts$city_lab_size, min.segment.length = 0, seed = 42
      )
  }
  if (isTRUE(opts$show_islands) && !is.null(opts$islands_df)) {
    p <- p + ggrepel::geom_text_repel(
      data = opts$islands_df, aes(lon, lat, label = name),
      size = opts$island_lab_size, min.segment.length = 0, box.padding = 0.3, seed = 42
    )
  }
  
  p <- p +
    ggspatial::annotation_scale(
      location = opts$scale_loc, width_hint = 0.3,
      text_cex = opts$scale_text_cex, line_width = opts$scale_linewidth
    ) +
    ggspatial::annotation_north_arrow(
      location = opts$north_loc, which_north = "true",
      style = ggspatial::north_arrow_fancy_orienteering
    ) +
    coord_sf(
      xlim = c(bbox$xmin, bbox$xmax),
      ylim = c(bbox$ymin, bbox$ymax),
      expand = FALSE, crs = sf::st_crs(4326)
    ) +
    scale_fill_gradientn(
      name   = opts$legend_title,
      colours = pal,
      trans   = opts$depth_trans,
      guide   = guide_colorbar(reverse = isTRUE(opts$legend_reverse))
    ) +
    labs(
      title    = opts$title,
      subtitle = opts$subtitle,
      caption  = opts$caption
    ) +
    theme_minimal(base_size = opts$base_size) +
    theme(
      legend.position   = opts$legend_pos,
      panel.grid.major  = element_line(color = opts$grid_color, linewidth = opts$grid_width),
      panel.grid.minor  = element_blank(),
      axis.title        = element_blank()
    )
  
  p
}

# 2) UI ----------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Santa Barbara Channel Bathymetry Map Builder"),
  sidebarLayout(
    sidebarPanel(
      h4("Extent & Data"),
      numericInput("center_lon", "Center longitude (°)", value = -119.75, min = -180, max = 180, step = 0.01),
      numericInput("center_lat", "Center latitude (°)",  value = 34.10,  min =  -90, max =  90, step = 0.01),
      sliderInput("width_km",  "Map width (km)",  min = 10, max = 800, value = 180, step = 10),
      sliderInput("height_km", "Map height (km)", min = 10, max = 800, value = 120, step = 10),
      
      selectInput("bathy_res", "Bathymetry resolution (arc-min; lower = finer, slower)",
                  choices = c("0.25","0.5","1.0","0.1"), selected = "0.25"),
      
      radioButtons("ne_scale", "Natural Earth scale (land & coast)",
                   choices = c("small (110m)" = "small", "medium (50m)" = "medium", "large (10m)" = "large"),
                   selected = "large", inline = TRUE),
      
      numericInput("densify_m", "Densify spacing (meters, after crop)", value = 200, min = 50, max = 1000, step = 50),
      
      tags$hr(),
      
      h4("Map theme"),
      selectInput("map_theme", "Choose a look",
                  choices = c("Modern", "Vintage nautical"),
                  selected = "Modern"),
      
      tags$hr(),
      
      h4("Isobaths (Modern theme)"),
      checkboxInput("show_contours", "Draw isobath contours", value = TRUE),
      radioButtons("contour_mode", "Breaks", choices = c("Auto by step", "Manual list"),
                   selected = "Auto by step", inline = TRUE),
      conditionalPanel(
        "input.contour_mode == 'Auto by step'",
        numericInput("iso_min_abs", "Start depth |m| (near shore)", value = 10,  min = 1,  max = 10000, step = 5),
        numericInput("iso_max_abs", "End depth |m| (offshore)",     value = 5000, min = 10, max = 10000, step = 10),
        numericInput("iso_step",    "Step (m)",                     value = 50,   min = 1,  max = 2000,  step = 1)
      ),
      conditionalPanel(
        "input.contour_mode == 'Manual list'",
        helpText("Enter comma-separated negative depths (m), e.g. -10, -20, -50, -100, -200, -500, -1000"),
        textInput("iso_manual", NULL, value = "-10,-20,-50,-100,-200,-500,-1000")
      ),
      
      tags$hr(),
      
      h4("Aesthetics (Modern theme)"),
      selectInput("palette", "Depth palette",
                  choices = c("Blues", "Viridis", "Magma", "CyanBlue", "Custom"), selected = "Blues"),
      sliderInput("pal_n", "Palette resolution (colors)", min = 3, max = 15, value = 9, step = 1),
      checkboxInput("pal_reverse", "Reverse palette", value = FALSE),
      uiOutput("custom_palette_ui"),
      selectInput("depth_trans", "Depth color transform", choices = c("sqrt", "identity", "log10"), selected = "sqrt"),
      checkboxInput("interpolate", "Interpolate raster shading", value = TRUE),
      
      colourInput("land_fill",  "Land fill",        value = "#ECECEC"),
      colourInput("coast_color","Coastline color",  value = "#5A5A5A"),
      numericInput("coast_width","Coastline linewidth", value = 0.3, min = 0.05, step = 0.05),
      
      colourInput("contour_color", "Isobath color",   value = "#FFFFFF"),
      numericInput("contour_width", "Isobath linewidth", value = 0.25, min = 0.05, step = 0.05),
      sliderInput("contour_alpha", "Isobath alpha",  min = 0, max = 1, value = 0.8, step = 0.05),
      
      selectInput("legend_pos", "Legend position", choices = c("right","left","bottom","top","none"), selected = "right"),
      textInput("legend_title", "Legend title", "Depth (m)"),
      checkboxInput("legend_reverse", "Reverse legend scale", value = TRUE),
      numericInput("base_size", "Base font size", value = 11, min = 6, max = 24, step = 1),
      colourInput("grid_color", "Gridline color", value = "grey80"),
      numericInput("grid_width", "Gridline width", value = 0.2, min = 0.05, max = 1, step = 0.05),
      
      tags$hr(),
      
      h4("Labels, Scale, Arrow"),
      checkboxInput("show_cities",  "Show city points & labels", value = TRUE),
      numericInput("city_pt_size",  "City point size", value = 2.2, min = 0.5, step = 0.1),
      numericInput("city_lab_size", "City label size", value = 3.2, min = 1, step = 0.1),
      checkboxInput("show_islands", "Show island labels", value = TRUE),
      numericInput("island_lab_size", "Island label size", value = 3.0, min = 1, step = 0.1),
      
      selectInput("scale_loc", "Scale bar location", choices = c("bl","br","tl","tr"), selected = "bl"),
      numericInput("scale_text_cex", "Scale bar text cex", value = 0.8, min = 0.5, max = 2, step = 0.1),
      numericInput("scale_linewidth", "Scale bar line width", value = 0.3, min = 0.1, max = 1, step = 0.05),
      
      selectInput("north_loc", "North arrow location", choices = c("bl","br","tl","tr"), selected = "tr"),
      
      tags$hr(),
      
      h4("Titles"),
      textInput("title",    "Title"),
      textInput("subtitle", "Subtitle", "Bathymetry: NOAA ETOPO 2022"),
      textInput("caption",  "Caption", "Your caption here"),
      
      tags$hr(),
      
      actionButton("render", "Render map", class = "btn-primary"),
      br(), br(),
      
      h4("Export PNG (high-res)"),
      textInput("file_base", "File name (no extension)", "custom_map"),
      numericInput("out_w", "Width (inches)", 10, min = 2, step = 0.5),
      numericInput("out_h", "Height (inches)", 7,  min = 2, step = 0.5),
      numericInput("out_dpi", "DPI", 300, min = 72, step = 10),
      downloadButton("download_png", "Download PNG")
    ),
    
    mainPanel(
      plotOutput("map_plot", height = "720px"),
      tags$hr(),
      verbatimTextOutput("status")
    )
  )
)

# 3) Server ------------------------------------------------------------

server <- function(input, output, session) {
  
  # Cache last good bathy-derived payload
  last_good_payload <- reactiveVal(NULL)
  last_used_res     <- reactiveVal(NA_real_)
  
  # Constant labels
  cities_default <- tibble(
    name = c("Santa Barbara","Goleta","Ventura"),
    lon  = c(-119.6982,-119.8276,-119.2290),
    lat  = c( 34.4208, 34.4358, 34.2746)
  )
  islands_default <- tibble(
    name = c("San Miguel I.","Santa Rosa I.","Santa Cruz I.","Anacapa I."),
    lon  = c(-120.37, -120.08, -119.77, -119.36),
    lat  = c(  34.04,   33.95,   34.03,   34.01)
  )
  
  # Custom palette UI
  output$custom_palette_ui <- renderUI({
    if (input$palette != "Custom") return(NULL)
    tagList(
      fluidRow(
        column(6, colourInput("c1", "Color 1", value = "#ffffff")),
        column(6, colourInput("c2", "Color 2", value = "#dff3ff"))
      ),
      fluidRow(
        column(6, colourInput("c3", "Color 3", value = "#8fc7ff")),
        column(6, colourInput("c4", "Color 4", value = "#5ea8ff"))
      ),
      fluidRow(
        column(6, colourInput("c5", "Color 5", value = "#0b5bd6")),
        column(6, colourInput("c6", "Color 6", value = "#03274d"))
      )
    )
  })
  
  # Reactive bbox
  bbox_reactive <- reactive({
    make_bbox_from_center_km(input$center_lon, input$center_lat, input$width_km, input$height_km)
  })
  
  # Heavy work behind the Render button
  map_data <- eventReactive(input$render, {
    withProgress(message = "Building map…", value = 0, {
      # 1) Extent
      bb <- bbox_reactive(); incProgress(0.05)
      lon_min <- bb$xmin; lon_max <- bb$xmax
      lat_min <- bb$ymin; lat_max <- bb$ymax
      
      # 2) NOAA bathymetry (robust)
      res_req <- as.numeric(input$bathy_res)
      incProgress(0.1, detail = "Fetching bathymetry (NOAA)…")
      
      payload <- tryCatch({
        bathy <- fetch_bathy_resilient(lon_min, lon_max, lat_min, lat_max, res_pref = res_req)
        res_used <- attr(bathy, "resolution_tried")
        incProgress(0.25, detail = sprintf("Processing bathymetry (%.2f arc-min)…", res_used))
        
        bathy_df <- marmap::fortify.bathy(bathy) |>
          dplyr::rename(lon = x, lat = y, z = z)
        
        sea_df <- bathy_df |>
          dplyr::filter(z <= 0) |>
          dplyr::mutate(depth_m = -z)  # depth_fm added in build_map when needed
        
        # Save cache & meta
        last_used_res(res_used)
        out <- list(bathy_df = bathy_df, sea_df = sea_df, res_used = res_used)
        last_good_payload(out)
        out
      }, error = function(e) {
        # Fall back to cache if available
        cached <- last_good_payload()
        if (!is.null(cached)) {
          showNotification(
            paste0("NOAA 503/unavailable. Using cached bathymetry (", cached$res_used,
                   " arc-min). Error: ", e$message),
            type = "warning", duration = 10
          )
          cached
        } else {
          showNotification(
            paste0("NOAA bathymetry unavailable: ", e$message,
                   " Try a coarser resolution (0.5–1.0 arc-min) and click Render again."),
            type = "error", duration = 10
          )
          validate("NOAA bathymetry temporarily unavailable.")
        }
      })
      
      # 3) Natural Earth land/coast (crop → transform → segmentize → back)
      incProgress(0.45, detail = "Fetching Natural Earth…")
      land_hi  <- rnaturalearth::ne_download(scale = input$ne_scale, type = "land",
                                             category = "physical", returnclass = "sf")
      coast_hi <- rnaturalearth::ne_coastline(scale = input$ne_scale, returnclass = "sf")
      land_hi  <- sf::st_make_valid(land_hi)
      
      incProgress(0.55, detail = "Cropping layers…")
      s2_prev <- sf::sf_use_s2()
      sf::sf_use_s2(TRUE)
      bbox_sf <- bb$bbox
      land_clip  <- suppressWarnings(sf::st_crop(land_hi,  bbox_sf))
      coast_clip <- suppressWarnings(sf::st_crop(coast_hi, bbox_sf))
      sf::sf_use_s2(s2_prev)
      
      land_clip  <- land_clip[!sf::st_is_empty(land_clip), , drop = FALSE]
      coast_clip <- coast_clip[!sf::st_is_empty(coast_clip), , drop = FALSE]
      
      incProgress(0.65, detail = "Densifying coastlines (meters)…")
      epsg <- utm_epsg_from_lonlat((lon_min + lon_max)/2, (lat_min + lat_max)/2)
      
      land_clip  <- land_clip |>
        sf::st_transform(epsg) |>
        sf::st_segmentize(units::set_units(input$densify_m, "m")) |>
        sf::st_transform(4326) |>
        sf::st_make_valid()
      
      coast_clip <- coast_clip |>
        sf::st_cast("MULTILINESTRING", warn = FALSE) |>
        sf::st_transform(epsg) |>
        sf::st_segmentize(units::set_units(input$densify_m, "m")) |>
        sf::st_transform(4326)
      
      # 4) Isobaths (Modern UI)
      incProgress(0.75, detail = "Computing isobaths…")
      isobaths <- if (input$contour_mode == "Manual list") {
        parse_manual_isobaths(input$iso_manual)
      } else {
        -as.numeric(seq(input$iso_min_abs, input$iso_max_abs, by = input$iso_step))
      }
      isobaths <- isobaths[is.finite(isobaths)]
      
      incProgress(0.9, detail = "Composing map…")
      list(
        bathy_df    = payload$bathy_df,
        sea_df      = payload$sea_df,
        land_clip   = land_clip,
        coast_clip  = coast_clip,
        isobaths    = isobaths,
        bbox        = bb,
        res_used    = payload$res_used
      )
    })
  })
  
  # Preview plot
  output$map_plot <- renderPlot({
    md <- map_data(); req(md)
    
    pal_custom <- NULL
    if (input$palette == "Custom") {
      pal_custom <- c(input$c1, input$c2, input$c3, input$c4, input$c5, input$c6)
      pal_custom <- pal_custom[nzchar(pal_custom)]
    }
    
    opts <- list(
      # theme selector
      map_theme = input$map_theme,
      
      # aesthetic controls (Modern)
      palette = input$palette, pal_n = input$pal_n, pal_reverse = input$pal_reverse,
      pal_custom = pal_custom,
      depth_trans = input$depth_trans, interpolate = input$interpolate,
      land_fill = input$land_fill, coast_color = input$coast_color, coast_width = input$coast_width,
      show_contours = input$show_contours, contour_color = input$contour_color,
      contour_width = input$contour_width, contour_alpha = input$contour_alpha,
      legend_pos = input$legend_pos, legend_title = input$legend_title,
      legend_reverse = input$legend_reverse, base_size = input$base_size,
      grid_color = input$grid_color, grid_width = input$grid_width,
      
      # labels & meta
      title = input$title, subtitle = input$subtitle, caption = input$caption,
      scale_loc = input$scale_loc, scale_text_cex = input$scale_text_cex, scale_linewidth = input$scale_linewidth,
      north_loc = input$north_loc,
      show_cities = input$show_cities, city_pt_size = input$city_pt_size, city_lab_size = input$city_lab_size,
      show_islands = input$show_islands, island_lab_size = input$island_lab_size,
      cities_df = cities_default, islands_df = islands_default
    )
    
    build_map(md$sea_df, md$bathy_df, md$land_clip, md$coast_clip, md$isobaths, md$bbox, opts)
  }, res = 120)
  
  # Status/debug
  output$status <- renderPrint({
    bb <- bbox_reactive()
    epsg <- utm_epsg_from_lonlat((bb$xmin + bb$xmax)/2, (bb$ymin + bb$ymax)/2)
    md <- isolate({ if (!is.null(map_data())) map_data() else NULL })
    list(
      bbox = list(xmin = bb$xmin, xmax = bb$xmax, ymin = bb$ymin, ymax = bb$ymax),
      utm_epsg = epsg,
      bathy_resolution_used_arcmin = if (!is.null(md$res_used)) md$res_used else last_used_res(),
      map_theme = input$map_theme,
      notes = "NOAA fetch: retries + fallback; crop → UTM → segmentize(m) → back to EPSG:4326"
    )
  })
  
  # Download: high-res PNG using ragg
  output$download_png <- downloadHandler(
    filename = function() paste0(ifelse(nzchar(input$file_base), input$file_base, "custom_map"), ".png"),
    content = function(file) {
      md <- map_data(); req(md)
      
      pal_custom <- NULL
      if (input$palette == "Custom") {
        pal_custom <- c(input$c1, input$c2, input$c3, input$c4, input$c5, input$c6)
        pal_custom <- pal_custom[nzchar(pal_custom)]
      }
      
      opts <- list(
        map_theme = input$map_theme,
        palette = input$palette, pal_n = input$pal_n, pal_reverse = input$pal_reverse,
        pal_custom = pal_custom,
        depth_trans = input$depth_trans, interpolate = input$interpolate,
        land_fill = input$land_fill, coast_color = input$coast_color, coast_width = input$coast_width,
        show_contours = input$show_contours, contour_color = input$contour_color,
        contour_width = input$contour_width, contour_alpha = input$contour_alpha,
        legend_pos = input$legend_pos, legend_title = input$legend_title,
        legend_reverse = input$legend_reverse, base_size = input$base_size,
        grid_color = input$grid_color, grid_width = input$grid_width,
        title = input$title, subtitle = input$subtitle, caption = input$caption,
        scale_loc = input$scale_loc, scale_text_cex = input$scale_text_cex, scale_linewidth = input$scale_linewidth,
        north_loc = input$north_loc,
        show_cities = input$show_cities, city_pt_size = input$city_pt_size, city_lab_size = input$city_lab_size,
        show_islands = input$show_islands, island_lab_size = input$island_lab_size,
        cities_df = cities_default, islands_df = islands_default
      )
      
      g <- build_map(md$sea_df, md$bathy_df, md$land_clip, md$coast_clip, md$isobaths, md$bbox, opts)
      
      ragg::agg_png(
        filename = file,
        width = input$out_w, height = input$out_h, units = "in",
        res = input$out_dpi, background = "white"
      )
      print(g); dev.off()
    }
  )
}

# 4) Run ---------------------------------------------------------------

shinyApp(ui, server)

