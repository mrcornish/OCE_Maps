# Pacific Map Lab — public web edition
# Lightweight Shiny + Leaflet edition for public hosting on Render.

options(
  shiny.maxRequestSize = 30 * 1024^2,
  scipen = 999
)

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(leaflet)
  library(htmltools)
  library(jsonlite)
})

APP_VERSION <- "3.0.0-web"

PACIFIC_LIMITS <- list(west = -140, east = -68, south = -60, north = 62)

PRESETS <- list(
  "Coastal California" = list(west = -125.5, east = -116.2, south = 31.3, north = 42.5),
  "Pacific Northwest" = list(west = -130.5, east = -121.0, south = 41.5, north = 51.2),
  "Baja California" = list(west = -119.5, east = -108.0, south = 22.0, north = 33.2),
  "Tropical Eastern Pacific" = list(west = -95.0, east = -76.0, south = -6.0, north = 16.0),
  "Peru margin" = list(west = -84.5, east = -69.0, south = -20.5, north = 1.0),
  "Central Chile" = list(west = -79.0, east = -68.0, south = -46.0, north = -27.0),
  "Pacific Americas overview" = list(west = -140.0, east = -68.0, south = -60.0, north = 62.0)
)

REFERENCE_POINTS <- data.frame(
  name = c(
    "Cape Mendocino", "San Francisco", "Monterey Canyon", "Point Conception",
    "Pismo Beach", "Los Angeles", "San Diego", "Bahía San Quintín",
    "Punta Eugenia", "Vancouver Island", "Gulf of Panama", "Galápagos",
    "Lima", "Pisco", "Valparaíso", "Concepción", "Chiloé Island"
  ),
  lon = c(
    -124.41, -122.4194, -122.08, -120.47,
    -120.6413, -118.2437, -117.1611, -115.98,
    -115.05, -125.30, -79.80, -90.60,
    -77.04, -76.22, -71.62, -73.05, -73.80
  ),
  lat = c(
    40.44, 37.7749, 36.72, 34.45,
    35.1428, 34.0522, 32.7157, 30.48,
    27.84, 49.60, 8.00, -0.70,
    -12.05, -13.71, -33.05, -36.83, -42.60
  ),
  kind = c(
    "Cape", "City", "Submarine canyon", "Cape",
    "Beach", "City", "City", "Bay",
    "Cape", "Island", "Gulf", "Archipelago",
    "City", "City", "City", "City", "Island"
  ),
  stringsAsFactors = FALSE
)

haversine_km <- function(lon1, lat1, lon2, lat2) {
  radius <- 6371.0088
  to_rad <- pi / 180
  lon1 <- lon1 * to_rad
  lon2 <- lon2 * to_rad
  lat1 <- lat1 * to_rad
  lat2 <- lat2 * to_rad
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  2 * radius * asin(pmin(1, sqrt(a)))
}

safe_hex <- function(x, fallback = "#38C8B5") {
  x <- trimws(as.character(x %||% ""))
  if (grepl("^#[0-9A-Fa-f]{6}$", x)) x else fallback
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

read_site_file <- function(path, filename) {
  extension <- tolower(tools::file_ext(filename))
  dat <- if (extension %in% c("tsv", "txt")) {
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  if (!nrow(dat)) stop("The uploaded table contains no rows.")

  original_names <- names(dat)
  normalized <- tolower(gsub("[^a-z0-9]+", "_", original_names))

  find_column <- function(candidates) {
    hit <- match(candidates, normalized, nomatch = 0)
    hit <- hit[hit > 0]
    if (length(hit)) original_names[hit[1]] else NULL
  }

  lon_col <- find_column(c("lon", "longitude", "long", "x", "decimal_longitude"))
  lat_col <- find_column(c("lat", "latitude", "y", "decimal_latitude"))
  name_col <- find_column(c("name", "site", "station", "label", "id", "location"))
  group_col <- find_column(c("group", "category", "type", "region", "project"))

  if (is.null(lon_col) || is.null(lat_col)) {
    stop("Include longitude and latitude columns, such as lon/lat or longitude/latitude.")
  }

  out <- data.frame(
    name = if (is.null(name_col)) paste("Site", seq_len(nrow(dat))) else as.character(dat[[name_col]]),
    lon = suppressWarnings(as.numeric(dat[[lon_col]])),
    lat = suppressWarnings(as.numeric(dat[[lat_col]])),
    group = if (is.null(group_col)) "Uploaded site" else as.character(dat[[group_col]]),
    stringsAsFactors = FALSE
  )

  out$name[is.na(out$name) | !nzchar(trimws(out$name))] <- paste("Site", which(is.na(out$name) | !nzchar(trimws(out$name))))
  out$group[is.na(out$group) | !nzchar(trimws(out$group))] <- "Uploaded site"
  out <- out[is.finite(out$lon) & is.finite(out$lat), , drop = FALSE]
  out <- out[
    out$lon >= PACIFIC_LIMITS$west & out$lon <= PACIFIC_LIMITS$east &
      out$lat >= PACIFIC_LIMITS$south & out$lat <= PACIFIC_LIMITS$north,
    , drop = FALSE
  ]

  if (!nrow(out)) stop("No valid points remained inside the Pacific Americas map window.")
  rownames(out) <- NULL
  out
}

app_theme <- bs_theme(
  version = 5,
  bg = "#06111C",
  fg = "#E8F4F8",
  primary = "#35BDD2",
  secondary = "#38C8B5",
  success = "#38C8B5",
  warning = "#F2B66D",
  danger = "#D16C7A",
  base_font = font_collection("Inter", "Avenir Next", "Helvetica Neue", "Arial", "sans-serif"),
  heading_font = font_collection("Inter", "Avenir Next", "Helvetica Neue", "Arial", "sans-serif"),
  "border-radius" = "0.9rem",
  "input-bg" = "#0B1C2A",
  "input-color" = "#E8F4F8",
  "input-border-color" = "rgba(140,210,227,.24)"
)

app_css <- tags$style(HTML("
  :root {
    --pm-bg: #06111c;
    --pm-panel: rgba(9, 27, 42, .88);
    --pm-panel-2: rgba(13, 39, 58, .86);
    --pm-line: rgba(140, 210, 227, .18);
    --pm-text: #e8f4f8;
    --pm-muted: #9fb8c6;
    --pm-cyan: #67e0ee;
    --pm-blue: #35bdd2;
    --pm-green: #38c8b5;
    --pm-lime: #b5f76f;
  }

  html, body { height: 100%; }

  body {
    margin: 0;
    background:
      radial-gradient(circle at 10% 0%, rgba(53,189,210,.14), transparent 28%),
      radial-gradient(circle at 92% 6%, rgba(181,247,111,.08), transparent 24%),
      linear-gradient(180deg, #06111c, #081827 48%, #06111c);
    letter-spacing: .006em;
  }

  .pm-shell { min-height: 100vh; }

  .pm-header {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 1rem;
    padding: .75rem 1rem;
    border: 1px solid var(--pm-line);
    border-radius: 1rem;
    background: rgba(6, 17, 28, .78);
    backdrop-filter: blur(14px);
    box-shadow: 0 16px 42px rgba(0,0,0,.25);
  }

  .pm-brand { display: flex; align-items: center; gap: .8rem; min-width: 0; }

  .pm-mark {
    width: 2.55rem;
    height: 2.55rem;
    border: 1px solid rgba(103,224,238,.65);
    border-radius: 50%;
    display: grid;
    place-items: center;
    color: var(--pm-cyan);
    font-size: .72rem;
    font-weight: 900;
    letter-spacing: .09em;
    box-shadow: inset 0 0 18px rgba(53,189,210,.14), 0 0 22px rgba(53,189,210,.12);
    flex: 0 0 auto;
  }

  .pm-title { font-weight: 850; letter-spacing: .08em; line-height: 1; white-space: nowrap; }
  .pm-subtitle { color: var(--pm-muted); font-size: .68rem; letter-spacing: .12em; margin-top: .28rem; }

  .pm-live {
    display: inline-flex;
    align-items: center;
    gap: .5rem;
    color: var(--pm-muted);
    font-size: .78rem;
    white-space: nowrap;
  }

  .pm-dot {
    width: .58rem;
    height: .58rem;
    border-radius: 50%;
    background: var(--pm-green);
    box-shadow: 0 0 0 .28rem rgba(56,200,181,.10);
  }

  .bslib-sidebar-layout > .sidebar {
    background: linear-gradient(180deg, var(--pm-panel), var(--pm-panel-2));
    border: 1px solid var(--pm-line) !important;
    border-radius: 1rem;
    box-shadow: 0 16px 42px rgba(0,0,0,.23);
    backdrop-filter: blur(12px);
  }

  .card {
    border: 1px solid var(--pm-line) !important;
    background: rgba(5, 15, 25, .82) !important;
    box-shadow: 0 18px 50px rgba(0,0,0,.28);
    overflow: hidden;
  }

  .pm-map-card { min-height: 720px; }
  .pm-map-card .card-body { padding: 0 !important; }

  .pm-map-toolbar {
    display: flex;
    justify-content: space-between;
    align-items: center;
    gap: 1rem;
    padding: .65rem .85rem;
    color: var(--pm-muted);
    font-size: .78rem;
    border-bottom: 1px solid var(--pm-line);
    background: rgba(6,17,28,.94);
  }

  .leaflet { background: #071522 !important; }
  .leaflet-control-layers, .leaflet-popup-content-wrapper, .leaflet-popup-tip {
    background: rgba(6,17,28,.94) !important;
    color: var(--pm-text) !important;
    border: 1px solid var(--pm-line);
  }
  .leaflet-control-layers label { color: var(--pm-text); }

  .form-label, .control-label {
    color: var(--pm-muted);
    font-size: .72rem;
    font-weight: 800;
    letter-spacing: .06em;
    text-transform: uppercase;
  }

  .form-control, .form-select, .selectize-input, .selectize-dropdown {
    background: rgba(5, 17, 28, .86) !important;
    color: var(--pm-text) !important;
    border-color: rgba(140,210,227,.20) !important;
  }

  .accordion-item { background: transparent !important; border-color: var(--pm-line) !important; }
  .accordion-button { background: transparent !important; color: var(--pm-text) !important; font-weight: 780; }

  .btn-pm {
    color: #06111c !important;
    background: linear-gradient(135deg, #67e0ee, #38c8b5) !important;
    border: 0 !important;
    font-weight: 850 !important;
    box-shadow: 0 10px 24px rgba(53,189,210,.18);
  }
  .btn-pm:hover { filter: brightness(1.05); transform: translateY(-1px); }

  .btn-pm-ghost {
    color: var(--pm-text) !important;
    background: rgba(103,224,238,.06) !important;
    border-color: rgba(103,224,238,.28) !important;
  }

  .pm-metric {
    padding: .75rem;
    border: 1px solid var(--pm-line);
    border-radius: .75rem;
    background: rgba(53,189,210,.05);
    margin-bottom: .7rem;
  }
  .pm-metric-label { color: var(--pm-muted); font-size: .66rem; letter-spacing: .1em; text-transform: uppercase; }
  .pm-metric-value { font-size: 1.1rem; font-weight: 780; margin-top: .12rem; }

  .pm-note { color: var(--pm-muted); font-size: .77rem; line-height: 1.5; }
  .pm-kicker { color: var(--pm-cyan); font-size: .66rem; letter-spacing: .14em; text-transform: uppercase; font-weight: 850; }
  .pm-footer { color: var(--pm-muted); font-size: .68rem; text-align: center; padding: .3rem 0 .1rem; }

  table { color: var(--pm-text) !important; font-size: .78rem; }

  @media (max-width: 900px) {
    .pm-subtitle { display: none; }
    .pm-live span:last-child { display: none; }
    .pm-map-card { min-height: 610px; }
  }
"))

sidebar_ui <- sidebar(
  width = 340,
  open = "always",
  div(class = "pm-kicker", "Map controls"),
  selectInput("preset", "Coastal extent", choices = names(PRESETS), selected = "Coastal California"),
  actionButton("fit_preset", "Frame selected coast", class = "btn-pm w-100", icon = icon("location-crosshairs")),
  br(), br(),
  accordion(
    id = "controls",
    open = c("Layers", "Field data"),
    accordion_panel(
      "Layers",
      checkboxInput("show_reference", "Reference places", TRUE),
      checkboxInput("show_labels", "Show point labels", TRUE),
      checkboxInput("cluster_sites", "Cluster uploaded sites", FALSE),
      selectInput(
        "marker_color", "Uploaded-site color",
        choices = c(
          "Sea glass" = "#38C8B5",
          "Electric cyan" = "#67E0EE",
          "Marine blue" = "#3383B2",
          "Signal lime" = "#B5F76F",
          "Coral" = "#F07C75",
          "Brass" = "#D2A85B",
          "White" = "#F4F8FA"
        ),
        selected = "#38C8B5"
      ),
      sliderInput("marker_size", "Marker size", min = 3, max = 14, value = 7, step = 1),
      sliderInput("marker_opacity", "Marker opacity", min = .2, max = 1, value = .9, step = .05)
    ),
    accordion_panel(
      "Field data",
      fileInput("site_file", "Upload CSV or TSV", accept = c(".csv", ".tsv", ".txt")),
      div(class = "pm-note", "Required: longitude and latitude. Optional: name and group."),
      br(),
      fluidRow(
        column(6, actionButton("fit_sites", "Fit sites", class = "btn-pm-ghost w-100")),
        column(6, actionButton("clear_sites", "Clear", class = "btn-pm-ghost w-100"))
      ),
      br(),
      downloadButton("download_sites", "Download mapped sites", class = "btn-pm-ghost w-100")
    ),
    accordion_panel(
      "Click survey",
      div(class = "pm-note", "Click the map to add survey points. Consecutive clicks report geodesic distance."),
      br(),
      actionButton("clear_clicks", "Clear clicked points", class = "btn-pm-ghost w-100"),
      br(), br(),
      downloadButton("download_clicks", "Download clicked points", class = "btn-pm-ghost w-100")
    ),
    accordion_panel(
      "Project",
      textInput("project_title", "Project title", "Coastal California"),
      textAreaInput("project_note", "Notes", "", rows = 3),
      downloadButton("download_settings", "Download settings JSON", class = "btn-pm-ghost w-100"),
      br(), br(),
      div(class = "pm-note", "Use the layer control on the map to switch between ocean bathymetry, satellite imagery, dark cartography, and terrain.")
    )
  ),
  hr(),
  uiOutput("point_metric"),
  uiOutput("click_metric")
)

ui <- page_fillable(
  theme = app_theme,
  app_css,
  div(
    class = "pm-shell",
    div(
      class = "pm-header",
      div(
        class = "pm-brand",
        div(class = "pm-mark", "PAC"),
        div(
          div(class = "pm-title", "PACIFIC MAP LAB"),
          div(class = "pm-subtitle", "WEST COAST AMERICAS · INTERACTIVE CARTOGRAPHY")
        )
      ),
      div(class = "pm-live", span(class = "pm-dot"), span("PUBLIC WEB EDITION"))
    ),
    br(),
    layout_sidebar(
      sidebar = sidebar_ui,
      card(
        class = "pm-map-card",
        full_screen = TRUE,
        div(
          class = "pm-map-toolbar",
          span(textOutput("extent_status", inline = TRUE)),
          span(textOutput("cursor_status", inline = TRUE))
        ),
        leafletOutput("map", height = "calc(100vh - 168px)")
      )
    ),
    div(class = "pm-footer", paste0("Pacific Map Lab v", APP_VERSION, " · Michael R. Cornish · Esri, OpenStreetMap, CARTO, and NOAA basemap services"))
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(
    uploaded = NULL,
    clicks = data.frame(
      name = character(), lon = numeric(), lat = numeric(),
      distance_from_previous_km = numeric(), stringsAsFactors = FALSE
    ),
    last_cursor = NULL
  )

  output$map <- renderLeaflet({
    ca <- PRESETS[["Coastal California"]]

    leaflet(options = leafletOptions(preferCanvas = TRUE, minZoom = 2, worldCopyJump = FALSE)) |>
      addProviderTiles(providers$Esri.OceanBasemap, group = "Ocean / bathymetry") |>
      addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") |>
      addProviderTiles(providers$CartoDB.DarkMatterNoLabels, group = "Dark cartography") |>
      addProviderTiles(providers$CartoDB.PositronNoLabels, group = "Light cartography") |>
      addProviderTiles(providers$Esri.WorldTopoMap, group = "Terrain") |>
      hideGroup(c("Satellite", "Dark cartography", "Light cartography", "Terrain")) |>
      addLayersControl(
        baseGroups = c("Ocean / bathymetry", "Satellite", "Dark cartography", "Light cartography", "Terrain"),
        overlayGroups = c("Reference places", "Uploaded sites", "Clicked points"),
        options = layersControlOptions(collapsed = FALSE)
      ) |>
      addScaleBar(position = "bottomleft", options = scaleBarOptions(metric = TRUE, imperial = TRUE)) |>
      fitBounds(ca$west, ca$south, ca$east, ca$north)
  })

  observe({
    proxy <- leafletProxy("map", session = session) |>
      clearGroup("Reference places")

    if (isTRUE(input$show_reference)) {
      labels <- if (isTRUE(input$show_labels)) REFERENCE_POINTS$name else NULL
      proxy |>
        addCircleMarkers(
          data = REFERENCE_POINTS,
          lng = ~lon, lat = ~lat,
          radius = 4.5,
          stroke = TRUE,
          weight = 1,
          color = "#E8F4F8",
          fillColor = "#35BDD2",
          fillOpacity = .82,
          label = labels,
          popup = ~paste0("<strong>", htmlEscape(name), "</strong><br>", htmlEscape(kind), "<br>", sprintf("%.4f, %.4f", lat, lon)),
          group = "Reference places"
        )
    }
  })

  observe({
    proxy <- leafletProxy("map", session = session) |>
      clearGroup("Uploaded sites")

    dat <- rv$uploaded
    if (!is.null(dat) && nrow(dat)) {
      labels <- if (isTRUE(input$show_labels)) dat$name else NULL
      cluster <- if (isTRUE(input$cluster_sites)) markerClusterOptions() else NULL
      proxy |>
        addCircleMarkers(
          data = dat,
          lng = ~lon, lat = ~lat,
          radius = input$marker_size %||% 7,
          stroke = TRUE,
          weight = 1.2,
          color = "#06111C",
          fillColor = safe_hex(input$marker_color),
          fillOpacity = input$marker_opacity %||% .9,
          label = labels,
          popup = ~paste0("<strong>", htmlEscape(name), "</strong><br>", htmlEscape(group), "<br>", sprintf("%.5f, %.5f", lat, lon)),
          clusterOptions = cluster,
          group = "Uploaded sites"
        )
    }
  })

  observe({
    proxy <- leafletProxy("map", session = session) |>
      clearGroup("Clicked points")

    dat <- rv$clicks
    if (nrow(dat)) {
      proxy |>
        addCircleMarkers(
          data = dat,
          lng = ~lon, lat = ~lat,
          radius = 6,
          stroke = TRUE,
          weight = 2,
          color = "#06111C",
          fillColor = "#B5F76F",
          fillOpacity = .95,
          label = ~name,
          popup = ~paste0(
            "<strong>", htmlEscape(name), "</strong><br>",
            sprintf("%.5f, %.5f", lat, lon),
            ifelse(is.finite(distance_from_previous_km), paste0("<br>", sprintf("%.2f km from previous point", distance_from_previous_km)), "")
          ),
          group = "Clicked points"
        )
    }
  })

  observeEvent(input$fit_preset, {
    preset <- PRESETS[[input$preset]] %||% PRESETS[["Coastal California"]]
    leafletProxy("map", session = session) |>
      fitBounds(preset$west, preset$south, preset$east, preset$north)
    updateTextInput(session, "project_title", value = input$preset)
  })

  observeEvent(input$preset, {
    preset <- PRESETS[[input$preset]] %||% PRESETS[["Coastal California"]]
    leafletProxy("map", session = session) |>
      fitBounds(preset$west, preset$south, preset$east, preset$north)
    updateTextInput(session, "project_title", value = input$preset)
  }, ignoreInit = TRUE)

  observeEvent(input$site_file, {
    req(input$site_file)
    tryCatch({
      rv$uploaded <- read_site_file(input$site_file$datapath, input$site_file$name)
      showNotification(paste(nrow(rv$uploaded), "sites added to the map."), type = "message")
    }, error = function(e) {
      rv$uploaded <- NULL
      showNotification(conditionMessage(e), type = "error", duration = 10)
    })
  })

  observeEvent(input$clear_sites, {
    rv$uploaded <- NULL
    showNotification("Uploaded sites cleared.", type = "message")
  })

  observeEvent(input$fit_sites, {
    dat <- rv$uploaded
    if (is.null(dat) || !nrow(dat)) {
      showNotification("Upload site data first.", type = "warning")
      return()
    }
    lon_pad <- max(.08, diff(range(dat$lon)) * .12)
    lat_pad <- max(.08, diff(range(dat$lat)) * .12)
    leafletProxy("map", session = session) |>
      fitBounds(
        min(dat$lon) - lon_pad, min(dat$lat) - lat_pad,
        max(dat$lon) + lon_pad, max(dat$lat) + lat_pad
      )
  })

  observeEvent(input$map_click, {
    click <- input$map_click
    req(click$lng, click$lat)

    old <- rv$clicks
    distance <- NA_real_
    if (nrow(old)) {
      previous <- old[nrow(old), ]
      distance <- haversine_km(previous$lon, previous$lat, click$lng, click$lat)
    }

    next_point <- data.frame(
      name = paste("Survey point", nrow(old) + 1L),
      lon = as.numeric(click$lng),
      lat = as.numeric(click$lat),
      distance_from_previous_km = distance,
      stringsAsFactors = FALSE
    )
    rv$clicks <- rbind(old, next_point)
    rv$last_cursor <- click

    if (is.finite(distance)) {
      showNotification(sprintf("Point added · %.2f km from the previous point", distance), type = "message", duration = 4)
    } else {
      showNotification("First survey point added.", type = "message", duration = 3)
    }
  })

  observeEvent(input$map_mouseout, {
    rv$last_cursor <- NULL
  })

  observeEvent(input$clear_clicks, {
    rv$clicks <- rv$clicks[0, , drop = FALSE]
    showNotification("Clicked survey points cleared.", type = "message")
  })

  output$extent_status <- renderText({
    paste0("EXTENT · ", toupper(input$preset %||% "Coastal California"))
  })

  output$cursor_status <- renderText({
    n <- nrow(rv$clicks)
    if (!n) "CLICK MAP TO START A SURVEY" else paste(n, if (n == 1) "SURVEY POINT" else "SURVEY POINTS")
  })

  output$point_metric <- renderUI({
    n <- if (is.null(rv$uploaded)) 0L else nrow(rv$uploaded)
    div(
      class = "pm-metric",
      div(class = "pm-metric-label", "Uploaded field sites"),
      div(class = "pm-metric-value", format(n, big.mark = ","))
    )
  })

  output$click_metric <- renderUI({
    dat <- rv$clicks
    total <- if (nrow(dat) > 1) sum(dat$distance_from_previous_km, na.rm = TRUE) else 0
    div(
      class = "pm-metric",
      div(class = "pm-metric-label", "Survey path"),
      div(class = "pm-metric-value", sprintf("%.1f km", total))
    )
  })

  output$download_sites <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_-]+", "_", input$project_title %||% "pacific_map"), "_sites.csv"),
    content = function(file) {
      dat <- rv$uploaded
      if (is.null(dat)) dat <- data.frame(name = character(), lon = numeric(), lat = numeric(), group = character())
      utils::write.csv(dat, file, row.names = FALSE)
    }
  )

  output$download_clicks <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_-]+", "_", input$project_title %||% "pacific_map"), "_survey_points.csv"),
    content = function(file) utils::write.csv(rv$clicks, file, row.names = FALSE)
  )

  output$download_settings <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_-]+", "_", input$project_title %||% "pacific_map"), "_settings.json"),
    content = function(file) {
      settings <- list(
        application = "Pacific Map Lab",
        version = APP_VERSION,
        created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
        project_title = input$project_title,
        notes = input$project_note,
        preset = input$preset,
        uploaded_site_count = if (is.null(rv$uploaded)) 0L else nrow(rv$uploaded),
        clicked_point_count = nrow(rv$clicks),
        marker_color = input$marker_color,
        marker_size = input$marker_size,
        marker_opacity = input$marker_opacity,
        show_reference = isTRUE(input$show_reference),
        show_labels = isTRUE(input$show_labels)
      )
      writeLines(toJSON(settings, pretty = TRUE, auto_unbox = TRUE, null = "null"), file)
    }
  )
}

shinyApp(ui, server)
