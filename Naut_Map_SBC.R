# Ultra‑detail vintage nautical chart (ggplot) — Santa Barbara Channel
# Robust to GDAL/sf differences; no wkt_filter; safe fallbacks.
# Produces: santa_barbara_channel_vintage_chart.png

# 0) Packages
req <- c("marmap","ggplot2","sf","ggspatial","ggrepel","dplyr","metR",
         "scales","grid","tibble","rnaturalearth","rnaturalearthdata")
new <- setdiff(req, rownames(installed.packages()))
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
lapply(req, library, character.only = TRUE)

# Turn off s2 for robustness with boundaries/crops on some builds
if (sf::sf_use_s2()) sf::sf_use_s2(FALSE)

# 1) Extent (Santa Barbara Channel)
set.seed(123)
lon_min <- -120.9; lon_max <- -119.0
lat_min <-  33.8;  lat_max <-  34.7
bbox_sf <- sf::st_bbox(c(xmin = lon_min, ymin = lat_min,
                         xmax = lon_max, ymax = lat_max),
                       crs = sf::st_crs(4326))

# 2) Bathymetry (15″ sampling via NOAA ETOPO 2022 15s)
options(timeout = max(600, getOption("timeout", 60)))
bathy <- marmap::getNOAA.bathy(lon1 = lon_min, lon2 = lon_max,
                               lat1 = lat_min, lat2 = lat_max,
                               resolution = 0.25)
bathy_df <- marmap::fortify.bathy(bathy) |>
  dplyr::rename(lon = x, lat = y, z = z)
sea_df <- bathy_df |>
  dplyr::filter(z <= 0) |>
  dplyr::mutate(depth_m = -z, depth_fm = depth_m / 1.8288)

# 3) HIGH‑RES shoreline: try GSHHG (full) safely; fall back to NE 1:10m
get_gshhg_safe <- function(bbox_sf) {
  urls <- c(
    "https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-2.3.7.zip",
    "https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/2.3.7/gshhg-shp-2.3.7.zip",
    "https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/2.3.6/gshhg-shp-2.3.6.zip"
  )
  dl <- tempfile(fileext = ".zip")
  ok <- FALSE
  for (u in urls) {
    try({
      utils::download.file(u, dl, mode = "wb", quiet = TRUE)
      if (file.exists(dl) && file.info(dl)$size > 1e6) { ok <- TRUE; break }
    }, silent = TRUE)
  }
  if (!ok) return(NULL)
  
  exdir <- file.path(tempdir(), "gshhg")
  if (!dir.exists(exdir)) dir.create(exdir, recursive = TRUE)
  utils::unzip(dl, exdir = exdir)
  
  shp_poly <- file.path(exdir, "GSHHS_shp", "f", "GSHHS_f_L1.shp")
  if (!file.exists(shp_poly)) return(NULL)
  
  # Plain read (no wkt_filter); then robust clean/cast/crop
  land_raw <- try(sf::st_read(shp_poly, quiet = TRUE), silent = TRUE)
  if (inherits(land_raw, "try-error")) return(NULL)
  
  # Make valid and extract polygons strictly (avoid GEOMETRYCOLLECTION issues)
  land_geom <- try({
    land_valid <- suppressWarnings(sf::st_make_valid(land_raw))
    # Extract only POLYGON parts (drop other geometry types)
    poly_only  <- suppressWarnings(sf::st_collection_extract(sf::st_geometry(land_valid), "POLYGON"))
    # Build sf; ensure CRS
    sf::st_sf(geometry = poly_only, crs = sf::st_crs(4326))
  }, silent = TRUE)
  if (inherits(land_geom, "try-error") || length(land_geom) == 0) return(NULL)
  
  land_clip <- try(suppressWarnings(sf::st_crop(land_geom, bbox_sf)), silent = TRUE)
  if (inherits(land_clip, "try-error") || is.null(land_clip) || nrow(land_clip) == 0) return(NULL)
  
  # Coastline as crisp lines
  coast <- try({
    # Use boundary of union to avoid internal borders
    bnd <- sf::st_boundary(sf::st_union(sf::st_geometry(land_clip)))
    # Ensure we have LINESTRING/MULTILINESTRING
    ln  <- suppressWarnings(sf::st_collection_extract(bnd, "LINESTRING"))
    sf::st_sf(geometry = ln, crs = sf::st_crs(4326))
  }, silent = TRUE)
  if (inherits(coast, "try-error") || is.null(coast) || length(coast) == 0) {
    # Fallback: simple boundary per feature
    bnd <- sf::st_boundary(sf::st_geometry(land_clip))
    ln  <- suppressWarnings(sf::st_collection_extract(bnd, "LINESTRING"))
    coast <- sf::st_sf(geometry = ln, crs = sf::st_crs(4326))
  }
  list(land = land_clip, coast = coast)
}

g <- try(get_gshhg_safe(bbox_sf), silent = TRUE)

if (inherits(g, "try-error") || is.null(g)) {
  message("GSHHG unavailable; using Natural Earth 1:10m as fallback.")
  land_hi <- rnaturalearth::ne_download(scale = "large", type = "land",
                                        category = "physical", returnclass = "sf")
  land <- suppressWarnings(sf::st_crop(land_hi, bbox_sf))
  coast <- suppressWarnings(sf::st_crop(rnaturalearth::ne_coastline(scale = "large", returnclass = "sf"), bbox_sf))
  # Clean just in case
  land  <- suppressWarnings(sf::st_make_valid(land))
  coast <- suppressWarnings(sf::st_make_valid(coast))
} else {
  land  <- g$land
  coast <- g$coast
}

# 4) Paper‑grain underlay (subtle)
make_paper_grain <- function(lon_min, lon_max, lat_min, lat_max, nx = 1000, ny = 1000) {
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
grain_df <- make_paper_grain(lon_min, lon_max, lat_min, lat_max)

# 5) Labels
cities <- tibble::tibble(
  name = c("Santa Barbara", "Goleta", "Ventura"),
  lon  = c(-119.6982, -119.8276, -119.2290),
  lat  = c( 34.4208,   34.4358,   34.2746)
)
islands <- tibble::tibble(
  name = c("San Miguel I.", "Santa Rosa I.", "Santa Cruz I.", "Anacapa I."),
  lon  = c(-120.37,      -120.08,        -119.77,        -119.36),
  lat  = c(  34.04,        33.95,          34.03,          34.01)
)

# 6) Fathom isobaths & palette
iso_fm_shallow <- c(3, 6, 10, 20)
iso_fm_deep    <- c(50, 100)
fm_to_m        <- function(fm) -(fm * 1.8288)
iso_shallow_m  <- fm_to_m(iso_fm_shallow)
iso_deep_m     <- fm_to_m(iso_fm_deep)
water_cols     <- c("#bfe3f1","#d4edf7","#ecf5fa","#ffffff")
max_fm         <- max(sea_df$depth_fm, na.rm = TRUE)

# 7) Nautical‑miles scale bar
make_nmi_scalebar <- function(lon0, lat0, nmi = 5, segments = 5,
                              height_deg = 0.025, label = "5 nmi") {
  lon_per_nmi <- (1/60) / cos(lat0 * pi/180)
  rects <- tibble::tibble(
    xmin = lon0 + (0:(segments-1)) * lon_per_nmi,
    xmax = lon0 + (1:segments)     * lon_per_nmi,
    ymin = lat0,
    ymax = lat0 + height_deg,
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
sb_lat <- lat_min + 0.08
sb_lon <- lon_min + 0.12
sb <- make_nmi_scalebar(sb_lon, sb_lat, nmi = 5, segments = 5, height_deg = 0.03,
                        label = "5 nautical miles")

# 8) Neatlines
frame <- tibble::tibble(xmin = lon_min, xmax = lon_max, ymin = lat_min, ymax = lat_max)
frame_inner <- tibble::tibble(xmin = lon_min + 0.01, xmax = lon_max - 0.01,
                              ymin = lat_min + 0.01, ymax = lat_max - 0.01)


# 9) Degree–minute axes  (ASCII-safe: use "'" instead of the Unicode prime)
to_dm_lon <- function(x){
  deg  <- floor(abs(x))
  mins <- round((abs(x) - deg) * 60)
  if (any(mins == 60)) { deg[mins == 60] <- deg[mins == 60] + 1; mins[mins == 60] <- 0 }
  paste0(deg, "\u00B0", sprintf("%02d", mins), "'", ifelse(x < 0, "W", "E"))
}
to_dm_lat <- function(x){
  deg  <- floor(abs(x))
  mins <- round((abs(x) - deg) * 60)
  if (any(mins == 60)) { deg[mins == 60] <- deg[mins == 60] + 1; mins[mins == 60] <- 0 }
  paste0(deg, "\u00B0", sprintf("%02d", mins), "'", ifelse(x < 0, "S", "N"))
}



# 10) Plot
p <- ggplot() +
  # Paper grain
  geom_raster(data = grain_df, aes(lon, lat, alpha = grain),
              inherit.aes = FALSE, fill = "#3b2f20") +
  scale_alpha(range = c(0.01, 0.05), guide = "none") +
  # Water
  geom_raster(data = sea_df, aes(lon, lat, fill = depth_fm), interpolate = TRUE) +
  # Isobaths + labels (fathoms)
  geom_contour(data = bathy_df, aes(lon, lat, z = z),
               breaks = iso_shallow_m, linewidth = 0.35, color = "#2c6486", alpha = 0.9) +
  geom_contour(data = bathy_df, aes(lon, lat, z = z),
               breaks = iso_deep_m,    linewidth = 0.25, color = "#6f8ea3", alpha = 0.8) +
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
  # Land & coastline (high‑res)
  geom_sf(data = land,  fill = "#eae2cf", color = NA) +
  geom_sf(data = coast, color = "#2f2c24", linewidth = 0.55,
          lineend = "round", linejoin = "round") +
  # Labels
  geom_point(data = cities, aes(lon, lat), size = 1.8, color = "#333333") +
  ggrepel::geom_text_repel(data = cities, aes(lon, lat, label = name),
                           size = 3.2, family = "serif", seed = 42, min.segment.length = 0) +
  ggrepel::geom_text_repel(data = islands, aes(lon, lat, label = name),
                           size = 3.0, family = "serif", seed = 42, min.segment.length = 0,
                           box.padding = 0.3, color = "#2b2b2b") +
  # Nautical‑miles scale bar
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
    location = "tr", which_north = "true",
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  # Neatlines
  annotate("rect", xmin = frame$xmin, xmax = frame$xmax, ymin = frame$ymin, ymax = frame$ymax,
           fill = NA, color = "#2e2c24", linewidth = 0.7) +
  annotate("rect", xmin = frame_inner$xmin, xmax = frame_inner$xmax,
           ymin = frame_inner$ymin, ymax = frame_inner$ymax,
           fill = NA, color = "#2e2c24", linewidth = 0.4) +
  # Axes, palette, theme
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max),
           expand = FALSE, crs = sf::st_crs(4326)) +
  scale_fill_gradientn(colours = water_cols,
                       values  = scales::rescale(c(0, 5, 10, 30, max_fm)),
                       guide   = "none") +
  scale_x_continuous(breaks = lon_breaks, labels = to_dm_lon) +
  scale_y_continuous(breaks = lat_breaks, labels = to_dm_lat) +
  labs(
    title    = "Santa Barbara Channel",
    subtitle = "Vintage-style nautical chart • High-resolution shoreline • Depths in fathoms",
    caption  = "Bathymetry: NOAA ETOPO 2022 • Code: https://github.com/mrcornish/OCE_Maps"
  ) +
  theme_minimal(base_size = 12) +
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

# 11) Save & show
outfile <- "santa_barbara_channel_vintage_chart.png"
print(p)
message("Saved: ", normalizePath(outfile))
