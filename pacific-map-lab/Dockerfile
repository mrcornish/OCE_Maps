FROM rocker/shiny:4.5.1

LABEL org.opencontainers.image.title="Pacific Map Lab"
LABEL org.opencontainers.image.description="Reproducible ggplot2 and sf coastal cartography for the Pacific Americas"

RUN install2.r --error --skipinstalled --ncpus -1 \
    bslib htmltools leaflet ggplot2 sf marmap rnaturalearth \
    rnaturalearthdata ggspatial ggrepel dplyr tibble readr \
    colourpicker ragg svglite scales jsonlite zip viridisLite units

WORKDIR /app
COPY app.R /app/app.R

ENV PORT=10000
ENV PACIFIC_MAP_CACHE_DIR=/tmp/PacificMapLab-cache
EXPOSE 10000

CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=as.numeric(Sys.getenv('PORT', '10000')), launch.browser=FALSE)"]
