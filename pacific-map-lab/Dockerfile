FROM rocker/geospatial:4.5.1

LABEL org.opencontainers.image.title="Pacific Map Lab"
LABEL org.opencontainers.image.description="Reproducible ggplot2 and sf coastal cartography for the Pacific Americas"

RUN install2.r --error --skipinstalled --ncpus -1 \
    bslib htmltools leaflet marmap rnaturalearth rnaturalearthdata \
    ggspatial ggrepel colourpicker ragg svglite jsonlite zip viridisLite

WORKDIR /app
COPY app4.R.gz.b64.part* /tmp/
RUN cat /tmp/app4.R.gz.b64.part* | base64 -d | gzip -dc > /app/app.R \
    && test -s /app/app.R \
    && rm -f /tmp/app4.R.gz.b64.part*

ENV PORT=10000
ENV PACIFIC_MAP_CACHE_DIR=/tmp/PacificMapLab-cache
EXPOSE 10000

CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=as.numeric(Sys.getenv('PORT', '10000')), launch.browser=FALSE)"]
