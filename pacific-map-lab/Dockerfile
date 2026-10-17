FROM rocker/geospatial:4.5.1

LABEL org.opencontainers.image.title="Pacific Map Lab"
LABEL org.opencontainers.image.description="Interactive coastal cartography for the Pacific margin of the Americas"

RUN install2.r --error --skipinstalled --ncpus -1 \
    shiny bslib htmltools marmap rnaturalearth rnaturalearthdata \
    ggspatial ggrepel colourpicker leaflet DT cowplot ragg svglite \
    jsonlite zip viridisLite scales readr tibble dplyr units terra sf ggplot2

WORKDIR /app

# The complete Shiny source is stored as numbered text chunks so it can be
# managed through GitHub's text-file API. Reassemble it during the image build.
COPY app.tar.gz.b64.part* /tmp/
RUN cat /tmp/app.tar.gz.b64.part* | base64 -d > /tmp/app.tar.gz \
    && tar -xzf /tmp/app.tar.gz -C /app \
    && test -s /app/app.R \
    && rm -f /tmp/app.tar.gz /tmp/app.tar.gz.b64.part*

ENV PORT=10000
ENV PACIFIC_MAP_CACHE_DIR=/tmp/PacificMapLab-cache
EXPOSE 10000

CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=as.numeric(Sys.getenv('PORT', '10000')), launch.browser=FALSE)"]
