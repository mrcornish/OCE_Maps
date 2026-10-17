FROM rocker/geospatial:4.5.1

LABEL org.opencontainers.image.title="Pacific Map Lab"
LABEL org.opencontainers.image.description="Reproducible ggplot2 and sf coastal cartography for the Pacific Americas"

RUN apt-get update \
    && apt-get install -y --no-install-recommends patch \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled --ncpus -1 \
    bslib rnaturalearth rnaturalearthdata ggspatial ggrepel \
    colourpicker ragg svglite jsonlite zip readr

WORKDIR /app
COPY app4.R.gz.b64.part* /tmp/
COPY optimize_app.patch.part* /tmp/
RUN cat /tmp/app4.R.gz.b64.part* | base64 -d | gzip -dc > /app/app.R \
    && cat /tmp/optimize_app.patch.part* > /tmp/optimize_app.patch \
    && cd /app \
    && patch --batch --silent -p0 < /tmp/optimize_app.patch \
    && sed -i 's/+ coord_cartesian(/+ ggplot2::coord_cartesian(/g; s/+ theme_void()/+ ggplot2::theme_void()/g' /app/app.R \
    && Rscript -e "invisible(parse(file='/app/app.R'))" \
    && test -s /app/app.R \
    && rm -f /tmp/app4.R.gz.b64.part* /tmp/optimize_app.patch.part* /tmp/optimize_app.patch

ENV PORT=10000
ENV PACIFIC_MAP_CACHE_DIR=/tmp/PacificMapLab-cache
EXPOSE 10000

CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=as.numeric(Sys.getenv('PORT', '10000')), launch.browser=FALSE)"]
