FROM rocker/shiny:4.5.1

LABEL org.opencontainers.image.title="Pacific Map Lab"
LABEL org.opencontainers.image.description="Interactive West Coast Americas coastal cartography tool"

RUN install2.r --error --skipinstalled --ncpus -1 \
    bslib leaflet htmltools jsonlite

WORKDIR /app
COPY app.R /app/app.R

ENV PORT=10000
EXPOSE 10000

CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=as.numeric(Sys.getenv('PORT', '10000')), launch.browser=FALSE)"]
