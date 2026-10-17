# Pacific Map Lab

Interactive R Shiny cartography for the Pacific coast of North and South America.

The public application opens with a preset Coastal California map, then lets visitors change the extent, bathymetric surface, palette, relief, contours, labels, and export settings.

## Deployment

- Host: Render
- Runtime: Docker
- Root directory: `pacific-map-lab`
- Health check: `/`
- Cache: `/tmp/PacificMapLab-cache`

The numbered `app.tar.gz.b64.part*` files are reassembled into `app.R` during the Docker image build.
