# Pacific Map Lab

Pacific Map Lab is a browser-based publication-map studio for the Pacific margin of North and South America. It builds real `ggplot2` maps from NOAA ETOPO bathymetry, Natural Earth `sf` layers, and user-supplied research data.

The public app opens on Coastal California and supports presets for the Pacific Northwest, Southern and Central California, Baja California, the tropical eastern Pacific, Peru, and Chile.

## R-style mapping workflow

Users can:

- query a protected-resolution NOAA ETOPO 2022 grid for a chosen coastal extent;
- style bathymetry, land, coastline, isobaths, coordinate grids, legends, scale bars, north arrows, titles, and captions;
- upload CSV/TSV coordinates, GeoJSON, GeoPackages, and zipped Shapefiles;
- display and label point, line, and polygon research layers;
- export PNG, TIFF, PDF, and SVG figures; and
- download the full editable R project behind the map.

The R project bundle contains:

- the rendered figure;
- an editable `ggplot` object (`.rds`);
- NOAA bathymetry as CSV;
- land, coastline, and uploaded research layers as a GeoPackage;
- `recreate_map.R`;
- exact settings and data provenance; and
- `sessionInfo()` output.

## Deployment

- Host: Render
- Runtime: Docker / R 4.5.1
- Service: `pacific-coast-map-lab`
- Branch: `render-lowmem`
- Root directory: `pacific-map-lab`
- Health check: `/`
- Public URL: `https://pacific-coast-map-lab.onrender.com`

The hosted build uses a protected grid-cell ceiling and a low-memory package-loading strategy so regional bathymetry, `sf` layers, and high-resolution figure export can coexist on the public service.
