# Pacific Map Lab

Interactive coastal cartography for the Pacific margin of North and South America.

The public web edition opens on Coastal California and is intentionally lightweight enough for public hosting. Visitors can:

- switch among ocean/bathymetry, satellite, terrain, dark, and light basemaps;
- jump to Pacific Northwest, Baja California, tropical eastern Pacific, Peru, and Chile presets;
- upload CSV or TSV field sites with longitude and latitude;
- style, label, cluster, and frame uploaded sites;
- click the map to build a geodesic survey path; and
- download mapped sites, survey points, and project settings.

## Deployment

- Host: Render
- Runtime: Docker
- Service: `pacific-coast-map-lab`
- Root directory: `pacific-map-lab`
- Health check: `/`
- Public URL: `https://pacific-coast-map-lab.onrender.com`

The desktop publication-map application remains separate from this public Leaflet edition so the hosted service stays responsive within Render Free's memory limit.
