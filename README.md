# celmap

Tools for working with **Computational Ecology Lab** map data


### Requirements

+ [R 3+](https://www.r-project.org/)

+ [sp](https://cran.r-project.org/package=sp)

+ [raster](https://cran.r-project.org/package=raster)

+ [nearth](https://github.com/BigelowLab/nearth)

### Installation

```r
if (!('nearth' %in% rownames(installed.packages())))
    devtools::install("BigelowLab/nearth")

devtools::install("BigelowLab/celmap")
```

### Coastal vectors

```r
# get a 50m interval coastline
Maine_coast = celmap::get_vectors(where = 'maine', what = 'coast50m')

# get a 50m interval state/province for Northeast North Atlantic region
NEAC = celmap::get_vectors(where = 'neac', what = 'boundary50m')
```

### Bounding Boxes

`celmap` provides a uniform repository for a number of bounding boxes by region.

    + `maine`  State of Maine  `[-71.1, -67, 43, 47.5]`

    + `gom`  Gulf of Maine  `[-72, -63, 39, 46]`

    + `nwa`  Northwest Atlantic  `[-77.0, -51.5, 37.9, 56.7]`

    + `neac`  New England Atlantic Canada  `[-74, -59.75, 41, 48.15]`

    + `world` the globe `[-180, 180, -90, 90]`

### Rasters

Use `read_grd()` to read raster `.grd` files.  This can be useful for comparing
a suite of rasters to each other - to look for anomalies in extent, dimensions,
projections, etc. - without needing to open the larger binary `.gri` files.

Use `closest_ocean_cell(x, lut)` to reassign the locations in x (with [lon, lat])
along shorelines to the closest cells.  The lut is a precomputed raster where each
cells values is the cell number of the closest ocean cell.  Obviously, for a
location already on the ocean it returns its own cell number.

Use `cellFromPts(R, pts)` to simply convert from lon/lat to cell numbers.

Use `layers_extractPoints(R, pts)` to extract points from a multi-layer raster
object.

Use `layers_randomPoints(R, pts, N)` to select N non-NA random points that avoid
pts locations in a multi-layer raster.

