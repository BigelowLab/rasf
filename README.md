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

### Usage

```r
# get a bounding box fro Gulf of Maine as SpatialPolygons
bb = celmap::celmap_bb("gom", form = 'sp')

# get a 50m interval coastline
Maine_coast = celmap::get_vectors(where = 'maine', what = 'coast50m')

# get a 50m interval state/province for Northeast North Atlantic region
NEAC = celmap::get_vectors(where = 'neac', what = 'boundary50m')
```