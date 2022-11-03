``` r
options(knitr.kable.NA = '')
```

# Climate futures for [Haskell Indian Nations University](https://www.haskell.edu)

## Prepare your environment

1.  Install Climate Futures Toolbox

``` r
install.packages("cft")
```

2.  Load other packages

``` r
library(tidyverse)
library(tidync)
library(cft)
library(sf)
library(ggplot2)
library(ggthemes)
library(ggpattern)
library(magick)
library(future)
library(tidytable)
library(janitor)
options(timeout = 600L)
```

## Finding basic layers in OpenStreetMap

Explain the basics of APIs and the premise of fast downloads.

### Use plain lnaguage to request a bounding box

1.  Find the general area on Open Street Map We use a function from the
    osmdata package to find a bounding box for our area of interest.
    This is a nice function for this purpose because it can use plan
    language declarations (e.g. “Lawrence, Kansas” or “Boulder,
    Colorado”) for the location. You do not need to use this function to
    define a bounding box. You can define your bounding box from any
    source. The benefit of this method is that is it is rather easy and
    reliable.

``` r
bb <- getbb("Lawrence, Kansas")
bb
```

    ##         min       max
    ## x -95.34454 -95.16662
    ## y  38.90447  39.03350

``` r
bbb <- bb
bbb[1,1] <- bbb[1,1] - 0.001
bbb[1,2] <- bbb[1,2] + 0.001
bbb[2,1] <- bbb[2,1] - 0.03
bbb[2,2] <- bbb[2,2] + 0.001
xlimit <- bbb[1,]
ylimit <- bbb[2,] 
xmid <- xlimit[1] + diff(xlimit) / 2 
ratio <- diff(xlimit) / diff(ylimit)
```

2.  Find any buildings associated with any University in our “Lawrence,
    Kansas” bounding box. This request first calls the opq() function,
    which mounts to the OpenStreetMap database and then queries the
    “building” key (i.e. all the building footprints) for any building
    types matching the value “university”. This is a representation of
    the “key” and “value” system that OSM uses to query data. The final
    step is to convert the OSM output into a spatial format that works
    well in R, called sf.

``` r
my_boundary <- opq(bb) %>% 
  add_osm_feature(key = "building", value = "university") %>% 
osmdata_sf() 

summary(my_boundary)
```

    ##                   Length Class  Mode     
    ## bbox               1     -none- character
    ## overpass_call      1     -none- character
    ## meta               3     -none- list     
    ## osm_points         5     sf     list     
    ## osm_lines         40     sf     list     
    ## osm_polygons      40     sf     list     
    ## osm_multilines     0     -none- NULL     
    ## osm_multipolygons 17     sf     list

3.  The output from this request shows a list of multipolygons,
    polygons, linestrings, and points. Each of these data types have a
    different storage structure, so we can’t look at them all at the
    same time. Instead, lets start with polygons, which likely represent
    a single building footprint. Printing ‘my_boundary$osm_polygons’
    shows that there are two Universities in Lawrence and we need to
    filter those results down to only include Haskell building.

``` r
boundaries <- my_boundary$osm_polygons %>% 
  filter(operator == "Haskell Indian Nations University") 
boundaries1 <- boundaries[1,] #take the first building (e.g. first row) of the returns
```

4.  Plot our discovered footprint to visually confirm It looks like we
    found [Winona Hall](https://www.kansasmemory.org/item/449914) in the
    OpenStreetMap database. This is how we plot the perimeters
    associated with it. ![Winona Hall](winona%20hall%20haskell.png)

``` r
basemap <- ggplot(data = st_as_sf(boundaries1)) +
  geom_sf(fill = "blueviolet") +
  geom_sf_text(aes(label = name), size=10) +
  theme_tufte()

basemap
```

![This seems to match.](Haskell_files/figure-gfm/unnamed-chunk-61-1.png)

## Build a basemap from OpenStreetMap

### Set your color palette

1.  Download the Haskell university logo

``` r
seamless_image_filenames <- c(
  'Haskell_logo copy.jpg'
)
```

<img src="/Users/ty/Documents/Haskell%20api%20demo/Haskell_logo.png"
style="width:50.0%" />

2.  Sample the colors on that logo to make a custom color palette for
    our basemap

``` r
our_blue <- "#3A4E8B"
our_yellow <- "#FFD60F"
our_beige <- "#EDEDF0"
our_purple <- "#3E1471"
our_purple_alpha <- "#3E1471999"
```

### Download all the layers you want to include in your basemap

1.  Download the Haskell University footprint

``` r
my_boundary <- opq(bb) %>% 
  add_osm_feature(key = "amenity", value = "university") %>% 
osmdata_sf() 

haskell_poly <- my_boundary$osm_multipolygons[1,]
```

![](Haskell_files/figure-gfm/unnamed-chunk-65-1.png)<!-- -->

2.  Download street vector layers The street vector is divided into two
    different downloads in order to create two different objects for
    coloring in the final figure. This first download will be in the
    foreground. It includes the larger and faster roadways.

``` r
# the big streets
streets <- 
  opq(bb) %>% 
  add_osm_feature(key = "highway", 
                  value = c("motorway", "trunk",  "primary", 
                            "secondary", "tertiary")) %>% 
  osmdata_sf() 

streets_crop <- streets$osm_lines %>%
  st_crop(y = c(ymin = bb[2,1], ymax = bb[2,2], xmin = bb[1,1], xmax = bb[1,2]))
```

![](Haskell_files/figure-gfm/unnamed-chunk-67-1.png)<!-- -->

The second street download is for the small side streets and footpaths.
These lines will be more faint and in the background.

``` r
small_streets <-  
  opq(bb) %>% 
  add_osm_feature(
    key = "highway", 
    value = c("residential", "living", "unclassified", 
              "service", "footway")) %>% 
  osmdata_sf()

small_streets_crop <- small_streets$osm_lines %>%
  st_crop(y = c(ymin = bb[2,1], ymax = bb[2,2], xmin = bb[1,1], xmax = bb[1,2]))
```

![](Haskell_files/figure-gfm/unnamed-chunk-69-1.png)<!-- -->

3.  Download water features. The water features are first divided into
    moving and stationary water. We will download the river layer from
    the waterway key.

``` r
water <- 
  opq(bb) %>% 
  add_osm_feature(key = "waterway", value = "river") %>% 
  osmdata_sf()
```

We divide the water into large and small waterways in the same way we
did with the road. We are interested in making the main river much
larger and the remaining waterways collectively smaller. The Kansas
river is the large feature in this map so, we pull it out first.

``` r
Kansas_river_multi <- water$osm_multilines %>%
  filter(name == "Kansas River")  %>% 
  st_as_sf() %>%
  st_crop(y = c(ymin = bb[2,1], ymax = bb[2,2], xmin = bb[1,1], xmax = bb[1,2]))
```

![](Haskell_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->

After removing the Kansas river, we are left with a number of remaining
waterways that are stored as both linestrings and multilinestrings. We
need to download each of those data types individually.

``` r
small_water_multi <- water$osm_multilines %>%
  filter(name != "Kansas River")%>%
  st_as_sf() %>%
  st_crop(y = c(ymin = bb[2,1], ymax = bb[2,2], xmin = bb[1,1], xmax = bb[1,2]))

small_water_lines <- water$osm_lines %>%
  filter(name != "Kansas River")%>%
  st_as_sf() %>%
  st_crop(y = c(ymin = bb[2,1], ymax = bb[2,2], xmin = bb[1,1], xmax = bb[1,2]))
```

![](Haskell_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->

The stationary water bodies are divided into two categories: lakes and
reservoirs. We need to pull each of these individually. First, the
reservoir.

``` r
reservoir <- 
  opq(bb) %>% 
  add_osm_feature(key = "water", value = "reservoir") %>% 
  osmdata_sf()

reservoir_crop_multi <- reservoir$osm_multipolygons %>%
  st_as_sf() %>%
  st_crop(y = c(ymin = bb[2,1], ymax = bb[2,2], xmin = bb[1,1], xmax = bb[1,2]))

reservoir_crop_poly <- reservoir$osm_multipolygons %>%
  st_as_sf() %>%
  st_crop(y = c(ymin = bb[2,1], ymax = bb[2,2], xmin = bb[1,1], xmax = bb[1,2]))

reservoir_crop <- cbind(reservoir_crop_multi, reservoir_crop_poly)
```

![](Haskell_files/figure-gfm/unnamed-chunk-76-1.png)<!-- -->

Now the lakes.

``` r
lake <- 
  opq(bb) %>% 
  add_osm_feature(key = "water", value = "lake") %>% 
  osmdata_sf()

lake_crop <- lake$osm_polygons %>%
  st_as_sf() %>%
  st_crop(y = c(ymin = bb[2,1], ymax = bb[2,2], xmin = bb[1,1], xmax = bb[1,2]))
```

![](Haskell_files/figure-gfm/unnamed-chunk-78-1.png)<!-- -->

### Stack downloaded OSM layers into a final basemap.

``` r
haskel_basemap <- ggplot() +
  # plot moving water layers first
  geom_sf(data = Kansas_river_multi, alpha = .8,
          size = 3, colour = our_blue) +
  geom_sf(data = small_water_multi, alpha = .8,
          size = 0.5, colour = our_blue) +
  geom_sf(data = small_water_lines, alpha = .8,
          size = 0.5, colour = our_blue) +
  # Layer bodies of water over the moving water layers
  geom_sf(data = lake_crop, alpha = 1, fill = our_blue, size=0) +
  geom_sf(data = reservoir_crop, alpha = 1, fill = our_blue, size=0) +
  
  # Plot small streets in the background with a light fade
  geom_sf(data = small_streets_crop, alpha = .6, 
          size = .1, colour = our_beige) +
  # Layer large streets over the top of the small streets with a bold color.
  geom_sf(data = streets_crop, alpha = .8, 
          size = .4, colour = our_yellow ) +
  # Layer Haskell university property polygon in the foreground
  geom_sf( data=haskell_poly, color=our_yellow, size=0.5) +
  # Fill Haskell property polygon with Haskell logo
   geom_sf_pattern( 
     data = haskell_poly,
     size=0,
    pattern       = 'image',
    pattern_type  = 'tile',
    pattern_scale = 1.1,
    pattern_filename = seamless_image_filenames
  ) +
  # set limits on final figure 
  coord_sf(ylim = ylimit, xlim = xlimit, expand = TRUE) +
  # adding labels
  annotate(geom = "text", y = bbb[2,1]+ 0.013, x = xmid, 
           label = "Haskell Indian Nations University", size = 12, colour = our_beige
           ) +
  annotate(geom = "errorbarh", xmin = xlimit[1], xmax = xlimit[2], y = bbb[2,1]+ 0.005,   
           height = 0, size = 0.5, colour = our_beige) +
  annotate(geom = "text", y = bbb[2,1]+ 0.0001, x =  xmid,
           label = "38°93'88\"N  95°23'29\"W",  size = 6,
           colour = our_beige) +
  # clean out unused elements and set background color
  theme_void() +
  theme(panel.background = element_rect(fill = our_purple),
        plot.background = element_rect(fill = NA))
```

![](All_roads_lead_to_Haskell.png)

### Save the basemap in high resolution print

``` r
ggsave(haskel_basemap, filename = "All_roads_lead_to_Haskell.png", height = 11, width=8.5, 
       units="in", dpi=600)
```

## Mount the climate dataset

This dataset is way too big to download to a particular machine.
Instead, you mount to the analysis ready data cube (i.e. netCDF) and
only download the subsetted data that you want to pull.

![](climate_grid.png)

1.  We calculate the center point for measuring to the nearest climate
    data point.

``` r
haskel_centroid <- st_coordinates(st_centroid(haskell_poly))
lat_pt <- haskel_centroid[1,2]
lon_pt <- haskel_centroid[1,1]
```

2.  Connect to the web server and activate the proper data dimensions.

``` r
web_link = "https://cida.usgs.gov/thredds/dodsC/macav2metdata_daily_future"

# Change to "https://cida.usgs.gov/thredds/catalog.html?dataset=cida.usgs.gov/macav2metdata_daily_historical" for historical data. 

src <- tidync::tidync(web_link)
lons <- src %>% activate("D2") %>% hyper_tibble()
lats <- src %>% activate("D1") %>% hyper_tibble()
```

3.  Search through the database of climate prediction points to find
    which one is closest to our centroid. We then spatially project that
    chosen pt into an sf object.

``` r
known_lon <- lons[which(abs(lons-lon_pt)==min(abs(lons-lon_pt))),]
known_lat <- lats[which(abs(lats-lat_pt)==min(abs(lats-lat_pt))),]

chosen_pt <- st_as_sf(cbind(known_lon,known_lat), coords = c("lon", "lat"), crs = "WGS84", agr = "constant")
```

![](grid_match.png)

## Downloading climate data for our bounding box.

### Connect to the downscaled dataset.

We cannot download the entire dataset. Instead, we mount to that dataset
by connecting to an api that route us to a particular part of the
dataset based on the bounding box specified. Mount to the USGS
downscaled dataset. The resultant object called ‘input’ includes three
elements The first is the full list of available data at each timestep.
The second is a list of possible time steps. The third is a list of
verbatim copy of the raw return from the server. This raw return shows
the dimensions of the data and how many elements are available in each
of those dimensions.

``` r
# Mount to the USGS downscaled dataset. 
inputs <- cft::available_data()
```

    ## Trying to connect to the USGS.gov API

    ## not a file: 
    ## ' https://cida.usgs.gov/thredds/dodsC/macav2metdata_daily_future '
    ## 
    ## ... attempting remote connection

    ## Connection succeeded.

    ## Reading results

    ## Converting into an R data.table

Element 1. Dataframe of availble data and descriptions of each of those
variables.

``` r
head(inputs[[1]])
```

Element 2. Short list of available data you can request

``` r
head(unique(inputs[[1]]$Variable))
```

    ## [1] "Specific Humidity"                   "Precipitation"                      
    ## [3] "Maximum Relative Humidity"           "Minimum Relative Humidity"          
    ## [5] "Surface Downswelling Shortwave Flux" "Maximum Temperature"

Element 3. Dataframe of available times

``` r
head(inputs[[2]])
```

## Compile your data request to send to the server.

Your data order needs to include three things: Variable, Scenario, and
Model. Your options can be found in the input elements we explored in
the previous section.

``` r
input_variables <- inputs$variable_names %>% 
  filter(Variable %in% c("Maximum Relative Humidity", 
                       "Minimum Relative Humidity",
                       "Maximum Temperature",
                       "Minimum Temperature",                 
                       "Precipitation",
                       "Eastward Wind",
                       "Northward Wind")) %>% 
  filter(Scenario %in% c( "RCP 8.5")) %>% 
  filter(Model %in% c(
    "Beijing Climate Center - Climate System Model 1.1",
    "Beijing Normal University - Earth System Model",
    "Canadian Earth System Model 2",                                                                
  "Centre National de Recherches Météorologiques - Climate Model 5",                              
  "Commonwealth Scientific and Industrial Research Organisation - Mk3.6.0",                       
  "Community Climate System Model 4",                                                             
  "Geophysical Fluid Dynamics Laboratory - Earth System Model 2 Generalized Ocean Layer Dynamics",
  "Geophysical Fluid Dynamics Laboratory - Earth System Model 2 Modular Ocean",                   
  "Hadley Global Environment Model 2 - Climate Chemistry 365 (day) ",                             
 "Hadley Global Environment Model 2 - Earth System 365 (day)",                                   
 "Institut Pierre Simon Laplace (IPSL) - Climate Model 5A - Low Resolution",                     
 "Institut Pierre Simon Laplace (IPSL) - Climate Model 5A - Medium Resolution",                  
 "Institut Pierre Simon Laplace (IPSL) - Climate Model 5B - Low Resolution",                     
 "Institute of Numerical Mathematics Climate Model 4",                                           
 "Meteorological Research Institute - Coupled Global Climate Model 3",                           
 "Model for Interdisciplinary Research On Climate - Earth System Model",                         
 "Model for Interdisciplinary Research On Climate - Earth System Model - Chemistry",             
 "Model for Interdisciplinary Research On Climate 5",                                            
 "Norwegian Earth System Model 1 - Medium Resolution"  )) %>%
  
  pull("Available variable")

head(as.data.frame(input_variables))
```

## Prepare for parallelization

``` r
n_cores <- availableCores() - 1
plan(multisession, workers = n_cores)
```

## Make parallel call to USGS server.

``` r
out <- single_point_firehose(input_variables, known_lat, known_lon )
head(out)
```

## Save output

``` r
haskell <- out
save(haskell, file = "haskell.RData")
```

## Load saved output

``` r
load(file = "haskell.RData")
```

``` r
df_precip <- haskell %>%
  st_drop_geometry() %>%
  select(time, which(startsWith(colnames(haskell), "pr"))) %>%
  gather(key = "variable", value = "value", -time)
head(df_precip)
```

``` r
df_min_temp <- haskell %>%
  st_drop_geometry() %>%
  select(time, which(startsWith(colnames(haskell), "tasmin"))) %>%
  gather(key = "variable", value = "value", -time)

df_min_temp <- df_min_temp[which(df_min_temp$value > 200), ]
df_min_temp$variable <- as.character(as.numeric(as.factor(df_min_temp$variable)))

head(df_min_temp)
```

``` r
df_RH <- haskell %>%
  st_drop_geometry() %>%
  select(time,which(startsWith(colnames(haskell), "rhs")) ) %>%
  gather(key = "variable", value = "value", -time)

df_RH_min <- haskell %>%
  st_drop_geometry() %>%
  select(time,which(startsWith(colnames(haskell), "rhsmin"))) %>%
  gather(key = "variable", value = "value", -time)

df_RH_max <- haskell %>%
  st_drop_geometry() %>%
  select(time,which(startsWith(colnames(haskell), "rhsmax"))) %>%
  gather(key = "variable", value = "value", -time)
head(df_RH_max)
```

``` r
all_climate_projections_RH <- ggplot(data= df_RH, aes(x = time, y = value, color = variable)) + 
  #geom_line()+
  
  scale_colour_viridis_d(option="A", alpha = 0.8) +
  geom_smooth() +
  theme_tufte() +
  theme(legend.position = "none") 
  

all_climate_projections_RH
```

![](Haskell_files/figure-gfm/unnamed-chunk-101-1.png)<!-- -->

``` r
all_climate_projections_RH <- ggplot(data= df_RH_min, aes(x = time, y = value) )+ 
  geom_smooth(data= df_RH_min, color=our_purple) +
  geom_smooth(data= df_RH_max, color=our_purple) +
  theme_tufte() +
  ylim(0,100)
  

all_climate_projections_RH
```

![](Haskell_files/figure-gfm/unnamed-chunk-102-1.png)<!-- -->

``` r
all_climate_projections <- ggplot(data= df_min_temp, aes(x = time, y = value, color = variable)) + 
  geom_line()+
  geom_smooth() +
  scale_colour_viridis_d(option="A")
  

all_climate_projections
```

![](Haskell_files/figure-gfm/unnamed-chunk-103-1.png)<!-- -->

``` r
ensemble_climate_projections <- ggplot(data= df_min_temp, aes(x = time, y = value)) + 
  geom_line(color=our_purple)+
  geom_smooth(color=our_yellow) +
  theme_tufte()
  

ensemble_climate_projections
```

![](Haskell_files/figure-gfm/unnamed-chunk-104-1.png)<!-- -->

``` r
ggplot(data= df_min_temp, aes(x = time, y = value)) + 
  geom_smooth(color=our_purple) +
  theme_tufte()
```

![](Haskell_files/figure-gfm/unnamed-chunk-104-2.png)<!-- -->

``` r
all_climate_projections <- ggplot(data= df_precip, aes(x = time, y = value, color = variable)) + 
  geom_line(color=our_purple)+
  theme_tufte()
  

all_climate_projections
```

![](Haskell_files/figure-gfm/unnamed-chunk-105-1.png)<!-- -->

``` r
df_precip_2 <- df_precip[which(df_precip$value != 0),]

ensemble_climate_projections <- ggplot(data= df_precip_2, aes(x = time, y = value)) + 
   geom_line(color=our_purple)+
  geom_smooth(color=our_yellow)+
  theme_tufte()
  

ensemble_climate_projections
```

![](Haskell_files/figure-gfm/unnamed-chunk-106-1.png)<!-- -->

``` r
ggplot(data= df_precip, aes(x = time, y = value)) + 
     geom_smooth(color=our_purple)+
  theme_tufte()
```

![](Haskell_files/figure-gfm/unnamed-chunk-106-2.png)<!-- -->

``` r
df_precip_quant_95 <- df_precip %>% 
  mutate(quant = ntile(value, 100)) %>%
  filter(quant == 95)

df_precip_quant_96 <- df_precip %>% 
  mutate(quant = ntile(value, 100)) %>%
  filter(quant == 96)

df_precip_quant_97 <- df_precip %>% 
  mutate(quant = ntile(value, 100)) %>%
  filter(quant == 97)

df_precip_quant_98 <- df_precip %>% 
  mutate(quant = ntile(value, 100)) %>%
  filter(quant == 98)

df_precip_quant_99 <- df_precip %>% 
  mutate(quant = ntile(value, 100)) %>%
  filter(quant == 99)

df_precip_quant_100 <- df_precip %>% 
  mutate(quant = ntile(value, 100)) %>%
  filter(quant == 100)

df_precip_quant_1000 <- df_precip %>% 
  mutate(quant = ntile(value, 1000)) %>%
  filter(quant == 1000)

df_precip_quant_100000 <- df_precip %>% 
  mutate(quant = ntile(value, 100000)) %>%
  filter(quant == 100000)
```

``` r
ggplot(data= df_precip, aes(x = time, y = value)) + 
     geom_smooth(data= df_precip, color=our_purple)+ 
     geom_smooth(data= df_precip_quant_95, color=our_purple) +
    geom_smooth(data= df_precip_quant_96, color=our_purple) +
    geom_smooth(data= df_precip_quant_97, color=our_purple) +
    geom_smooth(data= df_precip_quant_98, color=our_purple) +
    geom_smooth(data= df_precip_quant_99, color=our_purple) +
  geom_smooth(data= df_precip_quant_100, color=our_purple) +
  geom_smooth(data= df_precip_quant_1000, color=our_purple) +
    geom_smooth(data= df_precip_quant_100000, color=our_purple) +
  theme_tufte()
```

![](Haskell_files/figure-gfm/unnamed-chunk-108-1.png)<!-- -->
