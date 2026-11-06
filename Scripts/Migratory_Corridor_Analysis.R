###### Setup ######


install.packages(c("tidyverse", "rnaturalearth", 
                   "rnaturalearthdata", "MASS"), type = "binary") #Basics
install.packages("M3")

install.packages("aniMotum", 
                 repos = c("https://cloud.r-project.org",
                           "https://ianjonsen.r-universe.dev"),
                 dependencies = TRUE) #For filtering
install.packages("spatstat")
install.packages("adehabitatHR") #For KUD
install.packages("sp")   # for Spatial Points and lines

install.packages("marmap") #For bathymetry
install.packages("ggnewscale") #For transparent grapsh
install.packages("move") #Brownie bridge
install.packages("rasterVis") #Fancy rasters for the DBB
install.packages("gridExtra")
install.packages("SDLfilter") #overlap metrics
install.packages("concaveman") #concaveman small area polygon
install.packages("openxlsx")

library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(MASS)
library(M3)
library(aniMotum)
library(spatstat)
library(adehabitatHR)
library(sp)
library(raster)
library(marmap)
library(terra)
library(ggnewscale)
library(move)
library(rasterVis)
library(gridExtra)
library(SDLfilter)
library(concaveman)
library(zoo)


#Set your own working directory
setwd("~/GitHub/honours_corridors/Data")

#Set up the map for your graphs
world <- ne_countries(scale = "medium", returnclass = "sf")



##### Choose your Dataset #####

###### Great Shearwater PTT ######

#Load data
GSPTT <- read.csv("GreatShearwaterPTT_Raw.csv")

#Clean it to fit with animotum package
GSPTT$POSIX <- paste(GSPTT$date_gmt, GSPTT$time_gmt)
GSPTT$POSIX <- as.POSIXct(GSPTT$POSIX, format = "%d/%m/%Y %H:%M", tz = "UTC")
GSPTT_model <- GSPTT %>% 
  #Selects necessary data
  dplyr::select(scientific_name, common_name, site_name, colony_name, lat_colony, lon_colony,
                track_id, sex, argos_quality, POSIX, longitude, latitude) %>% 
  #Group by track_id
  group_by(track_id) %>% 
  #Filters out individuals with low number of points
  filter(n()>=10) %>% 
  #Filters out movements from breeding times
  filter(POSIX > as.POSIXct("2010-01-01 00:00:00", tz = "UTC")) 



df_sf_raw_dat <- GSPTT %>% 
  st_as_sf(coords = c('longitude', 'latitude')) %>% 
  st_set_crs(4326)

#Pre filtered Data
ggplot() +
  theme_bw() +
  geom_sf(aes(), data = world) +
  geom_sf(aes(colour = track_id), data = df_sf_raw_dat, show.legend = "point") +
  coord_sf(xlim = c(-80, 10), ylim = c(-60, 60)) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme(legend.position = "none")




#Format the data for the aniMotum package
AniMotum <- format_data(GSPTT_model, id = "track_id", date = "POSIX", 
                        lc = "argos_quality", coord = c("longitude", "latitude"))





##### Generate the SSM and MP #####

#Creates an SSM - random walk with time step of 12hrs
Fit <- fit_ssm(AniMotum, model = "rw", time.step = 12,
               control = ssm_control(verbose = 0))


#Reroutes the SSM to be over the ocean
Reroute <- route_path(Fit, what = "predicted", 
                                map_scale = 10, dist = 5000)


#Fits a Move persist model to the rerouted function
MPM <- fit_mpm(Reroute, what = "rerouted", model = "mpm", 
                     control = mpm_control(verbose = 0))



## Setup for map



laea_centre <- function(i) {
  range <- range(i)
  dif <- diff(range, na.rm = TRUE)
  mid <- min(i) + (dif/2)
  return(mid)
}





##### Extract the Rerouted locations and their corresponding Move Persist Values #####

df_locs <- Reroute |> grab("rerouted")
df_mpm <- MPM |> grab("fitted")

SA_shelf <- getNOAA.bathy(lon1=10,lon2=-70,lat1=10,lat2=-60, resolution=10)



### Map the total Points fitted

GSPTT_all <- merge(df_locs, df_mpm)

df_sf_all_filtered <- GSPTT_all %>% 
  st_as_sf(coords = c('lon', 'lat')) %>% 
  st_set_crs(4326)


#Get the coordinates from the bathymetry object
lon_vals_bathy <- as.numeric(rownames(SA_shelf))
lat_vals_bathy <- as.numeric(colnames(SA_shelf))

# Function to find the nearest depth value
get_nearest_depth <- function(lon, lat, bathy, lon_vals_bathy, lat_vals_bathy) {
  # Find nearest longitude and latitude indices
  lon_idx <- which.min(abs(lon_vals_bathy - lon))
  lat_idx <- which.min(abs(lat_vals_bathy - lat))
  # Extract depth value
  depth <- SA_shelf[lon_idx, lat_idx]
  return(depth)
}


# mapply nearest depth to the locations
GSPTT_all$depth <- mapply(get_nearest_depth,
                             lon = GSPTT_all$lon,
                             lat = GSPTT_all$lat,
                             MoreArgs = list(bathy = SA_shelf, lon_vals_bathy = lon_vals_bathy, lat_vals_bathy = lat_vals_bathy))


ggplot() +
  theme_bw() +
  geom_sf(aes(), data = world, fill = "grey") +
  geom_sf(aes(colour = g), data = df_sf_all_filtered, show.legend = "point") +
  coord_sf(xlim = c(-80, 10), ylim = c(-60, 60)) +
  labs(x = "Longitude",
       y = "Latitude") 




##### GSPTT Methods #####
###### Premigration Filtering ######


#Make a buffer around the Colony

# Colony Location
inaccessible_island_coords <- st_sfc(
   st_point(c(-12.6733, -37.3004)),
   crs = 4326
 )
  
# Create laea for colony

colony_lat_centre <- laea_centre(-37.3004)
colony_lon_centre <- laea_centre(-12.6733)

colony_laea_crs <- sprintf("+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m", 
                           colony_lat_centre, 
                           colony_lon_centre)


# Transform into laea to make buffer

colony_proj <- st_transform(inaccessible_island_coords, colony_laea_crs)

# Buffer in metres

colony_buffer <- st_buffer(colony_proj, 25000)
  

GSPTT_premig <- GSPTT_all %>% 
  group_by(id) %>% 
  # First Filter out the individuals that do not have valid data for their premigration
  filter(id != "663_96803a_7510") %>%  #filters out the bird that goes the wrong direction
  filter(id != "663_96806a_7504") %>%  #Bird that started on the shelf
  filter(id != "663_96808a_7518") %>%  #filters bird whose tag does not move
  filter(id != "663_96810a_7523") %>%  #Filters out the bird that starts on the shelf
  filter(id != "663_96811a_7524") %>%  #Doesn't migrate (starts away from colony, flies back, tag fails)
  filter(id != "663_96812a_7526") %>%  #Doesn't migrate (starts away from colony, flies back, tag fails)
  filter(id != "663_96814a_7528")   #filters bird whose tag does not move



## Filter the end by movement value and by depth

filtered_premig <- GSPTT_premig %>%
  group_by(id) %>%  #Groups by id 
  group_modify(~{
    
    dat <- .x %>% arrange(date)  # Orders by date
    
    # find first row where the move value and depth reach the values
    stop_idx <- which(dat$g <= 0.4 &
                        dat$depth >= -500)[1]
    
    # if ther is no match keeps entire track
    if (is.na(stop_idx)) return(dat)
    
    # otherwise stops at that row
    dat[1:stop_idx, ]
  }) %>%
  ungroup()



# Filter out the buffer area
# Make tracks

GSPTT_premig_tracks <- st_as_sf(filtered_premig, coords = c("lon", "lat"), crs = 4326) %>% 
  st_transform(colony_laea_crs)

# filter outside buffer
outside_buffer_premig <- GSPTT_premig_tracks[!st_intersects(GSPTT_premig_tracks, colony_buffer, sparse = FALSE), ]

GSPTT_premig <- outside_buffer_premig



# Orders for the DBBMM
GSPTT_premig <- GSPTT_premig[order(GSPTT_premig$id, GSPTT_premig$date), ] #%>%  #Orders it by id

GSPTT_premig <- st_transform(GSPTT_premig, crs = 4326) %>% 
  st_coordinates(GSPTT_premig) %>%
  as.data.frame() %>%
  rename(lon = X, lat = Y) %>%
  bind_cols(st_drop_geometry(GSPTT_premig))

GSPTT_premig <- GSPTT_premig %>% 
  mutate(SE = ifelse(x.se > y.se, y.se, x.se), #Adds a column of Standard Error, which is the smallest value of the error columns
         citation = NA,                        #Adds a column of NA citation values for DBBMM
         species = "bird") #%>%                 #Adds a column of "species" for DBBMM

#Set as sf
GSPTT_premig_sf_points <- st_as_sf(GSPTT_premig, coords = c("lon", "lat"), crs = 4326)



ggplot() +
  theme_bw() +
  geom_sf(aes(), data = world, fill = "grey") +
  geom_sf(aes(colour = id), data = GSPTT_premig_sf_points, show.legend = "point") +
  coord_sf(xlim = c(-80, 0), ylim = c(-60, 0)) +
  labs(x = "Longitude",
       y = "Latitude")








##### Premigration Methods #####

# Set the CRS for the laea for later


premig_lat_centre <- laea_centre(GSPTT_premig$lat)

premig_lon_centre <- laea_centre(GSPTT_premig$lon)

premig_laea_crs <- sprintf("+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m", 
                        premig_lat_centre, 
                        premig_lon_centre)


# World Map for Premig

world_premig <- st_transform(world, crs = premig_laea_crs)

# Set the extent for the maps

GSPTT_premig_ext <- SpatialPoints(coords = cbind(GSPTT_premig$lon, GSPTT_premig$lat), proj4string = CRS("WGS84")) 
GSPTT_premig_ext <- spTransform(GSPTT_premig_ext, premig_laea_crs)


GSPTT_premig_xmin <- min(GSPTT_premig_ext@coords[,1])* 1.5
GSPTT_premig_xmax <- max(GSPTT_premig_ext@coords[,1])* 1.5
GSPTT_premig_ymin <- min(GSPTT_premig_ext@coords[,2])* 1.5
GSPTT_premig_ymax <- max(GSPTT_premig_ext@coords[,2])* 1.5


# Average error for the points
mean_SE_premig <- mean(GSPTT_premig$SE)
print(mean_SE_premig)

# Number of individuals used for Premig
n_distinct(GSPTT_premig$id)



###### Kernel Utiilisation Distribution ######

# Change the data to be spatial points
GSPTT_premigKUD <- SpatialPoints(coords = cbind(GSPTT_premig$lon, GSPTT_premig$lat), 
                                 proj4string = CRS("+proj=longlat +datum=WGS84")) 
proj4string(GSPTT_premigKUD) <- CRS("+proj=longlat +datum=WGS84")
GSPTT_premigKUD <- spTransform(GSPTT_premigKUD, premig_laea_crs)


## Calc the grid size needed

# Extract bounding box
premig_bbox <- sp::bbox(GSPTT_premigKUD)

# Calc max range 
premig_x_range <- premig_bbox[1, "max"] - premig_bbox[1, "min"]
premig_y_range <- premig_bbox[2, "max"] - premig_bbox[2, "min"]

# Calc the max of the two
premig_max_range <- max(premig_x_range, premig_y_range)

# Calc the grid value from this
# the number represents the grid cell size in m that you want
# Due to the kernelUD taking an integer, it will round this value to the nearest
# whole number but it will be close enough
# if the grid number ends up super low, use the ceiling function to round it up
premig_grid_kud <- premig_max_range / 20000


## Apply the KUD

# The default kernel (bivariate normal kernel) is used by default in the adehabitatHR package
# smoothing parameter h is the width of the kernels and can be calculated
# Sigma (σ) is the standard deviation of points
# σ can be calculated by the mean of the standard error of the x or y multiplied 
# by the square root of the number of locations

premig_kud_h <- mean_SE_premig*1000

premig_kud <- kernelUD(GSPTT_premigKUD, h = premig_kud_h, grid = premig_grid_kud)


# Extract the Contours you want
premig_contour_50 <- getverticeshr(premig_kud, percent = 50)
premig_contour_95 <- getverticeshr(premig_kud, percent = 95)

# Convert contours to sf objects
premig_contour_50_sf <- st_as_sf(premig_contour_50)
premig_contour_95_sf <- st_as_sf(premig_contour_95)

#Turn the kud into a dataframe that ggplot can read
premig_Kuddf <- as.data.frame(premig_kud, xy = TRUE)


###Plot
GSPTT_premig_kud <- ggplot() +
  geom_sf(aes(), data = world_premig, fill = "grey") +
  geom_sf(data = premig_contour_95_sf, aes(fill = "Full Corridor", colour = "Full Corridor"), alpha = 0.5) +
  geom_sf(data = premig_contour_50_sf, aes(fill = "Core Corridor", colour = "Core Corridor"), alpha = 0.5) + 
  geom_sf(aes(colour = "Locations"), data = GSPTT_premig_sf_points, show.legend = "point") +
  coord_sf(xlim = c(GSPTT_premig_xmin, GSPTT_premig_xmax), 
           ylim = c(GSPTT_premig_ymin, GSPTT_premig_ymax), 
           expand = FALSE) +
  theme_bw() +
  labs(
    title = "A)",
    x = "Longitude",
    y = "Latitude") +
  scale_fill_manual(name = "Probability Contours", values = c("Locations" = "black",
                                                              "Full Corridor" = "red",
                                                              "Core Corridor" = "blue"),
                    guide = guide_legend(override.aes = list(alpha = 0.5))) +
  scale_colour_manual(name = "Probability Contours", values = c("Locations" = "black",
                                                                "Full Corridor" = "red", 
                                                                "Core Corridor" = "blue")) +
  theme(legend.position = "none")


plot(GSPTT_premig_kud)










###### Density Raster ######

#Make sf points in crs 
GSPTT_premig_sf_laea <- st_transform(GSPTT_premig_sf_points, premig_laea_crs)

#Choose grid cell size (in m)
premig_grid_size <- mean_SE_premig*1000

# Convert points to SpatVector
premig_points_vect <- vect(GSPTT_premig_sf_laea)
premig_points_vect <- project(premig_points_vect, premig_laea_crs)

#Create an empty raster in laea
premig_density_raster <- rast(ext(premig_points_vect), resolution = premig_grid_size, crs = premig_laea_crs)

# Rasterize with a density count
premig_r_density <- rasterize(premig_points_vect, premig_density_raster, 
                       fun = "length", background = 0)

# Convert to a dataframe for plotting in ggplot
premig_df_raster <- as.data.frame(premig_r_density, xy = TRUE)
colnames(premig_df_raster) <- c("x", "y", "Density")
premig_df_raster <- filter(premig_df_raster, Density > 0) %>% 
  filter(!is.na(Density))

#Calculate the top 10% bin values
premig_top_10 <- ceiling(0.10 * (nrow(premig_df_raster))) #Ceiling makes sure it never goes above 10%

premig_df_raster$Density <- as.numeric(premig_df_raster$Density)

#Order in descending order and bin
premig_df_raster <- premig_df_raster %>% 
  arrange(desc(Density)) %>% 
  mutate(rank = row_number()) %>% 
  mutate(bin = case_when(
    rank <= premig_top_10 ~ "Core Use",
    TRUE ~ "Corridor")) %>% 
  mutate(bin_num = case_when(
    rank <= premig_top_10 ~ 2,
    TRUE ~ 1))



## Removing noise

# re raster it
premig_cat_raster <- rast(premig_df_raster[, c("x", "y", "bin_num")], type = "xyz")

# remove the isolated points
premig_raster_mask <- premig_cat_raster == 1   #Just pulls the corridor points
premig_raster_patches <- patches(premig_raster_mask, directions = 8) # Creates patches based on whats touching in 8 directions
premig_raster_sizes <- freq(premig_raster_patches) #Counts the number of cells in each patch and puts into a df
premig_raster_ids_keep <- premig_raster_sizes$value[premig_raster_sizes$count > 2]  #Filters premig any patches with number of cells less than value
premig_clean_corridor <- mask(premig_raster_patches, premig_raster_patches %in% premig_raster_ids_keep, maskvalues = FALSE) #Makes a raster that only keeps the clean corridor
premig_clean_corridor[!is.na(premig_clean_corridor)] <- 1

# Extract the core use values to merge
premig_core_use <- premig_cat_raster
premig_core_use[premig_core_use != 2] <- NA

# Merge the rasters together for plotting
premig_clean_raster <- mosaic(premig_clean_corridor, premig_core_use, fun = "max")

# Turn into a df for plotting
premig_clean_df <- as.data.frame(premig_clean_raster, xy = TRUE)
colnames(premig_clean_df) <- c("x", "y", "bin_num")
premig_clean_df <- premig_clean_df %>% 
  mutate(bin = factor(bin_num,
                      levels = c(1, 2),
                      labels = c("Corridor", "Core Use")))



# Plot in ggplot
GSPTT_premig_density_raster <- ggplot() +
  geom_sf(aes(), data = world_premig, fill = "grey") +
  geom_tile(data = premig_df_raster, aes(x = x, y = y, fill = bin), colour = "black", alpha = 0.5) +
  geom_sf(aes(colour = "Locations"), data = GSPTT_premig_sf_points, show.legend = "point") +
  coord_sf(xlim = c(GSPTT_premig_xmin, GSPTT_premig_xmax), 
           ylim = c(GSPTT_premig_ymin, GSPTT_premig_ymax),
           expand = FALSE) +
  scale_fill_manual(values = c("Locations" = "black",
                              "Core Use" = "blue",
                               "Corridor" = "red"),
                     na.value = "transparent") +  # Transparent for missing values
  scale_colour_manual(values = c("Locations" = "black")) +
  labs(title = "B)",
    x = "Longitude",
    y = "Latitude") +
  theme_bw() +
  theme(legend.position = "none")


plot(GSPTT_premig_density_raster)










###### Line Density Estimation ######
# Turn the data into a format that can be exported into ARC



# Create line strings
GSPTT_premig_sf_lines <- GSPTT_premig_sf_points %>%
  group_by(id) %>%
  summarise(do_union = FALSE) %>%  # Avoid merging separate tracks
  st_cast("LINESTRING")  # Convert grouped points into lines

# Reproject to laea to be in metres
proj_premig_lines <- st_transform(GSPTT_premig_sf_lines, crs = premig_laea_crs)

# Save as shapefile
st_write(obj = proj_premig_lines,
         dsn = "~/GitHub/honours_corridors/Outputs/GSPTT_premig_LDE_lines.shp",
         driver = "ESRI Shapefile",
         append = FALSE)

# Load this in ARCgis in the Line Density fucntion from the Spatial Analyst Toolbox
# For search radius use the average or median Standard Error, whichever is more appropriate
# Use a grid cell size appropriate based on this error size and the distance of the journey
print(mean_SE_premig)

# The Grid Cell Size used:
# 20 000m

## Load the newly created LDE

premig_LDE_raster <- rast("~/GitHub/honours_corridors/Outputs/LDE_PTT_premig1.tif")

# Create the Contours

# Polygon 95
premig_LDE_95 <- premig_LDE_raster
premig_LDE_95[premig_LDE_95 < 0.05] <- NA
premig_LDE_95[premig_LDE_95 >= 0.05] <- 1

premig_LDE_95_poly <- terra::as.polygons(premig_LDE_95) %>% 
  sf::st_as_sf(crs = premig_laea_crs) %>%
  sf::st_buffer(dist = res(premig_LDE_95)[1]/4) %>%
  sf::st_simplify(dTolerance = res(premig_LDE_95)[1]) 



# Polygon 50
premig_LDE_50 <- premig_LDE_raster
premig_LDE_50[premig_LDE_50 < 0.50] <- NA
premig_LDE_50[premig_LDE_50 >= 0.50] <- 1

premig_LDE_50_poly <- terra::as.polygons(premig_LDE_50) %>% 
  sf::st_as_sf(crs = premig_laea_crs) %>%
  sf::st_buffer(dist = res(premig_LDE_50)[1]/4) %>%
  sf::st_simplify(dTolerance = res(premig_LDE_50)[1])



#Plot as a ggplot
GSPTT_premig_LDE <- ggplot() +
  geom_sf(aes(), data = world_premig, fill = "grey") +
  geom_sf(data = premig_LDE_95_poly, aes(fill = "Full Corridor", colour = "Full Corridor"), alpha = 0.5) +
  geom_sf(data = premig_LDE_50_poly, aes(fill = "Core Corridor", colour = "Core Corridor"), alpha = 0.5) +
  geom_sf(aes(colour = "Locations"), data = GSPTT_premig_sf_points, show.legend = "point") +
  coord_sf(xlim = c(GSPTT_premig_xmin, GSPTT_premig_xmax), 
           ylim = c(GSPTT_premig_ymin, GSPTT_premig_ymax),
           expand = FALSE) +
  labs(title = "C)",
       x = "Longitude",
       y = "Latitude") +
  scale_fill_manual(name = "Contours", values = c("Locations" = "black",
                                                  "Full Corridor" = "red", 
                                                  "Core Corridor" = "blue"),
                    guide = guide_legend(override.aes = list(alpha = 0.5))) +
  scale_colour_manual(name = "Contours", values = c("Locations" = "black",
                                                    "Full Corridor" = "red",
                                                    "Core Corridor" = "blue")) +
  theme_bw() +
  theme(legend.position = "none")



plot(GSPTT_premig_LDE)








###### Dynamic Brownian Bridge Movement Model ######
# Turn it into a move data for the DBBMM
GSPTT_premig_BBM <- move(x = GSPTT_premig$lon, y = GSPTT_premig$lat,
                         time = as.POSIXct(GSPTT_premig$date, format="%Y-%m-%d %H:%M:%OS", tz="UTC"), 
                         proj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
                         data = GSPTT_premig, animal = GSPTT_premig$id)

# Reproject into a metric measurment CRS for the DBBMM
GSPTT_premig_proj <- spTransform(GSPTT_premig_BBM, premig_laea_crs)


## Create the DBBMM

# Increasing window size increased the reliability of the brownian motion values
# but at the cost of missing rapid changes (not a problem here because of migration)
# Increasing the margin size increases the likelihood of identifying relatively weak breakpoints
# but at cost of potentially not detecting breakpoints due to limited observations in w
# Choices should be made based on the time interval you are looking at
DBB_GSPTT_premig <- brownian.bridge.dyn(GSPTT_premig_proj, 
                                        ext = 1.5, # Needs to be large enough to encompass the data. Increase until does
                                        raster = 20000, # raster box size in m
                                        location.error = mean_SE_premig * 1000, # error of location data in m
                                        window.size = 13, # Should be >2m, window of the breakpoints for behaviour
                                        margin = 5) # Holds out observations that wont be included in breakpoints from the window, 2m should be at most 1/2 of w


## Obtain the Utilisation Distributions of Each individual
premig_UD <- getVolumeUD(DBB_GSPTT_premig)
# plot(premig_UD, main = "UD")


# Stacks them all up
UD_premig <- stack(premig_UD)
premig_combined_UD <- calc(UD_premig, sum)


# Standardise the raster
premig_r_min <- cellStats(premig_combined_UD, "min")
premig_r_max <- cellStats(premig_combined_UD, "max")

premig_r_standardised <- (premig_combined_UD - premig_r_min) / (premig_r_max - premig_r_min)

# Reproject the raster
premig_r_reprojected <- projectRaster(premig_r_standardised, crs = premig_laea_crs)

# Turn it into a data frame for ggplot
premig_new_raster <- as.data.frame(premig_r_reprojected, xy = TRUE)
colnames(premig_new_raster) <- c("x", "y", "value")


# Set data frame as simple features and so ggplot can read it
premig_new_raster_sf <- st_as_sf(premig_new_raster, coords = c("x", "y"))
premig_new_raster_sf <- premig_new_raster_sf %>% st_set_crs(premig_laea_crs)
premig_new_raster_sf <- subset(premig_new_raster_sf, value <= 0.9995)


## Creating Contour layers
# 95 Full corridor
premig_c95 <- rasterToContour(premig_r_reprojected, levels = 0.95)
premig_c95_sf <- st_as_sf(premig_c95)
premig_c95_sf <- st_transform(premig_c95_sf, premig_laea_crs)
premig_c95_poly <- st_cast(st_union(premig_c95_sf), "POLYGON")


# 50 Core Corridor
premig_c50 <- rasterToContour(premig_r_reprojected, levels = 0.50)
premig_c50_sf <- st_as_sf(premig_c50)
premig_c50_sf <- st_transform(premig_c50_sf, premig_laea_crs)
premig_c50_poly <- st_cast(st_union(premig_c50_sf), "POLYGON")


## Plot

GSPTT_premig_DBBMM <- ggplot() +
  geom_sf(aes(), data = world_premig, fill = "grey") +
  geom_sf(data = premig_c95_poly, aes(fill = "Full Corridor", colour = "Full Corridor"), alpha = 0.5) +
  geom_sf(data = premig_c50_poly, aes(fill = "Core Corridor", colour = "Core Corridor"), alpha = 0.5) +
  geom_sf(aes(colour = "Locations"), data = GSPTT_premig_sf_points, show.legend = "point") +
  coord_sf(xlim = c(GSPTT_premig_xmin, GSPTT_premig_xmax), 
           ylim = c(GSPTT_premig_ymin, GSPTT_premig_ymax),
           expand = FALSE) +
  theme_bw() +
  labs(title = "D)",
       x = "Longitude",
       y = "Latitude") +
  scale_fill_manual(name = "Contours", values = c("Locations" = "black",
                                                  "Full Corridor" = "red", 
                                                  "Core Corridor" = "blue"),
                    guide = guide_legend(override.aes = list(alpha = 0.5))) +
  scale_colour_manual(name = "Contours", values = c("Locations" = "black",
                                                    "Full Corridor" = "red",
                                                    "Core Corridor" = "blue")) +
  theme(legend.position = "none")

plot(GSPTT_premig_DBBMM)










###### Concave Hull Shape ######

GSPTT_premig_hull_points <- st_transform(GSPTT_premig_sf_points, crs = premig_laea_crs)

# Concave Hull function (Concavity controls how sharp it is)
GSPTT_premig_hull <- concaveman(GSPTT_premig_hull_points, concavity = 4)

GSPTT_premig_hull_plot <- ggplot() +
  geom_sf(aes(), data = world_premig, fill = "grey") +
  geom_sf(data = GSPTT_premig_hull, aes(fill = "Concave Hull", colour = "Concave Hull"), alpha = 0.5) +
  geom_sf(aes(colour = "Locations"), data = GSPTT_premig_sf_points, show.legend = "point") +
  coord_sf(xlim = c(GSPTT_premig_xmin, GSPTT_premig_xmax), 
           ylim = c(GSPTT_premig_ymin, GSPTT_premig_ymax),
           expand = FALSE) +
  theme_bw() +
  scale_fill_manual(name = "Contours", values = c("Locations" = "black",
                                                 "Concave Hull" = "red")) +
  scale_colour_manual(name = "Contours", values = c("Locations" = "black",
                                                    "Concave Hull" = "red")) +
  labs(title = "E)",
       x = "Longitude",
       y = "Latitude") + 
  theme(legend.position = "none")

plot(GSPTT_premig_hull_plot)  





  
###### All plots together ######

### Areas of each of the polygons

## KUD Areas
premig_KUD_50_poly <- st_transform(premig_contour_50_sf, crs = premig_laea_crs)
premig_KUD_95_poly <- st_transform(premig_contour_95_sf, crs = premig_laea_crs)


## Density Raster
## Full density raster

# Turns raster binary
DR_premig_density_95 <- premig_clean_raster
crs(DR_premig_density_95) <- premig_laea_crs
DR_premig_shape <- classify(DR_premig_density_95, rcl = matrix(c(0, Inf, 1), ncol = 3, byrow = TRUE)) 

# Masks and removes the 0 values
DR_premig_masked <- mask(DR_premig_shape, DR_premig_shape, maskvalue = 0)

# Turns into polygon
DR_premig_polygons <- as.polygons(DR_premig_masked)

# Turns into sf
DR_premig_poly_sf <- st_as_sf(DR_premig_polygons)
DR_premig_poly <- st_transform(DR_premig_poly_sf, crs = premig_laea_crs)

# Checks the area
DR_premig_poly$area <- st_area(DR_premig_poly)
DR_premig_area <- sum(as.numeric(DR_premig_poly$area), na.RM = TRUE)


## Core Corridor of Density Raster
premig_r_density_filt <- premig_core_use
crs(premig_r_density_filt) <- premig_laea_crs


# Turns raster binary
DR_premig_filt_shape <- classify(premig_r_density_filt, rcl = matrix(c(0, Inf, 1), ncol = 3, byrow = TRUE)) 

# Masks and removes the 0 values
DR_premig_filt_masked <- mask(DR_premig_filt_shape, DR_premig_filt_shape, maskvalue = 0)

# Turns into polygon
DR_premig_filt_polygons <- as.polygons(DR_premig_filt_masked, na.rm = TRUE)

# Turns into sf
DR_premig_filt_poly_sf <- st_as_sf(DR_premig_filt_polygons)
DR_premig_filt_poly <- st_transform(DR_premig_filt_poly_sf, crs = premig_laea_crs)

# Checks the area
DR_premig_filt_poly$area <- st_area(DR_premig_filt_poly)
DR_premig_filt_area <- sum(as.numeric(DR_premig_filt_poly$area), na.RM = TRUE)



## LDE Areas
premig_LDE_50_poly <- st_transform(premig_LDE_50_poly, crs = premig_laea_crs)
premig_LDE_95_poly <- st_transform(premig_LDE_95_poly, crs = premig_laea_crs)



## DBBMM Areas
premig_c95_sf <- st_transform(premig_c95_sf , crs = premig_laea_crs)
premig_c50_sf <- st_transform(premig_c50_sf, crs = premig_laea_crs)


# Hull
premig_hull_poly <- st_transform(GSPTT_premig_hull, crs = premig_laea_crs)








# All 5 together
grid.arrange(GSPTT_premig_kud, 
             GSPTT_premig_density_raster, 
             GSPTT_premig_LDE,
             GSPTT_premig_DBBMM,
             GSPTT_premig_hull_plot,
             ncol = 3,
             nrow = 2)


# Overlapping the 95% confidence interval


GSPTT_95_contours_plot <- ggplot() +
  geom_sf(aes(), data = world_premig, fill = "grey") +
  geom_sf(data = GSPTT_premig_hull, aes(fill = "Concave Hull"), colour = "orange", alpha = 0.1) +
  geom_sf(data = premig_contour_95_sf, aes(fill = "Kernel Utilisation Distribution"), colour = "red", alpha = 0.1) +
  geom_sf(data = premig_LDE_95_poly, aes(fill = "Line Density Estimation"), colour = "blue", alpha = 0.1) +
  geom_sf(data = premig_c95_poly, aes(fill = "Dynamic Brownian Bridge Movement Model"), colour = "green", alpha = 0.1) +
  geom_sf(data = DR_premig_poly, aes(fill = "Density Raster"), colour = "black", alpha = 0.1) +
  coord_sf(xlim = c(GSPTT_premig_xmin, GSPTT_premig_xmax), 
           ylim = c(GSPTT_premig_ymin, GSPTT_premig_ymax), 
           expand = FALSE) +
  scale_fill_manual(name = "Corridor Delineation Method",
                    values = c("Concave Hull" = "orange",
                               "Kernel Utilisation Distribution" = "red",
                               "Line Density Estimation" = "blue",
                               "Dynamic Brownian Bridge Movement Model" = "green",
                               "Density Raster" = "black"),
                    guide = guide_legend(override.aes = list(alpha = 0.1))) +
  theme_bw() +
  labs(title = "F)",
    x = "Longitude",
    y = "Latitude") +
  theme(legend.position = "none")

plot(GSPTT_95_contours_plot)



## Overlapping 50%

GSPTT_50_contours_plot <- ggplot() +
  geom_sf(aes(), data = world_premig, fill = "grey") +
  geom_sf(data = premig_contour_50_sf, aes(fill = "Kernel Utilisation Distribution"),  colour = "red", alpha = 0.1) +
  geom_sf(data = premig_LDE_50_poly, aes(fill = "Line Density Estimation"),  colour = "blue", alpha = 0.1) +
  geom_sf(data = premig_c50_poly, aes(fill = "Dynamic Brownian Bridge Movement Model"),  colour = "green", alpha = 0.1) +
  geom_sf(data = DR_premig_filt_poly, aes(fill = "Density Raster"),  colour = "black", alpha = 0.1) +
  coord_sf(xlim = c(GSPTT_premig_xmin, GSPTT_premig_xmax), 
           ylim = c(GSPTT_premig_ymin, GSPTT_premig_ymax), 
           expand = FALSE) +
  scale_fill_manual(name = "Corridor Delineation Method",
                    values = c("Kernel Utilisation Distribution" = "red",
                               "Line Density Estimation" = "blue",
                               "Dynamic Brownian Bridge Movement Model" = "green",
                               "Density Raster" = "black"),
                    guide = guide_legend(override.aes = list(alpha = 0.1))) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude") 
plot(GSPTT_50_contours_plot)







##### 6 plots together #####

grid.arrange(GSPTT_premig_kud, 
             GSPTT_premig_density_raster, 
             GSPTT_premig_LDE,
             GSPTT_premig_DBBMM,
             GSPTT_premig_hull_plot,
             GSPTT_95_contours_plot,
             ncol = 2,
             nrow = 3)







###### Areas and Comparisons ######

## function for basic overlaps


basic_overlap_info <- function(polygons_sf, id_col = "ID") {
  n <- nrow(polygons_sf)
  
  # Ensure geometries are valid
  polygons_sf <- st_make_valid(polygons_sf)
  
  # Calculate area (in km²) for each polygon once
  polygon_areas_km2 <- data.frame(
    ID = as.character(polygons_sf[[id_col]]),
    area_km2 = as.numeric(st_area(polygons_sf)) / 1e6,
    stringsAsFactors = FALSE
  )
  
  # Output table
  results <- data.frame(
    from = character(),
    to = character(),
    overlap_area_km2 = numeric(),
    percent_overlap = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n) {
    for (j in 1:n) {
      
      # Skip self-comparison
      if (i == j) next
      
      poly_i <- polygons_sf[i, ]
      poly_j <- polygons_sf[j, ]
      
      # Check for valid geometries before intersecting
      if (st_is_empty(poly_i) || st_is_empty(poly_j)) next
      
      # Try intersection
      intersection <- tryCatch({
        st_intersection(st_geometry(poly_i), st_geometry(poly_j))
      }, error = function(e) {
        message(sprintf("Intersection failed between %s and %s", poly_i[[id_col]], poly_j[[id_col]]))
        return(NULL)
      })
      
      # If the intersection failed or is empty
      if (is.null(intersection) || length(intersection) == 0 || st_is_empty(intersection)) {
        overlap_area <- 0
        percent <- 0
      } else {
        overlap_area <- as.numeric(st_area(intersection))
        area_j <- as.numeric(st_area(poly_j))
        percent <- round(100 * overlap_area / area_j, 2)
      }
      
      # Add to results
      results <- rbind(results, data.frame(
        from = as.character(poly_i[[id_col]]),
        to = as.character(poly_j[[id_col]]),
        overlap_area_km2 = overlap_area/1000000,
        percent_overlap = percent,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Return both the overlap results and the area per polygon
  return(list(
    overlap_table = results,
    polygon_areas = polygon_areas_km2))
}


# 95% comparisons
premig_95_area_polys <- list(premig_KUD_95_poly, 
                            DR_premig_poly, 
                            premig_LDE_95_poly, 
                            premig_c95_poly,
                            premig_hull_poly)
premig_95_area_geoms <- lapply(premig_95_area_polys, st_geometry)
premig_95_area_geoms_sfc <- do.call(c, premig_95_area_geoms)
premig_95_area_combined <- st_sf(ID = paste0("poly_", seq_along(premig_95_area_geoms_sfc)), 
                                 geometry = premig_95_area_geoms_sfc)


# Print info
premig_95_overlap_info <- basic_overlap_info(premig_95_area_combined)
print(premig_95_overlap_info)



# 50% comparisons
premig_50_area_polys <- list(premig_KUD_50_poly,
                             DR_premig_filt_poly,
                             premig_LDE_50_poly, 
                             premig_c50_poly)
premig_50_area_geoms <- lapply(premig_50_area_polys, st_geometry)
premig_50_area_geoms_sfc <- do.call(c, premig_50_area_geoms)
premig_50_area_combined <- st_sf(ID = paste0("poly_", seq_along(premig_50_area_geoms_sfc)), 
                                 geometry = premig_50_area_geoms_sfc)


# Print info
premig_50_overlap_info <- basic_overlap_info(premig_50_area_combined)
print(premig_50_overlap_info)





## KUD 95 to 50

premig_KUD_area_polys <- list(premig_KUD_50_poly, 
                             premig_KUD_95_poly)
premig_KUD_area_geoms <- lapply(premig_KUD_area_polys, st_geometry)
premig_KUD_area_geoms_sfc <- do.call(c, premig_KUD_area_geoms)
premig_KUD_area_combined <- st_sf(ID = paste0("poly_", seq_along(premig_KUD_area_geoms_sfc)), 
                                 geometry = premig_KUD_area_geoms_sfc)


# Print info
premig_KUD_overlap_info <- basic_overlap_info(premig_KUD_area_combined)
print(premig_KUD_overlap_info)





## LDE 95 to 50

premig_LDE_area_polys <- list(premig_LDE_50_poly, 
                              premig_LDE_95_poly)
premig_LDE_area_geoms <- lapply(premig_LDE_area_polys, st_geometry)
premig_LDE_area_geoms_sfc <- do.call(c, premig_LDE_area_geoms)
premig_LDE_area_combined <- st_sf(ID = paste0("poly_", seq_along(premig_LDE_area_geoms_sfc)), 
                                  geometry = premig_LDE_area_geoms_sfc)


# Print info
premig_LDE_overlap_info <- basic_overlap_info(premig_LDE_area_combined)
print(premig_LDE_overlap_info)





## DBBMM 95 to 50

premig_DBBMM_area_polys <- list(premig_c95_poly, 
                              premig_c50_poly)
premig_DBBMM_area_geoms <- lapply(premig_DBBMM_area_polys, st_geometry)
premig_DBBMM_area_geoms_sfc <- do.call(c, premig_DBBMM_area_geoms)
premig_DBBMM_area_combined <- st_sf(ID = paste0("poly_", seq_along(premig_DBBMM_area_geoms_sfc)), 
                                  geometry = premig_DBBMM_area_geoms_sfc)


# Print info
premig_DBBMM_overlap_info <- basic_overlap_info(premig_DBBMM_area_combined)
print(premig_DBBMM_overlap_info)





## Density Raster full to part

premig_DR_area_polys <- list(DR_premig_poly, 
                                DR_premig_filt_poly)
premig_DR_area_geoms <- lapply(premig_DR_area_polys, st_geometry)
premig_DR_area_geoms_sfc <- do.call(c, premig_DR_area_geoms)
premig_DR_area_combined <- st_sf(ID = paste0("poly_", seq_along(premig_DR_area_geoms_sfc)), 
                                    geometry = premig_DR_area_geoms_sfc)


# Print info
premig_DR_overlap_info <- basic_overlap_info(premig_DR_area_combined)
print(premig_DR_overlap_info)



###### Spatial Comparison ######
#  We need to turn them all into the same format to be able to compare.

# Empty raster template
premig_comp_rast_grid <- rast(ext(GSPTT_premig_xmin, GSPTT_premig_xmax, 
                           GSPTT_premig_ymin, GSPTT_premig_ymax),
                           res = 5000,
                           crs = premig_laea_crs)




# resample the LDE 
premig_LDE_rast_re <- terra::resample(premig_LDE_raster, premig_comp_rast_grid, method = "bilinear")
premig_LDE_rast_re[is.na(premig_LDE_rast_re)] <- 0 #Ensure that NA cells are no resampled wrong

# resample the KUD
premig_KUD_rast <- rast(premig_kud)
premig_KUD_rast_re <- terra::resample(premig_KUD_rast, premig_comp_rast_grid, metho = "bilinear")
premig_KUD_rast_re[is.na(premig_KUD_rast_re)] <- 0 #Ensure that NA cells are no resampled wrong

# resample the DBBMM 
premig_DBBMM_rast <- rast(premig_r_reprojected)
premig_DBBMM_rast_re <- terra::resample(premig_DBBMM_rast, premig_comp_rast_grid, method = "bilinear")
premig_DBBMM_rast_re[is.na(premig_DBBMM_rast_re)] <- 0 #Ensure that NA cells are no resampled wrong


###### Boot overlap ######

premig_UD_list <- list(premig_KUD_rast_re, premig_LDE_rast_re, premig_DBBMM_rast_re)


# BA
# 50%

premig_BA_Stats_50 <- boot_overlap(
  premig_UD_list,
  R = 20,
  method = "BA",
  percent = 50,
  quantiles = c(0.25, 0.5, 0.75)
)


print(premig_BA_Stats_50)


# 95%
premig_BA_Stats_95 <- boot_overlap(
  premig_UD_list,
  R = 20,
  method = "BA",
  percent = 95,
  quantiles = c(0.25, 0.5, 0.75)
)

print(premig_BA_Stats_95)
















##### GSPTT Migration #####


###### Filtering ######  
  

GSPTT_mig <- merge(df_locs, df_mpm) %>% 
    filter(date > as.POSIXct("2010-05-01 00:00:00", tz = "UTC")) %>%   #Filters dates before migration
    filter(date < as.POSIXct("2010-08-01 00:00:00", tz = "UTC")) %>%  #Filters that one bird that returns
    filter(lat < 40) #Filters so max is 40 lat

  
# Filters to start when move value is >0.9 and latitude is increasing for more than 5 occurances

GSPTT_mig <- GSPTT_mig %>%
  group_by(id) %>%
  group_modify(~{
      
    dat <- .x %>% arrange(date)   # make sure ordered properly
      
    # Check lat increasing
    lat_increase <- c(FALSE, diff(dat$lat) > 0)
    
    # Move Value
    move_condition <- dat$g > 0.9   #0.9 is the move value
    
    # Must validate both
    condition <- lat_increase & move_condition
    
    # run length encoding
    r <- rle(condition)
    
    # find when run length is true with conditions
    idx <- which(r$values & r$lengths >= 5)[1]    # Change 5 to whatever for consecutives
    
    # if none found return nothing
    if (is.na(idx)) return(dat[0, ])
    
    # starting row of that run
    start_row <- sum(r$lengths[1:(idx - 1)]) + 1
    
    # keep everything from that row onward (Since most end)
    dat[start_row:nrow(dat), ]
  }) %>%
  ungroup()
  

#Orders for the DBBMM
GSPTT_mig <-GSPTT_mig[order(GSPTT_mig$id, GSPTT_mig$date), ] #%>%  #Orders it by id

GSPTT_mig <- GSPTT_mig %>% 
  mutate(SE = ifelse(x.se > y.se, y.se, x.se), #Adds a column of Standard Error, which is the smallest value of the error columns
         citation = NA,                        #Adds a column of NA citation values for DBBMM
         species = "bird")                #Adds a column of "species" for DBBMM


GSPTT_mig <- as.data.frame(GSPTT_mig)

#Set as sf
GSPTT_mig_sf_points <- st_as_sf(GSPTT_mig, coords = c("lon", "lat"), crs = 4326)


#Plot
ggplot() +
  theme_bw() +
  geom_sf(aes(), data = world, fill = "grey") +
  geom_sf(aes(colour = id), data = GSPTT_mig_sf_points, show.legend = "point") +
  coord_sf(xlim = c(-70, -20), ylim = c(-50, 45)) +
  theme(legend.position = "none") +
  labs(x = "Longitude",
       y = "Latitude") 





##### Migration Methods #####


# Set the CRS for the laea for later

mig_lat_centre <- laea_centre(GSPTT_mig$lat)
mig_lon_centre <- laea_centre(GSPTT_mig$lon)

mig_laea_crs <- sprintf("+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m", 
                           mig_lat_centre, 
                           mig_lon_centre)

# World Map for mig

world_mig <- st_transform(world, crs = mig_laea_crs)

# Set the extent for the maps


GSPTT_mig_ext <- SpatialPoints(coords = cbind(GSPTT_mig$lon, GSPTT_mig$lat), proj4string = CRS("WGS84")) 
GSPTT_mig_ext <- spTransform(GSPTT_mig_ext, mig_laea_crs)


GSPTT_mig_xmin <- min(GSPTT_mig_ext@coords[,1])* 1.5
GSPTT_mig_xmax <- max(GSPTT_mig_ext@coords[,1])* 1.5
GSPTT_mig_ymin <- min(GSPTT_mig_ext@coords[,2])* 1.2
GSPTT_mig_ymax <- max(GSPTT_mig_ext@coords[,2])* 1.2


# Average error for the points
mean_SE_mig <- mean(GSPTT_mig$SE)
print(mean_SE_mig)


#Number of individuals used
n_distinct(GSPTT_mig$id)




###### Kernel Utiilisation Distribution ######

# Change the data to be spatial points
GSPTT_migKUD <- SpatialPoints(coords = cbind(GSPTT_mig$lon, GSPTT_mig$lat), 
                                 proj4string = CRS("+proj=longlat +datum=WGS84")) 
proj4string(GSPTT_migKUD) <- CRS("+proj=longlat +datum=WGS84")
GSPTT_migKUD <- spTransform(GSPTT_migKUD, mig_laea_crs)


## Calc the grid size needed

# Extract bounding box
mig_bbox <- sp::bbox(GSPTT_migKUD)

# Calc max range 
mig_x_range <- mig_bbox[1, "max"] - mig_bbox[1, "min"]
mig_y_range <- mig_bbox[2, "max"] - mig_bbox[2, "min"]

# Calc the max of the two
mig_max_range <- max(mig_x_range, mig_y_range)

# Calc the grid value from this
# the number represents the grid cell size in m that you want
# Due to the kernelUD taking an integer, it will round this value to the nearest
# whole number but it will be close enough
# if the grid number ends up super low, use the ceiling function to round it up
mig_grid_kud <- mig_max_range / 20000



## Apply the KUD

# The default kernel (bivariate normal kernel) is used by default in the adehabitatHR package
# smoothing parameter h is the width of the kernels and can be calculated
# Sigma (σ) is the standard deviation of points
# σ can be calculated by the mean of the standard error of the x or y multiplied 
# by the square root of the number of locations

mig_kud_h = mean_SE_mig*1000

mig_kud <- kernelUD(GSPTT_migKUD, h = mig_kud_h, grid = mig_grid_kud)


# Extract the Contours you want
mig_contour_50 <- getverticeshr(mig_kud, percent = 50)
mig_contour_75 <- getverticeshr(mig_kud, percent = 75)
mig_contour_95 <- getverticeshr(mig_kud, percent = 95)

# Convert contours to sf objects
mig_contour_50_sf <- st_as_sf(mig_contour_50)
mig_contour_75_sf <- st_as_sf(mig_contour_75)
mig_contour_95_sf <- st_as_sf(mig_contour_95)


# Turn the kud into a dataframe that ggplot can read
mig_Kuddf <- as.data.frame(mig_kud, xy = TRUE)




### Plot
GSPTT_mig_kud <- ggplot() +
  geom_sf(aes(), data = world_mig, fill = "grey") +
  geom_sf(data = mig_contour_95_sf, aes(fill = "Full Corridor", colour = "Full Corridor"), alpha = 0.5) +
  geom_sf(data = mig_contour_50_sf, aes(fill = "Core Corridor", colour = "Core Corridor"), alpha = 0.5) + # Need a new geom_sf for each contour you add
  geom_sf(aes(colour = "Locations"), data = GSPTT_mig_sf_points, show.legend = "point") +
  coord_sf(xlim = c(GSPTT_mig_xmin, GSPTT_mig_xmax), 
           ylim = c(GSPTT_mig_ymin, GSPTT_mig_ymax), 
           expand = FALSE) +
  theme_bw() +
  labs(
    title = "A)",
    x = "Longitude",
    y = "Latitude") +
  scale_fill_manual(name = "Probability Contours", values = c("Locations" = "black",
                                                              "Full Corridor" = "red",
                                                              "Core Corridor" = "blue"),
                    guide = guide_legend(override.aes = list(alpha = 0.5))) +
  scale_colour_manual(name = "Probability Contours", values = c("Locations" = "black",
                                                                "Full Corridor" = "red", 
                                                                "Core Corridor" = "blue")) +
  theme(legend.position = "none")

plot(GSPTT_mig_kud)










###### Density Raster ######

# Make sf points in crs 
GSPTT_mig_sf_laea <- st_transform(GSPTT_mig_sf_points, mig_laea_crs)

# Choose grid cell size (in m)
mig_grid_size <- mean_SE_mig * 1000


# Convert points to SpatVector
mig_points_vect <- vect(GSPTT_mig_sf_laea)
mig_points_vect <- project(mig_points_vect, mig_laea_crs)

# Create an empty raster in laea
mig_density_raster <- rast(ext(mig_points_vect), resolution = mig_grid_size, crs = mig_laea_crs)

# Rasterize with a density count
mig_r_density <- rasterize(mig_points_vect, mig_density_raster, 
                              fun = "length", background = 0)

# Convert to a dataframe for plotting in ggplot
mig_df_raster <- as.data.frame(mig_r_density, xy = TRUE)
colnames(mig_df_raster) <- c("x", "y", "Density")
mig_df_raster <- filter(mig_df_raster, Density > 0) %>% 
  filter(!is.na(Density))

# Calculate the top 10% bin values
mig_top_10 <- ceiling(0.10 * (nrow(mig_df_raster)))

mig_df_raster$Density <- as.numeric(mig_df_raster$Density)

# Order in descending order and bin
mig_df_raster <- mig_df_raster %>% 
  arrange(desc(Density)) %>% 
  mutate(rank = row_number()) %>% 
  mutate(bin = case_when(
    rank <= mig_top_10 ~ "Core Use",
    TRUE ~ "Corridor")) %>% 
  mutate(bin_num = case_when(
    rank <= mig_top_10 ~ 2,
    TRUE ~ 1))


## Removing noise

# re raster it
mig_cat_raster <- rast(mig_df_raster[, c("x", "y", "bin_num")], type = "xyz")

# remove the isolated points
mig_raster_mask <- mig_cat_raster == 1   # Just pulls the corridor points
mig_raster_patches <- patches(mig_raster_mask, directions = 8) #  Creates patches based on whats touching in 8 directions
mig_raster_sizes <- freq(mig_raster_patches) # Counts the number of cells in each patch and puts into a df
mig_raster_ids_keep <- mig_raster_sizes$value[mig_raster_sizes$count > 2]  # Filters mig any patches with number of cells less than value
mig_clean_corridor <- mask(mig_raster_patches, mig_raster_patches %in% mig_raster_ids_keep, maskvalues = FALSE) # Makes a raster that only keeps the clean corridor
mig_clean_corridor[!is.na(mig_clean_corridor)] <- 1

# Extract the core use values to merge
mig_core_use <- mig_cat_raster
mig_core_use[mig_core_use != 2] <- NA

# Merge the rasters together for plotting
mig_clean_raster <- mosaic(mig_clean_corridor, mig_core_use, fun = "max")

# Turn into a df for plotting
mig_clean_df <- as.data.frame(mig_clean_raster, xy = TRUE)
colnames(mig_clean_df) <- c("x", "y", "bin_num")
mig_clean_df <- mig_clean_df %>% 
  mutate(bin = factor(bin_num,
                      levels = c(1, 2),
                      labels = c("Corridor", "Core Use")))



# Plot in ggplot
GSPTT_mig_density_raster <- ggplot() +
  geom_sf(aes(), data = world_mig, fill = "grey") +
  geom_tile(data = mig_df_raster, aes(x = x, y = y, fill = bin), colour = "black", alpha = 0.5) +
  geom_sf(aes(colour = "Locations"), data = GSPTT_mig_sf_points, show.legend = "point") +
  coord_sf(xlim = c(GSPTT_mig_xmin, GSPTT_mig_xmax), 
           ylim = c(GSPTT_mig_ymin, GSPTT_mig_ymax),
           expand = FALSE) +
  scale_fill_manual(values = c("Locations" = "black",
                               "Core Use" = "blue",
                               "Corridor" = "red"),
                    na.value = "transparent") +  # Transparent for missing values
  scale_colour_manual(values = c("Locations" = "black")) +
  labs(title = "B)",
       x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  theme(legend.position = "none")


plot(GSPTT_mig_density_raster)







###### Line Density Estimation ######
# Turn the data into a format that can be exported into ARC



# Create line strings
GSPTT_mig_sf_lines <- GSPTT_mig_sf_points %>%
  group_by(id) %>%
  summarise(do_union = FALSE) %>%  #  Avoid merging separate tracks
  st_cast("LINESTRING")  #  Convert grouped points into lines

# Reproject to aeqd to be in km
proj_mig_lines <- st_transform(GSPTT_mig_sf_lines, crs = mig_laea_crs)

# Finding bbox
proj_mig_lines_bbox <- st_bbox(proj_mig_lines)
print(proj_mig_lines_bbox)

# Save as shapefile
st_write(obj = proj_mig_lines,
         dsn = "~/GitHub/honours_corridors/Outputs/LDE_PTT_mig.shp",
         driver = "ESRI Shapefile",
         append = FALSE)

# Load this in ARC in the Line Density fucntion from the Spatial Analyst Toolbox
# For search radius use the average or median Standard Error, whichever is more appropriate
# Use a grid cell size appropriate based on this error size and the distance of the journey
print(mean_SE_mig)

# Write the cell size you used here:
#              20 000m

## Load the newly created LDE

mig_LDE_raster <- rast("~/GitHub/honours_corridors/Outputs/LDE_PTT_Outward_mig.tif")

# Create the Contours



# Polygon 95
mig_LDE_95 <- mig_LDE_raster
mig_LDE_95[mig_LDE_95 < 0.05] <- NA
mig_LDE_95[mig_LDE_95 >= 0.05] <- 1

mig_LDE_95_poly <- terra::as.polygons(mig_LDE_95) %>% 
  sf::st_as_sf(crs = mig_laea_crs) %>%
  sf::st_buffer(dist = res(mig_LDE_95)[1]/4) %>%
  sf::st_simplify(dTolerance = res(mig_LDE_95)[1]) 


# Polygon 50
mig_LDE_50 <- mig_LDE_raster
mig_LDE_50[mig_LDE_50 < 0.50] <- NA
mig_LDE_50[mig_LDE_50 >= 0.50] <- 1

mig_LDE_50_poly <- terra::as.polygons(mig_LDE_50) %>% 
  sf::st_as_sf(crs = mig_laea_crs) %>%
  sf::st_buffer(dist = res(mig_LDE_50)[1]/4) %>%
  sf::st_simplify(dTolerance = res(mig_LDE_50)[1])



# Plot as a ggplot
GSPTT_mig_LDE <- ggplot() +
  geom_sf(aes(), data = world_mig, fill = "grey") +
  geom_sf(data = mig_LDE_95_poly, aes(fill = "Full Corridor", colour = "Full Corridor"), alpha = 0.5) +
  geom_sf(data = mig_LDE_50_poly, aes(fill = "Core Corridor", colour = "Core Corridor"), alpha = 0.5) +
  geom_sf(aes(colour = "Locations"), data = GSPTT_mig_sf_points, show.legend = "point") +
  coord_sf(xlim = c(GSPTT_mig_xmin, GSPTT_mig_xmax), 
           ylim = c(GSPTT_mig_ymin, GSPTT_mig_ymax),
           expand = FALSE) +
  labs(title = "C)",
       x = "Longitude",
       y = "Latitude") +
  scale_fill_manual(name = "Contours", values = c("Locations" = "black",
                                                  "Full Corridor" = "red", 
                                                  "Core Corridor" = "blue"),
                    guide = guide_legend(override.aes = list(alpha = 0.5))) +
  scale_colour_manual(name = "Contours", values = c("Locations" = "black",
                                                    "Full Corridor" = "red",
                                                    "Core Corridor" = "blue")) +
  theme_bw() +
  theme(legend.position = "none")



plot(GSPTT_mig_LDE)








###### Dynamic Brownian Bridge Movement Model ######
# Turn it into a move data for the DBBMM
GSPTT_mig_BBM <- move(x = GSPTT_mig$lon, y = GSPTT_mig$lat,
                         time = as.POSIXct(GSPTT_mig$date, format="%Y-%m-%d %H:%M:%OS", tz="UTC"), 
                         proj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
                         data = GSPTT_mig, animal = GSPTT_mig$id)

# Reproject into a metric measurment CRS for the DBBMM
GSPTT_mig_proj <- spTransform(GSPTT_mig_BBM, mig_laea_crs)


## Create the DBBMM

# Increasing window size increased the reliability of the brownian motion values
# but at the cost of missing rapid changes
# Increasing the margin size increases the likelihood of identifying relatively weak breakpoints
# but at cost of potentially not detecting breakpoints due to limited observations in w
# Choices should be made based on the time interval you are looking at
DBB_GSPTT_mig <- brownian.bridge.dyn(GSPTT_mig_proj, 
                                        ext = 1.5, # Needs to be large enough to encompass the data. Increase until does
                                        raster = 20000, # raster box size in m
                                        location.error = mean_SE_mig * 1000, # error of location data in m
                                        window.size = 13, # Should be >2m, window of the breakpoints for behaviour
                                        margin = 5) # Holds out observations that wont be included in breakpoints from the window, 2m should be at most 1/2 of w

## Obtain the Utilisation Distributions of Each individual
mig_UD <- getVolumeUD(DBB_GSPTT_mig)
# plot(mig_UD, main = "UD")


# Stacks them all up
UD_mig <- stack(mig_UD)
mig_combined_UD <- calc(UD_mig, sum)


## Standardise the raster
mig_r_min <- cellStats(mig_combined_UD, "min")
mig_r_max <- cellStats(mig_combined_UD, "max")

mig_r_standardised <- (mig_combined_UD - mig_r_min) / (mig_r_max - mig_r_min)

# Reproject the raster
mig_r_reprojected <- projectRaster(mig_r_standardised, crs = mig_laea_crs)

# Turn it into a data frame for ggplot
mig_new_raster <- as.data.frame(mig_r_reprojected, xy = TRUE)
colnames(mig_new_raster) <- c("x", "y", "value")


# Set data frame as simple features and so ggplot can read it
mig_new_raster_sf <- st_as_sf(mig_new_raster, coords = c("x", "y"))
mig_new_raster_sf <- mig_new_raster_sf %>% st_set_crs(mig_laea_crs)
mig_new_raster_sf <- subset(mig_new_raster_sf, value <= 0.9995)


## Creating Contour layers
mig_c95 <- rasterToContour(mig_r_reprojected, levels = 0.95)
mig_c95_sf <- st_as_sf(mig_c95)
mig_c95_sf <- st_transform(mig_c95_sf, mig_laea_crs)
mig_c95_poly <- st_cast(st_union(mig_c95_sf), "POLYGON")

mig_c50 <- rasterToContour(mig_r_reprojected, levels = 0.50)
mig_c50_sf <- st_as_sf(mig_c50)
mig_c50_sf <- st_transform(mig_c50_sf, mig_laea_crs)
mig_c50_poly <- st_cast(st_union(mig_c50_sf), "POLYGON")


## Plot

GSPTT_mig_DBBMM <- ggplot() +
  geom_sf(aes(), data = world_mig, fill = "grey") +
  geom_sf(data = mig_c95_poly, aes(fill = "Full Corridor", colour = "Full Corridor"), alpha = 0.5) +
  geom_sf(data = mig_c50_poly, aes(fill = "Core Corridor", colour = "Core Corridor"), alpha = 0.5) +
  geom_sf(aes(colour = "Locations"), data = GSPTT_mig_sf_points, show.legend = "point") +
  coord_sf(xlim = c(GSPTT_mig_xmin, GSPTT_mig_xmax), 
           ylim = c(GSPTT_mig_ymin, GSPTT_mig_ymax),
           expand = FALSE) +
  theme_bw() +
  labs(title = "D)",
       x = "Longitude",
       y = "Latitude") +
  scale_fill_manual(name = "Contours", values = c("Locations" = "black",
                                                  "Full Corridor" = "red", 
                                                  "Core Corridor" = "blue"),
                    guide = guide_legend(override.aes = list(alpha = 0.5))) +
  scale_colour_manual(name = "Contours", values = c("Locations" = "black",
                                                    "Full Corridor" = "red",
                                                    "Core Corridor" = "blue")) +
  theme(legend.position = "none")

plot(GSPTT_mig_DBBMM)

















###### Concave Hull Shape ######

GSPTT_mig_hull_points <- st_transform(GSPTT_mig_sf_points, crs = mig_laea_crs)

GSPTT_mig_hull <- concaveman(GSPTT_mig_hull_points, concavity = 4)

GSPTT_mig_hull_plot <- ggplot() +
  geom_sf(aes(), data = world_mig, fill = "grey") +
  geom_sf(data = GSPTT_mig_hull, aes(colour = "Concave Hull", fill = "Concave Hull"), alpha =0.5) +
  geom_sf(aes(colour = "Locations"), data = GSPTT_mig_sf_points, show.legend = "point") +
  coord_sf(xlim = c(GSPTT_mig_xmin, GSPTT_mig_xmax), 
           ylim = c(GSPTT_mig_ymin, GSPTT_mig_ymax),
           expand = FALSE) +
  theme_bw() +
  scale_fill_manual(name = "Contours", values = c("Locations" = "black",
                                                  "Concave Hull" = "red")) +
  scale_colour_manual(name = "Contours", values = c("Locations" = "black",
                                                    "Concave Hull" = "red")) +
  labs(title = "E)",
       x = "Longitude",
       y = "Latitude") + 
  theme(legend.position = "none")


plot(GSPTT_mig_hull_plot)  





##### All plots together #####

# All together
grid.arrange(GSPTT_mig_kud, 
             GSPTT_mig_density_raster, 
             GSPTT_mig_LDE,
             GSPTT_mig_DBBMM,
             GSPTT_mig_hull_plot,
             ncol = 3,
             nrow = 2)












###### Areas and Comparisons ######

### Areas of each of the polygons

## KUD Areas
mig_KUD_50_poly <- st_transform(mig_contour_50_sf, crs = mig_laea_crs)
mig_KUD_95_poly <- st_transform(mig_contour_95_sf, crs = mig_laea_crs)



## Density Raster

# Turns raster binary
DR_mig_density_95 <- mig_clean_raster
crs(DR_mig_density_95) <- mig_laea_crs
DR_mig_shape <- classify(DR_mig_density_95, rcl = matrix(c(0, Inf, 1), ncol = 3, byrow = TRUE)) 
# Masks and removes the 0 values
DR_mig_masked <- mask(DR_mig_shape, DR_mig_shape, maskvalue = 0)
# Turns into polygon
DR_mig_polygons <- as.polygons(DR_mig_masked)
# Turns into sf
DR_mig_poly_sf <- st_as_sf(DR_mig_polygons)
DR_mig_poly <- st_transform(DR_mig_poly_sf, crs = mig_laea_crs)
# Checks the area
DR_mig_poly$area <- st_area(DR_mig_poly)
DR_mig_area <- sum(as.numeric(DR_mig_poly$area), na.RM = TRUE)



## Part density Raster
mig_r_density_filt <- mig_core_use
crs(mig_r_density_filt) <- mig_laea_crs

# Turns raster binary
DR_mig_filt_shape <- classify(mig_r_density_filt, rcl = matrix(c(0, Inf, 1), ncol = 3, byrow = TRUE)) 
# Masks and removes the 0 values
DR_mig_filt_masked <- mask(DR_mig_filt_shape, DR_mig_filt_shape, maskvalue = 0)
# Turns into polygon
DR_mig_filt_polygons <- as.polygons(DR_mig_filt_masked, na.rm = TRUE)
# Turns into sf
DR_mig_filt_poly_sf <- st_as_sf(DR_mig_filt_polygons)
DR_mig_filt_poly <- st_transform(DR_mig_filt_poly_sf, crs = mig_laea_crs)
# Checks the area
DR_mig_filt_poly$area <- st_area(DR_mig_filt_poly)
DR_mig_filt_area <- sum(as.numeric(DR_mig_filt_poly$area), na.RM = TRUE)



## LDE Areas
mig_LDE_50_poly <- st_transform(mig_LDE_50_poly, crs = mig_laea_crs)
mig_LDE_95_poly <- st_transform(mig_LDE_95_poly, crs = mig_laea_crs)



## DBBMM Areas
mig_c95_sf
mig_c50_sf



# hull poly
mig_hull_poly <- st_transform(GSPTT_mig_hull, crs = mig_laea_crs)


# Overlapping the 95% confidence interval

GSPTT_95_mig_contours_plot <- ggplot() +
  geom_sf(aes(), data = world_mig, fill = "grey") +
  geom_sf(data = GSPTT_mig_hull, aes(fill = "Concave Hull"), colour = "orange", alpha = 0.1) +
  geom_sf(data = mig_contour_95_sf, aes(fill = "Kernel Utilisation Distribution"), colour = "red", alpha = 0.1) +
  geom_sf(data = mig_LDE_95_poly, aes(fill = "Line Density Estimation"), colour = "blue", alpha = 0.1) +
  geom_sf(data = mig_c95_poly, aes(fill = "Dynamic Brownian Bridge Movement Model"), colour = "green", alpha = 0.1) +
  geom_sf(data = DR_mig_poly, aes(fill = "Density Raster"), colour = "black", alpha = 0.1) +
  coord_sf(xlim = c(GSPTT_mig_xmin, GSPTT_mig_xmax), 
           ylim = c(GSPTT_mig_ymin, GSPTT_mig_ymax), 
           expand = FALSE) +
  scale_fill_manual(name = "Corridor Delineation Method",
                    values = c("Concave Hull" = "orange",
                               "Kernel Utilisation Distribution" = "red",
                               "Line Density Estimation" = "blue",
                               "Dynamic Brownian Bridge Movement Model" = "green",
                               "Density Raster" = "black"),
                    guide = guide_legend(override.aes = list(alpha = 0.1))) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude") +
  theme(legend.position = "none")

plot(GSPTT_95_mig_contours_plot)




## Overlapping 50%

GSPTT_50_mig_contours_plot <- ggplot() +
  geom_sf(aes(), data = world_mig, fill = "grey") +
  geom_sf(data = mig_contour_50_sf, aes(fill = "Kernel Utilisation Distribution"), colour = "red", alpha = 0.1) +
  geom_sf(data = mig_LDE_50_poly, aes(fill = "Line Density Estimation"), colour = "blue", alpha = 0.1) +
  geom_sf(data = mig_c50_poly, aes(fill = "Dynamic Brownian Bridge Movement Model"), colour = "green", alpha = 0.1) +
  geom_sf(data = DR_mig_filt_poly, aes(fill = "Density Raster"),  colour = "black", alpha = 0.1) +
  coord_sf(xlim = c(GSPTT_mig_xmin, GSPTT_mig_xmax), 
           ylim = c(GSPTT_mig_ymin, GSPTT_mig_ymax), 
           expand = FALSE) +
  scale_fill_manual(name = "Corridor Delineation Method",
                      values = c("Kernel Utilisation Distribution" = "red",
                                 "Line Density Estimation" = "blue",
                                 "Dynamic Brownian Bridge Movement Model" = "green",
                                 "Density Raster" = "black"),
                      guide = guide_legend(override.aes = list(alpha = 0.1))) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude") 
plot(GSPTT_50_contours_plot)





##### 6 plots together #####

# All together
grid.arrange(GSPTT_mig_kud, 
             GSPTT_mig_density_raster, 
             GSPTT_mig_LDE,
             GSPTT_mig_DBBMM,
             GSPTT_mig_hull_plot,
             GSPTT_95_mig_contours_plot,
             ncol = 3,
             nrow = 2)












## function for basic overlaps


basic_overlap_info <- function(polygons_sf, id_col = "ID") {
  n <- nrow(polygons_sf)
  
  # Ensure geometries are valid
  polygons_sf <- st_make_valid(polygons_sf)
  
  # Calculate area (in km²) for each polygon once
  polygon_areas_km2 <- data.frame(
    ID = as.character(polygons_sf[[id_col]]),
    area_km2 = as.numeric(st_area(polygons_sf)) / 1e6,
    stringsAsFactors = FALSE
  )
  
  # Output table
  results <- data.frame(
    from = character(),
    to = character(),
    overlap_area_km2 = numeric(),
    percent_overlap = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n) {
    for (j in 1:n) {
      
      # Skip self-comparison
      if (i == j) next
      
      poly_i <- polygons_sf[i, ]
      poly_j <- polygons_sf[j, ]
      
      # Check for valid geometries before intersecting
      if (st_is_empty(poly_i) || st_is_empty(poly_j)) next
      
      # Try intersection
      intersection <- tryCatch({
        st_intersection(st_geometry(poly_i), st_geometry(poly_j))
      }, error = function(e) {
        message(sprintf("Intersection failed between %s and %s", poly_i[[id_col]], poly_j[[id_col]]))
        return(NULL)
      })
      
      # If the intersection failed or is empty
      if (is.null(intersection) || length(intersection) == 0 || st_is_empty(intersection)) {
        overlap_area <- 0
        percent <- 0
      } else {
        overlap_area <- as.numeric(st_area(intersection))
        area_j <- as.numeric(st_area(poly_j))
        percent <- round(100 * overlap_area / area_j, 2)
      }
      
      # Add to results
      results <- rbind(results, data.frame(
        from = as.character(poly_i[[id_col]]),
        to = as.character(poly_j[[id_col]]),
        overlap_area_km2 = overlap_area/1000000,
        percent_overlap = percent,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Return both the overlap results and the area per polygon
  return(list(
    overlap_table = results,
    polygon_areas = polygon_areas_km2))
}

## the results of these can be interpretted as the "from" polygon is being overlapped
## over the "to" polygon, so the overlap of 1 on 2, then the km2 and the percent

# 95% comparisons
# List geoms
mig_95_area_polys <- list(mig_KUD_95_poly, # KUD
                             DR_mig_poly, # Density Raster
                             mig_LDE_95_poly, # LDE
                             mig_c95_poly, # DBBMM
                             mig_hull_poly) # Concave Hull
mig_95_area_geoms <- lapply(mig_95_area_polys, st_geometry)
mig_95_area_geoms_sfc <- do.call(c, mig_95_area_geoms)
mig_95_area_combined <- st_sf(ID = paste0("poly_", seq_along(mig_95_area_geoms_sfc)), geometry = mig_95_area_geoms_sfc)


# Print info
mig_basic_overlap_info <- basic_overlap_info(mig_95_area_combined)

print(mig_basic_overlap_info)




# 50% comparisons
# List geoms
mig_50_area_polys <- list(mig_KUD_50_poly,
                          DR_mig_filt_poly,
                             mig_LDE_50_poly, 
                             mig_c50_poly)
mig_50_area_geoms <- lapply(mig_50_area_polys, st_geometry)
mig_50_area_geoms_sfc <- do.call(c, mig_50_area_geoms)
mig_50_area_combined <- st_sf(ID = paste0("poly_", seq_along(mig_50_area_geoms_sfc)), 
                                 geometry = mig_50_area_geoms_sfc)


# Print info
mig_50_overlap_info <- basic_overlap_info(mig_50_area_combined)
print(mig_50_overlap_info)



## KUD 95 to 50

mig_KUD_area_polys <- list(mig_KUD_50_poly, 
                              mig_KUD_95_poly)
mig_KUD_area_geoms <- lapply(mig_KUD_area_polys, st_geometry)
mig_KUD_area_geoms_sfc <- do.call(c, mig_KUD_area_geoms)
mig_KUD_area_combined <- st_sf(ID = paste0("poly_", seq_along(mig_KUD_area_geoms_sfc)), 
                                  geometry = mig_KUD_area_geoms_sfc)

# Print info
mig_KUD_overlap_info <- basic_overlap_info(mig_KUD_area_combined)
print(mig_KUD_overlap_info)





## LDE 95 to 50

mig_LDE_area_polys <- list(mig_LDE_50_poly, 
                              mig_LDE_95_poly)
mig_LDE_area_geoms <- lapply(mig_LDE_area_polys, st_geometry)
mig_LDE_area_geoms_sfc <- do.call(c, mig_LDE_area_geoms)
mig_LDE_area_combined <- st_sf(ID = paste0("poly_", seq_along(mig_LDE_area_geoms_sfc)), 
                                  geometry = mig_LDE_area_geoms_sfc)

# Print info
mig_LDE_overlap_info <- basic_overlap_info(mig_LDE_area_combined)
print(mig_LDE_overlap_info)





## DBBMM 95 to 50

mig_DBBMM_area_polys <- list(mig_c95_poly, 
                                mig_c50_poly)
mig_DBBMM_area_geoms <- lapply(mig_DBBMM_area_polys, st_geometry)
mig_DBBMM_area_geoms_sfc <- do.call(c, mig_DBBMM_area_geoms)
mig_DBBMM_area_combined <- st_sf(ID = paste0("poly_", seq_along(mig_DBBMM_area_geoms_sfc)), 
                                    geometry = mig_DBBMM_area_geoms_sfc)


# Print info
mig_DBBMM_overlap_info <- basic_overlap_info(mig_DBBMM_area_combined)
print(mig_DBBMM_overlap_info)





## Density Raster full to part

mig_DR_area_polys <- list(DR_mig_poly, 
                             DR_mig_filt_poly)
mig_DR_area_geoms <- lapply(mig_DR_area_polys, st_geometry)
mig_DR_area_geoms_sfc <- do.call(c, mig_DR_area_geoms)
mig_DR_area_combined <- st_sf(ID = paste0("poly_", seq_along(mig_DR_area_geoms_sfc)), 
                                 geometry = mig_DR_area_geoms_sfc)


# Print info
mig_DR_overlap_info <- basic_overlap_info(mig_DR_area_combined)
print(mig_DR_overlap_info)





###### Spatial Comparison ######
#    We need to turn them all into the same format to be able to compare.
#    preferably all with the same grid cell size

# Empty raster template
mig_comp_rast_grid <- rast(ext(GSPTT_mig_xmin, GSPTT_mig_xmax, 
                                  GSPTT_mig_ymin, GSPTT_mig_ymax),
                              res = 5000,
                              crs = mig_laea_crs)


# resample the LDE 
mig_LDE_rast_re <- terra::resample(mig_LDE_raster, mig_comp_rast_grid, method = "bilinear")
mig_LDE_rast_re[is.na(mig_LDE_rast_re)] <- 0

# resample the KUD
mig_KUD_rast <- rast(mig_kud)
mig_KUD_rast_re <- terra::resample(mig_KUD_rast, mig_comp_rast_grid, metho = "bilinear")
mig_KUD_rast_re[is.na(mig_KUD_rast_re)] <- 0

# resample the DBBMM 
mig_DBBMM_rast <- rast(mig_r_reprojected)
mig_DBBMM_rast_re <- terra::resample(mig_DBBMM_rast, mig_comp_rast_grid, method = "bilinear")
mig_DBBMM_rast_re[is.na(mig_DBBMM_rast_re)] <- 0



###### Boot overlap ######

mig_UD_list <- list(mig_KUD_rast_re, mig_LDE_rast_re, mig_DBBMM_rast_re)



# BA
# 50%

mig_BA_Stats_50 <- boot_overlap(
  mig_UD_list,
  R = 20,
  method = "BA",
  percent = 50,
  quantiles = c(0.25, 0.5, 0.75)
)
print(mig_BA_Stats_50)


# 95%

mig_BA_Stats_95 <- boot_overlap(
  mig_UD_list,
  R = 20,
  method = "BA",
  percent = 95,
  quantiles = c(0.25, 0.5, 0.75)
)

print(mig_BA_Stats_95)





















