##### Setup #####
 

library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(MASS)
library(M3)
library(aniMotum)
library(spatstat)
library(adehabitatHR) #KUD
library(sp)
library(raster)
library(marmap) #Bathymetry data
library(terra)
library(ggnewscale)
library(move) #DBBMM
library(rasterVis) #Rasters for the DBBMM
library(gridExtra)
library(SDLfilter) #overlap metrics
library(concaveman) #concaveman small area polygon
library(zoo) #For rollapply



# Set your own working directory
setwd("~/GitHub/honours_corridors/Data")

# Set up the map for your graphs
world <- ne_countries(scale = "medium", returnclass = "sf")

## Preset the GLS Error 
# in metres
GLSerror <- 200000

# Function for the centre of the  laea

laea_centre <- function(i) {
  range <- range(i)
  dif <- diff(range, na.rm = TRUE)
  mid <- min(i) + (dif/2)
  return(mid)
}



##### Great Shearwater GLS #####

# This data has been preprocessed and so cannot fit the ssm model to it


GSGLS <- read.csv("GreatShearwaterGLS.csv")
GSGLS$POSIX <- as.POSIXct(GSGLS$POSIX, format = "%d/%m/%Y %H:%M", tz = "UTC")

GSGLS_model <- GSGLS %>% 
  dplyr::select(track_id, POSIX, longitude, latitude, mig.stage2) %>% 
  group_by(track_id) %>% 
  # Filters ret individuals with low number of points
  filter(n()>=10)


colnames(GSGLS_model) <- c("id", "POSIX", "lon", "lat", "mig.stage")

# mig.stage is split into 5:
#  0premig
#  1outward
#  2winter
#  3return
#  4postmig









###### All Data #####

GSGLS_all <- GSGLS_model %>% 
  group_by(id) 


# Orders for the DBBMM
GSGLS_all <- GSGLS_all[order(GSGLS_all$id, GSGLS_all$POSIX), ] # Orders it by id

GSGLS_all <- as.data.frame(GSGLS_all)





df_sf_all <- GSGLS_all %>% 
  st_as_sf(coords = c('lon', 'lat')) %>% 
  st_set_crs(4326)


ggplot() +
  theme_bw() +
  geom_sf(aes(), data = world, colour = "grey") +
  geom_sf(aes(colour = id), data = df_sf_all, show.legend = "point") +
  coord_sf(xlim = c(-80, 40), ylim = c(-80, 80)) +
  labs(title = "GLS Location Data",
       x = "Longitude",
       y = "Latitude") +
  theme(legend.position = "none")

# Blank map
blank_plot <- ggplot(data = world) +
  geom_sf(fill = "darkgrey", color = "black", size = 0.2) +
  coord_sf(xlim = c(-140, 60), ylim = c(-80, 80)) +
  theme_void()
plot(blank_plot)
ggsave("blank_world_map.png", plot = blank_plot, width = 10, height = 6, dpi = 300, bg = "transparent")







##### Outward Migration Filtering #####


##### New Outward Migration Filtering #####



GSGLS_out <- GSGLS_model %>% 
  filter(mig.stage == "1outward") %>% 
  group_by(id) 



# Filters the by finding start at -20 lat and end at 40 lat
# then adds 5 days of occurences to either side
migration_filter <- function(df, start_lat, end_lat, x, y) {
  
  df <- df %>% arrange(id, POSIX)
  
  df %>% group_by(id) %>% group_modify(~{
    
    dat <- .x %>% arrange(POSIX)
    dat <- dat %>% mutate(date = as.Date(POSIX))
    
    # Find start at start_lat
    start_idx <- which(dat$lat >= start_lat)[1]
    if (is.na(start_idx)) return(dat[0, ])
    
    start_date <- dat$date[start_idx]
    start_date_expanded <- start_date - x # Extends by x days
    
    # Find end at end_lat
    end_idx <- which(dat$lat >= end_lat &
                       row_number(dat) > start_idx)[1]
    
    if (is.na(end_idx)) {
      end_date_expanded <- max(dat$date)
    } else {
      end_date <- dat$date[end_idx]
      end_date_expanded <- end_date + y #Extends by y days
    }
    
    # Filters data from above
    dat %>%
      filter(date >= start_date_expanded &
               date <= end_date_expanded)
  })
}



GSGLS_out <- migration_filter(
  GSGLS_out,       # Data
  start_lat = -20, # Starting latitude
  end_lat   = 40,  # Ending latitude
  x = 5,   # include x points before start
  y = 5    # include y points after end
)

nrow(GSGLS_out)

GSGLS_out_map <- GSGLS_out %>% 
   st_as_sf(coords = c('lon', 'lat')) %>% 
  st_set_crs(4326)
 

ggplot() +
   theme_bw() +
   geom_sf(aes(), data = world) +
   geom_sf(aes(colour = POSIX), data = GSGLS_out_map, show.legend = "point") +
   coord_sf(xlim = c(-80, 10), ylim = c(-60, 60)) +
   labs(title = "Outward Migration",
        x = "Longitude",
        y = "Latitude") 




###### SSM for GLS out ######


# Format the data for the aniMotum package
# AniMotum needs this column for some reason to work. So just put it in
GSGLS_out$lc <- "GL"
GSGLS_out$se.x <- as.numeric(200000)
GSGLS_out$se.y <- as.numeric(200000)
GSGLS_out <- GSGLS_out %>% 
  dplyr::select(-mig.stage, -date)



AniMotum_out <- format_data(GSGLS_out, id = "id", date = "POSIX", 
                            lc = "lc", coord = c("lon", "lat"), sderr = "se.x")


##### Generate the SSM and MP #####

## Same as PTT data

# Creates an SSM - random walk with time step of 12hrs
Fit_out <- fit_ssm(AniMotum_out, model = "rw", time.step = 12,
                   control = ssm_control(verbose = 0))


# Reroutes the SSM to be over the ocean
Reroute_out <- route_path(Fit_out, what = "predicted", 
                          map_scale = 10, dist = 5000)


# Fits a Move persist model to the rerouted function
MPM_out <- fit_mpm(Reroute_out, what = "rerouted", model = "mpm", 
                   control = mpm_control(verbose = 0))



df_locs_out <- Reroute_out |> grab("rerouted")
df_mpm_out <- MPM_out |> grab("fitted")


GSGLS_out <- merge(df_locs_out, df_mpm_out)
GSGLS_out$POSIX <- as.POSIXct(GSGLS_out$date, format = "%d/%m/%Y %H:%M", tz = "UTC")
GSGLS_out <- GSGLS_out[order(GSGLS_out$id, GSGLS_out$POSIX), ]   # Orders it by id

# Set the CRS for the laea for later


out_lat_centre <- laea_centre(GSGLS_out$lat)
out_lon_centre <- laea_centre(GSGLS_out$lon)

out_laea_crs <- sprintf("+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m", out_lat_centre, out_lon_centre)



## outward World Map

world_out <- st_transform(world, crs = out_laea_crs)


# Set the extent for the ggplots

GSGLS_outext <- SpatialPoints(coords = cbind(GSGLS_out$lon, GSGLS_out$lat), proj4string = CRS("WGS84")) 
GSGLS_outext <- spTransform(GSGLS_outext, out_laea_crs)


out_xmin <- min(GSGLS_outext@coords[,1])*1.1
out_xmax <- max(GSGLS_outext@coords[,1])*1.1
out_ymin <- min(GSGLS_outext@coords[,2])*1.05
out_ymax <- max(GSGLS_outext@coords[,2])*1.15









# Set as sf
GLS_out <- as.data.frame(GSGLS_out)

GLS_out_sf_points <- st_as_sf(GLS_out, coords = c("lon", "lat"), crs = 4326)


df_sf_out <- GSGLS_out %>% 
  st_as_sf(coords = c('lon', 'lat')) %>% 
  st_set_crs(4326)


ggplot() +
  theme_bw() +
  geom_sf(aes(), data = world, fill = "grey") +
  geom_sf(aes(colour = id), data = df_sf_out, show.legend = "point") +
  coord_sf(xlim = c(-80, 10), ylim = c(-60, 60)) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme(legend.position = "none")


#Number of unique individuals

# Number of individuals used for migration
n_distinct(GSGLS_out$id)


###### Line Density Estimation ######


# Turn the data into a format that can be exported into ARC

# Set as sf
GLS_out_sf <- st_as_sf(GSGLS_out, coords = c("lon", "lat"), crs = 4326)

# Create line strings
GLS_out_sf_lines <- GLS_out_sf %>%
  group_by(id) %>%
  summarise(do_union = FALSE) %>%  #  Avoid merging separate tracks
  st_cast("LINESTRING")  #  Convert grouped points into lines

# Reproject to aeqd to be in m
proj_out_lines <- st_transform(GLS_out_sf_lines, crs = out_laea_crs)

# Save as shapefile
st_write(obj = proj_out_lines,
         dsn = "~/GitHub/honours_corridors/Outputs/GSGLS_LDE_out.shp",
         driver = "ESRI Shapefile",
         append = FALSE)

# Load this in ARC in the Line Density function from the Spatial Analyst Toolbox

# Use a grid cell size appropriate based on this error size and the distance of the journey
## Grid cell used was 20000m
## error of 200000m for GLS


## Load the newly created LDE

out_LDE_raster <- rast("~/GitHub/honours_corridors/Outputs/GLS_Outward_mig_LDE2.tif")

# Create the Contours



# Polygon 95
out_LDE_95 <- out_LDE_raster
out_LDE_95[out_LDE_95 < 0.05] <- NA
out_LDE_95[out_LDE_95 >= 0.05] <- 1

out_LDE_95_poly <- as.polygons(out_LDE_95) %>% 
  st_as_sf(crs = out_laea_crs) %>%
  st_buffer(dist = res(out_LDE_95)[1]/4) %>%
  st_simplify(dTolerance = res(out_LDE_95)[1]) 


# Polygon 50
out_LDE_50 <- out_LDE_raster
out_LDE_50[out_LDE_50 < 0.50] <- NA
out_LDE_50[out_LDE_50 >= 0.50] <- 1

out_LDE_50_poly <- terra::as.polygons(out_LDE_50) %>% 
  sf::st_as_sf(crs = out_laea_crs) %>%
  sf::st_buffer(dist = res(out_LDE_50)[1]/4) %>%
  sf::st_simplify(dTolerance = res(out_LDE_50)[1])



# Plot as a ggplot
out_mig_LDE <- ggplot() +
  geom_sf(aes(), data = world_out, fill = "grey") +
  geom_sf(data = out_LDE_95_poly, aes(fill = "Full Corridor", colour = "Full Corridor"), alpha = 0.5) +
  geom_sf(data = out_LDE_50_poly, aes(fill = "Core Corridor", colour = "Core Corridor"), alpha = 0.5) +
  geom_sf(aes(colour = "Locations"), data = df_sf_out, show.legend = "point") +
  coord_sf(xlim = c(out_xmin, out_xmax), ylim = c(out_ymin, out_ymax), expand = FALSE) +
  labs(x = "Longitude",
       y = "Latitude") +
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



plot(out_mig_LDE)



###### Density Raster ######

GLS_out_sf_rast <- st_transform(GLS_out_sf, out_laea_crs)

grid_size <- 200000

# Convert points to SpatVector
out_points_vect <- vect(GLS_out_sf_rast)
out_points_vect <- project(out_points_vect, out_laea_crs)


# Create an empty raster in laea
out_density_raster <- rast(ext(out_points_vect), resolution = grid_size, crs = out_laea_crs) 

# Rasterize with a density count
out_r_density <- rasterize(out_points_vect, out_density_raster, 
                       fun = "length", background = 0)



# Convert to a dataframe for plotting in ggplot
out_df_raster <- as.data.frame(out_r_density, xy = TRUE)
colnames(out_df_raster) <- c("x", "y", "Density")
out_df_raster <- filter(out_df_raster, Density > 1) %>%  # Filters out the points with 0 or 1
  filter(!is.na(Density))

# Calculate the top 10% bin values
out_top_10 <- ceiling(0.10 * (nrow(out_df_raster)))

out_df_raster$Density <- as.numeric(out_df_raster$Density)

# Order in descending order and bin
out_df_raster <- out_df_raster %>% 
  arrange(desc(Density)) %>% 
  mutate(rank = row_number()) %>% 
  mutate(bin = case_when(
    rank <= out_top_10 ~ "Core Use",
    TRUE ~ "Corridor")) %>% 
  mutate(bin_num = case_when(
    rank <= out_top_10 ~ 2,
    TRUE ~ 1))


## Removing noise

# re raster it
out_cat_raster <- rast(out_df_raster[, c("x", "y", "bin_num")], type = "xyz")

# remove the isolated points
out_raster_mask <- out_cat_raster == 1   # Just pulls the corridor points
out_raster_patches <- patches(out_raster_mask, directions = 8) #  Creates patches based on whats touching in 8 directions
out_raster_sizes <- freq(out_raster_patches) # Counts the number of cells in each patch and puts into a df
out_raster_ids_keep <- out_raster_sizes$value[out_raster_sizes$count > 2]  # Filters out any patches with number of cells less than value
out_clean_corridor <- mask(out_raster_patches, out_raster_patches %in% out_raster_ids_keep, maskvalues = FALSE) # Makes a raster that only keeps the clean corridor
out_clean_corridor[!is.na(out_clean_corridor)] <- 1

# Extract the core use values to merge
out_core_use <- out_cat_raster
out_core_use[out_core_use != 2] <- NA

# Merge the rasters together for plotting
out_clean_raster <- mosaic(out_clean_corridor, out_core_use, fun = "max")

# Turn into a df for plotting
out_clean_df <- as.data.frame(out_clean_raster, xy = TRUE)
colnames(out_clean_df) <- c("x", "y", "bin_num")
out_clean_df <- out_clean_df %>% 
  mutate(bin = factor(bin_num,
                      levels = c(1, 2),
                      labels = c("Corridor", "Core Use")))

# Plot in ggplot
GLS_out_density_raster <- ggplot() +
  geom_sf(aes(), data = world_out, fill = "grey") +
  geom_tile(data = out_clean_df, aes(x = x, y = y, fill = bin), colour = "black", alpha = 0.5) +
  geom_sf(aes(colour = "Locations"), data = df_sf_out, show.legend = "point") +
  coord_sf(xlim = c(out_xmin, out_xmax), ylim = c(out_ymin, out_ymax), expand = FALSE) +
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


plot(GLS_out_density_raster)



###### Kernel Utilisation Distribution ######

# Change the data to be spatial points
GSGLS_outKUD <- SpatialPoints(coords = cbind(GSGLS_out$lon, GSGLS_out$lat), proj4string = CRS("WGS84")) 
proj4string(GSGLS_outKUD) <- CRS("+proj=longlat +datum=WGS84")
GSGLS_outKUD <- spTransform(GSGLS_outKUD, out_laea_crs)


## Calc the grid size needed

# Extract bounding box
out_bbox <- sp::bbox(GSGLS_outKUD)

# Calc max range 
out_x_range <- out_bbox[1, "max"] - out_bbox[1, "min"]
out_y_range <- out_bbox[2, "max"] - out_bbox[2, "min"]

# Calc the max of the two
out_max_range <- max(out_x_range, out_y_range)

# Calc the grid value from this
# the number represents the grid cell size in m that you want
# Due to the kernelUD taking an integer, it will round this value to the nearest
# whole number but it will be close enough
# if the grid number ends up super low, use the ceiling function to round it up
out_grid_kud <- out_max_range / 20000


## Apply the KUD

# The default kernel (bivariate normal kernel) is used by default in the adehabitatHR package
# smoothing parameter h is the width of the kernels and can be calculated
# Sigma (σ) is the standard deviation of points
# σ can be calculated by the mean of the standard error of the x or y multiplied 
# by the square root of the number of locations


out_kud <- kernelUD(GSGLS_outKUD, h = GLSerror, grid = out_grid_kud)


# Extract the Contours you want
out_contour_50 <- getverticeshr(out_kud, percent = 50)
out_contour_95 <- getverticeshr(out_kud, percent = 95)


# Convert contours to sf objects
out_contour_50_sf <- st_as_sf(out_contour_50)
out_contour_95_sf <- st_as_sf(out_contour_95)


#Turn the kud into a dataframe that ggplot can readI
out_Kuddf <- as.data.frame(out_kud, xy = TRUE)


###Plot
GSGLS_out_kud <- ggplot() +
  theme_bw() +
  geom_sf(aes(), data = world_out, fill = "grey") +
  geom_sf(data = out_contour_95_sf, aes(fill = "Full Corridor", colour = "Full Corridor"),alpha = 0.5) +
  geom_sf(data = out_contour_50_sf, aes(fill = "Core Corridor", colour = "Core Corridor"), alpha = 0.5) + 
  geom_sf(aes(colour = "Locations"), data = df_sf_out, show.legend = "point") +
  coord_sf(xlim = c(out_xmin, out_xmax), ylim = c(out_ymin, out_ymax), expand = FALSE) +
  labs(x = "Longitude",
       y = "Latitude") +
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
plot(GSGLS_out_kud)
















###### Dynamic Brownian Bridge Movement Model ######
# Turn it into a move data for the DBBMM
GSGLS_out_BBM <- move(x = GSGLS_out$lon, y = GSGLS_out$lat,
                      time = as.POSIXct(GSGLS_out$POSIX, format="%Y-%m-%d %H:%M:%OS", tz="UTC"), 
                      proj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
                      data = GSGLS_out, animal = GSGLS_out$id)

# Reproject into a metric measurment CRS for the DBBMM
GSGLS_out_proj <- spTransform(GSGLS_out_BBM, out_laea_crs)

## Create the DBBMM

# Increasing window size increased the reliability of the brownian motion values
# but at the cost of missing rapid changes
# Increasing the margin size increases the likelihood of identifying relatively weak breakpoints
# but at cost of potentially not detecting breakpoints due to limited observations in w
# Choices should be made based on the time interval you are looking at
DBB_GSGLS_out <- brownian.bridge.dyn(GSGLS_out_proj, 
                                     ext = 1.5, # Needs to be large enough to encompass the data. Increase until does
                                     raster = 20000, # raster box size in m
                                     location.error = GLSerror, # error of location data in m
                                     window.size = 13, # Should be >2m, window of the breakpoints for behaviour
                                     margin = 5) # Holds out observations that wont be included in breakpoints from the window, 2m should be at most 1/2 of w

##Obtain the Utilisation Distributions of Each individual
gls_out_UD <- getVolumeUD(DBB_GSGLS_out)
# plot(out_UD, main = "UD")


# Stacks them all up
UD_out <- stack(gls_out_UD)
out_combined_UD <- calc(UD_out, sum)


##Standardise the raster
out_r_min <- cellStats(out_combined_UD, "min")
out_r_max <- cellStats(out_combined_UD, "max")

out_r_standardised <- (out_combined_UD - out_r_min) / (out_r_max - out_r_min)



# Define the target CRS (WGS 84 Lat/Lon)
target_crs <- "+proj=longlat +datum=WGS84 +no_defs"

# Reproject the raster
out_r_reprojected <- projectRaster(out_r_standardised, crs = out_laea_crs)

# Turn it into a data frame for ggplot
out_raster <- as.data.frame(out_r_reprojected, xy = TRUE)
colnames(out_raster) <- c("x", "y", "value")

# Set data frame as simple features and so ggplot can read it
out_raster_sf <- st_as_sf(out_raster, coords = c("x", "y"))
out_raster_sf <- subset(out_raster_sf, value <= 0.9995)


## Creating Contour layers
out_c95 <- rasterToContour(out_r_reprojected, levels = 0.95)
out_c95_sf <- st_as_sf(out_c95)
out_c95_poly <- st_cast(st_union(out_c95_sf), "POLYGON")
st_area(out_c95_poly)/10^6


out_c50 <- rasterToContour(out_r_reprojected, levels = 0.50)
out_c50_sf <- st_as_sf(out_c50)
out_c50_poly <- st_cast(st_union(out_c50_sf), "POLYGON")
st_area(out_c50_poly)/10^6

## Plot

GSGLS_out_DBBMM <- ggplot() +
  geom_sf(aes(), data = world_out, fill = "grey") +
  geom_sf(data = out_c95_poly, aes(fill = "Full Corridor", colour = "Full Corridor"), alpha = 0.5) +
  geom_sf(data = out_c50_poly, aes(fill = "Core Corridor", colour = "Core Corridor"), alpha = 0.5) +
  geom_sf(aes(colour = "Locations"), data = df_sf_out, show.legend = "point") +
  coord_sf(xlim = c(out_xmin, out_xmax), ylim = c(out_ymin, out_ymax), expand = FALSE) +
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

plot(GSGLS_out_DBBMM)








###### Concave Hull Shape ######

GLS_out_hull_points <- st_transform(GLS_out_sf_points, crs = out_laea_crs)

GLS_out_hull <- concaveman(GLS_out_hull_points, concavity = 4)

GLS_out_hull_plot <- ggplot() +
  geom_sf(aes(), data = world_out, fill = "grey") +
  geom_sf(data = GLS_out_hull, aes(colour = "Concave Hull", fill = "Concave Hull"), alpha = 0.5) +
  geom_sf(aes(colour = "Locations"), data = df_sf_out, show.legend = "point") +
  coord_sf(xlim = c(out_xmin, out_xmax), 
           ylim = c(out_ymin, out_ymax),
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

plot(GLS_out_hull_plot) 












###### Spatial Comparison ######

###### Areas and Comparisons ######

### Areas of each of the polygons

## KUD Areas
out_KUD_50_poly <- st_transform(out_contour_50_sf, crs = out_laea_crs)
out_KUD_95_poly <- st_transform(out_contour_95_sf, crs = out_laea_crs)



## Density Raster

# Turns raster binary
DR_out_density_filt_95 <- out_clean_raster
crs(DR_out_density_filt_95) <- out_laea_crs
DR_out_shape <- classify(DR_out_density_filt_95, rcl = matrix(c(0, Inf, 1), ncol = 3, byrow = TRUE)) 
# Masks and removes the 0 values
DR_out_masked <- mask(DR_out_shape, DR_out_shape, maskvalue = 0)
# Turns into polygon
DR_out_polygons <- as.polygons(DR_out_masked)
# Turns into sf
DR_out_poly_sf <- st_as_sf(DR_out_polygons)
DR_out_poly <- st_transform(DR_out_poly_sf, crs = out_laea_crs)
# Checks the area
DR_out_poly$area <- st_area(DR_out_poly)
DR_out_area <- sum(as.numeric(DR_out_poly$area), na.RM = TRUE)



##  Part density Raster
out_r_density_filt <- out_core_use
crs(out_r_density_filt) <- out_laea_crs

# Turns raster binary
DR_out_filt_shape <- classify(out_r_density_filt, rcl = matrix(c(0, Inf, 1), ncol = 3, byrow = TRUE)) 
# Masks and removes the 0 values
DR_out_filt_masked <- mask(DR_out_filt_shape, DR_out_filt_shape, maskvalue = 0)
# Turns into polygon
DR_out_filt_polygons <- as.polygons(DR_out_filt_masked, na.rm = TRUE)
# Turns into sf
DR_out_filt_poly_sf <- st_as_sf(DR_out_filt_polygons)
DR_out_filt_poly <- st_transform(DR_out_filt_poly_sf, crs = out_laea_crs)
# Checks the area
DR_out_filt_poly$area <- st_area(DR_out_filt_poly)
DR_out_filt_area <- sum(as.numeric(DR_out_filt_poly$area), na.RM = TRUE)

plot(DR_out_filt_poly)


## LDE Areas
out_LDE_50_poly
out_LDE_95_poly



## DBBMM Areas
out_c95_sf
out_c50_sf



# hull poly
out_hull_poly <- st_transform(GLS_out_hull, crs = out_laea_crs)


# Overlapping the 95% confidence interval

out_95_contours_plot <- ggplot() +
  geom_sf(aes(), data = world_out, fill = "grey") +
  geom_sf(data = GLS_out_hull, aes(fill = "Concave Hull"), colour = "orange", alpha = 0.1) +
  geom_sf(data = out_contour_95_sf, aes(fill = "Kernel Utilisation Distribution"), colour = "red", alpha = 0.1) +
  geom_sf(data = out_LDE_95_poly, aes(fill = "Line Density Estimation"), colour = "blue", alpha = 0.1) +
  geom_sf(data = out_c95_poly, aes(fill = "Dynamic Brownian Bridge Movement Model"), colour = "green", alpha = 0.1) +
  geom_sf(data = DR_out_poly, aes(fill = "Density Raster"), colour = "black", alpha = 0.1) +
  coord_sf(xlim = c(out_xmin, out_xmax), 
           ylim = c(out_ymin, out_ymax), 
           expand = FALSE) +
  scale_fill_manual(name = "Corridor Delineation Method",
                    values = c(c("Concave Hull" = "orange",
                                 "Kernel Utilisation Distribution" = "red",
                                 "Line Density Estimation" = "blue",
                                 "Dynamic Brownian Bridge Movement Model" = "green",
                                 "Density Raster" = "black"),
                               guide = guide_legend(override.aes = list(alpha = 0.1)))) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude") +
  theme(legend.position = "none")

plot(out_95_contours_plot)



## Overlapping 50%

out_50_contours_plot <- ggplot() +
  geom_sf(aes(), data = world_out, fill = "grey") +
  geom_sf(data = out_contour_50_sf, aes(fill = "Kernel Utilisation Distribution"), colour = "red", alpha = 0.1) +
  geom_sf(data = out_LDE_50_poly, aes(fill = "Line Density Estimation"), colour = "blue", alpha = 0.1) +
  geom_sf(data = out_c50_poly, aes(fill = "Dynamic Brownian Bridge Movement Model"), colour = "green", alpha = 0.1) +
  geom_sf(data = DR_out_filt_poly, aes(fill = "Density Raster"),  colour = "black", alpha = 0.1) +
  coord_sf(xlim = c(out_xmin, out_xmax), 
           ylim = c(out_ymin, out_ymax), 
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
plot(out_50_contours_plot)




##### 6 Plots together #####


grid.arrange(GSGLS_out_kud, 
             GLS_out_density_raster, 
             out_mig_LDE,
             GSGLS_out_DBBMM ,
             GLS_out_hull_plot,
             out_95_contours_plot,
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
out_95_area_polys <- list(out_KUD_95_poly, # KUD
                          DR_out_poly, # Density Raster
                          out_LDE_95_poly, # LDE
                          out_c95_poly, # DBBMM
                          out_hull_poly) # Concave Hull
out_95_area_geoms <- lapply(out_95_area_polys, st_geometry)
out_95_area_geoms_sfc <- do.call(c, out_95_area_geoms)
out_95_area_combined <- st_sf(ID = paste0("poly_", seq_along(out_95_area_geoms_sfc)), geometry = out_95_area_geoms_sfc)


# Print Info
out_basic_overlap_info <- basic_overlap_info(out_95_area_combined)
print(out_basic_overlap_info)




# 50% comparisons
# List geoms
out_50_area_polys <- list(out_KUD_50_poly,
                          DR_out_filt_poly,
                          out_LDE_50_poly, 
                          out_c50_poly)
out_50_area_geoms <- lapply(out_50_area_polys, st_geometry)
out_50_area_geoms_sfc <- do.call(c, out_50_area_geoms)
out_50_area_combined <- st_sf(ID = paste0("poly_", seq_along(out_50_area_geoms_sfc)), 
                              geometry = out_50_area_geoms_sfc)


# Print Info
out_50_overlap_info <- basic_overlap_info(out_50_area_combined)
print(out_50_overlap_info)



## KUD 95 to 50

out_KUD_area_polys <- list(out_KUD_50_poly, 
                           out_KUD_95_poly)
out_KUD_area_geoms <- lapply(out_KUD_area_polys, st_geometry)
out_KUD_area_geoms_sfc <- do.call(c, out_KUD_area_geoms)
out_KUD_area_combined <- st_sf(ID = paste0("poly_", seq_along(out_KUD_area_geoms_sfc)), 
                               geometry = out_KUD_area_geoms_sfc)


# Print Info
out_KUD_overlap_info <- basic_overlap_info(out_KUD_area_combined)
print(out_KUD_overlap_info)





## LDE 95 to 50

out_LDE_area_polys <- list(out_LDE_50_poly, 
                           out_LDE_95_poly)
out_LDE_area_geoms <- lapply(out_LDE_area_polys, st_geometry)
out_LDE_area_geoms_sfc <- do.call(c, out_LDE_area_geoms)
out_LDE_area_combined <- st_sf(ID = paste0("poly_", seq_along(out_LDE_area_geoms_sfc)), 
                               geometry = out_LDE_area_geoms_sfc)


# Print Info
out_LDE_overlap_info <- basic_overlap_info(out_LDE_area_combined)
print(out_LDE_overlap_info)





## DBBMM 95 to 50

out_DBBMM_area_polys <- list(out_c95_poly, 
                             out_c50_poly)
out_DBBMM_area_geoms <- lapply(out_DBBMM_area_polys, st_geometry)
out_DBBMM_area_geoms_sfc <- do.call(c, out_DBBMM_area_geoms)
out_DBBMM_area_combined <- st_sf(ID = paste0("poly_", seq_along(out_DBBMM_area_geoms_sfc)), 
                                 geometry = out_DBBMM_area_geoms_sfc)


# Print Info
out_DBBMM_overlap_info <- basic_overlap_info(out_DBBMM_area_combined)
print(out_DBBMM_overlap_info)





## Density Raster full to part

out_DR_area_polys <- list(DR_out_poly, 
                          DR_out_filt_poly)
out_DR_area_geoms <- lapply(out_DR_area_polys, st_geometry)
out_DR_area_geoms_sfc <- do.call(c, out_DR_area_geoms)
out_DR_area_combined <- st_sf(ID = paste0("poly_", seq_along(out_DR_area_geoms_sfc)), 
                              geometry = out_DR_area_geoms_sfc)


# Print Info
out_DR_overlap_info <- basic_overlap_info(out_DR_area_combined)
print(out_DR_overlap_info)















#    We need to turn them all into the same format to be able to compare.
#    preferably all with the same grid cell size

# Empty raster template
comp_rast_grid <- rast(ext(out_xmin, out_xmax, out_ymin, out_ymax),
                       res = 5000,
                       crs = out_laea_crs)






# resample the LDE 
out_LDE_rast_re <- terra::resample(out_LDE_raster, comp_rast_grid, method = "bilinear")
out_LDE_rast_re[is.na(out_LDE_rast_re)] <- 0

# resample the KUD
out_KUD_rast <- rast(out_kud)
out_KUD_rast_re <- terra::resample(out_KUD_rast, comp_rast_grid, metho = "bilinear")
out_KUD_rast_re[is.na(out_KUD_rast_re)] <- 0


# resample the DBBMM 
out_DBBMM_rast <- rast(out_combined_UD)
out_DBBMM_rast_re <- terra::resample(out_DBBMM_rast, comp_rast_grid, method = "bilinear")
out_DBBMM_rast_re[is.na(out_DBBMM_rast_re)] <- 0




## Overlap comparison






###### Boot overlap ######

out_UD_list <- list(out_KUD_rast_re, out_LDE_rast_re, out_DBBMM_rast_re)

# BA
# 50%

out_BA_Stats_50 <- boot_overlap(
  out_UD_list,
  R = 20,
  method = "BA",
  percent = 50,
  quantiles = c(0.25, 0.5, 0.75)
)

print(out_BA_Stats_50)


# 95%
out_BA_Stats_95 <- boot_overlap(
  out_UD_list,
  R = 20,
  method = "BA",
  percent = 95,
  quantiles = c(0.25, 0.5, 0.75)
)

print(out_BA_Stats_95)

