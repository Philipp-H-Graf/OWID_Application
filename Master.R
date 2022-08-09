library(ncdf4)
library(ncdf4.helpers)
library(PCICt)

library(raster)
library(rgdal)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)

library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)

rm(list = ls());gc()

# Read Data and preparation to extract the data for January to February to 2100

projected_filepath <- paste0("./Data/evspsbl_Amon_CNRM-CM6-1-HR_ssp585_r1i1p1f2_gr_201501-210012_v20191202.nc")
projected_output <- nc_open(projected_filepath)

lon <- ncvar_get(projected_output, "lon")
lat <- ncvar_get(projected_output, "lat", verbose = F)
time <- ncvar_get(projected_output, "time")

projected.array <- ncvar_get(projected_output, "evspsbl")

finalprojected_tmp <- data.frame(matrix(ncol = 4, nrow = 0))
names(finalprojected_tmp) <- c("x", "y", "layer", "year")

# Extraction of the data of January to December to 2100

for (i in seq(1021, 1032, 1)) {
  
  #i = 1021
  
  print(i)
  tmp <- projected.array[, , i]
  tmp <- raster(t(tmp), xmn = min(lon), xmx = max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  tmp <- flip(tmp, direction='y')
  tmp <- rotate(tmp)
  tmp <- as.data.frame(rasterToPoints(tmp))
  tmp$year <- 2100
  finalprojected_tmp <- rbind(finalprojected_tmp, tmp)
  
}

# Calculate the mean over the twelve months and prepare dataset for assigning them to the corresponding countries

finalprojected_tmp_final <- finalprojected_tmp %>%
  mutate(x = as.character(x)) %>%
  mutate(y = as.character(y))
finalprojected_tmp_final$coords <- paste0(finalprojected_tmp_final$x, ",", finalprojected_tmp_final$y)
finalprojected_tmp_final <- finalprojected_tmp_final %>% 
  group_by(coords) %>% 
  mutate(mean = mean(layer)) %>% 
  ungroup()
finalprojected_tmp_final <- finalprojected_tmp_final %>% 
  dplyr::select(year, x, y, mean) %>% 
  unique()

rm(list = setdiff(ls(), "finalprojected_tmp_final"))

sf::sf_use_s2(FALSE)

finalprojected_tmp_final <- sf::st_as_sf(finalprojected_tmp_final, 
                                         coords = c("x", "y"),
                                         crs = 4326)

# Country Data: Transformation from Multipolygon to Polygon

map <- ne_countries(scale = "medium", type = "map_units", returnclass = "sf")
map <- sf::st_cast(map, "POLYGON")

final <- data.frame(Geounit = NA, Avg_Evaporation = NA); final <- final[-1,]

# merge Country- to Evaporation-Data and preparation for the visualization of
# Evaporation-data in Austria

for (country in unique(map$admin)){
  
  #country = "Austria"
  
  print(country)
  
  if(country %in% c("Antarctica", "Fiji")) {next}
  
  else {
  
  tmp <- map[map$admin == country,]
  
  tmp_merged <- spatialEco::point.in.poly(finalprojected_tmp_final, tmp)
  tmp_merged <- as.data.frame(tmp_merged)
  tmp_merged <- tmp_merged %>% 
    filter(admin == country) %>% 
    dplyr::select(admin, mean, coords.x1, coords.x2)
  tmp_plot <- sf::st_as_sf(tmp_merged, coords = c("coords.x1", "coords.x2"), crs = 4326)
  
  tmp_plot$geometry <- NULL
  
  tmp_plot <- tmp_plot%>%
    dplyr:: select(Geounit = admin, Avg_Evaporation = mean) %>% 
    mutate(Avg_Evaporation = mean(Avg_Evaporation)) %>% 
    unique()
  
  final <- rbind(final, tmp_plot)
  
  print(paste0(country, ": ", nrow(final)))
  
  }
  
}

# Visualization of Evaporation-data in Austria

ggplot() +
  geom_sf(data = map %>% filter(admin == "Austria"), fill = NA) +
  geom_sf(data = tmp_plot, mapping = aes(geometry = geometry, colour = mean), show.legend = "polygon") +
  scale_colour_viridis(option = "turbo", breaks = c(min(tmp_plot$mean), max(tmp_plot$mean))) +
  labs(title = "Visual example of data points on evaporation",
       subtitle = "Austria in 2100",
       caption = "Source: Copernicus") +
  theme(axis.text = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, colour = "darkgray"))

# Saving of data to visualize Austria and final two plots in the R Markdown

save(tmp_plot, file = "Austria.rda")
save(final, file = "~/Desktop/OWID_FILE.rda")
