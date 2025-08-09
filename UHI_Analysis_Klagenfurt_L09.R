# =====================================================
#1. LIBRARY IMPORT
# =====================================================
library(terra)
library(raster)
library(sf)
library(tidyverse)

# =====================================================
# 2. PREPARATION & DATA LOADING 
# =====================================================
#Defining Bounding Box for Klagenfurt
klagenfurt_bbox <- st_bbox(c(xmin = 438965.046613823, ymin = 5157669.01115719, xmax = 455197.537933982, ymax = 5173262.20815434), crs = 32633)
klagenfurt_sf <- st_as_sfc(klagenfurt_bbox)

#Load relevant bands (B4=Red, B5=NIR, B6=SWIR, B10=LST)
bands <- c(
  "Landsat09/LC09_L2SP_190028_20240811_20240812_02_T1_SR_B4.TIF",  
  "Landsat09/LC09_L2SP_190028_20240811_20240812_02_T1_SR_B5.TIF",  
  "Landsat09/LC09_L2SP_190028_20240811_20240812_02_T1_SR_B6.TIF",  
  "Landsat09/LC09_L2SP_190028_20240811_20240812_02_T1_ST_B10.TIF"  
)

# Load bands as a raster stack
raster_stack <- rast(bands)

#Load Quality Assessment Band
qa <- rast("Landsat09/LC09_L2SP_190028_20240811_20240812_02_T1_QA_PIXEL.TIF")

#Masking raster stack with Quality Assessment Band
# Function for extracting one bit per cell 
extract_bit <- function(x, bit) {
  bitwAnd(x, 2^bit) != 0
}

# Combines several bits into a “good pixel” mask
valid_mask <- app(qa, fun = function(x) {
  extract_bit(x, 6) &                 # Clear
    !extract_bit(x, 0) &              # Not Fill
    !extract_bit(x, 1) &              # Not Dilated Cloud
    !extract_bit(x, 2) &              # Not Cirrus
    !extract_bit(x, 3) &              # Not Cloud
    !extract_bit(x, 4) &              # Not Cloud Shadow
    !extract_bit(x, 5) &              # Not Snow
    !extract_bit(x, 7)                # Not Water
})
# Apply mask
raster_stack <- mask(raster_stack, valid_mask, maskvalue = FALSE)

# Assign band names for clarity (important step before plotting)
names(raster_stack) <- c("Red", "NIR", "SWIR", "Thermal")

# Convert LST band to Kelvin (Apply Scaling Factor)
# Formula: ST = (DN * 0.00341802) + 149.0
raster_stack$Thermal <- (raster_stack$Thermal * 0.00341802) + 149.0

# Convert LST band from Kelvin to Celsius
raster_stack$Thermal <- raster_stack$Thermal - 273.15

#Applying Scaling Factor to Surface Reflectance data
#Formula: SR = (DN * 0.0000275) - 0.2
raster_stack$Red<-(raster_stack$Red*0.0000275) - 0.2
raster_stack$NIR<-(raster_stack$NIR*0.0000275) - 0.2
raster_stack$SWIR<-(raster_stack$SWIR*0.0000275) - 0.2

# Clip raster stack using the Klagenfurt bounding box
raster_cropped <- crop(raster_stack, klagenfurt_sf)

# Plot each band individually
par(mfrow = c(2, 2))  # Set layout to 2x2 for better overview

plot(raster_cropped$Red, main = "Red Band (B4)", col = terrain.colors(10))
plot(raster_cropped$NIR, main = "NIR Band (B5)", col = terrain.colors(10))
plot(raster_cropped$SWIR, main = "SWIR Band (B6)", col = terrain.colors(10))
plot(raster_cropped$Thermal, main = "Thermal Band (B10)", col = terrain.colors(10))

# Reset plotting layout
par(mfrow = c(1, 1))

#Save the clipped raster stack
writeRaster(raster_cropped, "Klagenfurt_Landsat_Cropped.tif", overwrite = TRUE)

# =====================================================
# 3. CALCULATE VEGETATION INDEX (NDVI)
# =====================================================
# NDVI Formula: (NIR - Red) / (NIR + Red)
ndvi <- (raster_cropped$NIR - raster_cropped$Red)/(raster_cropped$NIR + raster_cropped$Red)
#ndvi[ndvi < -1 | ndvi > 1] <- NA
# Plot NDVI
plot(ndvi, main = "NDVI - Klagenfurt")

# Save NDVI
writeRaster(ndvi, "NDVI_Klagenfurt.tif", overwrite = TRUE)

# =====================================================
# 4. CALCULATE URBAN INDEX (NDBI)
# =====================================================
# NDBI Formula: (SWIR - NIR) / (SWIR + NIR)
ndbi <- (raster_cropped$SWIR - raster_cropped$NIR) / (raster_cropped$SWIR + raster_cropped$NIR)
#ndbi[ndbi < -1 | ndbi > 1] <- NA
# Plot NDBI
plot(ndbi, main = "NDBI - Klagenfurt")

# Save NDBI
writeRaster(ndbi, "NDBI_Klagenfurt.tif", overwrite = TRUE)

# =====================================================
# 5. PLOTTING LST
# =====================================================

# Extract the thermal band
thermal <- raster_cropped$Thermal

# Plot LST
plot(thermal, main = "LST (°C) - Klagenfurt")

# Save LST in Celsius
writeRaster(lst_celsius, "LST_Celsius_Klagenfurt.tif", overwrite = TRUE)

# =====================================================
# 6. APPLYING THRESHOLDS FOR LULC CLASSIFICATION
# =====================================================
# Creating masks for each category
# Vegetation
veg_mask <- (ndvi > 0.2) & (ndbi < 0)
veg_mask[veg_mask == 0] <- NA  # Set non-category values to NA

plot(veg_mask, main = "Vegetation", col = "#4f674f", legend = FALSE)

# Water Bodies
water_mask <- (ndvi < 0) & (ndbi < 0)
water_mask[water_mask == 0] <- NA

plot(water_mask, main = "Water Bodies", col = "blue", legend = FALSE)

# Built-up & Bare Land
built_mask <- (ndvi < 0.2) & (ndbi >= 0)
built_mask[built_mask == 0] <- NA

plot(built_mask, main = "Built-up & Bare Land", col = "#b7c", legend = FALSE)

#All plots together
# Plot Vegetation Mask
plot(veg_mask, main = "Land Use Classification", col = "#4f674f", legend = FALSE, axes = FALSE, box = FALSE)

# Add Water Mask to the plot
plot(water_mask, add = TRUE, col = "blue", legend = FALSE)

# Add Built-up & Bare Land Mask to the plot
plot(built_mask, add = TRUE, col = "#b7c", legend = FALSE)

# =====================================================
# 7. IDENTIFYING URBAN HEAT ISLANDS
# =====================================================
# Calculate mean and standard deviation of LST in Vienna
mean_LST <- mean(values(raster_cropped$Thermal), na.rm = TRUE)
std_LST <- sd(values(raster_cropped$Thermal), na.rm = TRUE)

# Apply the UHI formula:
# UHI: LST > mean + 0.5 * std

# Create a new raster for UHI
UHI_raster <- raster_cropped$Thermal

# Apply the condition for UHI (LST > mean + 0.5 * std)
UHI_raster[UHI_raster > (mean_LST + 0.5 * std_LST)] <- 1  # UHI marked as 1

# Apply the condition for Non-UHI (0 < LST <= mean + 0.5 * std)
UHI_raster[UHI_raster != 1] <- 0 # Non-UHI marked as 0

UHI_values<-values(UHI_raster)

# Plot the UHI map
plot(UHI_raster, main = "Urban Heat Island (UHI) Map for Klagenfurt",
     col = c("#4f674f", "#b7cdb7"), legend = FALSE)

# =====================================================
# 7. IDENTIFYING URBAN HOT SPOTS
# =====================================================
# Apply the UHI formula:
# UHS: LST > mean + 0.2 * std

# Create a new raster for UHS
UHS_raster <- raster_cropped$Thermal

# Apply the condition for UHI (LST > mean + 0.5 * std)
UHS_raster[UHS_raster > (mean_LST + 2 * std_LST)] <- 1  # UHI marked as 1

# Apply the condition for Non-UHI (0 < LST <= mean + 0.5 * std)
UHS_raster[UHS_raster != 1] <- 0 # Non-UHI marked as 0

UHS_values<-values(UHS_raster)

# Plot the UHI map
plot(UHI_raster, main = "Urban Hot Spots (UHS) Map for Klagenfurt",
     col = c("#4f674f", "#b7cdb7"))

# =====================================================
# 8. REGRESSION ANALYSIS BETWEEN LST AND INDICES
# =====================================================

# Extract LST values
lst_values <- values(raster_cropped$Thermal)

# Extract NDVI and NDBI values
ndvi_values <- values(ndvi)
ndbi_values <- values(ndbi)

# Create a dataframe with LST, NDVI, and NDBI
df <- data.frame(
  LST = lst_values,
  NDVI = ndvi_values,
  NDBI = ndbi_values
)
names(df) <- c("LST", "NDVI", "NDBI")
# Remove NA values to avoid errors in regression
df <- na.omit(df)

# Simple Linear Regression: NDVI vs. LST
ndvi_lm <- lm(LST ~ NDVI, data = df)
summary(ndvi_lm)

# Simple Linear Regression: NDBI vs. LST
ndbi_lm <- lm(LST ~ NDBI, data = df)
summary(ndbi_lm)

# Sample 5000 random points of the data for faster visualization 
set.seed(42)  
num_points <- 5000
sample_indices <- sample(nrow(df), size = num_points)

# Plot sampled data
plot(df$NDVI[sample_indices], df$LST[sample_indices], 
     main = "Regression: NDVI vs. LST ", 
     xlab = "NDVI", 
     ylab = "LST (°C)",
     pch = 20, col = "#b7cdb7")
abline(ndvi_lm, col = "#4f674f", lwd = 2)

# Plot NDBI vs. LST
plot(df$NDBI[sample_indices], df$LST[sample_indices], 
     main = "Regression: NDBI vs. LST", 
     xlab = "NDBI", 
     ylab = "LST (°C)",
     pch = 20, col = "#b7cdb7")
abline(ndbi_lm, col = "#4f674f", lwd = 2)


