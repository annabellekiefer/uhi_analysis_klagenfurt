# =====================================================
#1. LIBRARY IMPORT
# =====================================================
library(terra)
library(raster)
library(sf)
library(tidyverse)
library(sp)
library(spdep)
library(FNN)
library(RColorBrewer)
library(ggplot2)
library(rasterVis)   
library(dplyr)

# =====================================================
# 2. PREPARATION & DATA LOADING 
# =====================================================
#Loading bounding box for Klagenfurt
klagenfurt_bbox <- st_read("BB_small.shp")

#Load relevant bands (B3=Green, B4=Red, B5=NIR, B6=SWIR, B10=LST)
bands <- c(
  "Landsat09/LC09_L2SP_190028_20240811_20240812_02_T1_SR_B3.TIF",
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
raster_stack <- terra::mask(raster_stack, valid_mask, maskvalue = FALSE)

# Assign band names for clarity (important step before plotting)
names(raster_stack) <- c("Green", "Red", "NIR", "SWIR", "Thermal")

# Convert LST band to Kelvin (Apply Scaling Factor)
# Formula: ST = (DN * 0.00341802) + 149.0
raster_stack$Thermal <- (raster_stack$Thermal * 0.00341802) + 149.0

# Convert LST band from Kelvin to Celsius
raster_stack$Thermal <- raster_stack$Thermal - 273.15

#Applying Scaling Factor to Surface Reflectance data
#Formula: SR = (DN * 0.0000275) - 0.2
raster_stack$Green<-(raster_stack$Green*0.0000275) - 0.2
raster_stack$Red<-(raster_stack$Red*0.0000275) - 0.2
raster_stack$NIR<-(raster_stack$NIR*0.0000275) - 0.2
raster_stack$SWIR<-(raster_stack$SWIR*0.0000275) - 0.2

# Clip raster stack using the Klagenfurt bounding box
raster_cropped_small <- crop(raster_stack, klagenfurt_bbox)

# Plot each band individually
par(mfrow = c(2, 2))  # Set layout to 2x2 for better overview

plot(raster_cropped_small$Green, main = "Green Band (B3)", col = terrain.colors(10))
plot(raster_cropped_small$Red, main = "Red Band (B4)", col = terrain.colors(10))
plot(raster_cropped_small$NIR, main = "NIR Band (B5)", col = terrain.colors(10))
plot(raster_cropped_small$SWIR, main = "SWIR Band (B6)", col = terrain.colors(10))
plot(raster_cropped_small$Thermal, main = "Thermal Band (B10)", col = terrain.colors(10))

# Reset plotting layout
par(mfrow = c(1, 1))

#Save the clipped raster stack
writeRaster(raster_cropped_small, "Klagenfurt_Landsat_Cropped.tif", overwrite = TRUE)

# =====================================================
# 3. CALCULATE VEGETATION INDICES
# =====================================================
# 3.1 NDVI
# =====================================================
# NDVI Formula: (NIR - Red) / (NIR + Red)
ndvi <- (raster_cropped_small$NIR - raster_cropped_small$Red)/(raster_cropped_small$NIR + raster_cropped_small$Red)
#ndvi[ndvi < -1 | ndvi > 1] <- NA
#ndvi[is.infinite(ndvi)] <- NA
#vals <- values(ndvi, na.rm = TRUE)
# Plot NDVI
plot(ndvi, main = "NDVI - Klagenfurt")

# Save NDVI
writeRaster(ndvi, "NDVI_Klagenfurt.tif", overwrite = TRUE)

# =====================================================
# 3.2 SAVI
# =====================================================
#SAVI Formula: ((NIR-Red)(1+l))/(NIR+Red+l)
#l soil-adjusted correction factor, frequently accepted as 0.5
savi<- ((raster_cropped_small$NIR - raster_cropped_small$Red)*(1+0.5))/(raster_cropped_small$NIR + raster_cropped_small$Red+0.5)

# Plot SAVI
plot(savi, main = "SAVI - Klagenfurt")

# Save SAVI
writeRaster(savi, "SAVI_Klagenfurt.tif", overwrite = TRUE)
# =====================================================
# 4. CALCULATE URBAN INDICES
# =====================================================
# 4.1 NDBI
# =====================================================
# NDBI Formula: (SWIR - NIR) / (SWIR + NIR)
ndbi <- (raster_cropped_small$SWIR - raster_cropped_small$NIR) / (raster_cropped_small$SWIR + raster_cropped_small$NIR)
#ndbi[ndbi < -1 | ndbi > 1] <- NA
# Plot NDBI
plot(ndbi, main = "NDBI - Klagenfurt")

# Save NDBI
writeRaster(ndbi, "NDBI_Klagenfurt.tif", overwrite = TRUE)

# =====================================================
# 4.2 IBI
# =====================================================
#IBI Formula: ((2*SWIR1/(SWIR1+NIR))-((NIR/(NIR+Red))+(Green/(Green+SWIR1))))/((2*SWIR1/(SWIR1+NIR))+((NIR/(NIR+Red))+(Green/(Green+SWIR1))))
ibi<- ((2*raster_cropped_small$SWIR/(raster_cropped_small$SWIR+raster_cropped_small$NIR))-((raster_cropped_small$NIR/(raster_cropped_small$NIR+raster_cropped_small$Red))+(raster_cropped_small$Green/(raster_cropped_small$Green+raster_cropped_small$SWIR))))/((2*raster_cropped_small$SWIR/(raster_cropped_small$SWIR+raster_cropped_small$NIR))+((raster_cropped_small$NIR/(raster_cropped_small$NIR+raster_cropped_small$Red))+(raster_cropped_small$Green/(raster_cropped_small$Green+raster_cropped_small$SWIR))))
#ibi[ibi < -1 | ibi > 1] <- NA
# Plot IBI
plot(ibi, main = "IBI - Klagenfurt")

# Save IBI
writeRaster(ibi, "IBI_Klagenfurt.tif", overwrite = TRUE)

# =====================================================
# 5. CALCULATE Modified Normalized Difference Water Index (MNDWI)
# =====================================================
#MNDWI Formula: (Green-SWIR)/(Green+SWIR)

#mndwi<-(raster_cropped_small$Green-raster_cropped_small$SWIR)/(raster_cropped_small$Green+raster_cropped_small$SWIR)
#mndwi[mndwi < -1 | mndwi > 1] <- NA
# Plot MNDWI
#plot(mndwi, main = "MNDWI - Klagenfurt")

# Save MNDWI
#writeRaster(mndwi, "MNDWI_Klagenfurt.tif", overwrite = TRUE)

# =====================================================
# 6. REGRESSION ANALYSIS BETWEEN LST AND INDICES
# =====================================================

# Extract LST values
lst_values <- values(raster_cropped_small$Thermal)

# Extract NDVI and NDBI values
ndvi_values <- values(ndvi)
savi_values <- values (savi)
ndbi_values <- values(ndbi)
ibi_values <- values(ibi)

# Create a dataframe with LST, NDVI, and NDBI
df <- data.frame(
  LST = lst_values,
  NDVI = ndvi_values,
  SAVI = savi_values,
  NDBI = ndbi_values,
  IBI= ibi_values
)
names(df) <- c("LST", "NDVI", "SAVI", "NDBI", "IBI")
# Remove NA values to avoid errors in regression
df <- na.omit(df)

# Simple Linear Regression: NDVI vs. LST
ndvi_lm <- lm(LST ~ NDVI, data = df)
summary(ndvi_lm)

# Simple Linear Regression: SAVI vs. LST
savi_lm <- lm(LST ~ SAVI, data = df)
summary(savi_lm)

# Simple Linear Regression: NDBI vs. LST
ndbi_lm <- lm(LST ~ NDBI, data = df)
summary(ndbi_lm)

# Simple Linear Regression: IBI vs. LST
ibi_lm <- lm(LST ~ IBI, data = df)
summary(ibi_lm)

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

# Plot SAVI vs. LST
plot(df$SAVI[sample_indices], df$LST[sample_indices], 
     main = "Regression: SAVI vs. LST", 
     xlab = "SAVI", 
     ylab = "LST (°C)",
     pch = 20, col = "#b7cdb7")
abline(savi_lm, col = "#4f674f", lwd = 2)

# Plot NDBI vs. LST
plot(df$NDBI[sample_indices], df$LST[sample_indices], 
     main = "Regression: NDBI vs. LST", 
     xlab = "NDBI", 
     ylab = "LST (°C)",
     pch = 20, col = "#b7cdb7")
abline(ndbi_lm, col = "#4f674f", lwd = 2)

# Plot IBI vs. LST
plot(df$IBI[sample_indices], df$LST[sample_indices], 
     main = "Regression: IBI vs. LST", 
     xlab = "IBI", 
     ylab = "LST (°C)",
     pch = 20, col = "#b7cdb7")
abline(ibi_lm, col = "#4f674f", lwd = 2)

# =====================================================
# 7. LULC EXTRACTION BASED ON SPECTRAL INDICES
# =====================================================
thr_ndvi <- 0.4
thr_savi <- 0.3
thr_ndbi <- -0.25
thr_ibi <- -0.2

# Creating masks
veg_mask <- ndvi > thr_ndvi & savi > thr_savi

urban_mask <- ndvi < thr_ndvi &
  savi < thr_savi &
  ndbi > thr_ndbi &
  ibi > thr_ibi 
lulc <- ndvi
values(lulc) <- 0  

lulc[veg_mask]   <- 1 
lulc[urban_mask] <- 2 

# Classification
# Creating data frame
lulc_df <- as.data.frame(lulc, xy = TRUE)
colnames(lulc_df) <- c("x", "y", "class")

lulc_df$class <- factor(lulc_df$class,
                        levels = c(1, 2, 0),
                        labels = c("Vegetation", "Urban", "Other"))

class_colors <- c("forestgreen", "purple", "lightgrey")

# Plot with ggplot2
ggplot(lulc_df, aes(x = x, y = y, fill = class)) +
  geom_raster() +
  scale_fill_manual(values = class_colors, name = "LULC Class") +
  coord_equal() +
  theme_minimal() +
  labs(title = "Land Use Classification") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_blank()
  )

# =====================================================
# 8. PLOTTING LST
# =====================================================
# Extract the thermal band
thermal <- raster_cropped_small$Thermal

# Plot LST
plot(thermal, main = "LST (°C) - Klagenfurt")

# Save LST in Celsius
writeRaster(thermal, "LST_Celsius_Klagenfurt.tif", overwrite = TRUE)

# ========================
# 9. UHI ANALYSIS
# ========================
# 9.1 OVERALL UHI ANALYSIS
# ========================

# Extract LST data
lst_vals <- values(raster_cropped_small$Thermal)
mean_LST <- mean(lst_vals, na.rm = TRUE)
std_LST  <- sd(lst_vals, na.rm = TRUE)
UHI_thresh <- mean_LST + 0.5 * std_LST

# Mask: 1 = UHI, 0 = kein UHI
uhi_mask <- raster_cropped_small$Thermal > UHI_thresh
uhi_raster <- classify(uhi_mask, matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE))

# DataFrame for ggplot
uhi_df <- as.data.frame(uhi_raster, xy = TRUE)
colnames(uhi_df) <- c("x", "y", "uhi_flag")
uhi_df$uhi_flag <- factor(uhi_df$uhi_flag, levels = c(0, 1), labels = c("Non-UHI", "UHI"))

# Plot
ggplot(uhi_df, aes(x = x, y = y, fill = uhi_flag)) +
  geom_raster() +
  scale_fill_manual(
    values = c("lightblue", "red"),
    name = "UHI Status"
  ) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Urban Heat Island Analysis") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# ========================
# 9.2 UHI IN URBAN AREAS
# ========================
# Mask urban areas
urban_mask <- classify(lulc, matrix(c(0, NA, 1, NA, 2, 1), ncol = 2, byrow = TRUE))
thermal_urban <- terra::mask(raster_cropped_small$Thermal, urban_mask)

# Extract LST sata
urban_vals <- values(thermal_urban)
urban_mean <- mean(urban_vals, na.rm = TRUE)
urban_sd   <- sd(urban_vals, na.rm = TRUE)
urban_thresh <- urban_mean + 0.5 * urban_sd

# Mask
urban_uhi <- thermal_urban > urban_thresh
urban_uhi_raster <- classify(urban_uhi, matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE))

# Dataframe for ggplot
urban_uhi_df <- as.data.frame(urban_uhi_raster, xy = TRUE)
colnames(urban_uhi_df) <- c("x", "y", "uhi_flag")

urban_uhi_df$uhi_flag <- factor(
  urban_uhi_df$uhi_flag,
  levels = c(0, 1),
  labels = c("Non-UHI", "UHI")
)

# Plot
ggplot(urban_uhi_df, aes(x = x, y = y, fill = uhi_flag)) +
  geom_raster() +
  scale_fill_manual(
    values = c("lightblue", "red"),
    name = "UHI Status"
  ) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Urban Heat Islands Analysis for urban areas") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# ========================
# 9.3 LST FOR VEGETATION
# ========================
# Mask vegetation areas
veg_mask <- classify(lulc, matrix(c(0, NA, 1, 1, 2, NA), ncol = 2, byrow = TRUE))

thermal_veg <- terra::mask(raster_cropped_small$Thermal, veg_mask)

# Extract LST data
veg_vals <- values(thermal_veg)
veg_mean <- mean(veg_vals, na.rm = TRUE)
veg_sd   <- sd(veg_vals, na.rm = TRUE)

# ========================
# 10 IDENTIFYING SUHI
# ========================
#Formula: (LST-LSTmean)/std

suhi_raster <- (raster_cropped_small$Thermal - mean_LST) / std_LST

# Step 2: Reclassify into 3 categories
# 0 = Non-SUHI (≤ 1.5)
# 1 = Moderate SUHI (> 1.5 & ≤ 2.0)
# 2 = Severe SUHI (> 2.0)
suhi_classes <- classify(
  suhi_raster,
  rcl = matrix(c(
    -Inf, 1.5, 0,
    1.5, 2.0, 1,
    2.0, Inf, 2
  ), ncol = 3, byrow = TRUE)
)

# Step 3: Convert to DataFrame
suhi_df <- as.data.frame(suhi_classes, xy = TRUE)
colnames(suhi_df) <- c("x", "y", "suhi_class")

suhi_df$suhi_class <- factor(suhi_df$suhi_class,
                             levels = c(0, 1, 2),
                             labels = c("Non-SUHI", "SUHI > 1.5", "SUHI > 2.0"))

# Step 4: Plot
library(ggplot2)
ggplot(suhi_df, aes(x = x, y = y, fill = suhi_class)) +
  geom_raster() +
  scale_fill_manual(
    values = c("lightblue", "red", "darkred"),
    name = "SUHI Intensity"
  ) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Standardized Surface Urban Heat Islands") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )