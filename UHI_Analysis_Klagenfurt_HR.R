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
#2. UHI DEFINITION HIGH-RESOLUTION LARGE AREA (DAY)
# =====================================================
klagenfurt_hr<-rast("32121_TH-TOP_GradC_DOM_1m_Tag.tif")
plot(klagenfurt_hr)

# ========================
#2.1 OVERALL UHI ANALYSIS
# ========================
# Extract LST data
lst_vals <- values(klagenfurt_hr)
mean_LST <- mean(lst_vals, na.rm = TRUE)
std_LST  <- sd(lst_vals, na.rm = TRUE)
UHI_thresh <- mean_LST + 0.5 * std_LST

# Mask: 1 = UHI, 0 = kein UHI
uhi_mask <- klagenfurt_hr > UHI_thresh
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
#2.2 IDENTIFYING SUHI
# ========================
#Formula: (LST-LSTmean)/std

suhi_raster <- (klagenfurt_hr - mean_LST) / std_LST

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

# =====================================================
#3. UHI DEFINITION HIGH-RESOLUTION SMALL AREA (DAY)
# =====================================================
klagenfurt_hrs<-rast("c32121_THTOP_GradC_DOM_1m_Tag_Crop.tif")
plot(klagenfurt_hrs)

# ========================
#3.1 OVERALL UHI ANALYSIS
# ========================
# Extract LST data
lst_vals <- values(klagenfurt_hrs)
mean_LST <- mean(lst_vals, na.rm = TRUE)
std_LST  <- sd(lst_vals, na.rm = TRUE)
UHI_thresh <- mean_LST + 0.5 * std_LST

# Mask: 1 = UHI, 0 = kein UHI
uhi_mask <- klagenfurt_hrs > UHI_thresh
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
#3.2 IDENTIFYING SUHI
# ========================
#Formula: (LST-LSTmean)/std

suhi_raster <- (klagenfurt_hrs - mean_LST) / std_LST

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

# =====================================================
#4. UHI DEFINITION HIGH-RESOLUTION LARGE AREA (NIGHT)
# =====================================================
klagenfurt_hrn<-rast("32121_TH-TOP_GradC_DOM_1m_Nacht.tif")
plot(klagenfurt_hrn)

# ========================
#4.1 OVERALL UHI ANALYSIS
# ========================
# Extract LST data
lst_vals <- values(klagenfurt_hrn)
mean_LST <- mean(lst_vals, na.rm = TRUE)
std_LST  <- sd(lst_vals, na.rm = TRUE)
UHI_thresh <- mean_LST + 0.5 * std_LST

# Mask: 1 = UHI, 0 = kein UHI
uhi_mask <- klagenfurt_hrn > UHI_thresh
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
#4.2 IDENTIFYING SUHI
# ========================
#Formula: (LST-LSTmean)/std

suhi_raster <- (klagenfurt_hrn - mean_LST) / std_LST

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

# =====================================================
#5. UHI DEFINITION HIGH-RESOLUTION SMALL AREA (NIGHT)
# =====================================================
klagenfurt_hrsn<-rast("c32121_THTOP_GradC_DOM_1m_Nacht_Crop.tif")
plot(klagenfurt_hrsn)

# ========================
#5.1 OVERALL UHI ANALYSIS
# ========================
# Extract LST data
lst_vals <- values(klagenfurt_hrsn)
mean_LST <- mean(lst_vals, na.rm = TRUE)
std_LST  <- sd(lst_vals, na.rm = TRUE)
UHI_thresh <- mean_LST + 0.5 * std_LST

# Mask: 1 = UHI, 0 = kein UHI
uhi_mask <- klagenfurt_hrsn > UHI_thresh
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
#5.2 IDENTIFYING SUHI
# ========================
#Formula: (LST-LSTmean)/std

suhi_raster <- (klagenfurt_hrsn - mean_LST) / std_LST

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
