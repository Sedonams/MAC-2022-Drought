#-----Loading Libraries and Data Sets----
# Load necessary libraries
library(readxl)   # For reading Excel files
library(janitor)  # For cleaning column names
library(dplyr)    # For using the pipe operator (%>%)

# Read and clean the data
dmgr <- read_xlsx("C:/Users/sedon/Downloads/MAC 2022 Mother.xlsx", sheet = "DMGR Data")

amf <- read_xlsx("C:/Users/sedon/Downloads/MAC 2022 Mother.xlsx", sheet = "AMF Score")
mac <- read_xlsx("C:/Users/sedon/Downloads/MAC 2022 Mother.xlsx", sheet = "MAC22 Weights")
emh<- read_xlsx("C:/Users/sedon/Downloads/MAC 2022 Mother.xlsx", sheet = "EMH")
nuts<- read_xlsx("C:/Users/sedon/Downloads/MAC 2022 Mother.xlsx", sheet = "MAC22 Nutrients")

g2c <- read_xlsx("C:/Users/sedon/Downloads/G2C All Data.xlsx", sheet = "Weights")
g2camf <- read_xlsx("C:/Users/sedon/Downloads/G2C All Data.xlsx", sheet = "AMF Score")

dsLegacy <- read_xlsx("C:/Users/sedon/Downloads/Legacy Study 2024.xlsx", sheet = "LP2")

dsdrought<- read_xlsx("C:/Users/sms689/Downloads/z6-19231(z6-19231)-1736291922.xlsx", sheet = "Config 1")
dswet<- read_xlsx("C:/Users/sms689/Downloads/z6-19230(z6-19230)-1736199516.xlsx", sheet = "Config 1")


# Load unnecessary libraries
library(tidyr)
library(ggplot2)
library(EnvStats)
library(lubridate)
library(gridExtra)

#-----------Combining Data----

mac$Rep <- gsub("[^A-Za-z]", "", mac$`Sample ID`)  # Remove everything that's not a letter
mac$Position <- gsub("[^0-9]", "", mac$`Sample ID`)  # Remove everything that's not a digit

nuts$Rep <- gsub("[^A-Za-z]", "", nuts$`Sample ID`)  # Remove everything that's not a letter
nuts$Position <- gsub("[^0-9]", "", nuts$`Sample ID`)  # Remove everything that's not a digit

colnames(amf)[colnames(amf) == "Treatment"] <- "Treat" 
colnames(nuts)[colnames(nuts) == "Treatment"] <- "Treat"
colnames(emh)[colnames(emh) == "Treatment"] <- "Treat"

colnames(amf)[colnames(amf) == "Sample_ID"] <- "Sample"
colnames(nuts)[colnames(nuts) == "Sample ID"] <- "Sample" 
colnames(emh)[colnames(emh) == "Sample ID"] <- "Sample"
colnames(mac)[colnames(mac) == "Sample ID"] <- "Sample"

mac$Genotype <- gsub("\\.0$", "", as.character(mac$Genotype))
amf$Genotype <- gsub("\\.0$", "", as.character(amf$Genotype))
emh$Genotype <- gsub("\\.0$", "", as.character(emh$Genotype))
nuts$Genotype <- gsub("\\.0$", "", as.character(nuts$Genotype))

mac$Position <- as.integer(mac$Position)
nuts$Position <- as.integer(nuts$Position)
emh$Position <- as.integer(emh$Position)
amf$Position <- as.integer(amf$Position)

macombo <- mac %>%
  left_join(amf, by = c("Genotype", "Treat", "Rep","Position", "Sample")) %>%
  left_join(emh, by = c("Genotype", "Treat", "Rep","Position", "Sample")) %>%
  left_join(nuts,by = c("Genotype", "Treat", "Rep","Position", "Sample"))


g2c$Rep <- as.integer(g2c$Rep)
g2c$Genotype <- gsub("\\.0$", "", as.character(g2c$Genotype))
g2camf$Genotype <- gsub("\\.0$", "", as.character(g2camf$Genotype))

g2combo <- g2c %>%  
  left_join(g2camf, by = c("Genotype", "Treat", "Rep"))

#------macombo plots----

# Filter out 0 values and NA values
macombo$AMFInDrySoil <- as.numeric(as.character(macombo$AMFInDrySoil))

#subset the data in macombo into just the genotypes scored for emh
macombosubset <- macombo %>%
  filter(!is.na(TOT) & 
           TOT != 0 & 
           TOT != "#DIV/0!")

# Fit linear models for each treatment
model_watered <- lm(AMFInDrySoil ~ RLC_P, data = macombosubset %>% filter(Treat == "Watered"))
model_droughted <- lm(AMFInDrySoil ~ RLC_P, data = macombosubset %>% filter(Treat == "Droughted"))

# Extract slope and R-squared values
slope_watered <- coef(model_watered)[2]  # Slope for "Watered"
r2_watered <- summary(model_watered)$r.squared  # R-squared for "Watered"

slope_droughted <- coef(model_droughted)[2]  # Slope for "Droughted"
r2_droughted <- summary(model_droughted)$r.squared  # R-squared for "Droughted"

# Plot
macombosubset %>%
  ggplot(aes(x = RLC_P, y = AMFInDrySoil, fill = Genotype, color = Treat)) +
  geom_jitter(width = 0.2, height = 0.2, size = 3, shape = 21, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  scale_shape_manual(values = c("Watered" = 16, "Droughted" = 17)) +
  scale_color_manual(values = c("Watered" = "lightblue", "Droughted" = "darkred")) +
  geom_text(
    aes(x = 0.1, y = max(AMFInDrySoil) * 0.8, label = paste("Watered: Slope =", round(slope_watered, 2),
                                                            "\nR² =", round(r2_watered, 2))),
    color = "lightblue", hjust = 0, size = 3
  ) +
  geom_text(
    aes(x = 0.1, y = max(AMFInDrySoil) * 0.9, label = paste("Droughted: Slope =", round(slope_droughted, 2),
                                                            "\nR² =", round(r2_droughted, 2))),
    color = "darkred", hjust = 0, size = 3
  )







# Create the AMF boxplot# Create the AMF boxplotGenotype
macombo %>%
  ggplot(aes(x = Genotype, y = AMFInDrySoil, fill = Treat)) +
  geom_boxplot() +
  stat_n_text()+
  labs(y = "AMF in Dry Soil (m/g)")

# Create the DSE boxplot
macombo %>%
  ggplot(aes(x = Genotype, y = DSEInDrySoil, fill = Treat)) +
  geom_boxplot() +
  stat_n_text()+
  labs(y = "DSE in Dry Soil (m/g)")

#violin of Treat
emh_filtered %>%
  ggplot(aes(x = Treat, y = AMFInDrySoil, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank())



#-------g2combo plots----

g2combo$TotPlantDryWt <- as.numeric(as.character(g2combo$TotPlantDryWt))

# Filter out NA and 0 values from amf
g2combosubset <- g2combo %>%
  filter(!is.na(TOT) & TOT != 0 & Treat != "Dead")

# Calculate the IQR and remove outliers
Q1 <- quantile(g2combosubset$RLC_Percent, 0.25)
Q3 <- quantile(g2combosubset$RLC_Percent, 0.75)
IQR <- Q3 - Q1

g2combosubset_no_outliers <- g2combosubset %>%
  filter(RLC_Percent >= (Q1 - 1.5 * IQR) & RLC_Percent <= (Q3 + 1.5 * IQR))

# Plot the filtered data
ggplot(data = g2combosubset_no_outliers, aes(x = Genotype, y = RLC_Percent)) +
  geom_boxplot(aes(fill = Genotype)) +  # Color boxplot by 'Genotype'
  theme(
    legend.position = "none",  # Remove the legend
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate Genotype labels by 45 degrees
  ) + 
  scale_x_discrete(labels = function(x) gsub("\\.0$", "", x)) +  # Remove ".0" from labels
  stat_n_text() +  # Add sample size text inside the boxplot
  labs(y = "% Root Length Colonized", x = "Genotype",title = "G2C Root Length Colonized by Genotype")  # Y label and X label





#violin of Treat
g2combosubset_clean <- g2combosubset %>%
  filter(!is.na(Genotype), !is.na(Arb), !is.na(TOT)) %>%  # Remove rows with NA in 'Genotype', 'Arb', or 'TOT'
  mutate(Arb_perc = (Arb / TOT) * 100)  # Calculate percentage

# Plot
g2combosubset_clean %>%
  ggplot(aes(x = Genotype, y = Arb_perc, fill = Genotype)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with transparency
  stat_n_text() +  # Adds count labels on violin plots
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  labs(y = "% Root Length of Arbuscules", x = "")



g2combo %>%
  ggplot(aes(x = AM_Hyphae/TOT, y = Arb/TOT, fill = Treat, color = Treat)) +  # Color points and trendlines based on 'Treat'
  geom_jitter(width = 0.2, height = 0.2, size = 3, shape = 21, alpha = 0.7) +  # Customize jittered points
  geom_smooth(method = "lm", se = FALSE) +  # Add trendline with color based on 'Treat'
  theme_minimal() +  # Apply a minimal theme
  scale_fill_manual(values = c("Live" = "forestgreen", "Dead" = "darkred")) +  # Custom colors for points
  scale_color_manual(values = c("Live" = "forestgreen", "Dead" = "darkred")) +  # Custom colors for trendlines
  geom_text(
    aes(x = .1, y = .7, label = paste("Live: Slope =", round(slope_watered, 2), 
                                      "\nR² =", round(r2_watered, 2))),
    color = "lightblue", hjust = 0, size = 3
  ) +  # Annotate Watered slope and R²
  geom_text(
    aes(x = .1, y = 1, label = paste("Dead: Slope =", round(slope_droughted, 2), 
                                     "\nR² =", round(r2_droughted, 2))),
    color = "darkred", hjust = 0, size = 3
  )  # Annotate Droughted slope and R²






# Create the 'mgr' variable with the correct calculation
ds <- g2combo %>%
  group_by(Genotype, Treat) %>% # Grouping by Genotype and Treat
  mutate(mgr = log((TotPlantDryWt) / mean((TotPlantDryWt)[Treat == "Dead"], na.rm = TRUE))) %>% 
  ungroup() # Ungroup after mutate to avoid issues with further operations

# Remove outliers using IQR method
Q1 <- quantile(ds$mgr, 0.25, na.rm = TRUE)
Q3 <- quantile(ds$mgr, 0.75, na.rm = TRUE)
IQR_value <- IQR(ds$mgr, na.rm = TRUE)

ds_no_outliers <- ds %>%
  filter(mgr >= (Q1 - 1.5 * IQR_value) & mgr <= (Q3 + 1.5 * IQR_value))

# Arrange bars by mean of 'mgr' per 'Genotype'
ds_no_outliers <- ds_no_outliers %>%
  group_by(Genotype) %>%
  mutate(mean_mgr = mean(mgr, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Genotype = factor(Genotype, levels = unique(Genotype[order(mean_mgr,decreasing= TRUE)])))

# Create the boxplot
ggplot(ds_no_outliers, aes(x = Genotype, y = mgr, fill = Genotype)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + # Horizontal line at y=0
  labs(title = "Total Weight Response Ratio", y = "Live Over Sterile Inoculation") +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), # Hide legend and remove x-axis label
        axis.text.x = element_text(angle = 45, hjust = 1)) # Angle x-axis labels by 45 degrees






#----------MAC 2022 Probe Data----

# Add a 'Source' column to differentiate between the two datasets
dsdrought$Source <- "Drought"
dswet$Source <- "Wet"


# Combine the datasets by stacking them together
combined_data <- bind_rows(dsdrought, dswet)

# Convert Timestamps to datetime (POSIXct) if not already done
# Assume Timestamps are currently in 'YYYY-MM-DD' format without times, add a dummy time if needed
dsdrought$Timestamps <- as.POSIXct(dsdrought$Timestamps, format = "%Y-%m-%d", tz = "UTC")
dswet$Timestamps <- as.POSIXct(dswet$Timestamps, format = "%Y-%m-%d", tz = "UTC")

# Reshape data from wide to long format for easy processing
combined_data_long <- combined_data %>%
  pivot_longer(cols = starts_with("Port"), 
               names_to = "Port", 
               values_to = "Value")

# Calculate the average and standard deviation for each Port, grouped by Source and Port
summary_data <- combined_data_long %>%
  group_by(Source, Port, Timestamps) %>%
  summarise(Average_Value = mean(Value, na.rm = TRUE),
            StdDev = sd(Value, na.rm = TRUE),
            .groups = "drop")

# Create the line graph with averaged values and error shadows (standard deviation)
ggplot(data = summary_data, 
       aes(x = Timestamps, y = Average_Value, color = Port, group = interaction(Source, Port))) +
  geom_line() +  # Create line plot
  geom_ribbon(aes(ymin = Average_Value - StdDev, ymax = Average_Value + StdDev, fill = Port), 
              alpha = 0.2) +  # Add error shadow (ribbon)
  facet_wrap(~ Source, scales = "free_y") +  # Facet by Port to have separate plots for each
  #grid.arrange( ncol = 1) +  # ncol = 1 arranges the plots in one column (top-bottom)
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +  # Remove X title
  labs(y = "m³/m³ Water Content at 25in")+  # Label y-axis and legend
  scale_x_datetime(
    breaks = scales::date_breaks("4 days"),  # Adjust for monthly breaks (change to "1 day", "2 weeks", etc.)
    labels = scales::date_format("%b %d"))  # Format with month-day-year hour:minute


# Create a new column to define whether a port is for Drought or Wet conditions
summary_data$Condition <- ifelse(as.numeric(gsub("Port", "", summary_data$Port)) <= 6, "Drought", "Wet")

# Define color scales for each condition
ggplot(data = summary_data, 
       aes(x = Timestamps, y = Average_Value, color = Port, group = interaction(Source, Port))) +
  geom_line(size = 1) +  # Make the lines thicker by setting the size (e.g., 1 or higher)
  facet_wrap(~ Source, scales = "free_y", ncol = 1) +  # Stack the graphs vertically, one on top of the other
  theme_minimal() +  # Minimal theme for cleaner visualization
  theme(
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "none",  # No legend displayed
    panel.grid.major = element_line(size = 0.2, color = "gray80"),  # Lighter grid lines
    panel.grid.minor = element_blank()  # No minor grid lines
  ) +  
  labs(
    y = "m³/m³ Water Content at 25in",  # Y-axis label
    x = "Timestamp"  # X-axis label
  ) +
  scale_x_datetime(
    breaks = scales::date_breaks("2 days"),  # Adjust for breaks every 4 days
    labels = scales::date_format("%b %d")  # Format as month-day
  ) +
  # Define two separate gradients for color (one for Drought, one for Wet)
  scale_color_manual(
    values = c(
      # Red to Green for Drought (Port1 to Port6)
      setNames(colorRampPalette(c("darkred", "chartreuse3"))(6), paste0("Port", 1:6)),
      # Green to Blue for Wet (Port7 to Port12)
      setNames(colorRampPalette(c("darkorange", "darkblue"))(6), paste0("Port", 7:12))
    )
  )


ggplot(data = summary_data, 
       aes(x = Timestamps, y = Average_Value, color = Source, group = Source)) +
  geom_line() +  # Create line plot
  theme(axis.title.x = element_blank()) +  # Remove X title
  labs(y = "m³/m³ Water Content at 25in", color = "Treat") +  # Label y-axis and legend
  scale_x_datetime(
    breaks = scales::date_breaks("2 days"),  # Adjust for monthly breaks (change to "1 day", "2 weeks", etc.)
    labels = scales::date_format("%b %d")  # Format with month-day-year hour:minute
  ) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
  )



#---------MAC 2022 EMH Graph-----

# Filter out 0 values and NA values
emh_filtered <- emh %>%
  filter(!is.na(AMFInDrySoil) & 
           AMFInDrySoil != 0 & 
           AMFInDrySoil != "#DIV/0!"&
         !is.na(DSEInDrySoil))

emh_filtered$AMFInDrySoil <- as.numeric(as.character(emh_filtered$AMFInDrySoil))

# Create the AMF boxplot
emh_filtered %>%
  ggplot(aes(x = Genotype, y = AMFInDrySoil, fill = Treat)) +
  geom_boxplot() +
  stat_n_text()+
  labs(y = "AMF in Dry Soil (m/g)")

# Create the DSE boxplot
emh_filtered %>%
  ggplot(aes(x = Genotype, y = DSEInDrySoil, fill = Treat)) +
  geom_boxplot() +
  stat_n_text()+
  labs(y = "DSE in Dry Soil (m/g)")

#violin of Treat
emh_filtered %>%
  ggplot(aes(x = Treat, y = AMFInDrySoil, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank())


#--------Mac 2022 Colonization Graphs------

# Filter out NA and 0 values from amf
amf_filtered <- amf %>%
  filter(!is.na(RLC_P) & RLC_P != 0) %>%
  group_by(Genotype) %>%
  filter(sum(!is.na(RLC_P)) >= 5) %>%  # Ensure at least 5 non-NA values for each Genotype
  ungroup()

# Plot the filtered data
ggplot(data = amf_filtered,
       aes(x = Treat, y = RLC_P)) +
  geom_boxplot(aes(fill = Treat)) +  # color boxplot
  theme(axis.title.x = element_blank(),  # remove X title
        axis.text.x = element_blank()) +  # remove X axis labels
  facet_grid(~Genotype) +  # facet by Genotype
  stat_n_text() +  # Add sample size text
  labs(y = "% Root Length Colonized", x = "")

#violin of Treat
amf_filtered %>%
  ggplot(aes(x = Treat, y = (Arb/TOT)*100, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + 
  labs(y = "% Root Length of Arbuscules", x = "")


#------t-Test------

t.test(x, y = NULL,
       alternative = c("two.sided", "less", "greater"),
       mu = 0, paired = FALSE, var.equal = FALSE,
       conf.level = 0.95, ...)


#-------Anova------

one.way <- aov(MHyphDrySoil ~ Genotype, data= emh_filtered)
summary(one.way)
plot(one.way)

two.way <- aov(MHyphDrySoil ~ Genotype * Treat, data= emh_filtered)
summary(two.way)
plot(two.way)

blocking <- aov(MHyphDrySoil ~ Genotype + Treat + Position, data= emh_filtered)
summary(blocking)
plot(blocking)


one.way <- aov(RLC_P ~ Genotype, data= amf_filtered)
summary(one.way)
plot(one.way)

two.way <- aov(RLC_P ~ Genotype * Treat, data= amf_filtered)
summary(two.way)
plot(two.way)

blocking <- aov(RLC_P ~ Genotype + Treat + Position, data= amf_filtered)
summary(blocking)
plot(blocking)


model.set <- list(one.way, two.way, interaction, blocking)
model.names <- c("one.way", "two.way", "interaction", "blocking")

aictab(model.set, modnames = model.names)

par(mfrow=c(2,2))
plot(two.way)
par(mfrow=c(1,1))


# 95% confidence interval

tukey.two.way<-TukeyHSD(two.way)

tukey.two.way

tukey.plot.aov<-aov(AP_AM_tot ~ TP:Genotype, data=MAC_P)
tukey.plot.test<-TukeyHSD(tukey.plot.aov)
plot(tukey.plot.test, las = 1)

Mean.AP.AM.tot.data <- MAC_P %>%
  group_by(TP, Genotype) %>%
  summarise(AM_tot = mean(AP_AM_tot))

Mean.AP.AM.tot.data$group <- c("a","b","b","c")

Mean.AP.AM.tot.data

two.way.plot <- ggplot(Mean.AP.AM.tot.data, 
                       aes(x = Genotype, y = AM_tot, group=TP)) +
  geom_point(cex = 1.5, pch = 1.0,position = position_jitter(w = 0.1, h = 0))

two.way.plot <- two.way.plot +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0.2) +
  stat_summary(fun.data = 'mean_se', geom = 'pointrange') +
  geom_point(data=Mean.AP.AM.tot.data, aes(x=Genotype, y=AM_tot))

two.way.plot

#-----Legacy Plot----

dsLegacyfiltered <- dsLegacy %>%
  filter(!is.na(ShootWt))

# Calculate IQR and filter out outliers
Q1 <- quantile(dsLegacyfiltered$ShootWt, 0.25)
Q3 <- quantile(dsLegacyfiltered$ShootWt, 0.75)
IQR_value <- IQR(dsLegacyfiltered$ShootWt)

# Filter out outliers
dsLegacy_no_outliers <- dsLegacyfiltered %>%
  filter(ShootWt >= (Q1 - 1.5 * IQR_value) & ShootWt <= (Q3 + 1.5 * IQR_value))

# Create the boxplot without outliers
dsLegacy_no_outliers %>%
  ggplot(aes(x = GenoLetter, y = ShootWt, fill = Treat)) +
  geom_boxplot() +
  stat_n_text() +  # Adds count labels on boxplots
  theme(axis.title.x = element_blank())  # Removes x-axis title

#violin of Treat
dsLegacy_no_outliers %>%
  ggplot(aes(x = Treat, y = ShootWt, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank())

# Boxplot with Treat as fill
dsLegacy %>%
  ggplot(aes(x = GenoLetter, y = ShootWt, fill = Treat)) +
  geom_boxplot() +
  stat_n_text()+
  theme(axis.title.x = element_blank())  # Adds count labels on boxplots

#-------MAC Plots----

# Remove NAs
macfiltered <- mac %>%
  filter(!is.na(ShootWt))  %>%
  filter(!is.na(Genotype)) %>%
  filter(!is.na(Treat)) %>%
  filter(!is.na(Florets))

# Convert ShootWt to numeric if necessary
macfiltered$ShootWt <- as.numeric(as.character(macfiltered$ShootWt))
macfiltered$Florets <- as.numeric(as.character(macfiltered$Florets))

# Boxplot with Treat as fill
macfiltered %>%
  ggplot(aes(x = Genotype, y = ShootWt, fill = Treat)) +
  geom_boxplot() +
  stat_n_text()+
  theme(axis.title.x = element_blank(),  # Removes x-axis title
        axis.text.x = element_text(angle = 45, hjust = 1))  # Angles the x-axis labels

macfiltered %>%
  ggplot(aes(x = Genotype, y = Florets, fill = Treat)) +
  geom_boxplot() +
  stat_n_text()+
  theme(axis.title.x = element_blank(),  # Removes x-axis title
        axis.text.x = element_text(angle = 45, hjust = 1))  # Angles the x-axis labels


macfiltered %>%
  # Calculate the average Florets within each Genotype-Treat combination
  group_by(Genotype, Treat) %>%
  summarise(avg_Florets = mean(Florets, na.rm = TRUE)) %>%
  ungroup() %>%
  # Reorder Genotype based on the average Florets (median for ordering)
  mutate(Genotype = reorder(Genotype, avg_Florets, FUN = median)) %>%
  # Now join back to the original dataset to keep `Florets` for plotting
  left_join(macfiltered, by = c("Genotype", "Treat")) %>%
  # Create the boxplot with the reordered Genotype
  ggplot(aes(x = Genotype, y = Florets, fill = Treat)) +
  geom_boxplot() +
  stat_n_text() +  # Adds count labels on boxplots
  theme(axis.title.x = element_blank(),  # Removes x-axis title
  axis.text.x = element_text(angle = 45, hjust = 1))  # Angles the x-axis labels


# Violin plot 
macfiltered %>%
  ggplot(aes(x = Treat, y = ShootWt, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank())

# Violin plot 
macfiltered %>%
  ggplot(aes(x = Treat, y = Florets, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank())

macfiltered %>%
  ggplot(aes(x = Treat, y = Florets, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank())


dmgrfiltered <- dmgr %>%
  filter(!is.na(ShootWt))  %>%
  filter(!is.na(Genotype)) %>%
  filter(!is.na(Treat)) %>%
  filter(!is.na(DryRootWt))

# Violin plot 
dmgrfiltered %>%
  ggplot(aes(x = Treat, y = DryRootWt+ShootWt, fill = Treat)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.3, outlier.shape = NA) +  # Boxplot with some transparency
  stat_n_text()+  # Adds count labels on violin plots
  theme(axis.title.x = element_blank())


dmgrfiltered %>%
  filter(NumOfPlants %in% c("1","2","3","4","5")) %>%
  ggplot(aes(x = GenoNum, y = DryRootWt+ShootWt, fill = Treat)) +
  geom_boxplot()+
  stat_n_text()


dmgrfiltered %>%
  filter(NumOfPlants %in% c("1")) %>%
  ggplot(aes(x = GenoNum, y = DryRootWt+ShootWt, fill = Treat)) +
  geom_boxplot()+
  stat_n_text()

dmgrfiltered$NumOfPlants <- as.numeric(as.character(dmgrfiltered$NumOfPlants))

# Compute median ShootWt for each NumOfPlants
order_df <- dmgrfiltered %>%
  filter(NumOfPlants %in% c("1", "2", "3")) %>%
  group_by(NumOfPlants) %>%
  summarize(median_totwt = median(DryRootWt+ShootWt, na.rm = TRUE)) %>%
  arrange(median_totwt) %>%
  mutate(NumOfPlants = factor(NumOfPlants, levels = NumOfPlants))

# Join the order_df with the original data to get the ordered levels
ds_ordered <- dmgrfiltered %>%
  filter(NumOfPlants %in% c("1", "2", "3")) %>%
  left_join(order_df, by = "NumOfPlants") %>%
  mutate(NumOfPlants = factor(NumOfPlants, levels = order_df$NumOfPlants))

# Plot
ggplot(dmgrfiltered, aes(x = GenoNum, y = DryRootWt+ShootWt/NumOfPlants, fill = NumOfPlants)) +
  geom_boxplot() +
  stat_summary(aes(label = ..count..), fun = "count", geom = "text", vjust = -0.5) +
  facet_grid(~NumOfPlants) +
  theme_minimal() +
  labs(title = "Boxplot of Shoot Weight by Genotype and Number of Plants",
       x = "Genotype Number",
       y = "Shoot Weight")


#------DMGR Plots------

ds <- dmgr %>%
  filter(NumOfPlants == "1")

# Create the 'mgr' variable with the correct calculation
ds <- dmgr %>%
  group_by(Genotype, Treat) %>% # Grouping by Genotype and Treat
  mutate(mgr = log((DryRootWt+ShootWt) / mean((DryRootWt+ShootWt)[Treat == "S"], na.rm = TRUE))) %>% 
  ungroup() # Ungroup after mutate to avoid issues with further operations

# Remove outliers using IQR method
Q1 <- quantile(ds$mgr, 0.25, na.rm = TRUE)
Q3 <- quantile(ds$mgr, 0.75, na.rm = TRUE)
IQR_value <- IQR(ds$mgr, na.rm = TRUE)

ds_no_outliers <- ds %>%
  filter(mgr >= (Q1 - 1.5 * IQR_value) & mgr <= (Q3 + 1.5 * IQR_value))

# Arrange bars by mean of 'mgr' per 'Genotype'
ds_no_outliers <- ds_no_outliers %>%
  group_by(Genotype) %>%
  mutate(mean_mgr = mean(mgr, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Genotype = factor(Genotype, levels = unique(Genotype[order(mean_mgr,decreasing= TRUE)])))

# Create the boxplot
ggplot(ds_no_outliers, aes(x = Genotype, y = mgr, fill = Genotype)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + # Horizontal line at y=0
  labs(title = "Total Weight Response Ratio", y = "Live Over Sterile Inoculation") +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), # Hide legend and remove x-axis label
        axis.text.x = element_text(angle = 45, hjust = 1)) # Angle x-axis labels by 45 degrees


#-------- MAC 2022 amf ------

ggplot(data = g2combo,
       aes(x = Genotype, y =Olpidium/TOT))+
  geom_boxplot(aes(fill = Treat))+ # color boxplot
  theme(axis.title.x=element_blank()) #remove X title)
  #facet_grid(dmgr#facet_grid(~Genotype)+
  #stat_n_text()


amf_filtered <- amf %>%
  group_by(Genotype) %>%
  filter(sum(!is.na(RLC_P)) >= 5) %>%  # Filter genotypes with 5 or more non-NA RLC_P values
  ungroup()

ggplot(data = amf_filtered,
  aes(x = Treat, y = RLC_P)) +
  geom_boxplot(aes(fill = Treat)) +  # color boxplot
  theme(axis.title.x = element_blank()) +  # remove X title
  facet_grid(~Genotype)+
  stat_n_text()



# Join the two datasets by a common column, e.g., 'Genotype'
amf_filtered <- amf %>%
  group_by(Genotype) %>%
  filter(sum(!is.na(RLC_P)) >= 5) %>%  # Filter genotypes with 5 or more non-NA RLC_P values
  ungroup() %>%
  left_join(dmgr, by = "Genotype")  # Assuming 'Genotype' is the common column

ds <- amf_filtered %>%
  group_by(Genotype, Treat)%>% #We want to be doing MGR in groups by Genotype and Treat.
  mutate(mgr =  log(RLC_P / (mean(RLC_P[Treat == "Watered"], na.rm = TRUE))), 
         NA_real_)

ds %>%
  group_by(Genotype) %>%
  ggplot(aes(x = Genotype, 
             y = mgr, 
             fill = Genotype)) +
  geom_boxplot() +
  stat_n_text() +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank())  # Hide legend and rotate x-axis text
#facet_grid(~ Treat, labeller = labeller(soil = soil_labels))

# Plot with faceting by NumOfPlants from the dmgr dataset
ggplot(data = ds,
       aes(x = Genotype, y = mgr)) +
  geom_boxplot(aes(fill = Genotype)) +  # Color boxplot
  theme(axis.title.x = element_blank()) +  # Remove X title
  facet_grid(~NumOfPlants)+ # Facet by NumOfPlants from the dmgr dataset
  stat_n_text() 


  
# Filter and prepare the data
amf_filtered <- amf %>%
  group_by(Genotype) %>%
  filter(sum(!is.na(RLC_P)) >= 0) %>%  # Filter genotypes with 5 or more non-NA RLC_P values
  ungroup() %>%
  left_join(mac, by = "Genotype")  # Join with mac dataset on Genotype

# Calculate mgr (magnesium response ratio) by Genotype and Treat
ds <- amf_filtered %>%
  group_by(Genotype, Treat) %>%
  mutate(mgr = log(RLC_P / mean(RLC_P[Treat == "Watered"], na.rm = TRUE))) %>%
  ungroup()  # Remove grouping for further operations

# Boxplot for ShootWtg with count labels (stat_n_text should work if it's from a custom package)
ggplot(ds, aes(x = Genotype, y = ShootWtg, fill = Treat)) +
  geom_boxplot() +
  stat_n_text() +  # Assuming stat_n_text() is part of a custom package or extension
  theme_bw() +  # Use the 'black and white' theme
  theme(legend.position = "none", axis.title.x = element_blank())  # Hide legend and remove x-axis title



