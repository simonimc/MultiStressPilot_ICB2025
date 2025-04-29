#This code is for preliminary analyses presented in Integrative and Comparative Biology
#Simonis et al. 2025 "A collaborative multiple stressor approach for identifying spatial
#heterogeneities in wildlife health and conservation priorities."

#Code Author: Molly Simonis
#Last Updated: 28 April 2025


#Make sure everything is clear in the workspace 
rm(list=ls())##clear environment
graphics.off()##clear graphics
gc()##free up some RAM

##run cleanup code first for the starting df
#source('MultiStressors_CleanupCode.R')
#
##clear workspace again (df will be saved as a csv in working directory)
rm(list=ls())##clear environment
graphics.off()##clear graphics
gc()##free up some RAM

library(lme4)
library(lmerTest)
library(ggplot2)
library(car)
library(MuMIn)
library(mgcv)
library(statmod)
library(terra)
library(sf)
library(dplyr)
library(lubridate)
library(factoextra)
library(ggbiplot)
library(ggfortify)
library(ggeffects)
library(patchwork)

#theme for plotting later
th=theme_bw()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))+
  theme(plot.title = element_text(hjust = 0.5))


#read in data 
ms_df<- read.csv('MultiStressors2023_PilotData.csv', header = T, sep = ',')

#read in species home range data
sp_hr<- read.csv('Species_HomeRange.csv', header = T, sep = ',')

#merge data and home range info by species code
ms_df<- merge(ms_df, sp_hr, by = 'species_code')

#read in 2019 National Land Cover tif file as a raster
#can download the most recent version in the future at 
#https://www.mrlc.gov/data
nlcd_raster <- rast("Annual_NLCD_LndCov_2019_CU_C1V0.tif")

#check info on the raster file
print(nlcd_raster)

#check out the map
plot(nlcd_raster, main = "NLCD 2019 Land Cover")



#make lat/lon in data spatial points
# CRS (Coordinate Reference System) for geographic coordinates (EPSG:4326 - WGS84)
ms_df_spatial <- st_as_sf(ms_df, coords = c("lon", "lat"), crs = 4326)

#create projection to eventually get buffers
# EPSG:3857 (Web Mercator) works globally for buffers in meters
ms_df_projected <- st_transform(ms_df_spatial, crs = 3857)

#create buffers 
#need to be in m and data is in km
buffers <- st_buffer(ms_df_projected, dist = ms_df$home_range_radius_km*1000)

#transform the buffers back to geographic coordinates (optional)
buffers_geographic <- st_transform(buffers, crs = 4326)

#convert buffers to terra's SpatVector format for compatibility with the raster
buffers_vect <- vect(buffers_geographic)

#ensure buffer projection matches the raster CRS
buffers_vect <- project(buffers_vect, crs(nlcd_raster))


#extract raster values for each buffer
# terra::extract extracts pixel values for each buffer
extracted_values <- terra::extract(nlcd_raster, buffers_vect, fun = NULL, df = TRUE)

#add buffer IDs to the extracted values
#`ID` is automatically assigned by `extract()` to correspond to each buffer
names(extracted_values)<- c('buffer_id', 'land_cover')

#calculate land use proportions for each buffer
land_use_proportions <- extracted_values %>%
  dplyr::group_by(buffer_id, land_cover) %>%
  dplyr::summarise(count = length(land_cover), .groups = "drop") %>%  # Count occurrences of each land cover type
  dplyr::group_by(buffer_id) %>%
  dplyr::mutate(proportion = count / sum(count)) %>% # Calculate proportions
  dplyr::arrange(buffer_id, desc(proportion))

#create landcover descriptions
nlcd_descrip <- data.frame(
  land_cover = c(11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95),
  description = c(
    "Open Water", "Developed, Open Space", "Developed, Low Intensity", 
    "Developed, Medium Intensity", "Developed, High Intensity", 
    "Barren Land", "Deciduous Forest", "Evergreen Forest", 
    "Mixed Forest", "Shrub/Scrub", "Herbaceous", "Hay/Pasture", 
    "Cultivated Crops", "Woody Wetlands", "Emergent Herbaceous Wetlands"
  )
)

#merge proportions and descriptions together by 'land_cover'
land_use_proportions<- merge(land_use_proportions, nlcd_descrip, by = 'land_cover')

#subset land us prop df
land_use_proportions <- land_use_proportions %>%
  select(buffer_id, proportion, description)

#make the landuse df wide by buffer_id
land_use_wide <- land_use_proportions %>%
  pivot_wider(
    names_from = description,          # Use land cover descriptions as column names
    values_from = proportion,          # Proportions will be the values for these columns
    values_fill = list(proportion = 0) # Fill missing proportions with 0 (if any buffer_id does not have that land cover)
  )


#sort buffer IDs in order so they match row numbers in ms_df
land_use_wide <- land_use_wide %>%
  arrange(buffer_id)


#merge land use with ms_df
ms_df<- as.data.frame(cbind(ms_df, land_use_wide))


#make all character variables factors
ms_df <- ms_df %>%
  mutate_if(is.character, as.factor)

#re-level rep_stat factor
ms_df$rep_stat<- factor(ms_df$rep_stat, levels = 
                          c('non-reproductive', 'pregnant', 'lactating', 'testes-descended'))

#make dates again
ms_df$date<- ymd(ms_df$date)
#make dates a julian day
ms_df$jday<- yday(ms_df$date)




#perform PCA for all the stressors
#make a PCA matrix first
pca_m<- ms_df[ms_df$sex == 'female',]
pca_m<- pca_m[, c(36, 6, 41:55)]

#the PCA (with scaling)
ext_stress_pca <- prcomp(pca_m[,3:17], scale = TRUE)
#calculate total variance explained by each principal component
var_explained <- ext_stress_pca$sdev^2 / sum(ext_stress_pca$sdev^2)
print(var_explained)
#create scree plot
fviz_eig(ext_stress_pca)
#get loadings/eigenvectors for the first 2 PCs
sort(ext_stress_pca$rotation[, 1])
sort(ext_stress_pca$rotation[, 2])



#Figure S1
sp_pca<- autoplot(ext_stress_pca, data = pca_m, color = 'species_name', 
         size = 5, loadings = TRUE, loadings.colour = 'gray25', 
         loadings.label = TRUE, loadings.label.size = 3, 
         loadings.label.color = 'gray25', scale = 0) + 
  th +
  labs(x = 'PC1 (57%)', y = 'PC2 (38%)', color = 'species')



state_pca<- autoplot(ext_stress_pca, data = pca_m, color = 'state', 
         size = 5, loadings = TRUE, loadings.colour = 'gray25', 
         loadings.label = TRUE, loadings.label.size = 3, 
         loadings.label.color = 'gray25', scale = 0) + 
  th +
  labs(x = 'PC1 (57%)', y = 'PC2 (38%)', color = 'state')

state_pca + sp_pca + plot_layout(axis_titles = 'collect')


#extract PC1 and PC2 axes values
axes <- as.data.frame(predict(ext_stress_pca, newdata = ms_df))

#add to your ms_df
ms_df$PC1<- axes$PC1
ms_df$PC2<- axes$PC2


#make Hg in mg/kg for models
ms_df$Hg_mg_kg<- ms_df$Hg_ng_g/1000
#make nl ratios <0.00001 to 0.00001 (still functionally >0, but won't get flagged as 0 on log scale)
#ms_df$nl[ms_df$nl < 0.0001]<- 0.0001
#make sure repstat and gltA levels are in order
ms_df$rep_stat<- factor(ms_df$rep_stat, levels = c('non-reproductive', 'pregnant', 'lactating'))
ms_df$gltA<- factor(ms_df$gltA, levels = c('Positive', 'Negative'))

#make foraging distance categorical
ms_df$forg_cat<- ifelse(ms_df$home_range_radius_km < 10, 'short', 'long')

#make ms_df female only
ms_df<- ms_df[ms_df$sex == 'female', ]

#check correlation of mass and day of capture
mass_jday_corr<- cor(ms_df[, c(10, 56)], method = 'pearson')
mass_jday_corr

#create models
#using the gam() function, but still creating glms
#gam() function will provide an R^2 with a tweedie family
mod1<- gam(nl ~ mn_intensity*rep_stat + 
             hold_time, data = ms_df,
           family = tw, method = 'ML')
summary(mod1)
gam.check(mod1)
anova.gam(mod1)

mod2<- gam(nl ~ mn_intensity*gltA + 
             hold_time, data = ms_df,
           family = tw, method = 'ML')
summary(mod2)
gam.check(mod2)
anova.gam(mod2)

mod3<- gam(nl ~ Hg_mg_kg*rep_stat + 
                     hold_time, data = ms_df,
           family = tw, method = 'ML')
summary(mod3)
gam.check(mod3)
anova.gam(mod3)

mod4<- gam(nl ~ Hg_mg_kg*gltA + 
             hold_time, data = ms_df,
           family = tw, method = 'ML')
summary(mod4)
gam.check(mod4)
anova.gam(mod4)

mod5<- gam(nl ~ PC1*rep_stat + 
             hold_time, data = ms_df,
           family = tw, method = 'ML')
summary(mod5)
gam.check(mod5)
anova.gam(mod5)

mod6<- gam(nl ~ PC1*gltA + 
             hold_time, data = ms_df,
           family = tw, method = 'ML')
summary(mod6)
gam.check(mod6)
anova.gam(mod6)

mod7<- gam(nl ~ PC2*rep_stat + 
             hold_time, data = ms_df,
           family = tw, method = 'ML')
summary(mod7)
gam.check(mod7)
anova.gam(mod7)

mod8<- gam(nl ~ PC2*gltA + 
             hold_time, data = ms_df,  
           family = tw, method = 'ML')
summary(mod8)
gam.check(mod8)
anova.gam(mod8)

#mod9<- gam(nl ~ home_range_sqkm*rep_stat + 
#             hold_time, data = ms_df,   
#           family = tw, method = 'ML')
mod9<- gam(nl ~ forg_cat*rep_stat + 
             hold_time, data = ms_df,   
           family = tw, method = 'ML')
summary(mod9)
gam.check(mod9)
anova.gam(mod9)

mod10<- gam(nl ~ forg_cat*gltA 
            + hold_time, data = ms_df,  
            family = tw, method = 'ML')
summary(mod10)
gam.check(mod10)
anova.gam(mod10)

mod11<- gam(nl ~ mn_intensity + rep_stat + hold_time, data = ms_df,  
            family = tw, method = 'ML')
summary(mod11)
gam.check(mod11)
anova.gam(mod11)


mod12<- gam(nl ~ mn_intensity + gltA + hold_time, data = ms_df,  
            family = tw, method = 'ML')
summary(mod12)
gam.check(mod12)
anova.gam(mod12)



mod13<- gam(nl ~ Hg_mg_kg + rep_stat + hold_time, data = ms_df,  
            family = tw, method = 'ML')
summary(mod13)
gam.check(mod13)
anova.gam(mod13)


mod14<- gam(nl ~ Hg_mg_kg + gltA + hold_time, data = ms_df,  
            family = tw, method = 'ML')
summary(mod14)
gam.check(mod14)
anova.gam(mod14)


mod15<- gam(nl ~ PC1 + rep_stat + hold_time, data = ms_df,  
            family = tw, method = 'ML')
summary(mod15)
gam.check(mod15)
anova.gam(mod15)


mod16<- gam(nl ~ PC1 + gltA + hold_time, data = ms_df,  
            family = tw, method = 'ML')
summary(mod16)
gam.check(mod16)
anova.gam(mod16)


mod17<- gam(nl ~ PC2 + rep_stat + hold_time, data = ms_df,  
            family = tw, method = 'ML')
summary(mod17)
gam.check(mod17)
anova.gam(mod17)

mod18<- gam(nl ~ PC2 + gltA + hold_time, data = ms_df,  
            family = tw, method = 'ML')
summary(mod18)
gam.check(mod18)
anova.gam(mod18)

mod19<- gam(nl ~ forg_cat + rep_stat + hold_time, data = ms_df,  
            family = tw, method = 'ML')
summary(mod19)
gam.check(mod19)
anova.gam(mod19)

mod20<- gam(nl ~ forg_cat + gltA + hold_time, data = ms_df,  
            family = tw, method = 'ML')
summary(mod20)
gam.check(mod20)
anova.gam(mod20)




AIC_df<- AICc(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9, mod10, 
              mod11, mod12, mod13, mod14, mod15, mod16, mod17, mod18, mod19, mod20)
AIC_df$AIC_delta<- round(AIC_df$AICc - min(AIC_df$AICc), 2)

AIC_df<- AIC_df[order(AIC_df[ ,3]),]
##get AICc weights
AIC_df$weights<- round(Weights(AIC_df$AICc), 3)
AIC_df
#second top model is uninformative
#reporting top model (mod11; deltaAICc = 0), and third model (mod15; deltaAICc = 3)

#make estimated relationship and CIs for micronuclei intensities
newdat_11<- data.frame(mn_intensity = rep(seq(0, 0.0008, by = 0.00002), 3),
                 rep_stat = "non-reproductive",
                 hold_time = mean(ms_df$hold_time))
mod11_pred<- data.frame(newdat_11,
                predict(mod11,newdata = newdat_11, se.fit = T, type= "link"))

mod11_pred$lower<- mod11_pred$fit - (1.96*mod11_pred$se.fit)
mod11_pred$upper<- mod11_pred$fit + (1.96*mod11_pred$se.fit)
mod11_pred$fit<- mod11$fam$linkinv(mod11_pred$fit)
mod11_pred$lower<- mod11$fam$linkinv(mod11_pred$lower)
mod11_pred$upper<- mod11$fam$linkinv(mod11_pred$upper)

colnames(mod11_pred)<- c('mn_intensity', 'rep_stat', 'hold_time', 'nl', 'SE', 'lower', 'upper')

#Box 1 Fig C
MNplot<- ggplot(data = ms_df, aes(x = mn_intensity, y = nl)) +
  geom_jitter(aes(color = species_name), alpha = 0.75, size = 3, color = 'darkolivegreen3')+
  geom_line(data = mod11_pred, aes(x = mn_intensity, y = nl), 
              method = 'glm', method.args = list(family = "tweedie"),
             lwd = 2, color = 'gray25', se = F) +
  geom_ribbon(data = mod11_pred, aes(ymin = lower, ymax = upper), 
              alpha = 0.5, fill = 'gray50', color = NA) +
  th+ theme(legend.position = "top") +
  labs(x = expression(paste("Micronuclei Intensity (micronuclei * 150,200 ", RBCs^{-1}, ')')), y = 'NL Ratio', color = 'Species') 

#get estimated means for repstat
mod11_ems<- as.data.frame(emmeans(mod11, ~rep_stat))
colnames(mod11_ems)<- c('rep_stat', 'nl', 'SE', 'df', 'lower', 'upper')
mod11_ems$nl<- 10^mod11_ems$nl
mod11_ems$lower<- 10^mod11_ems$lower
mod11_ems$upper<- 10^mod11_ems$upper

pairs(emmeans(mod11, ~rep_stat), simple = 'each')

#Box 1 Fig B
repstat_plot<- ggplot(data = ms_df, aes(x = rep_stat, y = nl)) +
  geom_jitter(aes(color = species_name), alpha = 0.75, size = 3, color = 'darkolivegreen3')+
  geom_point(data = mod11_ems, aes(x = rep_stat, y = nl), 
              color = 'gray25', size = 5) +
  geom_errorbar(data = mod11_ems, aes(ymin = nl - lower, ymax = nl + upper),
                  color = 'gray25', lwd = 2, width = 0.1)+
  th+ theme(legend.position = "top") +
  labs(x = 'Reproductive Status', y = 'NL Ratio', color = 'Species') 



#get predictions for 2nd model (3rd top model, but 2nd top was uninformative)
#for mod15
#get relationship estimate and CIs for PC1
newdat_15<- data.frame(PC1 = seq(-10.25, 1.5, l = 125),
                       rep_stat = "non-reproductive",
                       hold_time = mean(ms_df$hold_time))
mod15_pred<- data.frame(newdat_15,
                        predict(mod15, newdata = newdat_15, se.fit = T, type= "link"))

mod15_pred$lower<- mod15_pred$fit - (1.96*mod15_pred$se.fit)
mod15_pred$upper<- mod15_pred$fit + (1.96*mod15_pred$se.fit)
mod15_pred$fit<- mod15$fam$linkinv(mod15_pred$fit)
mod15_pred$lower<- mod15$fam$linkinv(mod15_pred$lower)
mod15_pred$upper<- mod15$fam$linkinv(mod15_pred$upper)

colnames(mod15_pred)<- c('PC1', 'rep_stat', 'hold_time', 'nl', 'SE', 'lower', 'upper')

#Box 1 Fig D
PC1plot<- ggplot(data = ms_df, aes(x = PC1, y = nl)) +
  geom_jitter(aes(color = species_name), alpha = 0.75, size = 3, color = 'darkolivegreen3' )+
  geom_line(data = mod15_pred, aes(x = PC1, y = nl), 
            method = glm, method.args = list(family = "Gamma"), 
            lwd = 2, color = 'gray25', se = F) +
  geom_ribbon(data = mod15_pred, aes(ymin = lower, ymax = upper), 
              alpha = 0.5, fill = 'gray50', color = NA) +
  th+ theme(legend.position = "top") +
  labs(x = "PC1", y = 'NL Ratio', color = 'Species') 


repstat_plot 
MNplot 
PC1plot

