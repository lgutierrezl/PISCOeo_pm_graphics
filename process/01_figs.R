# Tmax : maximum temperature (°C)
# Tmin : minimum temperature (°C)
# Tmean : mean temperature (°C)
# Sd : sunshine duration (hours)
# Rs: solar radiation 
# Td :  dew point temperature (°C)
# Ws : wind speed (m/s)
# Hr: Relative humidity (%)
# Eo : reference evaporation/evapotranspiration (mm/day)
"%>%" = magrittr::`%>%`

library(xts)
library(lattice)
library(raster)
library(ggplot2)
library(geosphere)
library(openair)
library(ggspatial)
library(parallel)
library(trend)
library(ggpubr)
library(cowplot)
source('src/fc2plot.R')
#
## 01: Available data ---------
data_avl <- read.csv('data/available_data.csv')
data_avl$dat <- as.Date(data_avl$dat)

plot_info <- ggplot(data = data_avl, aes(x = dat, y = val, group=var)) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y")+
  geom_line(aes(col=var), lwd=0.5) + 
  scale_y_continuous(breaks = seq(0,400, 30))+
  labs(x='Year', y = 'Number of stations')+
  guides(color = guide_legend(title.position = "top",
                              reverse=F,
                              title = 'Meteorological subvariables',
                              title.hjust = 0.5,
                              ncol = 1))+
  theme_bw()+
  theme(legend.position = c(0.15,0.8),
        axis.title.x = element_blank())

ggsave(plot = plot_info, "img/Figura_S02.pdf", units = "mm", device = 'pdf', width = 200, height = 100, dpi = 900)
#
## 02: Correlation & Elevation & Distance -----
Td <- readRDS('data/raw/qc_td_obs.RDS')
MATRIX_D <- CORR_DIST_MATRIX(data.frame(Td$xyz))
MATRIX_H <- CORR_ALT_MATRIX(data.frame(Td$xyz))
MATRIX_R <- CORREL_FUNC(data.frame(data.frame(Td$values)), data.frame(Td$xyz))
GRAF_td <- COMP_DHR(MATRIX_D, MATRIX_H, MATRIX_R, Li = 0.7, Ls = 0.97)

Sd <- readRDS('data/raw/qc_sd_obs.RDS')
MATRIX_D <- CORR_DIST_MATRIX(data.frame(Sd$xyz))
MATRIX_H <- CORR_ALT_MATRIX(data.frame(Sd$xyz))
MATRIX_R <- CORREL_FUNC(data.frame(data.frame(Sd$values)), data.frame(Sd$xyz))
GRAF_sd <- COMP_DHR(MATRIX_D, MATRIX_H, MATRIX_R, Li = 0.5, Ls = 0.9)

Ws <- readRDS('data/raw/qc_ws_obs.RDS')
MATRIX_D <- CORR_DIST_MATRIX(data.frame(Ws$xyz))
MATRIX_H <- CORR_ALT_MATRIX(data.frame(Ws$xyz))
MATRIX_R <- CORREL_FUNC(data.frame(data.frame(Ws$values)), data.frame(Ws$xyz))
GRAF_ws <- COMP_DHR(MATRIX_D, MATRIX_H, MATRIX_R, Li = 0.05, Ls = 0.4)

df_dhr1 <- rbind(GRAF_td$tab1, GRAF_sd$tab1, GRAF_ws$tab1)
df_dhr1$Var <- c(rep('Td', nrow(GRAF_td$tab1)), rep('Sd', nrow(GRAF_sd$tab1)),
                 rep('Ws', nrow(GRAF_ws$tab1)))
df_dhr1$income <- factor(df_dhr1$income, levels = c('>1000', '1000', '500', '100'),
                         labels = c('> 1,000', '1,000', '500', '100'))

ext_data_tsw <- rbind(data.frame(data_td[,1:2], val = data_td[,3], var = 'Td'), 
                      data.frame(data_sd[,1:2], val = data_sd[,3], var = 'Sd'), 
                      data.frame(data_ws[,1:2], val = data_ws[,3], var = 'Ws'))
names(ext_data_tsw) <- c('alt', 'dis', 'val', 'var')
ext_data_tsw$dis <- as.numeric(as.character(ext_data_tsw$dis))
ext_data_tsw$var <- factor(ext_data_tsw$var, levels = unique(ext_data_tsw$var))
ext_data_tsw$alt <- factor(ext_data_tsw$alt, levels = c('>1000', '1000', '500', '100'),
                           labels = c('> 1,000', '1,000', '500', '100'))

dcor <- ggplot(ext_data_tsw, aes(x = dis, y = val)) + 
  geom_line(aes(color =alt), size = 1)+
  facet_wrap(vars(var), ncol = 3, scales = "fixed")+
  scale_y_continuous(breaks=seq(0, 1, 0.2), limits = c(0, 1))+
  scale_color_manual(values = c("100" = "#13AEE0","500" = "#CB12D7", "1,000" = "#478B04",'> 1,000' = '#D32626')) +
  labs(x = "Distance (km)", y = 'Correlation')+
  guides(color = guide_legend(title.position = "top",
                              reverse=T,title = 'Elevation \n Difference (m)',
                              title.vjust = 0.5))+
  theme_bw() +
  theme(legend.position='none')

df_dhr2 <- rbind(GRAF_td$tab2, GRAF_sd$tab2, GRAF_ws$tab2)
df_dhr2$Var <- c(rep('Td', nrow(GRAF_td$tab2)), rep('Sd', nrow(GRAF_sd$tab2)),
                 rep('Ws', nrow(GRAF_ws$tab2)))
df_dhr2$Var <- factor(df_dhr2$Var, levels = levels(ext_data_tsw$var))
df_dhr2$income <- factor(df_dhr2$income, levels = c('>1000', '1000', '500', '100'),
                         labels = c('> 1,000', '1,000', '500', '100'))

bar <- ggplot(df_dhr2, aes(x=DIST, y=count)) +
  geom_bar(aes(fill = income),stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = c("100" = "#13AEE0","500" = "#CB12D7", "1,000" = "#478B04",'> 1,000' = '#D32626')) +
  labs(x = "Distance (km)", y = 'Number of stations')+
  facet_wrap(vars(Var), ncol = 3, scales = "fixed")+
  guides(fill = guide_legend(title.position = "top",
                             reverse=T,title = 'Elevation \n Difference (m)',
                             title.vjust = 0.5))+
  theme_bw() +
  theme(legend.position=c(0.07,0.6))

plotcomb <- ggpubr::ggarrange(bar, dcor, labels = c("a", "b"),vjust=4,ncol=1, font.label=list(size=10,face="bold"))

ggsave(plot= plotcomb, 'img/Figura_S03.pdf', units="mm", device = 'pdf', width=200, height=140, dpi=900)
#
## 03: Monthly climatology correlation of sub-variables vs spatial co variable-----
dem <- raster('data/raw/spatial/Spatial_co_variables/DEM.nc')
lon <- raster('data/raw/spatial/Spatial_co_variables/X.nc')
lat <- raster('data/raw/spatial/Spatial_co_variables/Y.nc')
cc  <- raster::stack('data/raw/spatial/Spatial_co_variables/CC.nc')
lsm <- raster::stack('data/raw/spatial/Spatial_co_variables/LST_mean.nc')
wsc <- raster::stack('data/raw/spatial/Spatial_co_variables/WS_worldclim.nc')

#Td vs LST_mean, DEM, X, Y
Td <- readRDS('data/processed/TDWP_obs_hmgf.RDS')$xyz
Td_val<- raster::stack('data/processed/td_mean_1981-2010.nc')

Td_Td <- data.frame(raster::extract(Td_val, Td[,2:3]))
Td_lsm <- data.frame(raster::extract(lsm, Td[,2:3]))
Td_dem <- data.frame(raster::extract(dem, Td[,2:3]))
Td_lon <- data.frame(raster::extract(lon, Td[,2:3]))
Td_lat <- data.frame(raster::extract(lat, Td[,2:3]))

Td_eval <- list()
for (i in 1:12) {
  Tdls <- cor(Td_Td[,i], Td_lsm[,i], use="pairwise.complete.obs")
  Tddm <- cor(Td_Td[,i], Td_dem, use="pairwise.complete.obs")
  Tdlo <- cor(Td_Td[,i], Td_lon, use="pairwise.complete.obs")
  Tdla <- cor(Td_Td[,i], Td_lat, use="pairwise.complete.obs")
  Td_eval[[i]] <- data.frame(Tdls, Tddm, Tdlo, Tdla)
}
Td_eval <- do.call('rbind', Td_eval)
Td_eval$mes <- factor(c('Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun','Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic'), 
                      levels = c('Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic'))

Td_eval <- tidyr::pivot_longer(Td_eval, !c('mes'), names_to = 'var', values_to = 'val')
Td_eval$mes <- factor(Td_eval$mes)
Td_eval$var <- factor(Td_eval$var, levels = c("raster..extract.dem..Td...2.3..",
                                              "raster..extract.lat..Td...2.3..",
                                              "raster..extract.lon..Td...2.3..",
                                              "Tdls"),
                      labels = c('DEM','X','Y','LST_mean'))

#Sd vs CC, DEM, X, Y
Td <- readRDS('data/processed/HSOL_obs_hmgf.RDS')$xyz
Td_val<- raster::stack('data/processed/sd_mean_1981-2010.nc')

Td_Td <- data.frame(raster::extract(Td_val, Td[,2:3]))
Td_lsm <- data.frame(raster::extract(cc, Td[,2:3]))
Td_dem <- data.frame(raster::extract(dem, Td[,2:3]))
Td_lon <- data.frame(raster::extract(lon, Td[,2:3]))
Td_lat <- data.frame(raster::extract(lat, Td[,2:3]))

sd_eval <- list()
for (i in 1:12) {
  Tdls <- cor(Td_Td[,i], Td_lsm[,i], use="pairwise.complete.obs")
  Tddm <- cor(Td_Td[,i], Td_dem, use="pairwise.complete.obs")
  Tdlo <- cor(Td_Td[,i], Td_lon, use="pairwise.complete.obs")
  Tdla <- cor(Td_Td[,i], Td_lat, use="pairwise.complete.obs")
  sd_eval[[i]] <- data.frame(Tdls, Tddm, Tdlo, Tdla)
}
sd_eval <- do.call('rbind', sd_eval)

sd_eval$mes <- factor(c('Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic'), 
                      levels = c('Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic'))


sd_eval <- tidyr::pivot_longer(sd_eval, !c('mes'), names_to = 'var', values_to = 'val')
sd_eval$mes <- factor(sd_eval$mes)
sd_eval$var <- factor(sd_eval$var, levels = c("raster..extract.dem..Td...2.3..",
                                              "raster..extract.lat..Td...2.3..",
                                              "raster..extract.lon..Td...2.3..",
                                              "Tdls"),
                      labels = c('DEM', 'X','Y','CC'))

#Ws vs WS_Worldclim, DEM, X, Y
Td <- readRDS('data/processed/qc_ws_obs.RDS')$xyz
Td_val<- raster::stack('data/processed/ws_mean.nc')

Td_Td <- data.frame(raster::extract(Td_val, Td[,2:3]))
Td_lsm <- data.frame(raster::extract(wsc, Td[,2:3]))
Td_dem <- data.frame(raster::extract(dem, Td[,2:3]))
Td_lon <- data.frame(raster::extract(lon, Td[,2:3]))
Td_lat <- data.frame(raster::extract(lat, Td[,2:3]))

Ws_eval <- list()
for (i in 1:12) {
  Tdls <- cor(Td_Td[,i], Td_lsm[,i], use="pairwise.complete.obs")
  Tddm <- cor(Td_Td[,i], Td_dem, use="pairwise.complete.obs")
  Tdlo <- cor(Td_Td[,i], Td_lon, use="pairwise.complete.obs")
  Tdla <- cor(Td_Td[,i], Td_lat, use="pairwise.complete.obs")
  Ws_eval[[i]] <- data.frame(Tdls, Tddm, Tdlo, Tdla)
}
Ws_eval <- do.call('rbind', Ws_eval)

Ws_eval$mes <- factor(c('Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic'), 
                      levels = c('Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic'))

Ws_eval <- tidyr::pivot_longer(Ws_eval, !c('mes'), names_to = 'var', values_to = 'val')
Ws_eval$mes <- factor(Ws_eval$mes)
Ws_eval$var <- factor(Ws_eval$var, levels = c("raster..extract.dem..Td...2.3..",
                                              "raster..extract.lat..Td...2.3..",
                                              "raster..extract.lon..Td...2.3..",
                                              "Tdls"),
                      labels = c('DEM', 'X','Y','WS_worldclim'))

climt <- rbind(Td_eval,sd_eval,Ws_eval)
climt$clim <- c(rep('Td',48), rep('Sd',48), rep('Ws',48))
climt$mes <- factor(climt$mes, levels = c('Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sep', 'Oct', 'Nov', 'Dic'),
                    labels = month.abb)

gf_td <- ggplot(data = climt, aes(x = mes, y = val, group=var)) +
  geom_line(aes(col=var), lwd=0.6) +
  labs(y='Correlation')+
  geom_abline(intercept=0, slope=0, color='gray60', size=0.8, alpha=0.7, linetype = 'dashed')+ 
  facet_wrap(vars(clim), nrow = 1, scales = "fixed")+
  geom_point(aes(col=var))+
  guides(color = guide_legend(title.position = "left",
                              reverse=F,
                              title = 'Covariable',
                              title.hjust = 0.5, nrow = 1))+
  theme_bw()+
  theme(legend.position='bottom',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45))

ggsave(plot = gf_td, "img/Figura_S04.pdf", device = 'pdf', units = "mm", width = 200, height = 100, dpi = 900)
#
## 04: Spatial metrics by stations ----
# daily
piscoeo <- read.csv('data/processed/piscoeo_daily.csv')
pm_eo <- read.csv('data/processed/etp_pm.csv')
Eoc_xyz <- read.csv('data/processed/PISCOeo_pm_xyz_for_cv.csv')

Eoc <- piscoeo[,-1]
Eo_c <- pm_eo[,-1]
Eoc_v <- cbind(Eoc, Eo_c)
names(Eoc_v) <- paste0('X_',1:ncol(Eoc_v))
Eoc_v <- data.frame(Eoc_v, date=seq(as.POSIXct("1981-01-01"),as.POSIXct("2016-12-31"),by ="day"))
Eoc_eval <- list()
for (i in 1:ncol(Eoc)) {
  Eoc_eval[[i]] <- openair::modStats(mydata = Eoc_v, mod=names(Eoc_v)[i], obs=names(Eoc_v)[i+ncol(Eoc)], 
                                     statistic = c('MB','MGE','IOA'))
}
Eoc_eval <- do.call('rbind', Eoc_eval)
Eoc_Eop_map <- fortify(data.frame(COD = Eoc_xyz$ID, 
                                  LON = Eoc_xyz$LON, LAT = Eoc_xyz$LAT, ALT= Eoc_xyz$ALT,
                                  Eoc_eval[,2:4]))

Eoc_Eop_map <- data.frame(tidyr::pivot_longer(Eoc_Eop_map, !c("COD", "LON", "LAT", 'ALT'),
                                              names_to = 'METRIC', values_to = 'VAL'))
peru <- fortify(shapefile('data/raw/spatial/SEC_CLIM.shp'))

map1 <- MAP_METRIC_EP_2(data=Eoc_Eop_map[Eoc_Eop_map$METRIC == 'MB',], lab_id = 'a',
                        PERU = peru, metrics = 'bias', punto_medio = 0,
                        cortes = c(-1,-0.5,0,0.5,1), limites = c(-1,1))
map2 <- MAP_METRIC_EP_2(data=Eoc_Eop_map[Eoc_Eop_map$METRIC == 'MGE',], lab_id = 'b',
                        PERU = peru, metrics = 'MAE', punto_medio = 0.75,
                        cortes = c(0.2,0.6,0.8,1.2,1.5), limites = c(0.2,1.5))
map3 <- MAP_METRIC_EP_3(data=Eoc_Eop_map[Eoc_Eop_map$METRIC == 'IOA',], lab_id = 'c',
                        PERU = peru, metrics = expression(italic(d[r])), punto_medio = 0.75,
                        cortes = c(0.45,0.55,0.65,0.75,0.85,0.95), limites = c(0.4,1))
plotcomb <- ggarrange(map1, map2, map3, ncol = 3, nrow = 1)

ggsave(plot = plotcomb, 'img/Figura_03.pdf', units = "mm", device = 'pdf', width = 200, height = 120, dpi = 600)

#mensual
Eoc <- piscoeo[,-1]
Eo_c <- pm_eo[,-1]
Eoc_v <- cbind(Eoc, Eo_c)
names(Eoc_v) <- paste0('X_',1:ncol(Eoc_v))
Eoc_v <- xts(Eoc_v, order.by =seq(from = as.Date("1981-01-01"),to=as.Date("2016-12-31"),by ="day"))
Eoc_v <- data.frame(apply.monthly(Eoc_v, FUN = apply, MARGIN=2, sum))
Eoc_eval <- list()
for (i in 1:ncol(Eoc)) {
  Eoc_eval[[i]] <- openair::modStats(mydata = Eoc_v, mod=names(Eoc_v)[i], obs=names(Eoc_v)[i+ncol(Eoc)], 
                                     statistic = c('MB','MGE','IOA'))
}
Eoc_eval <- do.call('rbind', Eoc_eval)
Eoc_Eop_map <- fortify(data.frame(COD = Eoc_xyz$ID, 
                                  LON = Eoc_xyz$LON, LAT = Eoc_xyz$LAT, ALT= Eoc_xyz$ALT,
                                  Eoc_eval[,2:4]))

Eoc_Eop_map <- data.frame(tidyr::pivot_longer(Eoc_Eop_map, !c("COD", "LON", "LAT"),
                                              names_to = 'METRIC', values_to = 'VAL'))

map1 <- MAP_METRIC_EP_2(data=Eoc_Eop_map[Eoc_Eop_map$METRIC == 'MB',], lab_id = 'a', 
                        PERU = peru, metrics = 'bias', punto_medio = 0,
                        cortes = c(-20,-10,0,10,20), limites = c(-20,20))
map2 <- MAP_METRIC_EP_2(data=Eoc_Eop_map[Eoc_Eop_map$METRIC == 'MGE',], lab_id = 'b', 
                        PERU = peru, metrics = 'MAE', punto_medio = 15,
                        cortes = c(0,10,20,30,40), limites = c(0,45))
map3 <- MAP_METRIC_EP_3(data=Eoc_Eop_map[Eoc_Eop_map$METRIC == 'IOA',], lab_id = 'c',
                        PERU = peru, metrics = expression(italic(d[r])), punto_medio = 0.75,
                        cortes = c(0.45,0.55,0.65,0.75,0.85,0.95), limites = c(0.4,1))
plotcomb <- ggarrange(map1, map2, map3, ncol = 3, nrow = 1)

ggsave(plot = plotcomb, 'img/Figure_S06.pdf', units="mm", device='pdf', width = 200, height = 120, dpi = 600)

#anual
Eoc <- piscoeo[,-1]
Eo_c <- pm_eo[,-1]
Eoc_v <- cbind(Eoc, Eo_c)
names(Eoc_v) <- paste0('X_',1:ncol(Eoc_v))
Eoc_v <- xts(Eoc_v, order.by =seq(from = as.Date("1981-01-01"),to=as.Date("2016-12-31"),by ="day"))
Eoc_v <- data.frame(apply.yearly(Eoc_v, FUN = apply, MARGIN=2, sum))
Eoc_eval <- list()
for (i in 1:ncol(Eoc)) {
  Eoc_eval[[i]] <- openair::modStats(mydata = Eoc_v, mod=names(Eoc_v)[i], obs=names(Eoc_v)[i+ncol(Eoc)], 
                                     statistic = c('MB','MGE','IOA'))
}
Eoc_eval <- do.call('rbind', Eoc_eval)
Eoc_Eop_map <- fortify(data.frame(COD = Eoc_xyz$ID, 
                                  LON = Eoc_xyz$LON, LAT = Eoc_xyz$LAT, ALT= Eoc_xyz$ALT,
                                  Eoc_eval[,2:4]))

Eoc_Eop_map <- data.frame(tidyr::pivot_longer(Eoc_Eop_map, !c("COD", "LON", "LAT"),
                                              names_to = 'METRIC', values_to = 'VAL'))

map1 <- MAP_METRIC_EP_2(data=Eoc_Eop_map[Eoc_Eop_map$METRIC == 'MB',], lab_id = 'a', 
                        PERU = peru, metrics = 'bias', punto_medio = 0,
                        cortes = c(-300,-150,0,150,300), limites = c(-300,300))
map2 <- MAP_METRIC_EP_2(data=Eoc_Eop_map[Eoc_Eop_map$METRIC == 'MGE',], lab_id = 'b', 
                        PERU = peru, metrics = 'MAE', punto_medio = 100,
                        cortes = c(0,70,100,130,160,190), limites = c(0,190))
map3 <- MAP_METRIC_EP_3(data=Eoc_Eop_map[Eoc_Eop_map$METRIC == 'IOA',], lab_id = 'c',
                        PERU = peru, metrics = 'dr', punto_medio = 0.75,
                        cortes = c(0.45,0.55,0.65,0.75,0.85,0.95), limites = c(0.4,1))

plotcomb <- ggarrange(map1, map2, map3, ncol = 3, nrow = 1)

ggsave(plot = plotcomb, 'img/Figure_S07.pdf', units="mm", device ='pdf', width = 200, height = 120, dpi = 600)
#
## 05: Spatial variability of Eo ------
Eo_raster <- readRDS('data/processed/Eo_yr_raster.RDS')
peru <- shapefile("data/raw/spatial/SEC_CLIM.shp")
Eo_raster[[1]] <- mask(sum(raster::brick('data/processed/eo_mean_1981-2010.nc')), peru)
map_eop <- map_eo_prod(Eo_raster)

ggsave(plot=map_eop, "img/Figura_04.pdf", units="mm", device = 'pdf', width=200, height=160, dpi = 600)
#
## 06: Monthly climatology by climate sector --------
Eo_cr_psc <- read.csv('data/processed/PISCOeo_pm_normal_by_mcoregion.csv')
Eo_cr_trc <- read.csv('data/processed/TERRACLIMATE_Eo_Peru_1981_2010_monthly_sum_mean_mcoregion.csv')
Eo_cr_cru <- read.csv('data/processed/CRU_Eo_Peru_1981_2010_monthly_sum_mean_mcoregion.csv')
Eo_cr_era <- read.csv('data/processed/ERA5land_Eo_Peru_1981_2010_monthly_sum_mean_mcoregion.csv')

Eo_cr <- rbind(Eo_cr_psc, Eo_cr_trc, Eo_cr_cru, Eo_cr_era)
Eo_cr$SRC <- rep(c('PISCOeo_pm', 'TerraClimate', 'CRU_TS', 'ERA5-Land'), each=12)

Eo_clm_r <- tidyr::pivot_longer(Eo_cr, !c('month', 'SRC'), names_to = 'REG', values_to = 'VAL')

Eo_clm_r$SRC <- factor(Eo_clm_r$SRC, levels = c('PISCOeo_pm','CRU_TS','TerraClimate','ERA5-Land'))

Eo_clm_r$REG <- factor(Eo_clm_r$REG, levels = c("SEA","SEB","SIOR","SIOC","CO"),
                       labels = c('High Amazon','Low Amazon', 'Eastern Andes', 'Western Andes', 'Pacific Coast'))

Eo_clm_r$month <- factor(Eo_clm_r$month, levels = 1:12, labels = month.abb)

sbs_clm <- ggplot(Eo_clm_r, aes(x=month, y=VAL, group = SRC, fill= REG))+
  geom_line(aes(color=SRC), size=0.6)+
  labs(x = "Meses", y = 'ETo (mm/month)')+
  facet_wrap(vars(REG), nrow = 3, scales = "free")+
  scale_y_continuous(breaks = seq(0,200,15))+
  guides(col = guide_legend(title.position = "top", title = 'Product', title.vjust = 0.5))+
  theme_bw()+
  theme(legend.position=c(0.75,0.1),
        axis.title.x = element_blank())

ggsave(plot = sbs_clm, 'img/Figura_05.pdf', units = "mm", device = 'pdf', width = 200, height = 130, dpi = 600)
#
## 07: Yearly and monthly Eo by climate sector -----
Eo_cr_psc <- read.csv('data/processed/PISCOeo_pm_yearly_by_mcoregion.csv')
Eo_cr_trc <- read.csv('data/processed/TERRACLIMATE_Eo_Peru_1981_2019_yearly_sum_ts_mcoregion.csv')
Eo_cr_cru <- read.csv('data/processed/CRU_Eo_Peru_1981_2019_yearly_sum_ts_mcoregion.csv')
Eo_cr_era <- read.csv('data/processed/ERA5land_Eo_Peru_1981_2019_yearly_sum_ts_mcoregion.csv')

Eo_psc_ts <- TREND_TIMES(datast = Eo_cr_psc[,2:6])
Eo_trc_ts <- TREND_TIMES(datast = Eo_cr_trc[,2:6])
Eo_cru_ts <- TREND_TIMES(datast = Eo_cr_cru[,2:6])
Eo_era_ts <- TREND_TIMES(datast = Eo_cr_era[,2:6])

Eo_cr <- rbind(Eo_cr_psc, Eo_cr_trc, Eo_cr_cru, Eo_cr_era)
Eo_cr$SRC <- rep(c('PISCOeo_pm', 'TerraClimate', 'CRU_TS', 'ERA5-Land'), each=36)

Eo_clm_r <- tidyr::pivot_longer(Eo_cr, !c('time', 'SRC'), names_to = 'REG', values_to = 'VAL')
Eo_clm_r$SRC <- factor(Eo_clm_r$SRC, levels = c('PISCOeo_pm','CRU_TS','TerraClimate','ERA5-Land'))

Eo_clm_r$REG <- factor(Eo_clm_r$REG, levels = c("SEA","SEB","SIOR","SIOC","CO"),
                       labels = c('High Amazon','Low Amazon', 'Eastern Andes', 'Western Andes', 'Pacific Coast'))

Eo_clm_r$SRC_REG <- paste0(Eo_clm_r$SRC, '_', Eo_clm_r$REG)

Eo_slp_ts <- rbind(Eo_psc_ts, Eo_trc_ts, Eo_cru_ts, Eo_era_ts)
Eo_slp_ts$SRC_REG <- unique(Eo_clm_r$SRC_REG)

Eo_year_ts <- merge(Eo_clm_r, Eo_slp_ts, by='SRC_REG')
Eo_year_ts$time <- as.Date(Eo_year_ts$time)
Eo_year_ts$NS <- Eo_year_ts$MK_paa
Eo_year_ts$NS[Eo_year_ts$NS > 1] <- '  trend >= 95%'
Eo_year_ts$NS[Eo_year_ts$NS == '0'] <- '  trend < 95%'
Eo_year_ts$NS[Eo_year_ts$NS == '1'] <- '  trend < 95%'
Eo_year_ts_na <- Eo_year_ts
Eo_year_ts_na[Eo_year_ts_na$NS=='  trend < 95%',] <- NA
Eo_year_ts_na <- na.omit(Eo_year_ts_na)
Eo_year_ts$SLP <- round(Eo_year_ts$sslp_paa,3)
Eo_year_ts$SLP <- as.factor(Eo_year_ts$SLP)

Eo_year_text <- Eo_year_ts[,c(2:5, 9, 10)]
Eo_year_text <- Eo_year_text[Eo_year_text$time =="2016-12-31",]
Eo_year_text$time <- as.Date('2023-12-31')
row.names(Eo_year_text) <- NULL
Eo_year_text <- dplyr::arrange(Eo_year_text, REG)
Eo_year_text$NSLP <- paste0(substr(Eo_year_text$NS,1,1), 
                            Eo_year_text$SLP, substr(Eo_year_text$NS,2,2))
Eo_year_text[Eo_year_text$NS=='  trend < 95%',] <- NA
Eo_year_text <- na.omit(Eo_year_text)

names(Eo_year_ts)[9] <- 'Confidence interval'

sbs_clm <- ggplot(Eo_year_ts, aes(x=time, y=VAL, group = SRC))+
  geom_line(aes(color=SRC), size=0.5, alpha =0.6, show.legend = T)+ 
  geom_smooth(data=Eo_year_ts_na, aes(color=SRC), linetype= 'dashed',method = 'lm', 
              formula = 'y~x',se=FALSE, size=0.5, show.legend = F)+
  facet_wrap(vars(REG), nrow = 3, scales = "free")+
  geom_text(data = Eo_year_text, aes(x = time, y = VAL, label =  NSLP, color=SRC), show.legend = F)+
  labs(x = "Meses", y = 'ETo (mm/year)')+
  guides(col = guide_legend(title.position = "top", title = 'Product', title.vjust = 0.5, title.hjust = 0.5, ncol = 2))+
  theme_bw()+
  theme(legend.position=c(0.75,0.15),
        axis.title.x = element_blank())

ggsave(plot = sbs_clm, 'img/Figura_06.pdf', units = "mm", device = 'pdf', width = 200, height = 140, dpi = 900)

## 08: Delta PISCOeo_pm ----------
eo_month <- read.csv('data/raw/spatial/Deo_monthly_mean.csv')
eo_month <- tidyr::pivot_longer(eo_month, !c('X'),names_to = 'sec', values_to = 'eo')

eo_month$sec <- factor(eo_month$sec, levels = c("SEA","SEB","SIOR","SIOC","CO"),
                       labels = c('High Amazon','Low Amazon', 'Eastern Andes', 'Western Andes', 'Pacific Coast'))

eo_month$month <- factor(eo_month$X, levels = 1:12, labels = month.abb)

sbs_clm <- ggplot(eo_month, aes(x=month, y=eo, group = sec))+
  geom_line(aes(color=sec), size=0.75)+
  labs(x = "Month", y = 'ETo (mm/day)')+
  annotate(geom = "text", x = 1, y = 2.15, label = 'b', hjust = "left", size = 4, fontface =2)+
  guides(col = guide_legend(title.position = "left", title = 'Climate \nRegion', title.vjust = 0.5, nrow = 3))+
  theme_bw()+
  theme(legend.position='bottom',
        axis.title.x = element_blank())

#multianual 
eo_month <- raster('data/raw/spatial/Deo_annual_mean_clim.nc')
SCLIM     <- shapefile("data/raw/spatial/SEC_CLIM.shp")
eo_month <- mask(crop(eo_month, SCLIM), SCLIM)
PERU     <- fortify(SCLIM)
DEM <- na.omit(as.data.frame(eo_month, xy=T))

mapa <- ggplot() +
  geom_raster(data = DEM, aes(x=x, y=y, fill= DEM[,3]))+
  geom_polygon(data=PERU, aes(x=long, y=lat, group=group), fill=NA, color="gray40",size = 0.3)+
  annotate(geom = "text", x = -81, y = 0.4, label = 'a', hjust = "left", size = 4, fontface =2)+
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(10, "Spectral")),
                       breaks = seq(0,3,0.25), guide = 'legend', name='delta ETo (mm/day)')+
  coord_equal()+
  coord_sf(crs=4326, xlim = c(-81,-69), ylim = c(-18,-0))+
  theme_bw()+
  theme(legend.position = 'right')

plotcomb <- ggdraw() +
  draw_plot(mapa, x = 0, y = 0, width = 0.45, height = 1) +
  draw_plot(sbs_clm, x = 0.45, y = 0, width = 0.55, height = 1)

ggsave(plot = plotcomb, 'img/Figura_07.pdf', units = "mm", device = 'pdf', width = 200, height = 110, dpi = 600)
#
## 09: Correlation | Number stations vs distance vs elevation difference -----------
Td <- readRDS('data/raw/qc_td_obs.RDS')
Sd <- readRDS('data/raw/qc_sd_obs.RDS')
Ws <- readRDS('data/raw/qc_ws_obs.RDS')

td <- Td$values
sd <- Sd$values
ws <- Ws$values

nam_est <- unique(c(names(td), names(sd), names(ws)))
xyz_est <- rbind(data.frame(Td$xyz, SRC = 'td'), data.frame(Sd$xyz, SRC = 'sd'), data.frame(Ws$xyz, SRC = 'ws'))
xyz_est <- xyz_est[!duplicated(xyz_est$ID), ]
st_na <- xts(rep(NA, nrow(td)), order.by = index(td))

qc_vl <- list()
qc_xy <- list()
for (i in seq_along(nam_est)) {
  if(!is.na(match(nam_est[i], names(td)))){td_sel = td[,nam_est[i]]}else{td_sel=st_na}
  if(!is.na(match(nam_est[i], names(sd)))){sd_sel = sd[,nam_est[i]]}else{sd_sel=st_na}
  if(!is.na(match(nam_est[i], names(ws)))){ws_sel = ws[,nam_est[i]]}else{ws_sel=st_na}
  st_all <- merge(td_sel, sd_sel, ws_sel)
  names(st_all) <- c('td', 'sd', 'ws')
  qc_vl[[nam_est[i]]] <- st_all
  qc_xy[[i]] <- xyz_est[xyz_est$ID == nam_est[i],]
}
qc_xy <- do.call(rbind, qc_xy)
qc_xy <- data.frame(ID = qc_xy$ID, NAM=qc_xy$ID, qc_xy[,c(2:4, 6)])

qc01 <- list(values = qc_vl, xyz = qc_xy)

# number of stations vs distance vs elevation difference
no_cores <- detectCores(logical = TRUE)-1
cl <- makeCluster(no_cores-2)
registerDoParallel(cl)

clusterExport(cl,list('qc01', '%>%', 'td', 'sd', 'ws'))

parallel::parLapply(cl, 1:length(qc01$xyz$ID),
                    function(x){
                      x_can <- qc01$xyz[x, ]
                      x_nei <- qc01$xyz[-x,]
                      
                      x_res <- geosphere::distVincentyEllipsoid(x_can[, c("LON", "LAT")],
                                                                x_nei[, c("LON", "LAT")])
                      x_nei[, "DIS"]  <- x_res/1000
                      x_nei[, "ALT_diff"] <- abs(x_can[, "ALT"] - x_nei[, "ALT"])
                      
                      seq_dist <- c(5, 25, 50, 70, 100)
                      seq_elv <- c(100, 500, 1000, 5500)
                      
                      sapply(seq_dist, function(z){
                        sapply(seq_elv, function(y){
                          
                          c(dim(x_nei[x_nei$DIS < z & x_nei$ALT_diff <= y, ])[1])
                          
                        })
                      }) %>% as.data.frame.table -> response
                      
                      response$Var1 <- factor(response$Var1, levels = c("A","B","C","D"), labels = c("100", "500", "1000", ">1000"))
                      response$Var2 <- factor(response$Var2, levels = c("A","B","C","D", "E"), labels = c("5", "25", "50", "70", "100"))
                      
                      response
                      
                    }) %>%
  do.call(rbind, .) -> response_nei

bw_theme <- trellis.par.get()
bw_theme$box.rectangle$col <- "black"
bw_theme$box.umbrella$col <- "black"
bw_theme$plot.symbol$col <- "grey80"

bwplot(Freq~Var2, group = Var1, data = response_nei,
       panel = panel.superpose,
       box.width = 1/6,
       par.settings = bw_theme,
       panel.groups = function(x, y,..., group.number)
         panel.bwplot(x + (group.number-2.5)/6, y, ...),
       auto.key = list(corner = c(0, 1),
                       title = "Elevation\ndifference (m)",
                       cex.title = 1,
                       cex = .8,
                       rectangles = TRUE,
                       points = FALSE),
       xlab = "Distance (km)",
       ylab = "Number of stations") -> p3

# Correlation vs distance vs Elevation difference
td_xyz <- qc01$xyz[qc01$xyz$SRC=='td',]

parallel::parLapply(cl, 1:length(td_xyz$ID),
                    function(x){
                      td_xyz <- qc01$xyz[qc01$xyz$SRC=='td',]
                      x_can <- td_xyz[x, ]
                      x_nei <- td_xyz[-x,]
                      
                      x_res <- geosphere::distVincentyEllipsoid(x_can[, c("LON", "LAT")],
                                                                x_nei[, c("LON", "LAT")])
                      x_nei[, "DIS"]  <- x_res/1000
                      x_nei[, "ALT_diff"] <- abs(x_can[, "ALT"] - x_nei[, "ALT"])
                      
                      seq_dist <- seq(0, 150, 15)
                      seq_dist[1] <- 1
                      seq_elv <- c(100, 500, 1000, 5500)
                      
                      sapply(seq_dist, function(z){
                        sapply(seq_elv, function(y){
                          
                          sapply(x_nei[x_nei$DIS < z & x_nei$ALT_diff <= y, ]$ID,
                                 function(j) cor(td[,td_xyz$ID[x]], td[,j], use = "pairwise.complete.obs")) %>%
                            mean(na.rm = TRUE)
                        })
                      }) %>% as.data.frame.table -> response0
                      
                      sapply(seq_dist, function(z){
                        sapply(seq_elv, function(y){
                          
                          x_nei[x_nei$DIS < z & x_nei$ALT_diff <= y, ]$ID %>%
                            length()
                        })
                      }) %>% as.data.frame.table -> response3
                      
                      response <- cbind(response0, response3[, 3])
                      colnames(response) <- c("Var1", "Var2","Td", "Dist")
                      
                      response$Var1 <- factor(response$Var1, labels = c("100", "500", "1000", ">1000"))
                      response$Var2 <- factor(response$Var2, labels = seq_dist)
                      
                      response
                      
                    }) %>%
  do.call(rbind, .) -> response_cor_td

td_xyz <- qc01$xyz[qc01$xyz$SRC=='sd',]

parallel::parLapply(cl, 1:length(td_xyz$ID),
                    function(x){
                      td_xyz <- qc01$xyz[qc01$xyz$SRC=='sd',]
                      x_can <- td_xyz[x, ]
                      x_nei <- td_xyz[-x,]
                      
                      x_res <- geosphere::distVincentyEllipsoid(x_can[, c("LON", "LAT")],
                                                                x_nei[, c("LON", "LAT")])
                      x_nei[, "DIS"]  <- x_res/1000
                      x_nei[, "ALT_diff"] <- abs(x_can[, "ALT"] - x_nei[, "ALT"])
                      
                      seq_dist <- seq(0, 150, 15)
                      seq_dist[1] <- 1
                      seq_elv <- c(100, 500, 1000, 5500)
                      
                      sapply(seq_dist, function(z){
                        sapply(seq_elv, function(y){
                          
                          sapply(x_nei[x_nei$DIS < z & x_nei$ALT_diff <= y, ]$ID,
                                 function(j) cor(sd[,td_xyz$ID[x]], sd[,j], use = "pairwise.complete.obs")) %>%
                            mean(na.rm = TRUE)
                        })
                      }) %>% as.data.frame.table -> response0
                      
                      sapply(seq_dist, function(z){
                        sapply(seq_elv, function(y){
                          
                          x_nei[x_nei$DIS < z & x_nei$ALT_diff <= y, ]$ID %>%
                            length()
                        })
                      }) %>% as.data.frame.table -> response3
                      
                      response <- cbind(response0, response3[, 3])
                      colnames(response) <- c("Var1", "Var2","Sd", "Dist")
                      
                      response$Var1 <- factor(response$Var1, labels = c("100", "500", "1000", ">1000"))
                      response$Var2 <- factor(response$Var2, labels = seq_dist)
                      
                      response
                      
                    }) %>%
  do.call(rbind, .) -> response_cor_sd

td_xyz <- qc01$xyz[qc01$xyz$SRC=='ws',]

parallel::parLapply(cl, 1:length(td_xyz$ID),
                    function(x){
                      td_xyz <- qc01$xyz[qc01$xyz$SRC=='ws',]
                      x_can <- td_xyz[x, ]
                      x_nei <- td_xyz[-x,]
                      
                      x_res <- geosphere::distVincentyEllipsoid(x_can[, c("LON", "LAT")],
                                                                x_nei[, c("LON", "LAT")])
                      x_nei[, "DIS"]  <- x_res/1000
                      x_nei[, "ALT_diff"] <- abs(x_can[, "ALT"] - x_nei[, "ALT"])
                      
                      seq_dist <- seq(0, 150, 15)
                      seq_dist[1] <- 1
                      seq_elv <- c(100, 500, 1000, 5500)
                      
                      sapply(seq_dist, function(z){
                        sapply(seq_elv, function(y){
                          
                          sapply(x_nei[x_nei$DIS < z & x_nei$ALT_diff <= y, ]$ID,
                                 function(j) cor(ws[,td_xyz$ID[x]], ws[,j], use = "pairwise.complete.obs")) %>%
                            mean(na.rm = TRUE)
                        })
                      }) %>% as.data.frame.table -> response0
                      
                      sapply(seq_dist, function(z){
                        sapply(seq_elv, function(y){
                          
                          x_nei[x_nei$DIS < z & x_nei$ALT_diff <= y, ]$ID %>%
                            length()
                        })
                      }) %>% as.data.frame.table -> response3
                      
                      response <- cbind(response0, response3[, 3])
                      colnames(response) <- c("Var1", "Var2","Ws", "Dist")
                      
                      response$Var1 <- factor(response$Var1, labels = c("100", "500", "1000", ">1000"))
                      response$Var2 <- factor(response$Var2, labels = seq_dist)
                      
                      response
                      
                    }) %>%
  do.call(rbind, .) -> response_cor_ws


xyplot(value  ~ Var2 | variable, groups = Var1, type = c("l"), lwd = 3,
       auto.key = list(lines = TRUE, points = FALSE,
                       corner = c(.975, .075),
                       title = "Elevation\ndifference (m)",
                       cex.title = 1,
                       cex = .8),
       par.settings = list(superpose.line = list(lwd = 3)),
       xlab = "Distance (km)", ylab = "Pearson Correlation", ylim = c(0.75, 0.85),
       data_td = aggregate(response_cor_td[c("Td")], 
                           by = response_cor_td[c("Var1", "Var2")],
                           FUN = function(x) median(x, na.rm = TRUE)) %>%
         reshape2::melt()) -> p0t

xyplot(value  ~ Var2 | variable, groups = Var1, type = c("l"), lwd = 3,
       auto.key = list(lines = TRUE, points = FALSE,
                       corner = c(.975, .075),
                       title = "Elevation\ndifference (m)",
                       cex.title = 1,
                       cex = .8),
       par.settings = list(superpose.line = list(lwd = 3)),
       xlab = "Distance (km)", ylab = "Pearson Correlation", ylim = c(0.65, 0.9),
       data_sd = aggregate(response_cor_sd[c("Sd")], 
                           by = response_cor_sd[c("Var1", "Var2")],
                           FUN = function(x) median(x, na.rm = TRUE)) %>%
         reshape2::melt()) -> p0s

xyplot(value  ~ Var2 | variable, groups = Var1, type = c("l"), lwd = 3,
       auto.key = list(lines = TRUE, points = FALSE,
                       corner = c(.975, .075),
                       title = "Elevation\ndifference (m)",
                       cex.title = 1,
                       cex = .8),
       par.settings = list(superpose.line = list(lwd = 3)),
       xlab = "Distance (km)", ylab = "Pearson Correlation", ylim = c(0, 0.5),
       data_ws = aggregate(response_cor_ws[c("Ws")], 
                           by = response_cor_ws[c("Var1", "Var2")],
                           FUN = function(x) median(x, na.rm = TRUE)) %>%
         reshape2::melt()) -> p0w

update(c(p0w, p0s, p0t), layout = c(1, 3)) -> p0

png(filename= file.path('img', "Figure_S03.png"), width = 10, height = 5.75, units = "in", res = 600)
print(cowplot::plot_grid(p3, p0,  labels = c("a)", "b)")))
dev.off()
#
## 10: Evaluation obs - PISCOeo_pm -----------
peru <- fortify(shapefile('data/raw/spatial/SEC_CLIM.shp'))
Eoc_xyz <- read.csv('data/processed/PISCOeo_pm_xyz_for_cv.csv')

#td
td_mod <- readRDS('data/processed/qc_gf_td_obs_mod.RDS')
td_mod_val <- td_mod$values

td_obs <- readRDS('data/processed/qc_td_obs.RDS')
td_obs_val <- td_obs$values
td_obs_val <- td_obs_val[, match(names(td_mod_val), names(td_obs_val))]

est_sel <- readRDS('data/processed/Normals_OBS_td.RDS')$values[[1]]
td_obs_val <- td_obs_val[, match(colnames(est_sel), colnames(td_obs_val))]
td_mod_val <- td_mod_val[, match(colnames(est_sel), colnames(td_mod_val))]

Eoc_v <- cbind(td_obs_val, td_mod_val)
names(Eoc_v) <- paste0('X_',1:ncol(Eoc_v))
Eoc_v <- data.frame(Eoc_v, date=seq(as.POSIXct("1981-01-01"),as.POSIXct("2016-12-31"),by ="day"))

Eoc_eval <- list()
for (i in seq_along(names(td_obs_val))) {
  Eoc_eval[[i]] <- openair::modStats(mydata = Eoc_v, obs=names(Eoc_v)[i], mod=names(Eoc_v)[i+ncol(td_obs_val)], 
                                     statistic = c('MB','MGE','IOA'))
}
Eoc_eval <- do.call('rbind', Eoc_eval)
Eoc_eval$ID <- names(td_obs_val)
Eoc_Eop_map <- plyr::join(data.frame(td_obs$xyz)[,1:4], Eoc_eval[,2:5], by = 'ID')
Eoc_Eop_map <- na.omit(Eoc_Eop_map)
row.names(Eoc_Eop_map) <- NULL

td_df2map <- data.frame(tidyr::pivot_longer(Eoc_Eop_map, !c("ID", "LON", "LAT", 'ALT'), names_to = 'METRIC', values_to = 'VAL'))
#sd
td_mod <- readRDS('data/processed/qc_gf_sd_model.RDS')
td_mod_val <- td_mod$values

td_obs <- readRDS('data/processed/qc_sd_obs.RDS')
td_obs_val <- td_obs$values
td_obs_val <- td_obs_val[, match(names(td_mod_val), names(td_obs_val))]

est_sel <- readRDS('data/processed/OBS_sd.RDS')$values$sd
td_obs_val <- td_obs_val[, match(colnames(est_sel), colnames(td_obs_val))]
td_mod_val <- td_mod_val[, match(colnames(est_sel), colnames(td_mod_val))]

Eoc_v <- cbind(td_obs_val, td_mod_val)
names(Eoc_v) <- paste0('X_',1:ncol(Eoc_v))
Eoc_v <- data.frame(Eoc_v, date=seq(as.POSIXct("1981-01-01"),as.POSIXct("2019-12-31"),by ="day"))

Eoc_eval <- list()
for (i in seq_along(names(td_obs_val))) {
  Eoc_eval[[i]] <- openair::modStats(mydata = Eoc_v, obs=names(Eoc_v)[i], mod=names(Eoc_v)[i+ncol(td_obs_val)], 
                                     statistic = c('MB','MGE','IOA'))
}
Eoc_eval <- do.call('rbind', Eoc_eval)
Eoc_eval$ID <- names(td_obs_val)
Eoc_Eop_map <- plyr::join(data.frame(td_obs$xyz)[,1:4], Eoc_eval[,2:5], by = 'ID')
Eoc_Eop_map <- na.omit(Eoc_Eop_map)
row.names(Eoc_Eop_map) <- NULL

sd_df2map <- data.frame(tidyr::pivot_longer(Eoc_Eop_map, !c("ID", "LON", "LAT", 'ALT'), names_to = 'METRIC', values_to = 'VAL'))

sdtd_mts <- rbind(data.frame(td_df2map, var = 'Td'), data.frame(sd_df2map, var = 'Sd'))
sdtd_mts <- sdtd_mts[sdtd_mts$METRIC == 'IOA',]

map_td <- MAP_METRIC_EP_3(data=sdtd_mts, lab_id = c('Sd', 'Td'),
                          PERU = peru, metrics = 'dr', punto_medio = 0.75,
                          cortes = c(0.45,0.55,0.65,0.75,0.85,0.95), limites = c(0.4,1))

ggsave(plot = map_td, 'img/Figure_S08.pdf', units="mm", device='pdf', width = 130, height = 120, dpi = 600)
