CORR_DIST_MATRIX <- function(METAD){
  MATRIX <- data.frame()
  for (i in 1:nrow(METAD)) {
    MATY <- data.frame(1)
    for (y in 1:nrow(METAD)) {
      MAT <- geosphere::distm(METAD[i,2:3], METAD[y,2:3], fun = distHaversine)/1000
      MATY[,y] <- MAT
    }
    MATRIX <- rbind(MATRIX, MATY)
  }
  MATRIX$EST <- 1:nrow(MATRIX)
  colnames(MATRIX) <- c(paste0('EST_', 1:nrow(MATRIX)), 'EST')
  return(MATRIX)
}

CORR_ALT_MATRIX <- function(METAD){
  MATRIX <- data.frame()
  for (i in 1:nrow(METAD)) {
    MATY <- data.frame(1)
    for (y in 1:nrow(METAD)) {
      MAT <- abs(METAD[i,4]-METAD[y,4])
      MATY[,y] <- MAT
    }
    MATRIX <- rbind(MATRIX, MATY)
  }
  MATRIX$EST <- 1:nrow(MATRIX)
  colnames(MATRIX) <- c(paste0('EST_', 1:nrow(MATRIX)), 'EST')
  return(MATRIX)
}

CORREL_FUNC <- function(TEMP, METAD){
  yrs <- as.data.frame(t(matrix(1:ncol(TEMP), ncol = 1)))
  for (w in 1:nrow(METAD)) {
    for (e in 1:nrow(METAD)) {
      vect_ob <- cor(TEMP[,w], TEMP[,e], use="pairwise.complete.obs")
      yrs[w,e] <- vect_ob
    }
  }
  yrs$EST <- 1:nrow(METAD)
  colnames(yrs) <- colnames(TEMP)
  return(yrs)
}

COMP_DHR <- function(MATRIX_D, MATRIX_H, MATRIX_R, Li=0.5, Ls=0.85){
  cor_wind_df <- MATRIX_D[,1:(ncol(MATRIX_D)-1)]
  correl_dst <- data.frame(row=rownames(cor_wind_df)[row(cor_wind_df)[upper.tri(cor_wind_df)]], 
                           col=colnames(cor_wind_df)[col(cor_wind_df)[upper.tri(cor_wind_df)]], 
                           DIST=cor_wind_df[upper.tri(cor_wind_df)])
  cor_wind_df <- MATRIX_H[,1:(ncol(MATRIX_D)-1)]
  correl_alt <- data.frame(row=rownames(cor_wind_df)[row(cor_wind_df)[upper.tri(cor_wind_df)]], 
                           col=colnames(cor_wind_df)[col(cor_wind_df)[upper.tri(cor_wind_df)]], 
                           HALT=cor_wind_df[upper.tri(cor_wind_df)])
  cor_wind_df <- MATRIX_R[,1:(ncol(MATRIX_D)-1)]
  correl_cor <- data.frame(row=rownames(cor_wind_df)[row(cor_wind_df)[upper.tri(cor_wind_df)]], 
                           col=colnames(cor_wind_df)[col(cor_wind_df)[upper.tri(cor_wind_df)]], 
                           COR=cor_wind_df[upper.tri(cor_wind_df)])
  COMPR <- data.frame(correl_dst, HLT=correl_alt[,3], COR=correl_cor[,3])
  breaks_h <- c(0,100, 500, 1000, 6000)
  COMPR$HLF <- cut(COMPR[,4],breaks = breaks_h, right = FALSE,
                   labels=c('100','500','1000','>1000'))
  breaks_d <- c(1,seq(15, 150, 15),2500)
  COMPR$DTF <- cut(COMPR[,3],breaks = breaks_d, right = FALSE,
                   labels=c(1,seq(15, 150, 15)))
  breaks_d2 <- c(0, 5,25, 50, 70, 100)
  COMPR$DTF2 <- cut(COMPR[,3],breaks = breaks_d2, right = FALSE,
                    labels=c(5, 25, 50, 70, 100))

  COMPR[,6] <- as.factor(COMPR[,6])
  COMPR[,7] <- as.factor(COMPR[,7])
  COMPR[,8] <- as.factor(COMPR[,8])
  #matrix comp
  matrix_tp <- tapply(COMPR$COR, list(COMPR$HLF,COMPR$DTF), mean,na.rm = TRUE)
  matrix_tp <- as.data.frame(t(matrix_tp))
  matrix_tp$DIST <- row.names(matrix_tp)
  #df dist-cor-alt
  df_dhr1 <- tidyr::pivot_longer(matrix_tp, !DIST, names_to = "income", values_to = "count")
  df_dhr1$income <- factor(df_dhr1$income, levels = c(">1000","1000","500","100"))
  df_dhr1$DIST <- as.numeric(df_dhr1$DIST)
  
  #matrix comp
  COMPR$EST <- 1
  matrix_tp <- tapply(COMPR$EST, list(COMPR$HLF,COMPR$DTF2), sum, na.rm = TRUE)
  matrix_tp <- as.data.frame(t(matrix_tp))
  matrix_tp$DIST <- row.names(matrix_tp)
  #df dist-cor-alt
  df_dhr2 <- tidyr::pivot_longer(matrix_tp, !DIST, names_to = "income", values_to = "count")
  df_dhr2$income <- factor(df_dhr2$income, levels = c(">1000","1000","500","100"))
  df_dhr2$DIST <- factor(df_dhr2$DIST, levels = c("5","25","50","70",'100'))
  
  list_dhr <- list(tab1 = df_dhr1, tab2 = df_dhr2)
  return(list_dhr)
}

MAP_METRIC_EP_2 <- function(data = Eop_Eoc, PERU=peru, metrics = METRIC, lab_id = 'a',
                            punto_medio = 0, 
                            cortes = vecortes, limites= c(-1,1)){
  data$VAL[data$VAL < limites[1]] <- limites[1]
  data$VAL[data$VAL > limites[2]] <- limites[2]
  
  cortes_b <- c(-Inf, (cortes[1]+cortes[2])/2,
                (cortes[2]+cortes[3])/2,
                (cortes[3]+cortes[4])/2,
                (cortes[4]+cortes[5])/2, Inf)
  
  cortes_lb <- c(paste('<',cortes_b[2]),
                 paste(cortes_b[2], '-', cortes_b[3]),
                 paste(cortes_b[3], '-', cortes_b[4]),
                 paste(cortes_b[4], '-', cortes_b[5]),
                 paste('>',cortes_b[5]))
  
  data$VAL <- cut(data$VAL,
                  breaks = cortes_b,
                  right = FALSE,
                  labels= cortes_lb)
  
  mapa <- ggplot(data = data, aes(x=LON, y=LAT, col = VAL)) +
    geom_polygon(data=PERU, aes(x=long, y=lat, group=group), fill=NA, color="gray40", size = 0.3)+
    geom_point(size=2) +
    coord_equal() +
    coord_sf(crs=4326, xlim = c(-81,-69), ylim = c(-18,-0))+
    scale_color_manual(values=c("#D50600",'#ED6A66','#B9CFF0','#699DEC','#007AD5'), 
                       na.value = 'transparent', drop = F)+
    guides(color = guide_legend(title.position = "top",
                                reverse=F,
                                title.hjust = 0.5, title = metrics, nrow = 2))+
    annotate(geom = "text", x = -81.5, y = 0.2, label = lab_id, hjust = "left", size = 4, fontface =2)+
    theme_bw()
  return(mapa)
}

MAP_METRIC_EP_3 <- function(data = Eop_Eoc, PERU=peru, metrics = METRIC, lab_id = 'a',
                            punto_medio = 0, nrw = 2,
                            cortes = vecortes, limites= c(-1,1)){
  data$VAL[data$VAL < limites[1]] <- limites[1]
  data$VAL[data$VAL > limites[2]] <- limites[2]
  
  cortes_b <- c(-Inf, (cortes[1]+cortes[2])/2,
                (cortes[2]+cortes[3])/2,
                (cortes[3]+cortes[4])/2,
                (cortes[4]+cortes[5])/2,
                (cortes[5]+cortes[6])/2, Inf)
  
  cortes_lb <- c(paste('<',cortes_b[2]),
                 paste(cortes_b[2], '-', cortes_b[3]),
                 paste(cortes_b[3], '-', cortes_b[4]),
                 paste(cortes_b[4], '-', cortes_b[5]),
                 paste(cortes_b[5], '-', cortes_b[6]),
                 paste('>',cortes_b[6]))
  
  data$VAL <- cut(data$VAL,
                  breaks = cortes_b,
                  right = FALSE,
                  labels= cortes_lb)
  
  mapa <- ggplot(data = data, aes(x=LON, y=LAT, col = VAL)) +
    geom_polygon(data=PERU, aes(x=long, y=lat, group=group), fill=NA, color="gray40", size = 0.3)+
    geom_point(size=1.75) +
    coord_equal() +
    coord_sf(crs=4326, xlim = c(-81,-69), ylim = c(-18,-0))+
    scale_color_manual(values=c("#D01611", "#DBE8EC",'#9DD8FB','#49ABE6','#137CBB','#007AD5'), 
                       drop = FALSE,
                       labels = cortes_lb,
                       na.value = 'transparent')+
    guides(color = guide_legend(title.position = "top",
                                reverse=F,
                                title.hjust = 0.5, title = metrics, nrow = nrw))+
    annotate(geom = "text", x = -81.5, y = 0.2, label = lab_id, hjust = "left", size = 4, fontface =2)+
    theme_bw()
  return(mapa)
}

map_eo_prod <- function(Eo_raster = Eo_raster){
  Eo_MN_PER_df <- as.data.frame(Eo_raster[[1]], xy=T)
  Eo_MN_PER_df <- na.omit(Eo_MN_PER_df)
  Eo_MN_PER_df$SRC <- 'PISCOeo_pm'
  Eo_TRC_PER_df <- as.data.frame(Eo_raster[[2]], xy=T)
  Eo_TRC_PER_df <- na.omit(Eo_TRC_PER_df)
  Eo_TRC_PER_df$SRC <- 'TerraClimate'
  Eo_CRU_PER_df <- as.data.frame(Eo_raster[[3]], xy=T)
  Eo_CRU_PER_df[Eo_CRU_PER_df==0] <- NA
  Eo_CRU_PER_df <- na.omit(Eo_CRU_PER_df)
  Eo_CRU_PER_df$SRC <- 'CRU_TS'
  Eo_ERA_PER_df <- as.data.frame(Eo_raster[[4]], xy=T)
  Eo_ERA_PER_df <- na.omit(Eo_ERA_PER_df)
  Eo_ERA_PER_df$SRC <- 'ERA5-Land'

  Eo_TRCdPER_df <- Eo_raster[[2]]
  Eo_TRCdPER_df <- na.omit(as.data.frame(Eo_raster[[1]]-Eo_TRCdPER_df, xy=T))
  Eo_TRCdPER_df$SRC <- 'PISCOeo_pm - TerraClimate'
  
  Eo_CRUdPER_df <- Eo_raster[[3]]
  Eo_CRUdPER_df[Eo_CRUdPER_df == 0] <- NA
  Eo_CRUdPER_df <- na.omit(as.data.frame(Eo_raster[[1]]-Eo_CRUdPER_df, xy=T))
  Eo_CRUdPER_df$SRC <- 'PISCOeo_pm - CRU_TS'
  
  Eo_ERAdPER_df <- Eo_raster[[4]]
  Eo_ERAdPER_df <- na.omit(as.data.frame(Eo_raster[[1]]-Eo_ERAdPER_df, xy=T))
  Eo_ERAdPER_df$SRC <- 'PISCOeo_pm - ERA5-Land'
  
  Eo_yr_tif <- rbind(Eo_MN_PER_df,Eo_TRC_PER_df,Eo_CRU_PER_df,Eo_ERA_PER_df)
  
  Eo_yr_tif[Eo_yr_tif==0] <- NA
  Eo_yr_tif$clasd <- cut(Eo_yr_tif[,c(3)],
                         breaks = c(-Inf,seq(1050, 1950, 100),Inf),
                         right = FALSE,
                         labels=c('< 1,000','1,100','1,200','1,300','1,400','1,500',
                                  '1,600','1,700','1,800','1,900', "> 2,000"))
  
  Eo_df_tif <- rbind(Eo_TRCdPER_df, Eo_CRUdPER_df, Eo_ERAdPER_df)
  
  Eo_df_tif$clasd <- cut(Eo_df_tif[,c(3)],
                         breaks = c(-Inf,seq(-600, -150, 150),-50,0,50 ,seq(150, 600, 150),Inf),
                         right = FALSE,
                         labels=c('< -600', '-600', '-450','-300','-150', '-50',
                                  '0', '50','150', '300','450',"> 600"))
  
  colors_eo <- RColorBrewer::brewer.pal(11, "Spectral")
  colors_df <- RColorBrewer::brewer.pal(10, "RdYlBu")
  colors_df <- c(colors_df[1:5], 'white' ,colors_df[6:10])
  
  Eo_yr_tif$SRC <- factor(Eo_yr_tif$SRC, levels = c('PISCOeo_pm','CRU_TS','TerraClimate','ERA5-Land'))
  Eo_yr_tif <- na.omit(Eo_yr_tif)
  PERU     <- fortify(shapefile("data/raw/spatial/SEC_CLIM.shp"))
  Eo_df_tif$SRC <- factor(Eo_df_tif$SRC, levels = c('PISCOeo_pm - CRU_TS','PISCOeo_pm - TerraClimate','PISCOeo_pm - ERA5-Land'))
  Eo_df_tif <- na.omit(Eo_df_tif)
  
  map1 <- ggplot() +
    geom_raster(data=Eo_yr_tif, aes(x=x,y=y, fill = clasd)) +
    facet_wrap(vars(SRC), ncol = 4)+
    coord_equal() +
    coord_sf(crs=4326, xlim = c(-81,-69), ylim = c(-18,-0))+
    scale_fill_manual(values=colors_eo)+
    guides(fill = guide_legend(title.position = "top",
                               reverse=F,title = expression(ET[o] (mm/year)),
                               title.vjust = 0.5, ncol = 1))+
    theme_bw()+
    theme(legend.position = 'right')
  
  map2 <- ggplot() +
    geom_raster(data=Eo_df_tif, aes(x=x,y=y, fill = clasd)) +
    facet_wrap(vars(SRC), ncol = 4)+
    coord_equal() +
    coord_sf(crs=4326, xlim = c(-81,-69), ylim = c(-18,-0))+
    scale_fill_manual(values=colors_df)+
    guides(fill = guide_legend(title.position = "top",
                               reverse=F,title = expression(ET[o] (mm/year)),
                               title.vjust = 0.5, ncol = 1))+
    theme_bw()+
    theme(legend.position = 'right')
  
  
  plotcomb <- ggdraw()+
    draw_plot(map1, x = 0, y = 0.5, height = 0.5, width = 1)+
    draw_plot(map2, x = 0, y = 0, height = 0.5, width = 1)
  
  return(plotcomb)
}

TREND_TIMES <- function(datast=estac_p){
  ATrend_paa  <- data.frame()
  
  for(i in 1:length(datast)){
    ss      <- sens.slope(datast[,i])
    sslp_paa   <- ss[["estimates"]][["Sen's slope"]]
    
    res_mk   <- mk.test(datast[,i], continuity = TRUE)
    pval_paa <- res_mk$p.value
    MK99_paa <- ifelse(res_mk$p.value < 0.01, 3, 0) 
    MK95_paa <- ifelse(res_mk$p.value > 0.01 & res_mk$p.value < 0.05 , 2, 0) 
    MK90_paa <- ifelse(res_mk$p.value > 0.05 & res_mk$p.value < 0.1 , 1, 0)
    
    df <- data.frame(sslp_paa,pval_paa,MK_paa = sum(MK99_paa, MK95_paa, MK90_paa))
    ATrend_paa <- rbind(ATrend_paa,df)
  }
  return(ATrend_paa)
}