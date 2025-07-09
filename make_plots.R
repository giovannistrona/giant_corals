library(geosphere)
library(sf)
library(tidyr)
library(dplyr)
library(viridis)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(ggridges)
library(dbscan)


dir.create('plots')
or_data<-read.csv('MTG_data_16_04.csv',header=T)
genera<-or_data$Genus
gts<-read.csv('giant_sites_time_series.csv',header=T)
head(gts)
obs_xy<-gts[gts$rep==0 & gts$year==1986 & gts$month==1,c('lon','lat')]

null_xy<-unique(gts[gts$rep>0,c('lon','lat')])


####make map

lat2<-null_xy[,2]
lon2<-null_xy[,1]
lat1<-obs_xy[,2]
lon1<-obs_xy[,1]


pts1 <- st_as_sf(data.frame(
  lon = obs_xy[,1],
  lat = obs_xy[,2],
  gruppo = "observed_sites"
), coords = c("lon", "lat"), crs = 4326)

pts2 <- st_as_sf(data.frame(
  lon = null_xy[,1],
  lat = null_xy[,2],
  gruppo = "null_sites"
), coords = c("lon", "lat"), crs = 4326)


coords <- st_coordinates(pts1)
clust <- dbscan(coords, eps = 3, minPts = 1)  # eps ~1 km
pts1$cluster_id <- clust$cluster

cluster_summary <- pts1 %>%
  group_by(cluster_id) %>%
  summarize(
    n = n(),
    geometry = st_centroid(st_union(geometry)),
    .groups = "drop"
  ) %>%
  filter(n > 0)

cluster_summary_coords <- cbind(
  cluster_summary,
  st_coordinates(cluster_summary)
)

label_df <- cluster_summary_coords %>%
  mutate(
    lon_label = if_else(X > 140, X - 3, X + 3),
    lat_label = Y,
    xend_line = if_else(lon_label < 140, lon_label - 1, lon_label + 1),
    yend_line = lat_label
  )


land <- ne_countries(scale = "medium", returnclass = "sf")

all_points <- c(pts1$geometry, pts2$geometry)
bbox <- st_bbox(st_sfc(all_points, crs = 4326))
bbox_expanded <- bbox + c(-0.2, -0.2, 0.2, 0.2)

pdf("./plots/map_with_clusters.pdf", width = 8, height = 6)

ggplot() +
  geom_sf(data = land, fill = "grey80", color = "grey60") +
  geom_sf(data = pts2, aes(color = gruppo), size = 1, alpha = 0.05, show.legend = FALSE) +
  geom_sf(data = pts1, aes(color = gruppo), size = 0.8, alpha = 0.8, show.legend = FALSE) +
  # geom_segment(
  #   data = label_df,
  #   aes(x = X, y = Y, xend = xend_line, yend = yend_line),
  #   color = "black",
  #   size = 0.15,
  #   alpha = 0.6
  # ) +
  geom_text(
    data = label_df,
    aes(x = lon_label, y = lat_label, label = n),
    size = 5,
    color = "black"
  ) +
  scale_color_manual(
    values = c("observed_sites" = "gold", "null_sites" = "darkblue")
  ) +
  coord_sf(
    xlim = c(bbox_expanded["xmin"], bbox_expanded["xmax"]),
    ylim = c(bbox_expanded["ymin"], bbox_expanded["ymax"]),
    expand = FALSE
  ) +
  labs(x = NULL, y = NULL)+
  theme_minimal()+
  theme(
    axis.text = element_text(size = 12)
  )
  

dev.off()





###################
obs<-gts[which(gts$rep==0),]
null<-gts[which(gts$rep>0),]

sink('number_for_text.txt')
#############count only thermal anomalies >dhw8
obs$dhw_max<-1*(obs$dhw_max>=8)
null$dhw_max<-1*(null$dhw_max>=8)
obs_dhw<-aggregate(obs$dhw_max~obs$site_id,FUN='sum')

paste('number of colonies that did not experience DHW>8:',sum(obs_dhw$`obs$dhw_max`==0))
paste('% of colonies that did not experienced DHW>8:',round(100*sum(obs_dhw$`obs$dhw_max`==0)/nrow(obs_dhw),2)) #abstract, line 35

paste('number of colonies that experienced DHW>8:',sum(obs_dhw$`obs$dhw_max`>0))
paste('mean number of DHW>8 events for colonies that experienced DHW>8:',mean(obs_dhw[which(obs_dhw$`obs$dhw_max`>0),]$`obs$dhw_max`))
paste('median number of DHW>8 events for colonies that experienced DHW>8:',median(obs_dhw[which(obs_dhw$`obs$dhw_max`>0),]$`obs$dhw_max`))
paste('SD number of DHW>8 events for colonies that experienced DHW>8:',sd(obs_dhw[which(obs_dhw$`obs$dhw_max`>0),]$`obs$dhw_max`))
paste('min number of DHW>8 events for colonies that experienced DHW>8:',min(obs_dhw[which(obs_dhw$`obs$dhw_max`>0),]$`obs$dhw_max`))
paste('max number of DHW>8 events for colonies that experienced DHW>8:',max(obs_dhw[which(obs_dhw$`obs$dhw_max`>0),]$`obs$dhw_max`))


null_dhw<-aggregate(null$dhw_max~null$site_id+null$rep,FUN='sum')
mean_null<-aggregate(null_dhw$`null$dhw_max`~null_dhw$`null$site_id`,FUN='mean')
sd_null<-aggregate(null_dhw$`null$dhw_max`~null_dhw$`null$site_id`,FUN='sd')


p_values <- data.frame()
for (id in obs_dhw$`obs$site_id`) {
  obs_val <- as.numeric(obs_dhw[obs_dhw$`obs$site_id` == id, ]$`obs$dhw_max`)
  null_vals <- null_dhw[null_dhw$`null$site_id` == id, ]$`null$dhw_max`
  p_neg <- sum(null_vals <= obs_val) / length(null_vals)
  p_pos <- sum(null_vals >= obs_val) / length(null_vals)
  z <- (obs_val - mean(null_vals)) / sd(null_vals)
  
  row_data <- or_data[id + 1, ] %>%
    mutate(
      obs_val = obs_val,
      mean_null = mean(null_vals),
      sd_null = sd(null_vals),
      z = z,
      p_pos = p_pos,
      p_neg = p_neg
    )
  
  p_values <- bind_rows(p_values, row_data)
}


paste('% of giant sites with sum DHW>8 less than reference:',100*sum((obs_dhw$`obs$dhw_max`<=mean_null$`null_dhw$\`null$dhw_max\``))/nrow(obs_dhw))
paste('% of giant sites with sum DHW>8 more than reference:',100*sum((obs_dhw$`obs$dhw_max`>mean_null$`null_dhw$\`null$dhw_max\``))/nrow(obs_dhw))
paste('% of giant sites with sum DHW>8 == reference:',100*sum((obs_dhw$`obs$dhw_max`==mean_null$`null_dhw$\`null$dhw_max\``))/nrow(obs_dhw))


paste('number of colonies more exposed than control:',
      sum(p_values$obs_val>p_values$mean_null))

paste('number of colonies less exposed than control:',
      sum(p_values$obs_val<=p_values$mean_null))

paste('number of colonies significantly more exposed than control:',
      nrow(p_values[p_values$p_pos<0.01,]))

paste('number of colonies significantly less exposed than control:',
      nrow(p_values[p_values$p_neg<0.01,]))


paste('genera significantly more exposed than control:',
      p_values[p_values$p_pos<0.01,]$Genus)




######
library(geosphere)

coords <- matrix(c(
  115.5240, -8.676252,
  115.5138, -8.674901,
  115.5143, -8.674934,
  115.5154, -8.675428,
  115.5185, -8.670468,
  115.4572, -8.715528,
  115.5240, -8.676252  # chiusura del poligono
), ncol=2, byrow=TRUE)

area_m2 <- areaPolygon(coords)
area_km2 <- area_m2 / 1e6

cat("Area:", round(area_km2, 3), "km²\n")



#####
###distribution of dhw across years
gts<-read.csv('giant_sites_time_series.csv',header=T)

obs<-gts[which(gts$rep==0),]

obs_y_max<-aggregate(obs$dhw_max~obs$site_id+obs$year,FUN='max')

df<-data.frame("valore"=obs_y_max$`obs$dhw_max`,"anno"=obs_y_max$`obs$year`)
med_y<-median(unique(df$anno))

             
df$periodo <- ifelse(df$anno <= 2005, "1986–2005", "2006–2024")
df$valore_log <- log1p(df$valore)  # log(valore + 1)
df$anno <- factor(df$anno, levels = rev(sort(unique(df$anno))))


pdf('./plots/joyplot.pdf',height=4,width=6)

ggplot(df, aes(x = valore_log, y = factor(anno), fill = factor(anno))) +
  stat_density(
    aes(height = ..density.., color = after_scale(fill)),
    geom = "ridgeline",
    position = "identity",
    adjust = 0.5,
    scale = 0.4,
    alpha = 0.5,
    size = 0.3
  ) +
  geom_segment(
    data = df %>% group_by(periodo) %>% summarise(n_anni = n_distinct(anno)),
    aes(x = -Inf, xend = -Inf,
        y = n_anni, yend = 1),
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    color = "black",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ periodo, ncol = 2, scales = "free_y") +
  geom_vline(xintercept = log1p(8), linetype = "dashed", color = "red", linewidth = 1) +
  scale_x_continuous(
    name = "max DHW",
    breaks = log1p(c(0,2,4,8,16,32)),
    labels = function(x) round(expm1(x), 1)
  ) +
    labs(
    y = "year"
  ) +
  coord_cartesian(clip = "off")+
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey80", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey90", size = 0.25),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.title.x = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10)  
  )

dev.off()


df$periodo <- ifelse(as.numeric(as.character(df$anno)) <= 2100, "1986–2005", "2006–2024")

pdf('./plots/joyplot_stacked.pdf',height=12,width=6)
ggplot(df, aes(x = valore_log, y = factor(anno), fill = factor(anno))) +
  stat_density(
    aes(height = ..density.., color = after_scale(fill)),  
    geom = "ridgeline",
    position = "identity",
    adjust = 0.5,
    scale = 0.4,
    alpha = 0.5,
    size = 0.3
  ) +
  geom_segment(
    data = df %>% group_by(periodo) %>% summarise(n_anni = n_distinct(anno)),
    aes(x = -Inf, xend = -Inf,
        y = n_anni, yend = 1),
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    color = "black",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  geom_vline(xintercept = log1p(8), linetype = "dashed", color = "red", linewidth = 1) +
  scale_x_continuous(
    name = "max DHW",
    breaks = log1p(c(0,2,4,8,16,32)),
    labels = function(x) round(expm1(x), 1)
  ) +
  labs(
    y = "year"
  ) +
  coord_cartesian(clip = "off")+
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey80", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey90", size = 0.25),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.text.x = element_text(color = "black", size = 10),
    axis.title.x = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10)  
  )

dev.off()



####comparison with future
gts_fut<-read.csv('giant_sites_time_series_future.csv',header=T)
gts_fut<-gts_fut[gts_fut$year>2024,]
obs_fut<-gts_fut[which(gts_fut$rep==0),]
null_fut<-gts_fut[which(gts_fut$rep>0),]
obs_fut$dhw_max<-1*(obs_fut$dhw_max>=8)
null_fut$dhw_max<-1*(null_fut$dhw_max>=8)



obs_dhw_fut<-aggregate(obs_fut$dhw_max~obs_fut$site_id,FUN='sum')


#######plot aggregate time series for observed and null
finite_mean<-function(x){
  x<-x[is.finite(x)]
  return (mean(x))
}

gts_fut<-read.csv('giant_sites_time_series_future.csv',header=T)

obs_fut<-gts_fut[which(gts_fut$rep==0),]

null_fut<-gts_fut[which(gts_fut$rep>0),]


ref<-obs_fut[obs_fut$year<=2024,]
ref<-ref[ref$year>2014,]
fut<-obs_fut[obs_fut$year>2024,]
mon_mean<-aggregate(ref$dhw_max~ref$site_id+ref$month,FUN='mean')

res<-c()
for (year in 2025:2100){
  data<-fut[fut$year==year,]
  diff<-100*log(data$dhw_max/mon_mean$`ref$dhw_max`)
  diff_mean<-aggregate(diff~data$site_id,FUN='finite_mean')
  res<-rbind(res,cbind(year,diff_mean[,2]))
}



df<-data.frame("valore"=res[,2],"anno"=res[,1])
df$decade <- paste0(floor(df$anno / 10) * 10, "s")
pdf('./plots/future_exposure.pdf',width=8,height=4)
ggplot(df, aes(x = decade, y = valore, fill = decade)) +
  geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.9) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.6, color = "gray20") +
  #scale_fill_viridis(discrete = TRUE, option = "D") +
  geom_hline(yintercept = c(-100, 100, 200, 300), linetype = "dashed", color = "red", linewidth = 0.3) +
  geom_hline(yintercept = c(0), linetype = "dashed", color = "black", linewidth = 0.3) +
  theme_minimal() +
  labs(x = "decade", y = "log % difference DHW/DHW reference period (2014-2024)") +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey80", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey90", size = 0.25),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.text.x = element_text(color = "black", size = 10), 
    axis.title.x = element_text(size = 12),
    axis.line.y = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(size = 12),
    plot.margin = margin(10, 10, 10, 10)  
  )


dev.off()


sink()