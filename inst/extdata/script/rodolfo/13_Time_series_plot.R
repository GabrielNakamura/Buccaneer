# Creates plots present in research article and ESM ####

rm(list = ls());gc()

library(ggplot2)
library(cowplot)
library(gridExtra)
library(corrplot)
library(RColorBrewer)

dir.create("./Figs/MAIN_TEXT", recursive = T)
# loading ####

# REGIONAL_DIVERSITY ####

replicates <- 1:50 # here the "50" is the number of replicas

reg_files <- paste0("./Spatio_temporal_analyses/Diversity_curves/diversity", replicates, ".txt")
reg_bind <- do.call(rbind, lapply(reg_files, read.table, header =T))

names(reg_bind)[2] <- "val"
reg_bind$id <- "reg"

reg_bind_qs <- data.frame(do.call(rbind, tapply(reg_bind$val, reg_bind$time, function (x) quantile(x, c(0,0.05,0.5,0.95,1), na.rm = T))), time = unique(reg_bind$time), id = "reg", row.names = NULL)
reg_bind_qs[is.na(reg_bind_qs)] <- 0 # changes NAs for zeros



# REACH_DIVERSITY ####

reach_files <- paste0("./Spatio_temporal_analyses/Time_space/reach/reach_diversity", replicates, ".txt")
reach_bind <- do.call(rbind, lapply(reach_files, read.table, header =T))

names(reach_bind)[2] <- "val"
reach_bind$id <- "reach"

reach_bind_qs <- data.frame(do.call(rbind, tapply(reach_bind$val, reach_bind$time, function (x) quantile(x, c(0,0.05,0.5,0.95,1), na.rm = T))), time = unique(reach_bind$time), id = "reach", row.names = NULL)
reach_bind_qs[is.na(reach_bind_qs)] <- 0 # changes NAs for zeros


# SITE_DIVERSITY ####


site_files <- paste0("./Spatio_temporal_analyses/Time_space/site/site_diversity", replicates, ".txt")
site_bind <- do.call(rbind, lapply(site_files, read.table, header =T))

names(site_bind)[2] <- "val"
site_bind$id <- "site"

site_bind_qs <- data.frame(do.call(rbind, tapply(site_bind$val, site_bind$time, function (x) quantile(x, c(0,0.05,0.5,0.95,1), na.rm = T))), time = unique(site_bind$time), id = "site", row.names = NULL)
site_bind_qs[is.na(site_bind_qs)] <- 0 # changes NAs for zeros

# REGIONAL_MPD ####

reg_mpd_files <- paste0("./Ecomorphological_data/competition/regional/mpd/regional_mpd_diet_", replicates, ".txt")
reg_mpd_bind <- do.call(rbind, lapply(reg_mpd_files, read.table, header =T))

names(reg_mpd_bind)[2] <- "val"
reg_mpd_bind$id <- "reg_mpd"

reg_mpd_bind_qs <- data.frame(do.call(rbind, tapply(reg_mpd_bind$val, reg_mpd_bind$time, function (x) quantile(x, c(0,0.05,0.5,0.95,1), na.rm = T))), time = unique(reg_mpd_bind$time), id = "reg_mpd", row.names = NULL)
reg_mpd_bind_qs[is.na(reg_mpd_bind_qs)] <- 0 # changes NAs for zeros

# REACH_MPD ####

reach_mpd_files <- paste0("./Ecomorphological_data/competition/reach/mpd/reach_mpd_diet_", replicates, ".txt")
reach_mpd_bind <- do.call(rbind, lapply(reach_mpd_files, read.table, header =T))

names(reach_mpd_bind)[2] <- "val"
reach_mpd_bind$id <- "reach_mpd"

reach_mpd_bind_qs <- data.frame(do.call(rbind, tapply(reach_mpd_bind$val, reach_mpd_bind$time, function (x) quantile(x, c(0,0.05,0.5,0.95,1), na.rm = T))), time = unique(reach_mpd_bind$time), id = "reach_mpd", row.names = NULL)
reach_mpd_bind_qs[is.na(reach_mpd_bind_qs)] <- 0 # changes NAs for zeros


# SITE_MPD ####

site_mpd_files <- paste0("./Ecomorphological_data/competition/site/mpd/site_mpd_diet_", replicates, ".txt")
site_mpd_bind <- do.call(rbind, lapply(site_mpd_files, read.table, header =T))

names(site_mpd_bind)[2] <- "val"
site_mpd_bind$id <- "site_mpd"

site_mpd_bind_qs <- data.frame(do.call(rbind, tapply(site_mpd_bind$val, site_mpd_bind$time, function (x) quantile(x, c(0,0.05,0.5,0.95,1), na.rm = T))), time = unique(site_mpd_bind$time), id = "site_mpd", row.names = NULL)
site_mpd_bind_qs[is.na(site_mpd_bind_qs)] <- 0 # changes NAs for zeros

# REGIONAL_MNND ####

reg_mnnd_files <- paste0("./Ecomorphological_data/competition/regional/mnnd/regional_mnnd_", replicates, ".txt")
reg_mnnd_bind <- do.call(rbind, lapply(reg_mnnd_files, read.table, header =T))

names(reg_mnnd_bind)[2] <- "val"
reg_mnnd_bind$id <- "reg_mnnd"

reg_mnnd_bind_qs <- data.frame(do.call(rbind, tapply(reg_mnnd_bind$val, reg_mnnd_bind$time, function (x) quantile(x, c(0,0.05,0.5,0.95,1), na.rm = T))), time = unique(reg_mnnd_bind$time), id = "reg_mnnd", row.names = NULL)
reg_mnnd_bind_qs[is.na(reg_mnnd_bind_qs)] <- 0 # changes NAs for zeros


# REACH_MNND ####

reach_mnnd_files <- paste0("./Ecomorphological_data/competition/reach/mnnd/reach_mnnd_", replicates, ".txt")
reach_mnnd_bind <- do.call(rbind, lapply(reach_mnnd_files, read.table, header =T))

names(reach_mnnd_bind)[2] <- "val"
reach_mnnd_bind$id <- "reach_mnnd"

reach_mnnd_bind_qs <- data.frame(do.call(rbind, tapply(reach_mnnd_bind$val, reach_mnnd_bind$time, function (x) quantile(x, c(0,0.05,0.5,0.95,1), na.rm = T))), time = unique(reach_mnnd_bind$time), id = "reach_mnnd", row.names = NULL)
reach_mnnd_bind_qs[is.na(reach_mnnd_bind_qs)] <- 0 # changes NAs for zeros


# SITE_MNND ####

site_mnnd_files <- paste0("./Ecomorphological_data/competition/site/mnnd/site_mnnd_", replicates, ".txt")
site_mnnd_bind <- do.call(rbind, lapply(site_mnnd_files, read.table, header =T))

names(site_mnnd_bind)[2] <- "val"
site_mnnd_bind$id <- "site_mnnd"

site_mnnd_bind_qs <- data.frame(do.call(rbind, tapply(site_mnnd_bind$val, site_mnnd_bind$time, function (x) quantile(x, c(0,0.05,0.5,0.95,1), na.rm = T))), time = unique(site_mnnd_bind$time), id = "site_mnnd", row.names = NULL)
site_mnnd_bind_qs[is.na(site_mnnd_bind_qs)] <- 0 # changes NAs for zeros

# Temperature ####

temperature <- read.table("./Temperature/temperature.txt", header = T)

# Combining and plotting ####

tot_bind_qs <- rbind(reg_bind_qs, reach_bind_qs, site_bind_qs, reg_mpd_bind_qs, reach_mpd_bind_qs, site_mpd_bind_qs, reg_mnnd_bind_qs, reach_mnnd_bind_qs, site_mnnd_bind_qs)

tot_bind <- rbind(reg_bind, reach_bind, site_bind, reg_mpd_bind, reach_mpd_bind, site_mpd_bind, reg_mnnd_bind, reach_mnnd_bind, site_mnnd_bind)

save(tot_bind, tot_bind_qs, file = "./Ecomorphological_data/time_series_for_plot.Rdata", compress = "xz")

load("./Ecomorphological_data/time_series_for_plot.Rdata")

cols <- c("#00ee18", "#dfe947", "#eb75b0", "#008c4a", "#e7d661", "#b30041", "#054005", "#e0ab0e", "#7d129c")

#define the boundaries 
frames <- c(34, 23, 5, 2, 0) # Barnosky 2014
bottom <- frames[-length(frames)] # define the beginning of each boundary
top <- frames[-1] # define the end of each boundary

frames.df <- data.frame(bottom, top, col = rep(x = c("azure4", "white"), length.out = length(frames)-1)) # Here is the length of the objects bottom and top (they must have the same lengths)

scaleFUN <- function(x) sprintf("%.1f", x)
g1 <-
  ggplot(data = tot_bind_qs[which(tot_bind_qs$id == "reg"),], mapping = aes(x = time), cex = 2) +
  geom_rect(data = frames.df,
            aes(NULL, NULL, xmin = bottom, xmax = top, ymin = 0, ymax = 35, fill = col),colour = "white", alpha = 0.4, show.legend = F) +
  scale_fill_manual(values = c("azure4" = "azure4", "white" = "white")) +
  geom_ribbon(mapping = aes(ymin = X0., ymax = X100., group = id, colour = NULL), alpha = 0.5, show.legend = F) +
  geom_line(data = tot_bind_qs[which(tot_bind_qs$id == "reg"),],
            mapping = aes(x = time, y = X50.), linewidth = 2, col = cols[1]) +
  scale_x_reverse() + theme_cowplot() +
  ggtitle("A)") +
  labs(tag = "Regional") +
  ylab("Species (n)") +
  theme(axis.title.y = element_text(
    size = 60, face = "bold", margin = margin(r = 20)),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text = element_text(size = 60),
    plot.title = element_text(size = 60),
    plot.tag.position = "top",
    plot.tag = element_text(size = 60))

g2 <-
  ggplot(data = tot_bind_qs[which(tot_bind_qs$id == "reach"),], mapping = aes(x = time), cex = 2) +
  geom_rect(data = frames.df,
            aes(NULL, NULL, xmin = bottom, xmax = top, ymin = 0, ymax = 35, fill = col),colour = "white", alpha = 0.4, show.legend = F) +
  scale_fill_manual(values = c("azure4" = "azure4", "white" = "white")) +
  geom_ribbon(mapping = aes(ymin = X0., ymax = X100., group = id, colour = NULL), alpha = 0.5, show.legend = F) +
  geom_line(data = tot_bind_qs[which(tot_bind_qs$id == "reach"),],
            mapping = aes(x = time, y = X50.), linewidth = 2, col = cols[2]) +
  scale_x_reverse() + theme_cowplot() +
  ggtitle("B)") +
  labs(tag = "Reach") +
  theme(axis.title.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60),
        plot.tag.position = "top",
        plot.tag = element_text(size = 60))

g3 <-
  ggplot(data = tot_bind_qs[which(tot_bind_qs$id == "site"),], mapping = aes(x = time), cex = 2) +
  geom_rect(data = frames.df,
            aes(NULL, NULL, xmin = bottom, xmax = top, ymin = 0, ymax = 10, fill = col),colour = "white", alpha = 0.4, show.legend = F) +
  scale_fill_manual(values = c("azure4" = "azure4", "white" = "white")) +
  geom_ribbon(mapping = aes(ymin = X0., ymax = X100., group = id, colour = NULL), alpha = 0.5, show.legend = F) +
  geom_line(data = tot_bind_qs[which(tot_bind_qs$id == "site"),],
            mapping = aes(x = time, y = X50.), linewidth = 2, col = cols[3]) +
  scale_x_reverse() + theme_cowplot() +
  ggtitle("C)") +
  labs(tag = "Site") +
  theme(axis.title.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60),
        plot.tag.position = "top",
        plot.tag = element_text(size = 60))

g4 <-
  ggplot(data = tot_bind_qs[which(tot_bind_qs$id == "reg_mpd"),], mapping = aes(x = time), cex = 2) +
  geom_rect(data = frames.df,
            aes(NULL, NULL, xmin = bottom, xmax = top, ymin = 0, ymax = 3.6, fill = col),colour = "white", alpha = 0.4, show.legend = F) +
  scale_fill_manual(values = c("azure4" = "azure4", "white" = "white")) +
  geom_ribbon(mapping = aes(ymin = X0., ymax = X100., group = id, colour = NULL), alpha = 0.5, show.legend = F) +
  geom_line(data = tot_bind_qs[which(tot_bind_qs$id == "reg_mpd"),],
            mapping = aes(x = time, y = X50.), linewidth = 2, col = cols[4]) +
  scale_x_reverse() + 
  scale_y_continuous(labels = scaleFUN) + 
  theme_cowplot() +
  ggtitle("D)") +
  ylab("MPD_diet") +
  theme(axis.title.y = element_text(
    size = 60, face = "bold", margin = margin(r = 20)),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text = element_text(size = 60),
    plot.title = element_text(size = 60))

g5 <-
  ggplot(data = tot_bind_qs[which(tot_bind_qs$id == "reach_mpd"),], mapping = aes(x = time), cex = 2) +
  geom_rect(data = frames.df,
            aes(NULL, NULL, xmin = bottom, xmax = top, ymin = 0, ymax = 3.6, fill = col),colour = "white", alpha = 0.4, show.legend = F) +
  scale_fill_manual(values = c("azure4" = "azure4", "white" = "white")) +
  geom_ribbon(mapping = aes(ymin = X0., ymax = X100., group = id, colour = NULL), alpha = 0.5, show.legend = F) +
  geom_line(data = tot_bind_qs[which(tot_bind_qs$id == "reach_mpd"),],
            mapping = aes(x = time, y = X50.), linewidth = 2, col = cols[5]) +
  scale_x_reverse() + 
  scale_y_continuous(labels = scaleFUN) + 
  theme_cowplot() +
  ggtitle("E)") +
  theme(axis.title.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60))

g6 <-
  ggplot(data = tot_bind_qs[which(tot_bind_qs$id == "site_mpd"),], mapping = aes(x = time), cex = 2) +
  geom_rect(data = frames.df,
            aes(NULL, NULL, xmin = bottom, xmax = top, ymin = 0, ymax = 3.6, fill = col),colour = "white", alpha = 0.4, show.legend = F) +
  scale_fill_manual(values = c("azure4" = "azure4", "white" = "white")) +
  geom_ribbon(mapping = aes(ymin = X0., ymax = X100., group = id, colour = NULL), alpha = 0.5, show.legend = F) +
  geom_line(data = tot_bind_qs[which(tot_bind_qs$id == "site_mpd"),],
            mapping = aes(x = time, y = X50.), linewidth = 2, col = cols[6]) +
  scale_x_reverse() + 
  scale_y_continuous(labels = scaleFUN) + 
  theme_cowplot() +
  ggtitle("F)") +
  theme(axis.title.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60))

g7 <-
  ggplot(data = tot_bind_qs[which(tot_bind_qs$id == "reg_mnnd"),], mapping = aes(x = time), cex = 2) +
  geom_rect(data = frames.df,
            aes(NULL, NULL, xmin = bottom, xmax = top, ymin = 0, ymax = 1.5, fill = col),colour = "white", alpha = 0.4, show.legend = F) +
  scale_fill_manual(values = c("azure4" = "azure4", "white" = "white")) +
  geom_ribbon(mapping = aes(ymin = X0., ymax = X100., group = id, colour = NULL), alpha = 0.5, show.legend = F) +
  geom_line(data = tot_bind_qs[which(tot_bind_qs$id == "reg_mnnd"),],
            mapping = aes(x = time, y = X50.), linewidth = 2, col = cols[7]) +
  scale_x_reverse() + theme_cowplot() +
  ggtitle("G)") +    
  ylab("MNND") +
  xlab("Time (Ma)") +
  theme(axis.title.y = element_text(
    size = 60, face = "bold", margin = margin(r = 20)), 
    axis.text = element_text(size = 60),
    axis.title.x = element_text(size = 60, face = "bold",
                                margin = margin(t = 30)),
    plot.title = element_text(size = 60))

g8 <-
  ggplot(data = tot_bind_qs[which(tot_bind_qs$id == "reach_mnnd"),], mapping = aes(x = time), cex = 2) +
  geom_rect(data = frames.df,
            aes(NULL, NULL, xmin = bottom, xmax = top, ymin = 0, ymax = 1.5, fill = col),colour = "white", alpha = 0.4, show.legend = F) +
  scale_fill_manual(values = c("azure4" = "azure4", "white" = "white")) +
  geom_ribbon(mapping = aes(ymin = X0., ymax = X100., group = id, colour = NULL), alpha = 0.5, show.legend = F) +
  geom_line(data = tot_bind_qs[which(tot_bind_qs$id == "reach_mnnd"),],
            mapping = aes(x = time, y = X50.), linewidth = 2, col = cols[8]) +
  scale_x_reverse() + theme_cowplot() +
  ggtitle("H)") +
  xlab("Time (Ma)") +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 60),
        axis.title.x = element_text(size = 60, face ="bold",
                                    margin = margin(t = 30)),
        plot.title = element_text(size = 60))

g9 <-
  ggplot(data = tot_bind_qs[which(tot_bind_qs$id == "site_mnnd"),], mapping = aes(x = time), cex = 2) +
  geom_rect(data = frames.df,
            aes(NULL, NULL, xmin = bottom, xmax = top, ymin = 0, ymax = 3.6, fill = col),colour = "white", alpha = 0.4, show.legend = F) +
  scale_fill_manual(values = c("azure4" = "azure4", "white" = "white")) +
  geom_ribbon(mapping = aes(ymin = X0., ymax = X100., group = id, colour = NULL), alpha = 0.5, show.legend = F) +
  geom_line(data = tot_bind_qs[which(tot_bind_qs$id == "site_mnnd"),],
            mapping = aes(x = time, y = X50.), linewidth = 2, col = cols[9]) +
  scale_x_reverse() + 
  scale_y_continuous(labels = scaleFUN) + 
  theme_cowplot() +
  ggtitle("I)") +
  xlab("Time (Ma)") +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 60),
        axis.title.x = element_text(size = 60, face ="bold",
                                    margin = margin(t = 30)),
        plot.title = element_text(size = 60))

pdf(file = "./Figs/MAIN_TEXT/fig2.pdf", width = 40, height = 25)
grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, g9)
dev.off()

# Temperature
g10 <-
  ggplot(data = temperature, mapping = aes(x = time, y = temperature), cex = 2) +
  geom_rect(data = frames.df,
            aes(NULL, NULL, xmin = bottom, xmax = top, ymin = 0, ymax = 12, fill = col),colour = "white", alpha = 0.2, show.legend = F) +
  scale_fill_manual(values = c("azure4" = "azure4", "white" = "white")) +
  geom_line(linewidth = 2) +
  scale_x_reverse() + 
  theme_cowplot() +
  ylab("Sea-level temperature (°C)") +
  xlab("Time (Ma)") +
  theme(axis.title.y = element_text(
    size = 40, face = "bold", margin = margin(r = 20)),
    axis.text = element_text(size = 40),
    axis.title.x = element_text(
    size = 40, face ="bold", margin = margin(t = 30)))


pdf("./Figs/Temperature_time_series.pdf", width = 35, height = 20)
g10
dev.off()  
  
## Note that the graph will not be exacty the same as the one in the paper because a new imputation process 
## generates sligthly different values and hence time series.

M <- cor(data.frame(reg_bind_qs$X50., reach_bind_qs$X50., site_bind_qs$X50., reg_mpd_bind_qs$X50., reach_mpd_bind_qs$X50., site_mpd_bind_qs$X50., reg_mnnd_bind_qs$X50., reach_mnnd_bind_qs$X50., site_mnnd_bind_qs$X50.), method = "kendall")

colnames(M) <- c("Regional", "Reach", "Site", "Regional_MPD", "Reach_MPD", "Site_MPD", "Regional_MNND", "Reach_MNND", "Site_MNND")
rownames(M) <- c("Regional", "Reach", "Site", "Regional_MPD", "Reach_MPD", "Site_MPD", "Regional_MNND", "Reach_MNND", "Site_MNND")

# figure S8
corrplot(M, type = "upper", col = brewer.pal(n = 10, name = "RdYlBu"), diag = F, method = "color", outline = T, tl.col = "black", tl.srt = 45)
