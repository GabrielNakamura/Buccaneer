rm(list = ls());gc()

library(tidyverse)
# ggplot2 must be at least >= 3.5
library(cowplot)
theme_set(theme_cowplot())
library(coda)
library(gridExtra)

load("./Continuous/logs/epoch_results_tables.RData")

# define the function that extracts the summary ####
data_summary <- function(x){
  range <- coda::HPDinterval(as.mcmc((x)))
  median <- median(x)
  return(c(y = median, ymin = range[1], ymax = range[2]))
}

# define the colors used in main figures
cols <- c("#00ee18", "#dfe947", "#eb75b0", "#008c4a", "#e7d661", "#b30041", "#054005", "#e0ab0e", "#7d129c")
scaleFUN <- function(x) sprintf("%.1f", x)

# If you set different time series in the object results_tables, watch the indexes of the object and rename accordingly.

# Time windows - Epochs ####
# Speciation ####
# Diversity - results[[1]]
g1 <- 
  ggplot(epoch_results_tables[[1]][[1]], aes(x = ind, y = values)) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = -11, color = NA, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[1]) +
  coord_cartesian(ylim = c(-6.5,6.5)) +
  ggtitle("A)") +
  labs(tag = "Regional") +
  ylab("Gl - Spp") +
  theme(axis.title.y = element_text(size = 60, face = "bold",
                                    margin = margin(r = 20)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60),
        plot.tag.position = "top",
        plot.tag = element_text(size = 60)) +
  panel_border(color = "black")

# reach_diversity - results[[2]]
g2 <- 
  ggplot(epoch_results_tables[[2]][[1]], aes(x = ind, y = values)) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = -11, color = NA, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[2]) +
  coord_cartesian(ylim = c(-6.5,6.5)) +
  ggtitle("B)") +
  labs(tag = "Reach") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60),
        plot.tag.position = "top",
        plot.tag = element_text(size = 60)) +
  panel_border(color = "black")

# site_diversity - results[[3]]
g3 <- 
  ggplot(epoch_results_tables[[3]][[1]], aes(x = ind, y = values)) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = -11, color = NA, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[3]) +
  coord_cartesian(ylim = c(-6.5,6.5)) +
  ggtitle("C)") +
  labs(tag = "Site") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60),
        plot.tag.position = "top",
        plot.tag = element_text(size = 60)) +
  panel_border(color = "black")

# regional mpd - results[[4]]
g4 <- 
  ggplot(epoch_results_tables[[4]][[1]], aes(x = ind, y = values)) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = 11, color = NA, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[4]) +
  coord_cartesian(ylim = c(-3,3)) +
  ggtitle("D)") +
  ylab("Gl - MPD") +
  theme(axis.title.y = element_text(size = 60, face = "bold",
                                    margin = margin(r = 20)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60)) +
  panel_border(color = "black")

# reach_mpd - results[[5]]
g5 <- 
  ggplot(epoch_results_tables[[5]][[1]], aes(x = ind, y = values)) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = 11, color = NA, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[5]) +
  coord_cartesian(ylim = c(-3,3)) +
  ggtitle("E)") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60)) +
  panel_border(color = "black")

# site_mpd - results[[6]]
g6 <- 
  ggplot(epoch_results_tables[[6]][[1]], aes(x = ind, y = values)) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = 11, color = NA, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[6]) +
  coord_cartesian(ylim = c(-3,3)) +
  ggtitle("F)") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60)) +
  panel_border(color = "black")

# regional_mnnd - results[[7]]
g7 <- 
  ggplot(epoch_results_tables[[7]][[1]], aes(x = ind, y = values)) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = 11, color = NA, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[7]) +
  coord_cartesian(ylim = c(-3,3)) +
  ggtitle("G)") +
  ylab("Gl - MNND") +
  scale_x_discrete(labels = c("LE", "Oli", "Mio", "Plio", "Qua")) +
  theme(axis.title.y = 
          element_text(size = 60, 
                       face = "bold", margin = margin(r = 20)),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60),
        axis.title.x = element_blank()) +
  panel_border(color = "black")

# reach_mmnd - results[[8]]
g8 <- 
  ggplot(epoch_results_tables[[8]][[1]], aes(x = ind, y = values)) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = 11, color = NA, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[8]) +
  coord_cartesian(ylim = c(-3,3)) +
  ggtitle("H)") +
  scale_x_discrete(labels = c("LE", "Oli", "Mio", "Plio", "Qua")) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60),
        axis.title.x = element_blank()) +
  panel_border(color = "black")

# site_mnnd - results[[9]]
g9 <- 
  ggplot(epoch_results_tables[[9]][[1]], aes(x = ind, y = values)) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = 11, color = NA, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[9]) +
  coord_cartesian(ylim = c(-3,3)) +
  ggtitle("I)") +
  scale_x_discrete(labels = c("LE", "Oli", "Mio", "Plio", "Qua")) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60),
        axis.title.x = element_blank()) +
  panel_border(color = "black")

pdf(file = "./Figs/MAIN_TEXT/Fig_3.pdf", width = 40, height = 25)
grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, g9)
dev.off()

# Extinction ####
# Diversity - results[[1]]
g1 <- 
  ggplot(epoch_results_tables[[1]][[2]], aes(x = ind, y = values)) +
geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = 11, color = NA, alpha = 0.4) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[1]) +
  coord_cartesian(ylim = c(-6.5,6.5)) +
  ggtitle("A)") +
  labs(tag = "Regional") +
  ylab("Gm - Spp") +
  theme(axis.title.y = element_text(size = 60, face = "bold",
                                    margin = margin(r = 20)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60),
        plot.tag.position = "top",
        plot.tag = element_text(size = 60)) +
  panel_border(color = "black")

# reach_diversity - results[[2]]
g2 <- 
  ggplot(epoch_results_tables[[2]][[2]], aes(x = ind, y = values)) +
geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = 11, color = NA, alpha = 0.4) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[2]) +
  coord_cartesian(ylim = c(-6.5,6.5)) +
  ggtitle("B)") +
  labs(tag = "Reach") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60),
        plot.tag.position = "top",
        plot.tag = element_text(size = 60)) +
  panel_border(color = "black")

# site_diversity - results[[3]]
g3 <- 
  ggplot(epoch_results_tables[[3]][[2]], aes(x = ind, y = values)) +
geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = 11, color = NA, alpha = 0.4) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[3]) +
  coord_cartesian(ylim = c(-6.5,6.5)) +
  ggtitle("C)") +
  labs(tag = "Site") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60),
        plot.tag.position = "top",
        plot.tag = element_text(size = 60)) +
  panel_border(color = "black")

# regional mpd - results[[4]]
g4 <- 
  ggplot(epoch_results_tables[[4]][[2]], aes(x = ind, y = values)) +
geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = -11, color = NA, alpha = 0.4) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[4]) +
  coord_cartesian(ylim = c(-3,3)) +
  ggtitle("D)") +
  ylab("Gm - MPD") +
  theme(axis.title.y = element_text(size = 60, face = "bold",
                                    margin = margin(r = 20)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60)) +
  panel_border(color = "black")

# reach_mpd - results[[5]]
g5 <- 
  ggplot(epoch_results_tables[[5]][[2]], aes(x = ind, y = values)) +
geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = -11, color = NA, alpha = 0.4) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[5]) +
  coord_cartesian(ylim = c(-3,3)) +
  ggtitle("E)") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60)) +
  panel_border(color = "black")

# site_mpd - results[[6]]
g6 <- 
  ggplot(epoch_results_tables[[6]][[2]], aes(x = ind, y = values)) +
geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = -11, color = NA, alpha = 0.4) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[6]) +
  coord_cartesian(ylim = c(-3,3)) +
  ggtitle("F)") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60)) +
  panel_border(color = "black")

# regional_mnnd - results[[7]]
g7 <- 
  ggplot(epoch_results_tables[[7]][[2]], aes(x = ind, y = values)) +
geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = -11, color = NA, alpha = 0.4) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[7]) +
  coord_cartesian(ylim = c(-3,3)) +
  ggtitle("G)") +
  ylab("Gm - MNND") +
  scale_x_discrete(labels = c("LE", "Oli", "Mio", "Plio", "Qua")) +
  theme(axis.title.y = 
          element_text(size = 60, 
                       face = "bold", margin = margin(r = 20)),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60),
        axis.title.x = element_blank()) +
  panel_border(color = "black")

# reach_mmnd - results[[8]]
g8 <- 
  ggplot(epoch_results_tables[[8]][[2]], aes(x = ind, y = values)) +
geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = -11, color = NA, alpha = 0.4) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[8]) +
  coord_cartesian(ylim = c(-3,3)) +
  ggtitle("H)") +
  scale_x_discrete(labels = c("LE", "Oli", "Mio", "Plio", "Qua")) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60),
        axis.title.x = element_blank()) +
  panel_border(color = "black")

# site_mnnd - results[[9]]
g9 <- 
  ggplot(epoch_results_tables[[9]][[2]], aes(x = ind, y = values)) +
geom_hline(yintercept = 0, linetype = 2, linewidth = 3) +
  annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax = -11, color = NA, alpha = 0.4) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 4, linewidth = 8, col = cols[9]) +
  coord_cartesian(ylim = c(-3,3)) +
  ggtitle("I)") +
  scale_x_discrete(labels = c("LE", "Oli", "Mio", "Plio", "Qua")) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 60),
        plot.title = element_text(size = 60),
        axis.title.x = element_blank()) +
  panel_border(color = "black")

pdf(file = "./Figs/MAIN_TEXT/fig_4.pdf", width = 40, height = 25)
grid.arrange(g1, g2, g3, g4, g5, g6, g7, g8, g9)
dev.off()

# Temperature ####
# Speciation ####

g10 <-
  ggplot(epoch_results_tables[[10]][[1]], aes(x = ind, y = values)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 2) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 2, linewidth = 4, col = '#4c4cec') +
  coord_cartesian(ylim = c(-4,4)) +
  scale_x_discrete(labels = c("LE", "Oli", "Mio", "Plio", "Qua")) +
  ylab("Gl") +
  theme(axis.title.y = element_text(size = 40, face = "bold",
                                    margin = margin(r = 20)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 40),
        plot.title = element_text(size = 40),
        plot.tag.position = "top",
        plot.tag = element_text(size = 40)) +
  panel_border(color = "black")

# Extinction ####
g11 <-
  ggplot(epoch_results_tables[[10]][[2]], aes(x = ind, y = values)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 2) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 2, linewidth = 4, col = '#e34a33') +
  coord_cartesian(ylim = c(-4,4)) +
  scale_x_discrete(labels = c("LE", "Oli", "Mio", "Plio", "Qua")) +
  ylab("Gm") +
  theme(axis.title.y = element_text(size = 40, face = "bold",
                                    margin = margin(r = 20)),
        axis.text = element_text(size = 40),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 40),
        plot.tag.position = "top",
        plot.tag = element_text(size = 40)) +
  panel_border(color = "black")

pdf(file = "./Figs/MAIN_TEXT/fig_5.pdf", width = 10, height = 8)
grid.arrange(g10, g11)
dev.off()

# Time windows - Cenozoic ####
# Temperature ####
# Speciation ####

load("./Continuous/logs/ceno_results_table.RData")

g1 <-
  ggplot(ceno_results_table[which(ceno_results_table$ind == "Gl_t38.0"),], aes(x = ind, y = values)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 2) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 2, linewidth = 5, col = '#4c4cec') +
  coord_cartesian(ylim = c(-2,2)) +
  ylab("Gl") +
  theme(axis.title.y = element_text(size = 30, face = "bold",
                                    margin = margin(r = 20)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 30),
        plot.title = element_text(size = 30),
        plot.tag.position = "top",
        plot.tag = element_text(size = 30)) +
  panel_border(color = "black")

# Extinction ####
g2 <-
  ggplot(ceno_results_table[which(ceno_results_table$ind == "Gm_t38.0"),], aes(x = ind, y = values)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 2) +
  stat_summary(fun.data = data_summary, show.legend = F,size = 2, linewidth = 5, col = '#e34a33') +
  coord_cartesian(ylim = c(-2,2)) +
  ylab("Gm") +
  theme(axis.title.y = element_text(size = 30, face = "bold",
                                    margin = margin(r = 20)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text = element_text(size = 30),
        plot.title = element_text(size = 30),
        plot.tag.position = "top",
        plot.tag = element_text(size = 30)) +
  panel_border(color = "black")

pdf(file = "./Figs/SUPP/FigS15.pdf", width = 12, height = 6)
grid.arrange(g1, g2)
dev.off()