# Extract temperature curve of Zachos 2008 from RPANDA datasets and apply kernel smoothing to 0.1 myr scale

rm(list = ls());gc()

library(RPANDA)
library(zoo)

data("InfTemp")

k_smooth <- ksmooth(x = InfTemp$Age, y = InfTemp$Temperature, x.points = seq(from = 0, to = max(InfTemp$Age), by = 0.1), bandwidth = 0.1)
k_smooth_df <- data.frame(k_smooth)
colnames(k_smooth_df) <- c("time", "temperature")
k_smooth_df$temperature <- na.approx(k_smooth_df$temperature)

plot(x = InfTemp$Age, y = InfTemp$Temperature, type = "l", xlim = c(37.2, 0), ylim = c(0, 15), xlab = c("Age (Myr)"), ylab = expression(paste("Temperature", " ", "(", delta, ")")))
lines(x = k_smooth_df$time, y = k_smooth_df$temperature, col = "blue")

dir.create("./Temperature/", recursive = T)
write.table(x = k_smooth_df[which(k_smooth_df$time <= 37.9),], file = c("./Temperature/temperature.txt"), row.names = F, quote = F)

