df_toy <- data.frame(
  species = paste0("Species", 1:6),
  trait1 = c(1, 1, 1, 5, 5, 4),
  trait2 = c(2, 1, 3, 2, 1, 3),
  TS = c(10, 8, 12, 5, 9, 7),
  TE = c(2,  1,  6, 0, 4, 3),
  site = c(1, 2, 1, 2, 3, 3),
  group = rep("focal", 6))
