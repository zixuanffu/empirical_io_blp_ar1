library(data.table)
dt <- fread("Data/In/blp.csv")
name <- colnames(dt)
name <- gsub("BLP.", "", name, fixed = TRUE)
name <- gsub(".", "_", name, fixed = TRUE)
colnames(dt) <- name
saveRDS(dt, "Data/Out/blp.rds")
