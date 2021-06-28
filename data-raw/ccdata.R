## Code to prepare `ccdata` dataset

# Get the raw data
raw_data <- utils::read.csv("https://raw.githubusercontent.com/Jo-Karl/ccpsyc/master/data/example.csv")

# Select countries of Lebanon and Taiwan
ccdata  <- subset(raw_data, country %in% c("LEB", "TAIW"))

# Center and scale the data
ccdata[,-1]  <- scale(ccdata[,-1], center = TRUE, scale = TRUE)

# Rename items for convenience
colnames(ccdata) <- c("country", "h1", "h2", "h3", "h4", "h5", "h6", "h7",
                      "v1", "v2", "v3", "v4", "v5")

utils::write.csv(ccdata, "data-raw/ccdata.csv")
usethis::use_data(ccdata, overwrite = TRUE)
