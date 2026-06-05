
this_file <- rstudioapi::getSourceEditorContext()$path
 setwd(dirname(this_file))
 setwd("..")

data <- read.csv("./data/combined_data.csv", sep = ";")
