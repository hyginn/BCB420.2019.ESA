if (requireNamespace("biogridr", quietly = TRUE)) {
  devtools::install_github("NElnour/biogridr", force = TRUE)
  library(biogridr)
} else {
  library(biogridr)
}

# Load all required packages.
require(xlsx, quietly = TRUE)
require(readxl, quietly = TRUE)
require(dplyr, quietly = TRUE)
require(biomaRt, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(biogridr, quietly = TRUE)
require(visNetwork, quietly = TRUE)
