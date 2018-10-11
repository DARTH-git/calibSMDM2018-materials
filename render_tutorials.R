
library(rmarkdown)

# tutorial1
render("tutorials/random_search_3state.Rmd",
       knit_root_dir = "..",
       output_file = "random_search_3state.html")

# this extracts the R code into a separate script
purl("tutorials/random_search_3state.Rmd",
     output = "tutorials_scripts/random_search_3state.R")
