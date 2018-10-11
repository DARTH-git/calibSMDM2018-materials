
library(rmarkdown)

# tutorial1
render("tutorials/random_search_3state.Rmd",
       knit_root_dir = "..",
       output_file = "random_search_3state.html")
