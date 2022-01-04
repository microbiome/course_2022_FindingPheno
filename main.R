# This code below builds the Rmd file into a book that is in html format
# Restart R to clean up the environment, make sure that you've saved everything
#.rs.restartR()
# Removes all objects
rm(list = ls())
authors <- "Leo Lahti, Tuomas Borman, Henrik Eckermann, Chouaib Benchraka"

library(bookdown)
# Render html files
render_book("index.Rmd", "bookdown::gitbook")

####################################
# Builds book in pdf format
# Restart R to clean up the environment, make sure that you've saved everything
#.rs.restartR()
# Removes all objects
rm(list = ls())
authors <- "Leo Lahti, Tuomas Borman, Henrik Eckermann, Chouaib Benchraka"

library(bookdown)
# Render pdf files
render_book("index.Rmd", "bookdown::pdf_book")


# After you have rendered the book locally, you can observe rendered
# html files that are located in docs subfolder.