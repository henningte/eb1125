# This is a helper script to run the pipeline.
# Choose how to execute the pipeline below.
# See https://books.ropensci.org/targets/hpc.html
# to learn about your options.

# renv
renv::hydrate()
renv::repair()

# takes about 4 hours and 25 GB of RAM
targets::tar_make()

system2(getOption("pdfviewer", default=''), args=c("leaching-paper.pdf"), wait = FALSE)
system2(getOption("pdfviewer", default=''), args=c("leaching-supporting-info.pdf"), wait = FALSE)



