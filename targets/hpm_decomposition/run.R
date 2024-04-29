# This is a helper script to run the pipeline.
# Choose how to execute the pipeline below.
# See https://books.ropensci.org/targets/hpc.html
# to learn about your options.

# renv
# ---todo: install hpmdpredict from github and then make a snapshot with renv
renv::hydrate()
renv::repair()

targets::tar_make()
# targets::tar_make_clustermq(workers = 2) # nolint
# targets::tar_make_future(workers = 2) # nolint

system2(getOption("pdfviewer", default=''), args=c("hpmd-paper.pdf"), wait = FALSE)
system2(getOption("pdfviewer", default=''), args=c("hpmd-supporting-info.pdf"), wait = FALSE)
