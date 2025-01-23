# renv
renv::hydrate()
renv::repair()

targets::tar_make()

system2(getOption("pdfviewer", default=''), args=c("hpmd-paper.pdf"), wait = FALSE)
system2(getOption("pdfviewer", default=''), args=c("hpmd-supporting-info.pdf"), wait = FALSE)
