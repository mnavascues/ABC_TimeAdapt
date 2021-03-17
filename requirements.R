CRAN_repo = "https://cloud.r-project.org/" # change if needed
if (!requireNamespace("remotes"))
  install.packages("remotes", repos=CRAN_repo)
remotes::install_version("argparser", "0.4", repos=CRAN_repo)
remotes::install_version("extraDistr", "1.8.11", INSTALL_opts = '--no-lock', repos=CRAN_repo)
remotes::install_version("rcarbon", "1.2.0", repos=CRAN_repo)
remotes::install_version("testthat", "2.1.1", repos=CRAN_repo)

