######################################################################
#                          Website Creation
######################################################################

#FMestre
#04/08/2026

#Load library
library(pkgdown)

# Run once to configure your package to use and deploy pkgdown
usethis::use_pkgdown_github_pages()

# Preview your site locally before publishing
pkgdown::build_site()

#usethis::use_pkgdown_github_pages()
