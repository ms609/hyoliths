#!/bin/sh

#tlmgr install collection-fontsrecommended
tlmgr install ebgaramond
Rscript -e "devtools::install_github('ropensci/rcrossref')" # Need version 0.8.1 or higher
Rscript -e " rmarkdown::render_site(encoding = 'UTF-8')"
#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"
#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book')"
#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::epub_book')"
