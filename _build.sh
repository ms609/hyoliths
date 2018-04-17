#!/bin/sh

tlmgr install collection-fontsrecommended
tlmgr install collection-fontsextra
Rscript -e "devtools::install_github('ropensci/rcrossref')"
Rscript -e "bookdown::render_book('index.Rmd')"
#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"
#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book')"
#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::epub_book')"
