#!/bin/sh

Rscript -e "devtools::install_github('ms609/TreeSearch')" # Need version > 0.1.2
#Rscript -e " rmarkdown::render_site(encoding = 'UTF-8')"
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"
#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book')"
#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::epub_book')"
