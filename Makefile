install:
	R CMD INSTALL .

check:
	R CMD check .

test:
	Rscript -e "devtools::test()"

docs:
	Rscript -e "devtools::document()"

examples:
	Rscript -e "devtools::build_rmd('vignettes/tractr-single-bundle.Rmd')"
	Rscript -e "devtools::build_rmd('vignettes/gratia-plots.Rmd')"

clean:
	rm vignettes/*.html
	rm vignettes/*.png