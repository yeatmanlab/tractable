install:
	R CMD INSTALL .

check:
	R CMD check .

test:
	Rscript -e "devtools::test()"

docs:
	Rscript -e "devtools::document()"

pkgdown: 

