install:
	R CMD INSTALL .

check:
	R CMD check .

<<<<<<< HEAD
test:
	Rscript -e "devtools::test()"

docs:
	Rscript -e "devtools::document()"
=======
docs:
	Rscript -e "devtools::document()"
>>>>>>> d6f1f11 (Update docs.)
