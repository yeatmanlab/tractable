install:
	R CMD INSTALL .

check:
	R CMD check .

<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> c2aa996 (Add test directive in makefile and fix up gam test.)
test:
	Rscript -e "devtools::test()"

docs:
	Rscript -e "devtools::document()"
=======
docs:
	Rscript -e "devtools::document()"
>>>>>>> d6f1f11 (Update docs.)
