all: partdiff-par

partdiff-seq:
	g++   -Ofast -fassociative-math  -o partdiff-seq.exe partdiff-par.c
partdiff-omp:
	g++   -Ofast -fassociative-math  -o partdiff-omp.exe -DVARIANTE_OPENMP -fopenmp partdiff-par.c
partdiff-mpi:
	mpic++ -Ofast -fassociative-math  -o partdiff-mpi.exe -DVARIANTE_MPI partdiff-par.c
partdiff-mpi-omp:
	mpic++ -Ofast -fassociative-math  -o partdiff-mpi-omp.exe -DVARIANTE_MPI -DVARIANTE_OPENMP -fopenmp partdiff-par.c
partdiff-par:partdiff-mpi
	mv partdiff-mpi.exe partdiff-par	
partdiff-par-hybrid:partdiff-mpi-omp
	mv partdiff-mpi-omp.exe partdiff-par-hybrid
clean:
	rm -f *.exe
	rm partdiff-par
	rm partdiff-par-hybrid
