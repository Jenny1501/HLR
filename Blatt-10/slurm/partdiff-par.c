#define _POSIX_C_SOURCE 200809L
#include <inttypes.h>
#include <math.h>
#include <math.h>
#include <stdint.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#ifndef PI
#define PI 3.141592653589793
#endif
#define TWO_PI_SQUARE (2 * PI * PI)
#define MAX_INTERLINES 10240
#define MAX_ITERATION 200000
#define MAX_THREADS 1024
#define METH_GAUSS_SEIDEL 1
#define METH_JACOBI 2
#define FUNC_F0 1
#define FUNC_FPISIN 2
#define TERM_PREC 1
#define TERM_ITER 2
struct options {
	uint64_t number;		 /* Number of threads                              */
	uint64_t method;		 /* Gauss Seidel or Jacobi method of iteration     */
	uint64_t interlines;	 /* matrix size = interlines*8+9                   */
	uint64_t inf_func;		 /* inference function                             */
	uint64_t termination;	/* termination condition                          */
	uint64_t term_iteration; /* terminate if iteration number reached          */
	double   term_precision; /* terminate if precision reached                 */
};
struct calculation_arguments {
	uint64_t  N;			/* number of spaces between lines (lines=N+1) */
	uint64_t  num_matrices; /* number of matrices */
	double	h;			/* length of a space between two lines */
	double ***Matrix;		/* index matrix used for addressing M */
	double *  M;			/* two matrices with real values */
};
struct calculation_results {
	uint64_t m;
	uint64_t stat_iteration; /* number of current iteration */
	double   stat_precision; /* actual precision of all slaves in iteration */
};
void		AskParams (struct options *, int, char **);
static void displayStatistics (struct calculation_results const *results, struct options const *options);
/* ************************************************************************ */
/* Global variables */
/* ************************************************************************ */
int					 world_size; // wie viele Prozesse gibt es?
int					 world_rank; // der wieviele Prozess ist dies?
signed long long int memoryUsed; // korrekte Speicherberechnung
struct timeval		 start_time; /* time when program started */
struct timeval		 comp_time;  /* time when calculation completed */
static void usage (char *name) {
	printf ("Usage: %s [num] [method] [lines] [func] [term] [prec/iter]\n", name);
	printf ("\n");
	printf ("  - num:       number of threads (1 .. %d)\n", MAX_THREADS);
	printf ("  - method:    calculation method (1 .. 2)\n");
	printf ("                 %1d: Gauß-Seidel\n", METH_GAUSS_SEIDEL);
	printf ("                 %1d: Jacobi\n", METH_JACOBI);
	printf ("  - lines:     number of interlines (0 .. %d)\n", MAX_INTERLINES);
	printf ("                 matrixsize = (interlines * 8) + 9\n");
	printf ("  - func:      interference function (1 .. 2)\n");
	printf ("                 %1d: f(x,y) = 0\n", FUNC_F0);
	printf ("                 %1d: f(x,y) = 2 * pi^2 * sin(pi * x) * sin(pi * y)\n", FUNC_FPISIN);
	printf ("  - term:      termination condition ( 1.. 2)\n");
	printf ("                 %1d: sufficient precision\n", TERM_PREC);
	printf ("                 %1d: number of iterations\n", TERM_ITER);
	printf ("  - prec/iter: depending on term:\n");
	printf ("                 precision:  1e-4 .. 1e-20\n");
	printf ("                 iterations:    1 .. %d\n", MAX_ITERATION);
	printf ("\n");
	printf ("Example: %s 1 2 100 1 2 100 \n", name);
}
static int check_number (struct options *options) {
	return (options->number >= 1 && options->number <= MAX_THREADS);
}
static int check_method (struct options *options) {
	return (options->method == METH_GAUSS_SEIDEL || options->method == METH_JACOBI);
}
static int check_interlines (struct options *options) {
	return (options->interlines <= MAX_INTERLINES);
}
static int check_inf_func (struct options *options) {
	return (options->inf_func == FUNC_F0 || options->inf_func == FUNC_FPISIN);
}
static int check_termination (struct options *options) {
	return (options->termination == TERM_PREC || options->termination == TERM_ITER);
}
static int check_term_precision (struct options *options) {
	return (options->term_precision >= 1e-20 && options->term_precision <= 1e-3);
}
static int check_term_iteration (struct options *options) {
	return (options->term_iteration >= 1 && options->term_iteration <= MAX_ITERATION);
}
void AskParams (struct options *options, int argc, char **argv) {
	int ret;
	/*printf ("============================================================\n");
	 printf ("Program for calculation of partial differential equations.  \n");
	 printf ("============================================================\n");
	 printf ("(c) Dr. Thomas Ludwig, TU München.\n");
	 printf ("    Thomas A. Zochler, TU München.\n");
	 printf ("    Andreas C. Schmidt, TU München.\n");
	 printf ("============================================================\n");
	 printf ("\n");*/
	if (argc < 2) {
		/* ----------------------------------------------- */
		/* Get input: method, interlines, func, precision. */
		/* ----------------------------------------------- */
		do {
			printf ("\n");
			printf ("Select number of threads:\n");
			printf ("Number> ");
			fflush (stdout);
			ret = scanf ("%" SCNu64, &(options->number));
			while (getchar () != '\n')
				;
		} while (ret != 1 || !check_number (options));
		do {
			printf ("\n");
			printf ("Select calculation method:\n");
			printf ("  %1d: Gauß-Seidel.\n", METH_GAUSS_SEIDEL);
			printf ("  %1d: Jacobi.\n", METH_JACOBI);
			printf ("method> ");
			fflush (stdout);
			ret = scanf ("%" SCNu64, &(options->method));
			while (getchar () != '\n')
				;
		} while (ret != 1 || !check_method (options));
		do {
			printf ("\n");
			printf ("Matrixsize = Interlines*8+9\n");
			printf ("Interlines> ");
			fflush (stdout);
			ret = scanf ("%" SCNu64, &(options->interlines));
			while (getchar () != '\n')
				;
		} while (ret != 1 || !check_interlines (options));
		do {
			printf ("\n");
			printf ("Select interference function:\n");
			printf (" %1d: f(x,y)=0.\n", FUNC_F0);
			printf (" %1d: f(x,y)=2pi^2*sin(pi*x)sin(pi*y).\n", FUNC_FPISIN);
			printf ("interference function> ");
			fflush (stdout);
			ret = scanf ("%" SCNu64, &(options->inf_func));
			while (getchar () != '\n')
				;
		} while (ret != 1 || !check_inf_func (options));
		do {
			printf ("\n");
			printf ("Select termination:\n");
			printf (" %1d: sufficient precision.\n", TERM_PREC);
			printf (" %1d: number of iterations.\n", TERM_ITER);
			printf ("termination> ");
			fflush (stdout);
			ret = scanf ("%" SCNu64, &(options->termination));
			while (getchar () != '\n')
				;
		} while (ret != 1 || !check_termination (options));
		if (options->termination == TERM_PREC) {
			do {
				printf ("\n");
				printf ("Select precision:\n");
				printf ("  Range: 1e-4 .. 1e-20.\n");
				printf ("precision> ");
				fflush (stdout);
				ret = scanf ("%lf", &(options->term_precision));
				while (getchar () != '\n')
					;
			} while (ret != 1 || !check_term_precision (options));
			options->term_iteration = MAX_ITERATION;
		} else if (options->termination == TERM_ITER) {
			do {
				printf ("\n");
				printf ("Select number of iterations:\n");
				printf ("  Range: 1 .. %d.\n", MAX_ITERATION);
				printf ("Iterations> ");
				fflush (stdout);
				ret = scanf ("%" SCNu64, &(options->term_iteration));
				while (getchar () != '\n')
					;
			} while (ret != 1 || !check_term_iteration (options));
			options->term_precision = 0;
		}
	} else {
		if (argc < 7 || strcmp (argv[1], "-h") == 0 || strcmp (argv[1], "-?") == 0) {
			usage (argv[0]);
			exit (0);
		}
		ret = sscanf (argv[1], "%" SCNu64, &(options->number));
		if (ret != 1 || !check_number (options)) {
			usage (argv[0]);
			exit (1);
		}
		ret = sscanf (argv[2], "%" SCNu64, &(options->method));
		if (ret != 1 || !check_method (options)) {
			usage (argv[0]);
			exit (1);
		}
		ret = sscanf (argv[3], "%" SCNu64, &(options->interlines));
		if (ret != 1 || !check_interlines (options)) {
			usage (argv[0]);
			exit (1);
		}
		ret = sscanf (argv[4], "%" SCNu64, &(options->inf_func));
		if (ret != 1 || !check_inf_func (options)) {
			usage (argv[0]);
			exit (1);
		}
		ret = sscanf (argv[5], "%" SCNu64, &(options->termination));
		if (ret != 1 || !check_termination (options)) {
			usage (argv[0]);
			exit (1);
		}
		if (options->termination == TERM_PREC) {
			ret						= sscanf (argv[6], "%lf", &(options->term_precision));
			options->term_iteration = MAX_ITERATION;
			if (ret != 1 || !check_term_precision (options)) {
				usage (argv[0]);
				exit (1);
			}
		} else {
			ret						= sscanf (argv[6], "%" SCNu64, &(options->term_iteration));
			options->term_precision = 0;
			if (ret != 1 || !check_term_iteration (options)) {
				usage (argv[0]);
				exit (1);
			}
		}
	}
}
static void freeMemory (void *ptr, size_t size) {
	free (ptr);
}
static void *allocateMemory (size_t size) {
	memoryUsed += size;
	void *p;
	if ((p = malloc (size)) == NULL) {
		printf ("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
		exit (1);
	}
	return p;
}
//////////////////////////////////////////////////////////////////////////////////////////
// Includes wenn benötigt
//////////////////////////////////////////////////////////////////////////////////////////
#ifdef VARIANTE_MPI
#include <mpi.h>
#endif
#ifdef VARIANTE_OPENMP
#include <omp.h>
#endif
//////////////////////////////////////////////////////////////////////////////////////////
// Makros für die berechnung der Matrixaufteilung
//////////////////////////////////////////////////////////////////////////////////////////
#define TOTAL_SIZE (N + 1)
#define AVERAGE_HEIGHT (TOTAL_SIZE / world_size + (TOTAL_SIZE % world_size > 0) + 1)
#define FIRST_ROW(rank) (rank * (AVERAGE_HEIGHT - 2))
#define LAST_HEIGHT(rank) (TOTAL_SIZE - FIRST_ROW (rank))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
static void calculateGauss (struct calculation_arguments *arguments, struct calculation_results *results, struct options const *options) {
//////////////////////////////////////////////////////////////////////////////////////////
// Variablen-Definitionen + einige Initialisierungen ->
//////////////////////////////////////////////////////////////////////////////////////////
#ifdef VARIANTE_MPI
	int	termination_send		= 0; // the last process should send termination to process 0 just once
	double should_terminate		= 1; // terminate if negative
	double backup				= 0;
	double backup2				= 0;
	double tmp_should_terminate = 1;
#endif
	results->m					= 0;
	results->stat_iteration		= 0;
	results->stat_precision		= 0;
	arguments->num_matrices		= 1;
	int const	N				= (options->interlines * 8) + 9 - 1;
	double const h				= 1.0 / N;
	int			 i				= 0; // schleifenvar
	int			 j				= 0; // schleifenvar
	double		 star			= 0;
	double		 residuum		= 0;
	double		 maxresiduum	= 0;
	double		 pih			= 0;
	double		 fpisin			= 0;
	double		 fpisin_i		= 0;
	int			 width			= TOTAL_SIZE;
	int			 height			= AVERAGE_HEIGHT;		  // wenn mpi verwendet wird, durchschnittliche Zeilenzahl berechnen
	int			 firstRowNumber = FIRST_ROW (world_rank); // da mpi die Zeilen verteilt muss die Randfunction sowie die ausgabe der Werte den offset kennen.
	int			 term_iteration = options->term_iteration;
	height						= ((world_rank == world_size - 1) ? (LAST_HEIGHT (world_rank)) : height); // der letzte rank hat möglicherweise weniger zeilen als alle anderen
	arguments->M				= (double *) allocateMemory (arguments->num_matrices * width * height * sizeof (double));
	double **Matrix				= (double **) allocateMemory (height * sizeof (double *));
	arguments->Matrix			= &Matrix;
	int height_lower			= 1;
	int height_upper			= height - 1;
	//////////////////////////////////////////////////////////////////////////////////////////
	// Matrix Initialisierung ->
	//////////////////////////////////////////////////////////////////////////////////////////
	for (i = 0; i < height; i++) {
		Matrix[i] = arguments->M + (i * width);
		for (j = 0; j < width; j++) {
			Matrix[i][j] = 0.0;
		}
		if (options->inf_func == FUNC_F0) {
			Matrix[i][0]		 = 1.0 - (h * (firstRowNumber + i)); // offset beachten
			Matrix[i][width - 1] = h * (firstRowNumber + i);
		}
	}
	if (options->inf_func == FUNC_F0) {
		if (world_rank == 0) {
			for (i = 0; i <= N; i++) {
				Matrix[0][i] = 1.0 - (h * i);
			}
			Matrix[0][N] = 0.0;
		}
		if (world_rank == world_size - 1) {
			for (i = 0; i <= N; i++) {
				Matrix[height - 1][i] = h * i;
			}
			Matrix[height - 1][0] = 0.0;
		}
	}
	if (options->inf_func == FUNC_FPISIN) {
		pih	= PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	// Begin der Zeitmessung ->
	//////////////////////////////////////////////////////////////////////////////////////////
	gettimeofday (&start_time, NULL);
	results->stat_iteration = -(world_rank > 0);
	//////////////////////////////////////////////////////////////////////////////////////////
	// Hauptschleife ->
	//////////////////////////////////////////////////////////////////////////////////////////
	while (term_iteration > 0) {
		if (0 == world_rank) {
			maxresiduum = 0;
		}
//////////////////////////////////////////////////////////////////////////////////////////
// MPI synchronisierung empfange von oben->
//////////////////////////////////////////////////////////////////////////////////////////
#ifdef VARIANTE_MPI
		else {
			MPI_Recv (Matrix[0], width, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			maxresiduum		 = Matrix[0][N];
			should_terminate = Matrix[0][0];
			if (options->termination == TERM_PREC) {
				if (should_terminate < 0) {
					term_iteration = 0;
				}
			}
		}
#endif
		//////////////////////////////////////////////////////////////////////////////////////////
		// erste Zeile berechnen ->
		//////////////////////////////////////////////////////////////////////////////////////////
		if (options->inf_func == FUNC_FPISIN) {
			fpisin_i = fpisin * sin (pih * (double) (firstRowNumber + height_lower)); // offset beachten
		}
		for (j = 1; j < N; j++) {
			star = 0.25 * (Matrix[height_lower - 1][j] + Matrix[height_lower][j - 1] + Matrix[height_lower][j + 1] + Matrix[height_lower + 1][j]);
			if (options->inf_func == FUNC_FPISIN) {
				star += fpisin_i * sin (pih * (double) j);
			}
			if (options->termination == TERM_PREC || term_iteration == 1) {
				residuum = Matrix[height_lower][j] - star;
				residuum = (residuum < 0) ? -residuum : residuum;
				if (residuum > maxresiduum) {
					maxresiduum = (residuum > maxresiduum) ? residuum : maxresiduum;
				}
			}
			Matrix[height_lower][j] = star;
		}
//////////////////////////////////////////////////////////////////////////////////////////
// MPI synchronisierung sende nach oben->
//////////////////////////////////////////////////////////////////////////////////////////
#ifdef VARIANTE_MPI
		// TODO :: nach unten ... letzte iteration nicht abspalten
		if ((world_rank > 0) && (should_terminate > 0)) {
			MPI_Send (Matrix[1], width, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD);
		}
#endif
		//////////////////////////////////////////////////////////////////////////////////////////
		// andere Zeilen berechnen ->
		//////////////////////////////////////////////////////////////////////////////////////////
		for (i = height_lower + 1; i < height_upper; i++) {
			if (options->inf_func == FUNC_FPISIN) {
				fpisin_i = fpisin * sin (pih * (double) (firstRowNumber + i)); // offset beachten
			}
			for (j = 1; j < N; j++) {
				star = 0.25 * (Matrix[i - 1][j] + Matrix[i][j - 1] + Matrix[i][j + 1] + Matrix[i + 1][j]);
				if (options->inf_func == FUNC_FPISIN) {
					star += fpisin_i * sin (pih * (double) j);
				}
				if (options->termination == TERM_PREC || term_iteration == 1) {
					residuum = Matrix[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					if (residuum > maxresiduum) {
						maxresiduum = (residuum > maxresiduum) ? residuum : maxresiduum;
					}
				}
				Matrix[i][j] = star;
			}
		}
#ifdef VARIANTE_MPI
		if ((world_size > 1) && (world_rank < world_size - 1)) {
			//////////////////////////////////////////////////////////////////////////////////////////
			// MPI synchronisierung sende nach unten->
			//////////////////////////////////////////////////////////////////////////////////////////
			backup				  = Matrix[height - 2][N];
			backup2				  = Matrix[height - 2][0];
			Matrix[height - 2][N] = maxresiduum;
			Matrix[height - 2][0] = should_terminate;
			MPI_Send (Matrix[height - 2], width, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD);
			Matrix[height - 2][N] = backup;
			Matrix[height - 2][0] = backup2;
			//////////////////////////////////////////////////////////////////////////////////////////
			// MPI synchronisierung empfange von unten->
			//////////////////////////////////////////////////////////////////////////////////////////
			if (should_terminate > 0) {
				MPI_Recv (Matrix[height - 1], width, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
#endif
		results->stat_iteration++;
		if (options->termination == TERM_PREC) {
#ifdef VARIANTE_MPI
			if (world_size > 1) {
				if (should_terminate < 0) {
					term_iteration = 0;
				} else if ((world_rank == 0) && (results->stat_iteration >= world_size - 1)) {
					MPI_Recv (&should_terminate, 1, MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				} else if ((world_rank == world_size - 1) && (termination_send == 0)) {
					if (maxresiduum < options->term_precision) {
						tmp_should_terminate = -1;
						termination_send	 = 1;
					}
					MPI_Send (&tmp_should_terminate, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
				}
			} else {
				if (maxresiduum < options->term_precision) {
					term_iteration = 0;
				}
			}
#else
			if (maxresiduum < options->term_precision) {
				term_iteration = 0;
			}
#endif
		} else if (options->termination == TERM_ITER) {
			term_iteration--;
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	// Matrix Re-Initialisierung an den oberen Ecken ->
	//////////////////////////////////////////////////////////////////////////////////////////
	if (options->inf_func == FUNC_F0) {
		Matrix[0][0]		 = 1.0 - (h * firstRowNumber); // offset beachten
		Matrix[0][width - 1] = h * firstRowNumber;
	} else {
		Matrix[0][0]		 = 0.0;
		Matrix[0][width - 1] = 0.0;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	// Senden des Maxresiduums zur korrekten ausgabe ->
	//////////////////////////////////////////////////////////////////////////////////////////
    #ifdef VARIANTE_MPI
	if (world_size > 1) {
		if (world_rank == 0) {
			MPI_Recv (&results->stat_precision, 1, MPI_DOUBLE, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		} else if (world_rank == world_size - 1) {
			MPI_Send (&maxresiduum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}
    #endif
	//////////////////////////////////////////////////////////////////////////////////////////
	// Ausgabe der Matrix ->
	//////////////////////////////////////////////////////////////////////////////////////////
	results->m = 0;
	gettimeofday (&comp_time, NULL);
	displayStatistics (results, options);
	int const interlines = options->interlines;
	int		  resultRanks[9]; // nummer des Prozesses der die reihe hat
	int		  resultRows[9];  // reihe innerhalb des Prozesses
	for (int i = 0; i < 9; i++) {
		int targetRow = i * (interlines + 1);
#ifdef VARIANTE_MPI
		resultRanks[i] = MIN ((targetRow) / (AVERAGE_HEIGHT - 2), world_size - 1);
		resultRows[i]  = targetRow - FIRST_ROW (resultRanks[i]);
		if ((resultRanks[i] > 0) && (resultRanks[i] == world_rank)) {
			MPI_Send (Matrix[resultRows[i]], TOTAL_SIZE, MPI_DOUBLE, 0, i * 10, MPI_COMM_WORLD);
		}
#else
		resultRanks[i] = 0;
		resultRows[i]  = targetRow;
#endif
	}
	if (world_rank == 0) {
		int x, y;
		printf ("Matrix:\n");
		for (y = 0; y < 9; y++) {
			if (resultRanks[y] == 0) {
				for (x = 0; x < 9; x++) {
					printf ("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
				}
			} else {
#ifdef VARIANTE_MPI
				MPI_Recv (Matrix[0], TOTAL_SIZE, MPI_DOUBLE, resultRanks[y], y * 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (x = 0; x < 9; x++) {
					printf ("%7.4f", Matrix[0][x * (interlines + 1)]);
				}
#endif
			}
			printf ("\n");
		}
		fflush (stdout);
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	// Freigeben der Matrix ->
	//////////////////////////////////////////////////////////////////////////////////////////
	freeMemory (Matrix, height * sizeof (double *));
	freeMemory (arguments->M, arguments->num_matrices * width * height * sizeof (double));
}
static void calculateJacobi (struct calculation_arguments *arguments, struct calculation_results *results, struct options const *options) {
#if defined(VARIANTE_OPENMP)
	omp_set_num_threads (options->number); // setzen der vorgegebenen Threadanzahl
#endif
	//////////////////////////////////////////////////////////////////////////////////////////
	// Variablen-Definitionen + einige Initialisierungen ->
	//////////////////////////////////////////////////////////////////////////////////////////
	results->m					= 0;
	results->stat_iteration		= 0;
	results->stat_precision		= 0;
	arguments->num_matrices		= 2;
	int const	N				= (options->interlines * 8) + 9 - 1;
	double const h				= 1.0 / N;
	int			 g				= 0; // schleifenvar
	int			 i				= 0; // schleifenvar
	int			 j				= 0; // schleifenvar
	int			 m1				= 0;
	int			 m2				= 1;
	double		 star			= 0;
	double		 residuum		= 0;
	double		 maxresiduum	= 0;
	double		 pih			= 0;
	double		 fpisin			= 0;
	double		 fpisin_i		= 0;
	int			 width			= TOTAL_SIZE;
	int			 height			= AVERAGE_HEIGHT;		  // wenn mpi verwendet wird, durchschnittliche Zeilenzahl berechnen
	int			 firstRowNumber = FIRST_ROW (world_rank); // da mpi die Zeilen verteilt muss die Randfunction sowie die ausgabe der Werte den offset kennen.
	int			 term_iteration = options->term_iteration;
	height						= ((world_rank == world_size - 1) ? (LAST_HEIGHT (world_rank)) : height); // der letzte rank hat möglicherweise weniger zeilen als alle anderen
	arguments->M				= (double *) allocateMemory (2 * width * height * sizeof (double));
	arguments->Matrix			= (double ***) allocateMemory (2 * sizeof (double **));
	double ***Matrix			= arguments->Matrix;
	int		  height_lower		= 1;
	int		  height_upper		= height - 1;
	//////////////////////////////////////////////////////////////////////////////////////////
	// Matrix Initialisierung ->
	//////////////////////////////////////////////////////////////////////////////////////////
	for (g = 0; g < 2; g++) {
		Matrix[g] = (double **) allocateMemory (height * sizeof (double *));
#if defined(VARIANTE_OPENMP)
#pragma omp parallel for
		for (i = 0; i < height; i++) // parallele initialisierung (first touch localität)
#else
		for (i = 0; i < height; i++)
#endif
		{
			arguments->Matrix[g][i] = arguments->M + (g * height * width) + (i * width);
			for (j = 0; j < width; j++) {
				Matrix[g][i][j] = 0.0;
			}
			if (options->inf_func == FUNC_F0) {
				Matrix[g][i][0]			= 1.0 - (h * (firstRowNumber + i)); // offset beachten
				Matrix[g][i][width - 1] = h * (firstRowNumber + i);
			}
		}
		if (options->inf_func == FUNC_F0) {
			if (world_rank == 0) {
				for (i = 0; i <= N; i++) {
					Matrix[g][0][i] = 1.0 - (h * i);
				}
				Matrix[g][0][N] = 0.0;
			}
			if (world_rank == world_size - 1) {
				for (i = 0; i <= N; i++) {
					Matrix[g][height - 1][i] = h * i;
				}
				Matrix[g][height - 1][0] = 0.0;
			}
		}
	}
	if (options->inf_func == FUNC_FPISIN) {
		pih	= PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	// Begin der Zeitmessung ->
	//////////////////////////////////////////////////////////////////////////////////////////
	gettimeofday (&start_time, NULL);
	//////////////////////////////////////////////////////////////////////////////////////////
	// Hauptschleife ->
	//////////////////////////////////////////////////////////////////////////////////////////
	while (term_iteration > 0) {
		double **Matrix_Out = arguments->Matrix[m1];
		double **Matrix_In  = arguments->Matrix[m2];
		maxresiduum			= 0;
		if (options->termination == TERM_PREC || term_iteration == 1) {
//////////////////////////////////////////////////////////////////////////////////////////
// Residuum nur berechnen wenn nötig ->
//////////////////////////////////////////////////////////////////////////////////////////
#if defined(VARIANTE_OPENMP)
#pragma omp parallel for reduction(max : maxresiduum) private(i, j, star, residuum, fpisin_i)
			for (i = height_lower; i < height_upper; i++)
#else
			for (i = height_lower; i < height_upper; i++)
#endif
			{
				if (options->inf_func == FUNC_FPISIN) {
					fpisin_i = fpisin * sin (pih * (double) (firstRowNumber + i)); // offset beachten
				}
				for (j = 1; j < N; j++) {
					star = 0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] + Matrix_In[i][j + 1] + Matrix_In[i + 1][j]);
					if (options->inf_func == FUNC_FPISIN) {
						star += fpisin_i * sin (pih * (double) j);
					}
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					if (residuum > maxresiduum) {
						maxresiduum = (residuum > maxresiduum) ? residuum : maxresiduum;
					}
					Matrix_Out[i][j] = star;
				}
			}
		} else {
//////////////////////////////////////////////////////////////////////////////////////////
// Kein Residuum berechnen ->
//////////////////////////////////////////////////////////////////////////////////////////
#if defined(VARIANTE_OPENMP)
#pragma omp parallel for private(i, j, star, residuum, fpisin_i)
			for (i = height_lower; i < height_upper; i++)
#else
			for (i	 = height_lower; i < height_upper; i++)
#endif
			{
				if (options->inf_func == FUNC_FPISIN) {
					fpisin_i = fpisin * sin (pih * (double) (firstRowNumber + i));
				}
				for (j = 1; j < N; j++) {
					star = 0.25 * (Matrix_In[i - 1][j] + Matrix_In[i][j - 1] + Matrix_In[i][j + 1] + Matrix_In[i + 1][j]);
					if (options->inf_func == FUNC_FPISIN) {
						star += fpisin_i * sin (pih * (double) j);
					}
					Matrix_Out[i][j] = star;
				}
			}
		}
		if (world_size > 1) {
//////////////////////////////////////////////////////////////////////////////////////////
// MPI synchronisierung ->
//////////////////////////////////////////////////////////////////////////////////////////
#ifdef VARIANTE_MPI
			if (world_rank % 2) {
				if (world_rank < world_size - 1) {
					MPI_Send (Matrix_Out[height - 2], width, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD);
					MPI_Recv (Matrix_Out[height - 1], width, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				if (world_rank > 0) {
					MPI_Send (Matrix_Out[1], width, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD);
					MPI_Recv (Matrix_Out[0], width, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
			} else {
				if (world_rank > 0) {
					MPI_Recv (Matrix_Out[0], width, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Send (Matrix_Out[1], width, MPI_DOUBLE, world_rank - 1, 0, MPI_COMM_WORLD);
				}
				if (world_rank < world_size - 1) {
					MPI_Recv (Matrix_Out[height - 1], width, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Send (Matrix_Out[height - 2], width, MPI_DOUBLE, world_rank + 1, 0, MPI_COMM_WORLD);
				}
			}
			if (options->termination == TERM_PREC || term_iteration == 1) {
				MPI_Allreduce (&maxresiduum, &results->stat_precision, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			}
#endif
		} else {
			//////////////////////////////////////////////////////////////////////////////////////////
			// Andernfalls muss Maxresiduum direkt gesetzt werden ->
			//////////////////////////////////////////////////////////////////////////////////////////
			results->stat_precision = maxresiduum;
		}
		results->stat_iteration++;
		i  = m1;
		m1 = m2;
		m2 = i;
		if (options->termination == TERM_PREC) {
			if (results->stat_precision < options->term_precision) {
				term_iteration = 0;
			}
		} else if (options->termination == TERM_ITER) {
			term_iteration--;
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	// Ausgabe der Matrix ->
	//////////////////////////////////////////////////////////////////////////////////////////
	results->m = m2;
	gettimeofday (&comp_time, NULL);
	displayStatistics (results, options);
	int const interlines = options->interlines;
	int		  resultRanks[9]; // nummer des Prozesses der die reihe hat
	int		  resultRows[9];  // reihe innerhalb des Prozesses
	for (int i = 0; i < 9; i++) {
		int targetRow = i * (interlines + 1);
#ifdef VARIANTE_MPI
		resultRanks[i] = MIN ((targetRow) / (AVERAGE_HEIGHT - 2), world_size - 1);
		resultRows[i]  = targetRow - FIRST_ROW (resultRanks[i]);
		if ((resultRanks[i] > 0) && (resultRanks[i] == world_rank)) {
			MPI_Send (Matrix[results->m][resultRows[i]], TOTAL_SIZE, MPI_DOUBLE, 0, i * 10, MPI_COMM_WORLD);
		}
#else
		resultRanks[i] = 0;
		resultRows[i]  = targetRow;
#endif
	}
	if (world_rank == 0) {
		int x, y;
		printf ("Matrix:\n");
		for (y = 0; y < 9; y++) {
			if (resultRanks[y] == 0) {
				for (x = 0; x < 9; x++) {
					printf ("%7.4f", Matrix[results->m][y * (interlines + 1)][x * (interlines + 1)]);
				}
			} else {
#ifdef VARIANTE_MPI
				MPI_Recv (Matrix[results->m][0], TOTAL_SIZE, MPI_DOUBLE, resultRanks[y], y * 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (x = 0; x < 9; x++) {
					printf ("%7.4f", Matrix[results->m][0][x * (interlines + 1)]);
				}
#endif
			}
			printf ("\n");
		}
		fflush (stdout);
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	// Freigeben der Matrix ->
	//////////////////////////////////////////////////////////////////////////////////////////
	for (g = 0; g < 2; g++) {
		freeMemory (Matrix[g], height * sizeof (double *));
	}
	freeMemory (arguments->Matrix, 2 * sizeof (double **));
	freeMemory (arguments->M, 2 * width * height * sizeof (double));
}
/* ************************************************************************ */
/* calculate: solves the equation */
/* ************************************************************************ */
static void calculate (struct calculation_arguments *arguments, struct calculation_results *results, struct options const *options) {
	if (options->method == METH_JACOBI) {
		calculateJacobi (arguments, results, options);
	} else {
		calculateGauss (arguments, results, options);
	}
}
/* ************************************************************************ */
/* displayStatistics: displays some statistics about the calculation */
/* ************************************************************************ */
static void displayStatistics (struct calculation_results const *results, struct options const *options) {
	signed long long int allMemoryUsed = memoryUsed;
	signed long long int maxMemoryUsed = memoryUsed;
#ifdef VARIANTE_MPI
	MPI_Allreduce (&memoryUsed, &allMemoryUsed, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce (&memoryUsed, &maxMemoryUsed, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
#endif
	if (world_rank == 0) {
		int	N	= ((options->interlines * 8) + 9 - 1);
		double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;
		printf ("Berechnungszeit:    %f s \n", time);
		int width  = TOTAL_SIZE;
		int height = AVERAGE_HEIGHT;
		if (allMemoryUsed < 1024) {
			printf ("Speicher Gesamt:    %f B\n", allMemoryUsed * 1.0);
			printf ("Speicher Prozess:   %f B\n", maxMemoryUsed * 1.0);
		} else if (allMemoryUsed < 1024 * 1024) {
			printf ("Speicher Gesamt:    %f KiB\n", allMemoryUsed / 1024.0);
			printf ("Speicher Prozess:   %f KiB\n", maxMemoryUsed / 1024.0);
		} else if (allMemoryUsed < 1024 * 1024 * 1024) {
			printf ("Speicher Gesamt:    %f MiB\n", allMemoryUsed / 1024.0 / 1024.0);
			printf ("Speicher Prozess:   %f MiB\n", maxMemoryUsed / 1024.0 / 1024.0);
		} else {
			printf ("Speicher Gesamt:    %f GiB\n", allMemoryUsed / 1024.0 / 1024.0 / 1024.0);
			printf ("Speicher Prozess:   %f GiB\n", maxMemoryUsed / 1024.0 / 1024.0 / 1024.0);
		}
		printf ("Processanzahl:      %d \n", world_size);
		printf ("Threadanzahl:       %d \n", options->number);
		printf ("Taskanzahl:         %d \n", options->number * world_size);
		printf ("Berechnungsmethode: ");
		if (options->method == METH_GAUSS_SEIDEL) {
			printf ("Gauß-Seidel");
		} else if (options->method == METH_JACOBI) {
			printf ("Jacobi");
		}
		printf ("\n");
		printf ("Interlines:         %" PRIu64 "\n", options->interlines);
		printf ("Stoerfunktion:      ");
		if (options->inf_func == FUNC_F0) {
			printf ("f(x,y)=0");
		} else if (options->inf_func == FUNC_FPISIN) {
			printf ("f(x,y)=2pi^2*sin(pi*x)sin(pi*y)");
		}
		printf ("\n");
		printf ("Terminierung:       ");
		if (options->termination == TERM_PREC) {
			printf ("Hinreichende Genaugkeit");
		} else if (options->termination == TERM_ITER) {
			printf ("Anzahl der Iterationen");
		}
		printf ("\n");
		printf ("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
		printf ("Norm des Fehlers:   %e\n", results->stat_precision);
		printf ("\n");
	}
}
/* ************************************************************************ */
/* main */
/* ************************************************************************ */
int main (int argc, char **argv) {
#ifdef VARIANTE_MPI
	MPI_Init (&argc, &argv);
	MPI_Comm_size (MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank (MPI_COMM_WORLD, &world_rank);
#else
	world_size		   = 1;
	world_rank		   = 0;
#endif
	memoryUsed = 0;
	struct options				 options;
	struct calculation_arguments arguments;
	struct calculation_results   results;
	AskParams (&options, argc, argv);
	calculate (&arguments, &results, &options);
#ifdef VARIANTE_MPI
	MPI_Finalize ();
#endif
	return 0;
}
