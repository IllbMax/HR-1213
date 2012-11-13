/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                TU Muenchen - Institut fuer Informatik                  **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**            JK und andere  besseres Timing, FLOP Berechnung             **/
/**                                                                        **/
/** File:      partdiff-seq.c                                              **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauss-Seidel and   **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <pthread.h>
#include "partdiff-posix.h"

struct calculation_arguments
{
	int     N;              /* number of spaces between lines (lines=N+1)     */
	int     num_matrices;   /* number of matrices                             */
	double  h;              /* length of a space between two lines            */
	double  ***Matrix;      /* index matrix used for addressing M             */
	double  *M;             /* two matrices with real values                  */
};

struct calculation_results
{
	int     m;
	int     stat_iteration; /* number of current iteration                    */
	double  stat_precision; /* actual precision of all slaves in iteration    */
};



/*
 * Argumente die dem pthread beim Erstellen übergeben werden
 * 
 */
struct thread_args
{
	struct calculation_arguments const* arguments;
	struct options const* options;
	
	int from, to; /* zeile die der Thread berechnen soll */
	
	double maxresiduum; /* das aktuelle maxresiduum des Zeilen blocks des threads */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */


/* pthread Variablen */
pthread_barrier_t barrier; 	/* Barrier, die zum synchronisieren benutzt wird */
bool run; 			/* bei false beenden die Threads */


/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("\n\nSpeicherprobleme!\n");
		/* exit program */
		exit(1);
	}

	return p;
}


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}
/*
 * Init. die pthread variablen (globale und threadargs-array) 
 *  */
static void initPthreadVariables(struct thread_args** thread_args, pthread_t** threads, struct calculation_arguments* arguments, struct options const* options)
{
	int i, nextline;
	struct thread_args* ta;
	int n = options->number; /*Anzahl der Threads*/
	int N = arguments->N; /* Dim der Matrix */
	int linesPerThread; /* Zeilen pro thread */
	int overhead;       /* Rest bei Division */
	
	run = true;
	pthread_barrier_init(&barrier, NULL, n);
	*threads = allocateMemory(sizeof(pthread_t) * (n));
	
	ta = allocateMemory(sizeof(struct thread_args) * n);
	
	
	/* Zeilen pro thread bestimmen mit N-1 da der Rand nicht mit brechnet werden muss*/
	linesPerThread = (N-1) / n;
	overhead = (N-1) % n;
	
	///printf("n = %d, N = %d, linesPerThread = %d, overhead = %d\n",n, N, linesPerThread, overhead);
	for(i = 0, nextline = 1; i < n; i++) /*nextline beginnt bei 1, da rand nicht berechnet wird + increment von nextline*/
	{
		ta[i].arguments = arguments;
		ta[i].options = options;
		
		ta[i].from = nextline;
		nextline += linesPerThread + (i < overhead ? 1 : 0); 	/* nextline += range*/
		ta[i].to = nextline; 					/* -1 da der letzte nicht mitgezählt wird */		
	///	printf(" i = %i: from = %d, to = %d\n", i, ta[i].from, ta[i].to);		
		
		ta[i].maxresiduum = 0.0;
		
	}
	
	*thread_args = ta;
}
static void freePthreadVariables(struct thread_args* thread_args, pthread_t* threads)
{
	/* nur die struct-array selber free, da member schon free sind */
	if(thread_args != NULL)
	{
		free(thread_args);
	}
	
	if(threads != NULL)
	{
		free(threads);
	}
}
/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
	int i;

	for (i = 0; i < arguments->num_matrices; i++)
	{
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);

}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	int i, j;

	int const N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

		for (j = 0; j <= N; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
	int g, i, j;                                /*  local variables for loops   */

	int const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{
			for (i = 0; i <= N; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * i);
				Matrix[g][i][N] = h * i;
				Matrix[g][0][i] = 1.0 - (h * i);
				Matrix[g][N][i] = h * i;
			}
		}

		for (g = 0; g < arguments->num_matrices; g++)
		{
			Matrix[g][N][0] = 0.0;
			Matrix[g][0][N] = 0.0;
		}
	}
}

/* Führt die Iteration von calculate auf den Zeilen von from bis to aus 
 * und gibt das maxresiduum dieser Zeilen zurück */
static double calcLoop(double**  Matrix_In, double** Matrix_Out, int const from, int const to, int const N, double const h, struct options const* options)
{
	double residuum = 0, maxresiduum = 0, star = 0;
	int i, j;
	
	/* over all rows */
		for (i = from; i < to; i++)
		{
			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);
	
				if (options->inf_func == FUNC_FPISIN)
				{
					star += (0.25 * TWO_PI_SQUARE * h * h) * sin((PI * h) * (double)i) * sin((PI * h) * (double)j);
				}

				if (options->termination == TERM_PREC)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}
		return maxresiduum;
}

/*
 * calculateThreadMain: solves the equation für den Hauptthread
 * @thread_args:
 * 		array der Argument der anderen Threads
 * */
static
void
calculateThreadMain (struct calculation_results *results, struct thread_args* thread_args)
{
	struct calculation_arguments const* arguments = thread_args[0].arguments;
	struct options const* options = thread_args[0].options;
	
	int const from = thread_args[0].from;		/* Grenzen (von row i) für einen Thread	( in [1,N] )*/
	int const to = thread_args[0].to;
	
	//int i, j;                                   /* local variables for loops  */
	int m1, m2, i;                                 /* used as indices for old and new matrices       */
	//double star;                                /* four times center value minus 4 neigh.b values */
	//double residuum;                            /* residuum of current iteration                  */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	int const n = options->number;			/* n = Anzahl der Threads */
	int const N = arguments->N;
	double const h = arguments->h;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	while (run)
	{
		
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxresiduum = calcLoop(Matrix_In, Matrix_Out, from, to, N, h, options);


		/* darauf warten, dass alle Threads ihre Zeilen berechnet haben => gemeinsames Enden der Berechnung*/
		pthread_barrier_wait(&barrier); 
		

		/* maxresiduum aus den anderen Threads beachten */
		for(i = 1; i < n; i++)
		{
			if(maxresiduum < thread_args[i].maxresiduum)
			{
				maxresiduum = thread_args[i].maxresiduum;
			}
		}

		results->stat_iteration++;
		results->stat_precision = maxresiduum;
		
		
		/* check for stopping calculation, depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxresiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
		
		
		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		
		/* Abbruchbedingungen des loops */
		if(term_iteration <= 0) /* Abbruch */
		{
			run = false; /* den threads sagen, dass sie fertig sind */
		}
		
		/* darauf warten, dass der Leiterthread (dieser hier) die synchro beendet hat => gemeinsames Ende des Loops */
		pthread_barrier_wait(&barrier); 
		
	}

	
	results->m = m2;
}


/* ************************************************************************ */
/* calculateThread: solves the equation für einen Thread                    */
/*    @args:   Typ = thread_args; pointer auf arg für diesen Thread         */
/* ************************************************************************ */
static
void *
calculateThread (void* args)
{
	struct thread_args* thread_args = (struct thread_args*) args;
	struct calculation_arguments const* arguments = thread_args->arguments;
	struct options const* options = thread_args->options;
	
	int const from = thread_args->from;	     /* Grenzen (von row i) für einen Thread	( in [1,N] )*/
	int const to = thread_args->to;
	
	//int i, j;                                   /* local variables for loops  */
	int m1, m2, i;                                 /* used as indices for old and new matrices       */
	//double star;                                /* four times center value minus 4 neigh.b values */
	//double residuum;                            /* residuum of current iteration                  */
	//double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;


	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	while (run)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		thread_args->maxresiduum = calcLoop(Matrix_In, Matrix_Out, from, to, N, h, options);


		pthread_barrier_wait(&barrier); /* darauf warten, dass alle Threads ihre Zeilen berechnet haben */

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;
		
		pthread_barrier_wait(&barrier); /* darauf warten, dass der Leiterthread die synchro beendet hat => dann weitermachen*/
	}
	
	return NULL;
}




/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;

	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	//Calculate Flops
	// star op = 5 ASM ops (+1 XOR) with -O3, matrix korrektur = 1
	double q = 6;
	double mflops;
	long dataPoints = (long) (N - 1) * (N - 1) * results->stat_iteration;
	printf("Berechnungszeit:    %f s \n", time);

	if (options->inf_func == FUNC_F0)
	{
		// residuum: checked 1 flop in ASM, verified on Nehalem architecture.
		q += 1.0;
	}
	else
	{
		// residuum: 11 with O0, but 10 with "gcc -O3", without counting sin & cos
		q += 10.0;
	}

	/* calculate flops  */
	mflops = (q * dataPoints) * 1e-6;
	printf("Executed float ops: %f MFlop\n", mflops);

	printf("Speed:              %f MFlop/s\n", mflops / time);

	printf("Memory footprint:   %f MiB\n",   (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Mem thp (read):     %f MiB/s\n", dataPoints * 4.0 / 1024 / 1024  / time);
	printf("Mem thp (write):    %f MiB/s\n", dataPoints * 1.0 / 1024 / 1024 / time);
	printf("Mem thp (total):    %f MiB/s\n", dataPoints * 5.0 / 1024 / 1024  / time);

	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauss-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %d\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0)
	{
		printf("f(x,y)=0");
	}
	else if (options->inf_func == FUNC_FPISIN)
	{
		printf("f(x,y)=2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:       ");

	if (options->termination == TERM_PREC)
	{
		printf("Hinreichende Genaugkeit");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %d\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;
	
	struct thread_args* thread_args; 	/* array mit den argumenten für die Threads */
	pthread_t* threads;			/* array mit den threads */
	int i; 					/* Index für loop */
	int rc;					/* returncode */
	void *status;
	pthread_attr_t attr_arg;
	pthread_attr_init(&attr_arg);
	pthread_attr_setdetachstate(&attr_arg, PTHREAD_CREATE_JOINABLE);
	/* get parameters */
	AskParams(&options, argc, argv);              /* ************************* */

	initVariables(&arguments, &results, &options);           /* ******************************************* */
	initPthreadVariables(&thread_args, &threads, &arguments, &options); /* Vars fürs Threadding init */

	allocateMatrices(&arguments);        /*  get and initialize variables and matrices  */
	initMatrices(&arguments, &options);            /* ******************************************* */

	gettimeofday(&start_time, NULL);                   /*  start timer         */
	
	for(i = 1; i < options.number; i++)
	{	
		rc = pthread_create(&threads[i], &attr_arg, &calculateThread, (void*) &(thread_args[i]));
		if (rc){
         		printf("ERROR; return code from pthread_create() is %d\n", rc);
         		exit(-1);
      		}
	}
	pthread_attr_destroy(&attr_arg);
	calculateThreadMain(&results, thread_args);
	  
	gettimeofday(&comp_time, NULL);                   /*  stop timer          */
	
	for(i=1; i<options.number; i++) {
      		rc = pthread_join(threads[i], &status);
      		if (rc) {
         		printf("ERROR; return code from pthread_join() is %d\n", rc);
         		exit(-1);
         	}
        }
	
	displayStatistics(&arguments, &results, &options);                                  /* **************** */
	DisplayMatrix("Matrix:",                              /*  display some    */
			arguments.Matrix[results.m][0], options.interlines);            /*  statistics and  */

	freeMatrices(&arguments);                                       /*  free memory     */
	freePthreadVariables(thread_args, threads);
	pthread_exit(NULL);	
}
