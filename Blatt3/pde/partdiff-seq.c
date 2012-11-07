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
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include "partdiff-seq.h"


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

struct calculation_arguments
{
	int     N;              /* number of spaces between lines (lines=N+1)     */
	int     num_matrices;   /* number of matrices                             */
	//double  ***Matrix;      /* index matrix used for addressing M           */
	struct matrix_small** Mat;
	//double  *M;             /* two matrices with real values                  */
	double  h;              /* length of a space between two lines            */
	double  h_2;		/* neu: length of a space between two lines squared    */
	double  pi_h;		/* neu: value of PI * arguments->h precomputed         */
	double *mysin;		/* neu: array to precalculate sin		  */
};

struct calculation_results
{
	int     m;
	int     stat_iteration; /* number of current iteration                    */
	double  stat_precision; /* actual precision of all slaves in iteration    */
};

struct matrix_small
{
  double* M;
  int N;
  int method; /* Gauss Seidel or Jacobi method of iteration     */
  int inf_func;
};

static struct matrix_small* Matrix_init(int N, int method, int inf_func)
{
  struct matrix_small* mat = allocateMemory(sizeof(*mat));
  mat->N = N;
  mat->method = method;
  mat->inf_func = inf_func;
  if(mat->method == METH_GAUSS_SEIDEL)
  {
    int N_small = N+1; 
    mat->M = allocateMemory( sizeof(*(mat->M)) * (((N_small*(N_small+1))/ 2)));
  }
  else // Jacobi
  {
    int N_small = mat->inf_func == FUNC_F0 ? (N+1) :(N + 2) / 2; // N + 1 = lines, und + 1 mit / 2 ist round
    mat->M = allocateMemory( sizeof(*(mat->M)) * (((N_small*(N_small+1))/ 2))); 
  }
  return mat;
}
void Matrix_free(struct matrix_small* mat)
{
  if(mat != NULL)
  {
    free(mat->M);
  }
  free(mat);
}

static double Matrix_getValue(struct matrix_small* mat, int i, int j)
{
  int N = mat->N; // letztmoegliche Index
  int tmp; // index zum Tauschen
  
  // Man kann Gauss Seidel nur 1/2 ,  Jacobi kann man 1/4 - 1/8
  if(mat->method == METH_GAUSS_SEIDEL)
  {
    if(i < j) // Hauptdiag: i und j vertauschen
    {
      tmp = i;
      i = j;
      j = tmp;
      
    }
    return mat->M[((i*(i+1))/ 2)  + j];
  }
  else if(mat->inf_func == FUNC_F0) // Jakobi (1/4)
  {
    if(i < j) // Hauptdiag: i und j vertauschen
    {
      tmp = i;
      i = j;
      j = tmp;
    }
    if(i + j > N) // Nebendiag: i und j vertauschen
    {
      tmp = i;
      i = N - j;
      j = N - tmp;
    }
    return mat->M[((i*(i+1))/ 2)  + j];   
  }                                     
  else // Jacobi (1/8)
  {
    i = i >= N/2 ? N - i : i;
    j = j >= N/2 ? N - j : j;
     
    if(i < j) // Hauptdiag: i und j vertauschen
    {
      tmp = i;
      i = j;
      j = tmp;
    }
    return mat->M[((i*(i+1))/ 2)  + j];      
  }
  
      
}
static void Matrix_setValue(struct matrix_small* mat, int i, int j, double v)
{ 
  int N = mat->N; // letztmoegliche Index
  int tmp; // index zum Tauschen
  // Man kann Gauss Seidel nur 1/2 ,  Jacobi kann man 1/4 - 1/8
  if(mat->method == METH_GAUSS_SEIDEL)
  {
    if(i < j) // Hauptdiag: i und j vertauschen
    {
      tmp = i;
      i = j;
      j = tmp;      
    }
    mat->M[((i*(i+1))/ 2)  + j] = v;
  }
  else if(mat->inf_func == FUNC_F0) // Jakobi (1/4)
  {
    if(i < j) // Hauptdiag: i und j vertauschen
    {
      tmp = i;
      i = j;
      j = tmp;
    }
    if(i + j > N) // Nebendiag: i und j vertauschen
    {
      tmp = i;
      i = N - j;
      j = N - tmp;
    }
    mat->M[((i*(i+1))/ 2)  + j] = v;   
  }
  else // Jacobi (1/8)
  {
    i = i >= N/2 ? N - i : i;
    j = j >= N/2 ? N - j : j;
     
    if(i < j) // Hauptdiag: i und j vertauschen
    {
      tmp = i;
      i = j;
      j = tmp;
    }
    mat->M[((i*(i+1))/ 2)  + j] = v;      
  }        
      
}

double Matrix_interface_getValue(Matrix_interface* mat,int i, int j)
{
  return Matrix_getValue(mat, i,j);
}

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */

<<<<<<< HEAD
=======
/* free array mysin */
static void freeMysin(struct calculation_arguments* arguments)
{
  if(arguments->mysin != NULL)
    free(arguments->mysin);
}

/* precalculate the sinus values for the residuum calculation */
static void initMysin(struct calculation_arguments* arguments)
{
  int i;
  int len = (arguments->N + 2) / 2; 
  arguments->mysin = malloc(sizeof(*(arguments->mysin)) *  len);
  for(i = 0;  i < len; i++)
  {
    arguments->mysin[i] = sin((double)(i) * arguments->pi_h);
  }
}


>>>>>>> upstream/master

/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	arguments->N = options->interlines * 8 + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = (float)( ( (float)(1) ) / (arguments->N));
	arguments->h_2 = arguments->h * arguments->h; /* neu */
	arguments->pi_h = PI * arguments->h;	/* neu */


	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
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
		//free(arguments->Matrix[i]);
		Matrix_free(arguments->Mat[i]);
	}
  free(arguments->Mat);
	//free(arguments->Matrix);
	//free(arguments->M);
}


/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments, struct options* options)
{
	int i;//, j;

	//int N = arguments->N;

	//arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	//arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));
  // neu
  arguments->Mat = allocateMemory(arguments->num_matrices * sizeof(**(arguments->Mat))); 
  
	for (i = 0; i < arguments->num_matrices; i++)
	{
		//arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));
    // neu
    arguments->Mat[i] = allocateMemory(arguments->num_matrices * sizeof(*(arguments->Mat)));
    arguments->Mat[i] = Matrix_init(arguments->N, options->method, options->inf_func);
    
		//for (j = 0; j <= N; j++)
		//{
		//	arguments->Matrix[i][j] = (double*)(arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1)));
		//}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options* options)
{
	int g, i, j;                                /*  local variables for loops   */

	int N = arguments->N;
	double h = arguments->h;
	//double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
		{
			for (j = 0; j <= N; j++)
			{
				// Matrix[g][i][j] = 0;
				Matrix_setValue(arguments->Mat[g], i, j, 0);
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for(i = 0; i <= N; i++)
		{
			for (j = 0; j < arguments->num_matrices; j++)
			{
				//Matrix[j][i][0] = 1 - (h * i);
				//Matrix[j][i][N] = h * i;
				//Matrix[j][0][i] = 1 - (h * i);
				//Matrix[j][N][i] = h * i;
				Matrix_setValue(arguments->Mat[j], i, 0,  1 - (h * i));
				Matrix_setValue(arguments->Mat[j], i, N,  h * i      );
				Matrix_setValue(arguments->Mat[j], 0, i,  1 - (h * i));
				Matrix_setValue(arguments->Mat[j], N, i,  h * i      );
			}
		}

		for (j = 0; j < arguments->num_matrices; j++)
		{
			//Matrix[j][N][0] = 0;
			//Matrix[j][0][N] = 0;
			Matrix_setValue(arguments->Mat[j], N, 0,  0);
			Matrix_setValue(arguments->Mat[j], 0, N,  0);
		}
	}
	
	//for(i = 0; i <=N ; i++)
	// for(j = 0; j <=N; j++)
	//   for (g = 0; g < arguments->num_matrices; g++)
	//     Matrix_setValue(arguments->Mat[g], i, j, Matrix[g][i][j]); 
}

/* ************************************************************************ */
/* getResiduum: calculates residuum                                         */
/* Input: x,y - actual column and row                                       */
/* ************************************************************************ */
/*
double
getResiduum (struct calculation_arguments* arguments, struct options* options, int x, int y, double star)
{
	if (options->inf_func == FUNC_F0)
	{
		return ((-star) / 4.0);
	}
	else
	{
		return ((TWO_PI_SQUARE * sin((double)(y) * arguments->pi_h) * sin((double)(x) * arguments->pi_h) * arguments->h_2 - star) / 4.0); // neu
	}
}
*/

/* free array mysin */
static void freeMysin(struct calculation_arguments* arguments)
{
  if(arguments->mysin != NULL)
    free(arguments->mysin);
}

/* precalculate the sinus values for the residuum calculation */
static void initMysin(struct calculation_arguments* arguments)
{
  int i;
  int len = arguments->N; 
  arguments->mysin = malloc(sizeof(*(arguments->mysin)) *  len);
  for(i = 0;  i < len; i++)
  {
    arguments->mysin[i] = sin((double)(i) * arguments->pi_h);
  }
}



/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments* arguments, struct calculation_results *results, struct options* options)
{
	int i, j;                                   /* local variables for loops  */
	int m1, m2;                                 /* used as indices for old and new matrices       */
	double star;                                /* four times center value minus 4 neigh.b values */
	
	double residuum;                            /* residuum of current iteration                  */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */
  double *mysin = arguments->mysin;	/* neu */
  
	int N = arguments->N;
<<<<<<< HEAD
	double*** Matrix = arguments->Matrix;
	double **Mat2;  /* neu */	
	double *mysin = arguments->mysin;	/* neu */

=======
	int Nbound = (options->method == METH_GAUSS_SEIDEL || options->inf_func == FUNC_F0) ? N : (N / 2 + 1); // Man kann Gauss Seidel nur 1/2 ,  Jacobi kann man 1/4 - 1/8
	///double*** Matrix = arguments->Matrix;
	///double **Mat2;  /* neu */	
  
  struct matrix_small** Matrix = arguments->Mat;
  struct matrix_small* Mat2;
  
>>>>>>> upstream/master
	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_GAUSS_SEIDEL)
	{
		m1=0; m2=0;
	}
	else
	{
		m1=0; m2=1;
	}

	
	while (options->term_iteration > 0)
	{
	  Mat2 = Matrix[m2]; /* neu */
		maxresiduum = 0;
    
		/* over all rows */
		/* alt: for (j = 1; j < N; j++) */
		// for (i = 1; i < N; i++)
		for (i = 1; i < Nbound; i++)
		{
		  int jbound = ((options->method == METH_GAUSS_SEIDEL || options->inf_func == FUNC_FPISIN) || i < N/2 +1)  
         ? i : N - i;
         
			/* over all columns */
			/* for (i = 1; i < N; i++) */
			//for (j = 1; j < N; j++)
			for (j = 1; j <= jbound; j++)
			{
<<<<<<< HEAD
				double mij = Mat2[i][j]; /* neu */
				/* alt: star = -Matrix[m2][i-1][j] - Matrix[m2][i][j-1] - Matrix[m2][i][j+1] - Matrix[m2][i+1][j] + 4.0 * Matrix[m2][i][j];
	
				alt: star = -Mat2[i-1][j] - Mat2[i][j-1] - Mat2[i][j+1] - Mat2[i+1][j] + 4.0 * mij; */
				star = (- Mat2[i][j-1] - Mat2[i][j+1] - Mat2[i-1][j] - Mat2[i+1][j]) + 4.0 * mij ;

				//residuum = getResiduum(arguments, options, i, j, star);
=======
				///double mij = Mat2[i][j]; /* neu */
				double mij = Matrix_getValue(Mat2, i, j); /* neu */
				
        /* alt: star = -Matrix[m2][i-1][j] - Matrix[m2][i][j-1] - Matrix[m2][i][j+1] - Matrix[m2][i+1][j] + 4.0 * Matrix[m2][i][j];
				alt: star = -Mat2[i-1][j] - Mat2[i][j-1] - Mat2[i][j+1] - Mat2[i+1][j] + 4.0 * mij; 
				star = (- Mat2[i][j-1] - Mat2[i][j+1] - Mat2[i-1][j] - Mat2[i+1][j]) + 4.0 * mij ;  */
        /*  alt: ...
        alt: residuum = getResiduum(arguments, options, i, j, star); 
>>>>>>> upstream/master
				if (options->inf_func == FUNC_F0)
				{
					residuum = ((-star)/4);
				}
				else
				{
<<<<<<< HEAD
					/*residuum = ((TWO_PI_SQUARE * sin((double)(j) * arguments->pi_h) * sin((double)(i) * arguments->pi_h) * arguments->h_2 - star)/4);  neu */
					residuum = ((TWO_PI_SQUARE * mysin[j] * mysin[i] * arguments->h_2 - star)/4); /* neu */
=======
				  /* alt: residuum = ((TWO_PI_SQUARE * sin((double)(j) * arguments->pi_h) * sin((double)(i) * arguments->pi_h) * arguments->h_2 - star)/4); 
					residuum = ((TWO_PI_SQUARE * mysin[i]) * mysin[j] * arguments->h_2 - star)/4; /* neu 
>>>>>>> upstream/master
				}
				korrektur = residuum;
				residuum = (residuum < 0) ? -residuum : residuum;
				maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;

<<<<<<< HEAD
				/* alt: Matrix[m1][i][j] = Mat2[i][j] + korrektur; */
				Matrix[m1][i][j] = mij + korrektur;
				
			}
=======
				/* alt: Matrix[m1][i][j] = Mat2[i][j] + korrektur; 
				Matrix[m1][i][j] = mij + korrektur; 
        
        */
				
				///star = 0.25 * (+Mat2[i-1][j] + Mat2[i][j-1] + Mat2[i][j+1] + Mat2[i+1][j] );
        star = 0.25 * (+Matrix_getValue(Mat2, i-1, j) +Matrix_getValue(Mat2, i, j-1) +Matrix_getValue(Mat2, i, j+1) +Matrix_getValue(Mat2, i+1, j));

				if (options->inf_func == FUNC_FPISIN)
				{
					star += 0.25 * TWO_PI_SQUARE * mysin[i] * mysin[j] * arguments->h_2;
					// ohne mysin: star += (0.25 * TWO_PI_SQUARE * arguments->h_2) * sin((PI * arguments->h) * (double)i) * sin((PI * arguments->h) * (double)j);
				}
        
        if (options->termination == TERM_PREC || options->term_iteration == 1)
				{
					//residuum = mij - star;
					residuum = mij - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}


				
				///Matrix[m1][i][j] = star;
				Matrix_setValue(Matrix[m1], i, j, star);

        }   
>>>>>>> upstream/master
		}

		results->stat_iteration++;
		results->stat_precision = maxresiduum;
    
		/* exchange m1 and m2 */
		i=m1; m1=m2; m2=i;
    
		/* check for stopping calculation, depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxresiduum < options->term_precision)
			{
				options->term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			options->term_iteration--;
		}
	}

	results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments* arguments, struct calculation_results *results, struct options* options)
{
	int N = arguments->N;

	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;
	printf("Berechnungszeit:    %f s \n", time);

	//Calculate Flops
	// star op = 5 ASM ops (+1 XOR) with -O3, matrix korrektur = 1
	double q = 6;
	double mflops;

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
	mflops = (q * (N - 1) * (N - 1) * results->stat_iteration) * 1e-6;
	printf("Executed float ops: %f MFlop\n", mflops);
	printf("Speed:              %f MFlop/s\n", mflops / time);

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

	/* get parameters */
	AskParams(&options, argc, argv);              /* ************************* */

	initVariables(&arguments, &results, &options);           /* ******************************************* */
	initMysin(&arguments);	

	allocateMatrices(&arguments, &options);        /*  get and initialize variables and matrices  */
	initMatrices(&arguments, &options);            /* ******************************************* */

<<<<<<< HEAD
=======
  if (options.inf_func == FUNC_FPISIN) /* Init Cache für Sinus berechnung */
    initMysin(&arguments);
  

>>>>>>> upstream/master
	gettimeofday(&start_time, NULL);                   /*  start timer         */
	calculate(&arguments, &results, &options);                                      /*  solve the equation  */
	gettimeofday(&comp_time, NULL);                   /*  stop timer          */

//  int i, j;
//	for(i = 0; i <=arguments.N ; i++)
//	 for(j = 0; j <=arguments.N; j++)
//	     arguments.Matrix[results.m][i][j] = Matrix_getValue(arguments.Mat[results.m], i, j); 

	displayStatistics(&arguments, &results, &options);                                  /* **************** */
	DisplayMatrix("Matrix:",                              /*  display some    */
			arguments.Mat[results.m], options.interlines);            /*  statistics and  */

<<<<<<< HEAD
	freeMysin(&arguments);
=======
  freeMysin(&arguments);                     /* Cache für Sinus freigeben */
>>>>>>> upstream/master
	freeMatrices(&arguments);                                       /*  free memory     */
	

	return 0;
}
