
#include <mpi.h>
#include <sys/time.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#define MASTER 0



int
main (int argc, char** argv)
{
  int rank, size;
	int i;										// Laufindex für Schleifen
		
	MPI_Datatype mpi_timeval_t;
  	
	const int hostlen = HOST_NAME_MAX + 1;  // +1 wegen \0-Byte
	char host[hostlen];
  struct timeval time;
  char** rchost;                // Buffer fuer MASTER (receive)
  struct timeval *rctime;	      // Buffer fuer MASTER (receive)
	long min;                     // minimale musec

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);   // Nummer des Prozesses holen, auf dem diese Instanz des Programms läuft
	MPI_Comm_size(MPI_COMM_WORLD, &size);     // Anzahl der Prozesse merken

 
  // init timeval type
  MPI_Type_contiguous(2,MPI_LONG,&mpi_timeval_t);
  MPI_Type_commit(&mpi_timeval_t);
  
  
  
	if(rank == MASTER) // Der Mainproc ließt daten ein 
	{
     rchost = malloc(sizeof(*rchost) * size);
     rchost[0] = malloc(sizeof(**rchost) * hostlen * size);
          
     for(i = 1; i < size; i++)
     {
        rchost[i] = rchost[i-1] + hostlen;
     }
     char *asd = rchost[0];
     
     rctime = malloc(sizeof(*rctime) * size);
    
  }
	// host und zeit auslesen
  gethostname(host, sizeof(host));
	gettimeofday(&time, NULL);

  
  // host und zeit senden
	MPI_Gather( &host, hostlen, MPI_CHAR, (rchost[0]), hostlen, MPI_CHAR, MASTER, MPI_COMM_WORLD);
  MPI_Gather( &time, 1, mpi_timeval_t, rctime, 1, mpi_timeval_t, MASTER, MPI_COMM_WORLD);

  // Minimum finden
  MPI_Reduce(&(time.tv_usec), &min, 1, MPI_LONG, MPI_MIN, MASTER, MPI_COMM_WORLD);
  
  if(rank == MASTER)
  {
    // Minimum ausgeben
    printf("min: %ld\n", min);
    for( i = 1; i < size; i++)
    { 
        // Hoste und Timestamp ausgeben
        printf("%s: %d\n", rchost[i], rctime[i].tv_usec);
    }	
  }
  
  // sync vorm Beenden
  MPI_Barrier( MPI_COMM_WORLD );


  // free type
	MPI_Type_free(&mpi_timeval_t);
  if(rank == MASTER)
  {
     free(rchost[0]);
     free(rchost);
     
     free(rctime);
  }  

  // Beenden bekannt geben.
  printf("Rang %d beendet jetzt!\n", rank);
	
	MPI_Finalize();
	return 0;
}
