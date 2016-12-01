#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define printtag 2000
#define circletag 2001

int*
init (int my_N)
{
  int* buf = malloc(sizeof(int) * my_N);

  srand(time(NULL));

  for (int i = 0; i < my_N; i++)
  {
    buf[i] = rand() % 25; //do not modify % 25
  }
  return buf;
}

int*
circle (int* buf)
{
  // Todo
  return buf;
}

int
main (int argc, char** argv)
{
  char arg[256];
  int N, my_N;
  int rank, amount_procs, root_process;
  int endcondition;
  int* buf, *new_buf;
  int i;

  MPI_Status status;

  if (argc < 2)
  {
    printf("Arguments error!\n");
    return EXIT_FAILURE;
  }

  sscanf(argv[1], "%s", arg);

  /* start processes*/
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &amount_procs);

  // Array length
  N = atoi(arg);
  my_N = N / amount_procs;
  int rest = N % amount_procs;
  if(rank < rest)
  {
    my_N++;
  }
  buf = init(my_N);

  /* get process rank and amount of total processes */
  root_process = 0;


  /* root process: fetching data and printing*/
  if (rank == root_process) {
    printf("\nBEFORE\n");
    printf("rank %d:", rank);
    for (int j = 0; j < my_N; j++)
    {
      printf (" %d",buf[j]);
    }
    printf("\n");
    /* After this printing, the root process will already be in the first iteration of circle
      It seems dirty because it is dirty.*/
    new_buf = buf;
    for(i = 1; i < amount_procs; i++)
    {
      MPI_Recv(buf, my_N, MPI_INT, i, printtag, MPI_COMM_WORLD, &status);
      printf("rank %d:", i);
      for (int j = 0; j < my_N; j++)
      {
        printf (" %d", buf[j]);
      }
      printf("\n");
    }
  }
  /* non root process: sending data to root for print */
  else
  {
    MPI_Send(buf, my_N, MPI_INT, root_process, printtag, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  circle(buf);

  //printf("\nAFTER\n");

  MPI_Finalize();
  return EXIT_SUCCESS;
}
