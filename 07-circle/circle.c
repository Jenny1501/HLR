#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define printtag 2000
#define circletag 2001
#define defaulttag 2002
#define finishedtag 2003

int*
init (int max_N, int my_N)
{
  int* buf = malloc(sizeof(int) * max_N);

  srand(time(NULL));

  for (int i = 0; i < my_N; i++)
  {
    /* use rand_r for thread safety */
    buf[i] = rand_r() % 25; //do not modify % 25
  }
  return buf;
}

/*int*
circle (int* buf, int* temp_buf)
{
  return buf;
}*/

int
main (int argc, char** argv)
{
  char arg[256];
  int N, my_N, max_N;
  int rank, amount_procs, root_process;
  int endcondition;
  int* buf, *temp_buf;
  int i;
  int length = 0;

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
  max_N = (rest > 0) ? my_N : my_N-1;

  buf = init(max_N, my_N);
  temp_buf = malloc(sizeof(int) * max_N);
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

    for(i = 1; i < amount_procs; i++)
    {
      MPI_Recv(&length, 1, MPI_INT, i, printtag, MPI_COMM_WORLD,&status);
      MPI_Recv(temp_buf, length, MPI_INT, i, printtag, MPI_COMM_WORLD, &status);
      printf("rank %d:", i);
      for (int j = 0; j < length; j++)
      {
        printf (" %d", temp_buf[j]);
      }
      printf("\n");
    }
    /* send endcondition to last process */
    MPI_Send(&endcondition, 1, MPI_INT, amount_procs - 1, defaulttag, MPI_COMM_WORLD);
  }
  /* non root process: sending data to root for print */
  else
  {
    MPI_Send(&my_N, 1, MPI_INT, root_process, printtag, MPI_COMM_WORLD);
    MPI_Send(buf, my_N, MPI_INT, root_process, printtag, MPI_COMM_WORLD);
  }
  /* last process: get endcondition int */
  if(rank == amount_procs - 1)
  {
    MPI_Recv(&endcondition, 1, MPI_INT, root_process, defaulttag, MPI_COMM_WORLD, &status);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /* circle(buf) is the following: */
  int finished = 0;
  while(!finished)
  {
    if(rank == root_process)
    {
      MPI_Send(&my_N, 1, MPI_INT, rank+1, circletag, MPI_COMM_WORLD);
      MPI_Send(buf, my_N, MPI_INT, rank+1, circletag, MPI_COMM_WORLD);

      MPI_Recv(&length, 1, MPI_INT, amount_procs-1, circletag, MPI_COMM_WORLD, &status);
      MPI_Recv(&temp_buf, length, MPI_INT, amount_procs-1, circletag, MPI_COMM_WORLD, &status);
    }
    else if(rank < amount_procs -1)
    {
      MPI_Recv(&length, 1, MPI_INT, rank-1, circletag, MPI_COMM_WORLD, &status);
      MPI_Recv(&temp_buf, length, MPI_INT, rank-1, circletag, MPI_COMM_WORLD, &status);

      MPI_Send(&my_N, 1, MPI_INT, rank+1, circletag, MPI_COMM_WORLD);
      MPI_Send(buf, my_N, MPI_INT, rank+1, circletag, MPI_COMM_WORLD);
    }
    else
    {
      MPI_Recv(&length, 1, MPI_INT, rank-1, circletag, MPI_COMM_WORLD, &status);
      MPI_Recv(&temp_buf, length, MPI_INT, rank-1, circletag, MPI_COMM_WORLD, &status);

      MPI_Send(&my_N, 1, MPI_INT, root_process, circletag, MPI_COMM_WORLD);
      MPI_Send(buf, my_N, MPI_INT, root_process, circletag, MPI_COMM_WORLD);

      finished = temp_buf[0] == endcondition;
    }

    /* swap buffers and lengths to take the received values into buf */
    int *x = buf;
    buf = temp_buf;
    temp_buf = x;
    int y = length;
    length = my_N;
    my_N = y;

    MPI_Bcast(&finished, 1, MPI_INT, amount_procs-1, MPI_COMM_WORLD);
  }

  if (rank == root_process) {
    printf("\nAFTER\n");
    printf("rank %d:", rank);
    for (int j = 0; j < my_N; j++)
    {
      printf (" %d",buf[j]);
    }
    printf("\n");

    for(i = 1; i < amount_procs; i++)
    {
      MPI_Recv(&length, 1, MPI_INT, i, printtag, MPI_COMM_WORLD,&status);
      MPI_Recv(temp_buf, length, MPI_INT, i, printtag, MPI_COMM_WORLD, &status);
      printf("rank %d:", i);
      for (int j = 0; j < length; j++)
      {
        printf (" %d", temp_buf[j]);
      }
      printf("\n");
    }
  }
  /* non root process: sending data to root for print */
  else
  {
    MPI_Send(&my_N, 1, MPI_INT, root_process, printtag, MPI_COMM_WORLD);
    MPI_Send(buf, my_N, MPI_INT, root_process, printtag, MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}
