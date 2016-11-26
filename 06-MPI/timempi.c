




#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>

#define send_data_tag 1001
#define return_data_tag 1002

main(int argc, char **argv)
{
  MPI_Status status;

  int ierr; /* mpi error code */
  int root_process; /* id of root process */
  int my_id; /* id of the current process */
  int i; /* loop iterator */
  int amount_procs; /* total amount of processes */

  char timestamp[30]; /* process individual timestamp */

  /* start processes*/
  ierr = MPI_Init(&argc, &argv);

  root_process = 0;

  /* get process id and amount of total processes */
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &amount_procs);

  /* root process: monitoring and printing */
  if (my_id == root_process) {
    for(i = 1; i < amount_procs; i++)
    {
      ierr = MPI_Recv(timestamp, timestamp.length(), MPI_CHAR, MPI_ANY_SOURCE, return_data_tag, MPI_COMM_WORLD, &status);
      printf("%s\n", timestamp);
    }
  }
  else
  {
    char* hostname; /* hostname of current process */
    char* curtime; /* current time */

    struct timeval *tv; /* time */

    int err; /* debug error code */

    /* gets hostname and time */
    err = gethostname(hostname, 255);
    err = gettimeofday(tv,NULL);

    
  }

}
