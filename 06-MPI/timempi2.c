#define _BSD_SOURCE

#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>

#define send_data_tag 1001
#define return_data_tag 1002

int main(int argc, char **argv)
{
  MPI_Status status;

  int root_process; /* id of root process */
  int my_id; /* id of the current process */
  int i; /* loop iterator */
  int amount_procs; /* total amount of processes */

  long local_milliseconds = 1000000; /* set local_milliseconds so that any valid value is lower for the minimum calculation. */
  long min_milliseconds;

  const int timestamplength = 64;
  char timestamp[64]; /* process individual timestamp */

  /* start processes*/
  MPI_Init(&argc, &argv);

  root_process = 0;

  /* get process id and amount of total processes */
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &amount_procs);

  /* root process: monitoring and printing */
  if (my_id == root_process) {
    for(i = 1; i < amount_procs; i++)
    {
      MPI_Recv(timestamp, timestamplength, MPI_CHAR, i, return_data_tag, MPI_COMM_WORLD, &status);
      printf("%s\n", timestamp);
    }
  }
  else
  {
    char hostname[10]; /* hostname of current process */

    struct timeval tv; /* time struct */
    time_t nowtime;
    struct tm *nowtm;
    char tmbuf[64];

    /* gets hostname and time */
    gethostname(hostname, 40);
    gettimeofday(&tv,NULL);
    nowtime = tv.tv_sec;
    nowtm = localtime(&nowtime);
    strftime(tmbuf, sizeof tmbuf, "%Y-%m-%d %H:%M:%S", nowtm);
    snprintf(timestamp, sizeof timestamp, "%s: %s.%06ld", hostname, tmbuf, tv.tv_usec);
    local_milliseconds = tv.tv_usec;

    MPI_Send(timestamp, timestamplength, MPI_CHAR, root_process, return_data_tag, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Reduce(&local_milliseconds, &min_milliseconds, 1, MPI_LONG, MPI_MIN, root_process, MPI_COMM_WORLD);
  if(my_id == root_process && amount_procs > 1)
  {
    printf("%ld\n", min_milliseconds);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  printf("Rang %d beendet jetzt!\n", my_id);
  MPI_Finalize();
}
