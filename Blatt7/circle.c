#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define MASTER_PROCESS 0
#define SEND_DATA_TAG 123

int*
init (int part, int bufLength, int rest)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  int* buf = malloc(sizeof(int) * bufLength);

  srand(time(NULL) * rank);

  for (int i = 0; i < part; i++)
  {
    buf[i] = rand() % 25; //do not modify % 25
  }
  
  if (rank >= rest ) {
	  buf[part] = -1;
  }

  return buf;
}

int*
circle (int* buf)
{
  int rank, processCount;
  MPI_Status status;
	
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &processCount);
	
  int sendID = (rank + 1) % processCount;
  int recvID = rank == 0 ? processCount - 1 : rank - 1;
	
  MPI_Send( &buf, 1, MPI_INT, sendID, SEND_DATA_TAG, MPI_COMM_WORLD);
  
  MPI_Recv( &buf, 1, MPI_INT, recvID, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  
  return buf;
}

int
main (int argc, char** argv)
{
  char arg[256];
  int N;
  int rank;
  int processCount;
  int* buf;
  
  //Initialisierung von mpi
  MPI_Init(&argc, &argv);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &processCount);

  if (argc < 2)
  {
	  if (rank == MASTER_PROCESS)
    {
      printf("Arguments error!\n");
    }
	MPI_Finalize();
    return EXIT_FAILURE;
  }

  sscanf(argv[1], "%s", arg);

  // Array length
  N = atoi(arg);
  
  if (processCount > N)
  {
	  if (rank == MASTER_PROCESS)
    {
      printf("Array size must be higher than Process count!\n");
    }
	MPI_Finalize();
    return EXIT_FAILURE;
  }
  
  int part = N / processCount;
  int rest = N % processCount;
  int bufLength = (rest != 0) ? part + 1 : part;
  part = rank < rest ? part + 1 : part;
  
  //printf ("\nrank %d: %d\n", rank, part);
  buf = init(part, bufLength, rest);

  printf("\nBEFORE\n");

  for (int i = 0; i < part; i++)
  {
    printf ("rank %d: %d\n", rank, buf[i]);
  }

  circle(buf);

  printf("\nAFTER\n");

  for (int j = 0; j < part; j++)
  {
    printf ("rank %d: %d\n", rank, buf[j]);
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}
