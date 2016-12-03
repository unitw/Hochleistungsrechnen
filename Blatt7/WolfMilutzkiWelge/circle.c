#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
   
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#define ROOT 0
#define TAG_INIT 0
#define TAG_CIRCLE 1
#define TAG_PRINT 2
#define TAG_END_VALUE 3

typedef struct 
{
  int rank;
  int num_procs;
  int moveSize;
  int* bufsize;
  MPI_Status status;
} CommParams;

// In params müssen rank und num_procs gesetzt sein
int* init (int completeSize, CommParams* params)
{
  int chunkSize = completeSize / params->num_procs;
  int extra = completeSize % params->num_procs;
  params->moveSize = MAX(1, chunkSize);

  params->bufsize = malloc(sizeof(int) * params->num_procs);

  for (int i = 0; i < params->num_procs; ++i)
  {
    params->bufsize[i] = chunkSize;
    if ((i < extra) || (chunkSize == 0))
    {
      params->bufsize[i] += 1;
    }
  }

  int N = params->bufsize[params->rank];
  int* buf = malloc(sizeof(int) * N);

  // mit -1 initialisieren
  memset(buf, -1, sizeof(int) * N);

  if (params->rank == ROOT)
  {
    int allData[completeSize];
    srand(time(NULL));
    for (int i = 0; i < completeSize; i++)
    {
      allData[i] = rand() % 25; //do not modify %25
    }

    int start = params->bufsize[ROOT];
    for (int i = 1; i < params->num_procs; ++i)
    {
      int len = chunkSize;
      if (i < extra)
      {
        len += 1;
      }

      MPI_Send(&allData[start], len, MPI_INT, i, TAG_INIT, MPI_COMM_WORLD);   
      start += len;
    }
    memcpy(buf, allData, sizeof(int) * params->bufsize[ROOT]);
  }
  else
  {
    int len = chunkSize;
    if (params->rank < extra)
    {
      len += 1;
    }
    MPI_Recv(buf, len, MPI_INT, ROOT, TAG_INIT, MPI_COMM_WORLD, &params->status);
  }

  return buf;
}

int completeSize(int argc, char** argv)
{
  if (argc >= 2)
  {
    char arg[256];
    sscanf(argv[1], "%s", arg);
    return atoi(arg);
  } 

  return 0;
}

void printArray(int rank, int* buf, int len)
{
  for (int i = 0; i < len; ++i)
  {
    printf("rank %d:%d\n", rank, buf[i]);
  }
}

void printComplete(int* buf, CommParams* params)
{
  if (params->rank == ROOT)
  {
    printArray(ROOT, buf, params->bufsize[ROOT]);
    int procData[params->moveSize + 1];
    for (int i = 1; i < params->num_procs; ++i)
    {
      MPI_Recv(procData, params->bufsize[i], MPI_INT, i, TAG_PRINT, MPI_COMM_WORLD, &params->status);
      printArray(i, procData, params->bufsize[i]);
    }
  }
  else
  {
    MPI_Send(buf, params->bufsize[params->rank], MPI_INT, ROOT, TAG_PRINT, MPI_COMM_WORLD);
  }
}

void circle (int* buf, int endValue, CommParams* params)
{
  if (params->num_procs < 2)
  {
    return;
  }

  int myBufSize = params->bufsize[params->rank];
  int moveSize = params->moveSize;

  int sendTo = (params->rank + 1) % params->num_procs;
  int recFrom = (params->rank -1) % params->num_procs;
  int sendData[params->moveSize];

  int endReached = 0;

  while (1)
  {
    if (params->rank == (params->num_procs - 1) && buf[0] == endValue)
    {
      endReached = 1;
    }
    MPI_Bcast(&endReached, 1, MPI_INT, params->num_procs - 1, MPI_COMM_WORLD);

    if (endReached)
    {
      break;
    }

    // die zu sendenden Werte aus dem Array kopieren
    memcpy(sendData, &buf[myBufSize - moveSize], sizeof(int) * moveSize);
    buf[myBufSize - 1] = buf[0];  // den Wert den wir eventuell behalten wollen ans Ende kopieren

    if ((params->rank % 2) == 0)
    {
      MPI_Send(sendData, moveSize, MPI_INT, sendTo, TAG_CIRCLE, MPI_COMM_WORLD);
      MPI_Recv(buf, moveSize, MPI_INT, recFrom, TAG_CIRCLE, MPI_COMM_WORLD, &params->status);
    }
    else
    {
      MPI_Recv(buf, moveSize, MPI_INT, recFrom, TAG_CIRCLE, MPI_COMM_WORLD, &params->status);
      MPI_Send(sendData, moveSize, MPI_INT, sendTo, TAG_CIRCLE, MPI_COMM_WORLD);
    }
  }
}

int main (int argc, char** argv)
{
  int* buf;
  int endValue;
  CommParams params;

  int size = completeSize(argc, argv);

  if (size <= 0)
  {
    printf("Arguments error\n");
    return EXIT_FAILURE;
  }

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &params.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &params.num_procs);

  buf = init(size, &params);

  // dem letzten prozess mitteilen wann zu stoppen ist
  if (params.rank == ROOT)
  {
    MPI_Send(&buf[0], 1, MPI_INT, params.num_procs - 1, TAG_END_VALUE, MPI_COMM_WORLD);  
  }
  else if (params.rank == (params.num_procs - 1))
  {
    MPI_Recv(&endValue, 1, MPI_INT, ROOT, TAG_END_VALUE, MPI_COMM_WORLD, &params.status);
  }

  if (params.rank == ROOT)
  {
    printf("\nBEFORE\n");
  }

  printComplete(buf, &params);
  MPI_Barrier(MPI_COMM_WORLD);

  circle(buf, endValue, &params);
  MPI_Barrier(MPI_COMM_WORLD);

  if (params.rank == ROOT)
  {
    printf("\nAFTER\n");
  }
  printComplete(buf, &params);

  MPI_Finalize();
  free(buf);
  free(params.bufsize);

  return EXIT_SUCCESS;
}
