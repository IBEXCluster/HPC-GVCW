#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define MPI_MAX_PROCESSOR_NAME 150
int main(int argc, char** argv) 
{
    char processor_name[MPI_MAX_PROCESSOR_NAME], rank[10];
    int world_size, world_rank, name_len;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Get_processor_name(processor_name, &name_len);
    printf("This Job exected from processor %s, rank %d out of %d processors\n", processor_name, world_rank, world_size);
    sprintf(rank, "sh ./Phase2.sh  %d", world_rank);
    system(rank);
    MPI_Finalize();
}
