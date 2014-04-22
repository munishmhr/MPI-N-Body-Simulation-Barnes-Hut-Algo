#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char** argv){
    int myid, size;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    //printf("NO of processors are : %d",size);
    int ar[12] = {1,2,3,4,5,6,7,8,9,10,11,12};
    int i ;

    int count = 12/size;
    int j = myid;
    int localdata[size];
    int g = 0;
    for(i = (count*myid); i<(count*myid)+count; i++){
        localdata[g] = ar[i] + myid;
        printf("\nprinting array by %d procesor id %d",localdata[g],myid);
        g = g+1;
    }
    MPI_Gather(&localdata, 4, MPI_INT, ar, 4, MPI_INT, 0, MPI_COMM_WORLD);
    i = 0;
    if(myid == 0){
        for(i = 0 ; i <12;i++){
            printf("\nfrom my id %d : %d",myid,ar[i]);
        }
    }
    MPI_Bcast(ar,12,MPI_INT,0,MPI_COMM_WORLD);
    for(i=0;i<12;i++){
        printf("\nfrom my id %d : %d",myid,ar[i]);
    }

    MPI_Finalize();
return 0;
}
