#include<mpi.h>
#include<stdio.h>

intmain(intargc,char**argv){
    MPI_Init(&argc,&argv);//initializetheMPIenvironment
    intsize;
    MPI_Comm_size(MPI_COMM_WORLD,&size); //getthe numberofranks
    intrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank); //getthe rankofthecaller
    printf("Helloworldfromrank%dof%d\n",rank,size);
    double random_value;
    srand ( time ( NULL));
    int gen,hits = 0;
    int MAX_GEN = 1000000000;
    for (;gen<MAX_GEN;gen++) {
        random_value = (double)rand()/RAND_MAX*sqrt(2); //float in range 0 to sqaure root of 2
        hits = (random_value > 1) ? hits : hits++;
    }
    double pi = hits/gen;
    // TODO Send to rank 1
    // TODO Make avarage of all
    MPI_Finalize();//cleanup
}