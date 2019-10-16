
#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(){
    srand ( time ( NULL));
    int gen,hits = 0;
    long int MAX_GEN = 1000000000;
    for (;gen<MAX_GEN;gen++) {
        // https://itp.tugraz.at/MML/MonteCarlo/MCIntro.pdf
        double random_x = (double)rand()/RAND_MAX; //float in range 0 to 1
        double random_y = (double)rand()/RAND_MAX;
        double random_point = pow(random_x,2) + pow(random_y,2);
        hits = (random_point < 1) ? ++hits : hits;
    }
    double pi = ((double)hits/(double)gen)*4;
    printf("%d %d %f",hits, gen, pi);
    // TODO Send to rank 1
    // TODO Make avarage of all
}