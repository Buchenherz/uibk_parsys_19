
#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(){
    double random_value;
    srand ( time ( NULL));
    int gen,hits = 0;
    long long int MAX_GEN = 10000000000;
    for (;gen<MAX_GEN;gen++) {
        random_value = (double)rand()/RAND_MAX*sqrt(2); //float in range 0 to sqaure root of 2
        hits = (random_value > 1) ? hits : ++hits;
    }
    double pi = ((double)hits/(double)gen) *4;
    printf("%d %d %f",hits, gen, pi);
    // TODO Send to rank 1
    // TODO Make avarage of all
}