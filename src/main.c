#include "criteria.h"
#include "de.h"


void printDesign(int * a, int n, int m){
  printf("\n");
  for(int j = 0; j < m; ++j){
    for (int i = 0; i < n; ++i)  printf("%d ", a[i + j*n]);
    printf("\n");
  }

  printf("\n");
}

#define reps 1
#define n 70
#define m 7
#define NP 200

int main(int argc, char const *argv[])
{

   double vals[reps], timeTaken;
   int X[n*m];
   unsigned int seed = 12;
   static criteria PHI_FUNS[2] = {unipro, maxpro};
   DE_CC(n, m, n, NP, 1000, 0.2, 0.3, 0.9,
     reps, seed, vals, &timeTaken, X, 15, PHI_FUNS[0], 15);
   printf("\n");

   for(int i = 0; i < (reps > 6?6:reps); ++i) printf("%8.5f ",vals[i]);
   //printf("\n\n");
   printDesign(X, n, m);
   printf(" ...\n\n");
}
// how to compile
// gcc -I"/usr/share/R/include" -DNDEBUG -fopenmp -lgomp -fpic  -g -O2  prgressbar.c -lR -fopenmp -lm -o progress
// ./progress
