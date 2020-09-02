//Random numbers
#include "random.h"
int seed[4];
Random rnd;
//data
double S0 = 100.;
double K = 100.;
double T = 1.;
double r = 0.1;
double sigma = 0.25;
double t=0;
double S=0;
//Input
int nblk,nstep,nInt;
//Averages
double blk_av[4],blk_norm,stima;
double glob_av[4],glob_av2[4],err;
//Functions
void Input(void);
void Averages(int);
double Error(double,double,int);
void Reset(int);
void Accumulate(void);
