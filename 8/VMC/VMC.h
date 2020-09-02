//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//Input 
int nblk,nstep;
double par;
double mu,sigma;

//Metropolis
double x;
double acc,att;

//Averages 
double glob_av, glob_av2;
double blk_av,blk_norm;
double stima_en,err_en;

//Function
void Input(void);
void Reset(int);
void Move(void);
double psi(double);
double LocalEn(double);
void Accumulate(void);
void Averages(int);
double Error(double,double,int);
void ConfFinal(void);
void PrintConf(void);
