//Random numbers
#include "random.h"
int seed[4];
Random rnd;
//Input 
int nblk,nstep;
double par;
int key1,key2;
//Metropolis
double rold[3],Rho=0; //Vecchia conf e il suo modulo
double rnew[3],rho1=0; //Nuova conf e il suo modulo
double accepted=0,attempted=0;
//Averages
double blk_av, blk_norm, stima;
double glob_av,glob_av2,err;
//Function
void Input(void);
void Reset(int);
void Metropolis(void);
void Uniform(void);
void Gaussian(void);
void Gs(void);
void Ex(void);
void Spherical_Coord(void);
void Accumulate(void);
void Averages(int);
double Error(double,double,int);
void ConfFinal(void);
