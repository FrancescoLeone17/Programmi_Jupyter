#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "random.h"

using namespace std;

//Random numbers
int seed[4];
Random rnd;
//Input
const int M=100;  //Numero di step nel RW 
int nblk;   //Numero di blocchi
int nstep=100;  //Numero di step nel blocco
int key;
//Move
double rw[3]={},rw_old[3]={};   //Posizione corrente, posizione precedente
//Averages
double blk_av[M],blk_norm,stima;
double glob_av[M],glob_av2[M],err;

//Functions
void Input();
void Averages(int);
void Accumulate(int);
void Reset(int);
double Error(double,double,int);
void Move(void);

int main() 
{
    Input();
    for(int iblk=1; iblk<=nblk;iblk++) //Ciclo sui blocchi
    {  
       Reset(iblk);
       for(int istep=1;istep<=nstep; istep++)  //Ciclo sugli step di blocco
         {  
            for(int k=0;k<3;k++)  rw_old[k] = 0;  //Cond iniziale
            for(int irw=0;irw<M;irw++) //Ciclo sugli step del random walk
            {
               Move();
	       Accumulate(irw);
	       for(int k=0;k<3;k++)  rw_old[k] = rw[k];
            }
         }
         Averages(iblk);
     }

  return 0;
} 

void Input(void)
{
  //Lettura input
  ifstream inFile;
  inFile.open("2.in");
  inFile >> nblk;
  inFile >> key; //if zero rw on a lattice, otherwise in a random direction
  inFile.close();

   //Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();

   return;
}

void Reset(int iblk)
{
   if(iblk == 1)
   {  
      for(int k=0;k<M;k++)
      {
         glob_av[k]  =  0;
         glob_av2[k] =  0;
      }
   }

   for(int k=0;k<M;k++) blk_av[k] = 0;

   return;
}

void Move(void)
{ 
    //Random walk on a lattice	
    if(key==0)
    {
       for(int k=0;k<3;k++)
       {
          if(rnd.Rannyu()<0.5)  rw[k] = rw_old[k] + 1.;
          else rw[k] = rw_old[k] - 1.;
       }
    }

    //Random walk on random direction
    else
    {
      double teta = rnd.Rannyu()* M_PI;
      double phi =  2 * M_PI * rnd.Rannyu();
      rw[0] = rw_old[0] + sin(teta) * cos(phi);
      rw[1] = rw_old[1] + sin(teta) * sin(phi);
      rw[2] = rw_old[2] + cos(teta);
    }
    return; 
}

void Accumulate(int irw)
{    
      double r=0;
      for(int k=0;k<3;k++) r = r + rw[k] * rw[k];

      blk_av[irw] = blk_av[irw] + r;

   return;
}


void Averages(int iblk) //Print results for current block
{
   ofstream RW;
   const int wd=12;

    cout << "Block number " << iblk << endl;

    //Accumulating global props
    for(int k=0;k<M;k++)
    {
       stima = sqrt( blk_av[k] / nstep);
       glob_av[k]  += stima;
       glob_av2[k] += stima * stima;
    }
    
    //Write final values
    if(iblk==nblk)
    {
       RW.open("rw.out",ios::app);
       RW << setw(wd) << 0 <<  setw(wd) << 0 << setw(wd) << 0 << endl;
       for(int k=0;k<M;k++)
       {
          err=Error(glob_av[k],glob_av2[k],iblk);
          RW << setw(wd) << k+1 <<  setw(wd) << glob_av[k]/(double)iblk << setw(wd) << err << endl;
       }
       RW.close();
    }

    cout << "----------------------------" << endl << endl;

    return;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

