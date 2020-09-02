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
int nblk,nstep;
//Sample
double x_u,x_is; //x_u sampled uniform , x_is sampled by imp of sampling
const double norm = pow( 1 - pow( M_PI , 2 ) / 24 + pow( M_PI , 4 ) / 1920, -1); //normalizazione dello sviluppo di Taylor cos
//Averages
double blk_av[2],blk_norm,stima;
double glob_av[2],glob_av2[2],err;

//Functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
double Error(double,double,int);
double Tcos(double);

int main() 
{
  Input();
  
  for(int iblk=1;iblk<=nblk;iblk++) 
  {
     Reset(iblk);
     for(int istep=1;istep<=nstep;istep++) Accumulate();
     Averages(iblk);
  }

  return 0;
}

void Input(void)
{
  //Read input
  ifstream inFile;
  inFile.open("1.in");
  inFile >> nblk;     //numero di blocchi
  inFile >> nstep;    //numero di step per blocco
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

void Accumulate(void)
{  
   //Uniform Sampling
   x_u = rnd.Rannyu() * M_PI  / 2;
   blk_av[0] = blk_av[0] + M_PI / 2 * cos(x_u);

   //Importance of sampling
   int continua=0;
   while(continua==0)
   {
      x_is = M_PI / 2 * rnd.Rannyu();
      if(rnd.Rannyu() < Tcos(x_is))
      { 
         blk_av[1] = blk_av[1] + M_PI/2 *cos(x_is) / ( norm * Tcos(x_is));
	 continua = 1;
      }
   } 

   //Accumulo
   blk_norm = blk_norm + 1;
   return;
}

void Reset(int iblk)
{
   if(iblk == 1)
   {
      glob_av[0]  =  0;
      glob_av2[0] =  0;
      glob_av[1]  =  0;
      glob_av2[1] =  0;
   }

   blk_av[0] = 0;
   blk_av[1] = 0;

   blk_norm = 0;

   return;
}

void Averages(int iblk) //Print results for current block
{

   ofstream Ave_u,Ave_is;
   const int wd=12;

    cout << "Block number " << iblk << endl;

    //Uniform Sampling
    Ave_u.open("1unif.out",ios::app);
    stima = blk_av[0] / blk_norm;
    glob_av[0]  += stima;
    glob_av2[0] += stima*stima;
    err=Error(glob_av[0],glob_av2[0],iblk);
    Ave_u << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[0]/(double)iblk << setw(wd) << err << endl;
    Ave_u.close();

    //Importance of Sampling
    Ave_is.open("1is.out",ios::app);
    stima = blk_av[1] / blk_norm;
    glob_av[1]  += stima;
    glob_av2[1] += stima*stima;
    err=Error(glob_av[1],glob_av2[1],iblk);
    Ave_is << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[1]/(double)iblk << setw(wd) << err << endl;
    Ave_is.close();

    
    cout << "----------------------------" << endl << endl;

    return;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

//Sviluppo del coseno al quart'ordine
double Tcos(double var)
{
  double cl;

  cl = 1 - 0.5 * pow( var, 2) + 1/(double)24 * pow( var, 4);
  return cl;
}

