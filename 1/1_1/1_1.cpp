#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "random.h"
using namespace std;

//Random numbers
int seed[4];
Random rnd;
//Input
int nblk,nstep;
const int M = 100; //Numero di intervalli
double L= 1. / double(M); //Lunghezza di intervallo
//Averages 
double blk_av[2],blk_norm,stima; 
double glob_av[2],glob_av2[2];
double err;
//Chi
double n[M];
//Functions
void Input();
void Reset(int);
void Accumulate();
void Averages(int);
double Error(double,double,int);

int main() {
 
  Input();	
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
     Reset(iblk);   //Reset block average
     for(int istep=1; istep <= nstep; ++istep)
     {
        Accumulate(); //Update block average 
     }
     Averages(iblk);   //Write on file
  }
  
  return 0;
}

void Input()
{
  //Lettura input
  ifstream ReadInput;
  ReadInput.open("1_1.in");
  ReadInput >> nblk; //numero di blocchi
  ReadInput >> nstep; //numero di step interni al blocco
  ReadInput.close();
  
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

void Reset(int iblk) //Reset block averages
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
   
   for(int i=0;i<M;i++) n[i] = 0;

   blk_norm = 0;

   return;
}

void Accumulate(void) //Update block averages
{  
   int s;

   blk_av[0] = blk_av[0] + rnd.Rannyu();   // <r>
   blk_av[1] = blk_av[1] + pow( (rnd.Rannyu() - 0.5), 2);  // <(r-0.5)^2>
   
   //Chi^2
   s = int(rnd.Rannyu()/L);
   n[s] += 1; 

   blk_norm = blk_norm + 1.0;

   return;
}


void Averages(int iblk) //Print results for current block
{

   ofstream Ave_r, Ave_var, Chi;
   const int wd=12;

    cout << "Block number " << iblk << endl;
    
    //<r>
    Ave_r.open("1_1r.out",ios::app);
    stima = blk_av[0] / blk_norm;
    glob_av[0]  += stima;
    glob_av2[0] += stima*stima;
    err=Error(glob_av[0],glob_av2[0],iblk);
    Ave_r << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[0]/(double)iblk << setw(wd) << err << endl;
    Ave_r.close();
    
    //<(r-0.5)**2>
    Ave_var.open("1_1var.out",ios::app);
    stima = blk_av[1]/blk_norm;
    glob_av[1]  += stima;
    glob_av2[1] += stima*stima;
    err=Error(glob_av[1],glob_av2[1],iblk);
    Ave_var << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[1]/(double)iblk << setw(wd) << err << endl;
    Ave_var.close();
    
    //Chi^2
    double sum=0, c = (double) nstep / (double) M;
    Chi.open("1_1chi.out",ios::app);
    for(int i=0;i<M;i++)
    {
       sum += pow( (n[i]-c) , 2) / c;
    }
    Chi << setw(wd) << iblk <<  setw(wd) << sum << endl;
    Chi.close();

    cout << "----------------------------" << endl << endl;

    return;
}


double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

