#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "1.h"

using namespace std;

int main() {

  Input();
  
  for(int iblk=1;iblk<=nblk;iblk++)
  { 
     Reset(iblk);	  
     for(int istep=1;istep<=nstep;istep++)
     {  
        Accumulate();
     }
     Averages(iblk);
  }
  
  return 0;
} 

void Input(void)
{
   //Lettura input
  ifstream inFile;
  inFile.open("./1.in");
  inFile >> nblk; //numero di blocchi
  inFile >> nstep; //numero di step per blocco
  inFile >> nInt; //numero di intervalli
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

void Reset(int iblk) //Reset block averages
{
   if(iblk == 1)
   {  
      for(int k=0;k<4;k++)
      {
         glob_av[k]  =  0;
         glob_av2[k] =  0;
      }
   }

   for(int k=0;k<4;k++) blk_av[k] = 0;

   blk_norm = 0;

   return;
}


void Accumulate(void)
{  
    double S,C,P;
    double S1=S0,S2;
    double z;
    double deltaT = T / (double)nInt;

    //direct sampling	
    z = rnd.Gauss(0.,1.);
    S = S0 * exp( ( r - 0.5 * sigma * sigma) * T + sigma * z * sqrt(T) );
    C = exp(-r * T ) * fmax(0., S-K);
    P = exp(-r * T) * fmax(0., K-S);
    //Accumulate
    blk_av[0] = blk_av[0] + C;
    blk_av[1] = blk_av[1] + P;

    //Discretized sampling
    for(int i=0;i<nInt;i++)
    {
       z = rnd.Gauss(0.,1.);
       S2 = S1 * exp( ( r - 0.5 * sigma * sigma) * deltaT + sigma * z * sqrt( deltaT ) );
       S1 = S2;
    }
    C = exp(-r * T) * fmax(0.,S1-K);
    P = exp(-r * T) * fmax(0.,K-S1);
    //Accumulate
    blk_av[2] = blk_av[2] + C;
    blk_av[3] = blk_av[3] + P;

    blk_norm = blk_norm + 1;
   return;
}

void Averages(int iblk) //Print results for current block
{

   ofstream C_dir, P_dir,C_discr, P_discr;
   const int wd=12;

    cout << "Block number " << iblk << endl;

    //C camp diretta
    C_dir.open("C_dir.out",ios::app);
    stima = blk_av[0] / blk_norm;
    glob_av[0]  += stima;
    glob_av2[0] += stima*stima;
    err=Error(glob_av[0],glob_av2[0],iblk);
    C_dir << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[0]/(double)iblk << setw(wd) << err << endl;
    C_dir.close();

    //P camp diretta
    P_dir.open("P_dir.out",ios::app);
    stima = blk_av[1]/blk_norm;
    glob_av[1]  += stima;
    glob_av2[1] += stima*stima;
    err=Error(glob_av[1],glob_av2[1],iblk);
    P_dir << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[1]/(double)iblk << setw(wd) << err << endl;
    P_dir.close();
 
    //C camp discr
    C_discr.open("C_discr.out",ios::app);
    stima = blk_av[2]/blk_norm;
    glob_av[2]  += stima;
    glob_av2[2] += stima*stima;
    err=Error(glob_av[2],glob_av2[2],iblk);
    C_discr << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[2]/(double)iblk << setw(wd) << err << endl;
    C_discr.close();

    //P camp discr
    P_discr.open("P_discr.out",ios::app);
    stima = blk_av[3]/blk_norm;
    glob_av[3]  += stima;
    glob_av2[3] += stima*stima;
    err=Error(glob_av[3],glob_av2[3],iblk);
    P_discr << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[3]/(double)iblk << setw(wd) << err << endl;
    P_discr.close();

    cout << "----------------------------" << endl << endl;

    return;
}


double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

