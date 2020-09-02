#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // Exp
#include <algorithm>    // min function
#include <iomanip>      // setw
#include "1.h"

using namespace std;

int main()
{
   Input();
   for(int iblk=1;iblk<=nblk;iblk++)
   {
      Reset(iblk);
      for(int istep=1;istep<=nstep;istep++)
      {
         Metropolis();
         Accumulate();
      }
      Averages(iblk);
   }

   ConfFinal();
  return 0;
}



void Input(void)
{
   ifstream InFile;
   cout << "Campionamento funzione d'onda di H" << endl;
   cout << "Stato n=1, m=0, l=0" << endl << endl;


   InFile.open("input.dat");
   //Input blocking
   InFile >> nblk;
   InFile >> nstep;
   cout << "Numero di blocchi   " << nblk << endl;
   cout << "Lunghezza di blocco   " << nstep << endl;
  
   //Input Metropolis
   InFile >> par;  //Regola accettazione
   InFile >> key1; //Se 0 conf proposta in modo uniforme, altrimenti conf proposta in modo gauss
   InFile >> key2; //Se 0 campiono il ground state, altrimenti lo stato eccitato
   InFile.close();

   //Config iniziale
   InFile.open("conf.0");
   InFile >> rold[0];
   InFile >> rold[1];
   InFile >> rold[2];
   InFile.close();

   cout << "Dato iniziale del Metropolis   " << endl;
   cout << "x   " << rold[0] << endl;
   cout << "y   " << rold[1] << endl;
   cout << "z   " << rold[2] << endl;
 
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


void Metropolis(void)
{
   //Campionamento da matrice di trasferimento
   if(key1==0) Uniform();
   else Gaussian();
   //Conversione
   Spherical_Coord();
   //accept/reject
   if(key2==0) Gs();
   else Ex();
  
   attempted = attempted + 1;
  return;
}

//Uniform Sampling
void Uniform(void)
{
   for(int i=0;i<3;i++) rnew[i] = rold[i] + ( rnd.Rannyu() - 0.5 ) * par;
   return;
}

//Gaussian Sampling
void Gaussian(void)
{
   for(int i=0;i<3;i++) rnew[i] = rold[i] + par * rnd.Gauss(0,1);
   return;
}

//Accept/reject with gs function
void Gs(void)
{
   double alfa;

   alfa= min(1. , exp(2 * Rho) / exp(2 * rho1) );
   if(rnd.Rannyu() < alfa)
   {
      Rho = rho1;
      for(int i=0;i<3;i++) rold[i] = rnew[i];
      accepted = accepted + 1;
   }
   return;
}

//accept/reject with excited state function
void Ex(void)
{  
  double alfa;
  
  alfa= min( 1. ,  exp(-rho1) * pow( rnew[2], 2) / ( exp( -Rho ) * pow( rold[2], 2 ) ) );
  if(rnd.Rannyu() < alfa)
  {
     Rho = rho1;
     for(int i=0;i<3;i++) rold[i] = rnew[i];
     accepted = accepted + 1;
  }
   return;	
}

//Conversione from x,y,z to sqrt(x^2+y^2+z^2)
void Spherical_Coord(void)
{
  rho1 = 0;
  for(int i=0;i<3;i++) rho1 = rho1 + rnew[i] * rnew[i];
  rho1 = sqrt( rho1 );

  Rho = 0;
  for(int i=0;i<3;i++) Rho = Rho + rold[i] * rold[i];
  Rho = sqrt( Rho );

  return;
}

void Reset(int iblk)
{
   if(iblk == 1)
   {
      glob_av  =  0;
      glob_av2 =  0;
   }

   blk_av = 0;
   accepted = 0;
   attempted = 0;

   blk_norm = 0;
   return;
}

void Accumulate(void)
{
   blk_av = blk_av + Rho;
   blk_norm = blk_norm + 1;

   return;
}

void Averages(int iblk) //Print results for current block
{

   ofstream Ave;
   const int wd=12;

    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate" << accepted/attempted << endl; 

    //Averages
    Ave.open("ave.out",ios::app);
    stima = blk_av / blk_norm;
    glob_av  += stima;
    glob_av2 += stima*stima;
    err=Error(glob_av,glob_av2,iblk);
    Ave << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av/(double)iblk << setw(wd) << err << endl;
    Ave.close();

    cout << "----------------------------" << endl << endl;

    return;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void ConfFinal(void)
{
   ofstream OutFile;
   OutFile.open("conf.final");

   for(int i=0;i<3;i++)  OutFile << rold[i] << endl;
   
   OutFile.close();   
   return;
}
