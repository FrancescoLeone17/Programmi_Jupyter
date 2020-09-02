#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // Exp
#include <iomanip>
#include "VMC.h"        // global variable

using namespace std;

int main()
{
  Input();

  for(int iblk=1;iblk <= nblk;iblk++)
  {
     Reset(iblk);
     for(int istep=1;istep <= nstep; istep++)
     {
        Move();
	Accumulate();
	if(istep % 10 == 0) PrintConf();
     }
     Averages(iblk);
  } 
  ConfFinal();

  return 0;
}


void Input(void){
  ifstream InFile,Conf;
  cout << "Calcolo Energia di stato fondamentale" << endl;
  cout << "1 particella soggetta a potenziale" << endl;
  cout << "V(x) = x^4 - 2.5*x^2" << endl << endl;
  
  //Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();

  cout << "Lettura input";
  InFile.open("input.dat");

  //Input blocking
  InFile >> nblk;
  InFile >> nstep;
  cout << "Numero di blocchi   " << nblk << endl;
  cout << "Lunghezza di blocco   " << nstep << endl;
  
  //Input Metropolis
  InFile >> par;  //Regola accettazione
  
  InFile >> mu;
  InFile >> sigma;
  cout << "parametro mu:  " << mu << endl;
  cout << "parametro sigma:  " << sigma << endl << endl;

  InFile.close();

  Conf.open("config.0");
  Conf >> x;
  cout << "Dato iniziale del Metropolis   " << endl;
  cout << "x:   " << x << endl;
  Conf.close();

  return;
}


void Move(void){
  double alfa;
  double xold,xnew;

  //Old
  xold = x;
  
  //New
  xnew = xold + ( rnd.Rannyu() - 0.5 ) * par;
  
  //Metropolis Test
  alfa = pow( psi(xnew) / psi(xold), 2) ;

  if(alfa > rnd.Rannyu())
  {
     x = xnew;
     acc += 1;
  }

  att += 1;

  return;
}

void Accumulate(void)
{
   blk_av = blk_av + LocalEn(x);
   blk_norm = blk_norm + 1;
   return;
}

void Averages(int iblk)
{
   ofstream En;
   const int wd=12;

    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << acc/att << endl << endl;

    En.open("output.ene.0",ios::app);

    stima_en = blk_av/blk_norm; 
    glob_av += stima_en;
    glob_av2 += stima_en*stima_en;
    err_en=Error(glob_av,glob_av2,iblk);

    En << setw(wd) << iblk <<  setw(wd) << stima_en << setw(wd) << glob_av/(double)iblk << setw(wd) << err_en << endl;

   return;
}

double LocalEn(double x)
{
   double Eloc;
   double Ekin,Epot;
   double Ekin1,Ekin2;
   double y_m,y_p;

   y_m = pow((x - mu) / sigma, 2);
   y_p = pow((x + mu) / sigma, 2);
 
   Ekin1 = exp(- y_m / 2 ) * 1 / pow(sigma,2) * ( y_m - 1 ); 
   Ekin2 = exp(- y_p / 2 ) * 1 / pow(sigma,2) * ( y_p - 1 );
   Ekin = - 0.5 * (Ekin1 + Ekin2);  // -P^2/(2M) psi

   Epot = (pow(x,4) - 2.5 * pow(x,2)) * psi(x); // V(x) psi

   Eloc = (Ekin + Epot) / psi(x);

   return Eloc;
}

void Reset(int iblk) //Reset block averages
{
   if(iblk == 1)
   {
      glob_av = 0;
      glob_av2 = 0;
   }

   blk_av = 0;
   blk_norm = 0;
   att = 0;
   acc = 0;

   return;
}

double psi(double crd)
{
  double fz=0;
  double y_m, y_p;

  y_m = pow((crd - mu) / sigma, 2);
  y_p = pow((crd + mu) / sigma, 2);

  fz = exp(- y_m / 2 ) + exp(- y_p / 2 );

  return fz;
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

void ConfFinal(void){
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;

  WriteConf.open("config.final");
  WriteConf << x << endl;
  WriteConf.close();

  //rnd.SaveSeed();
}

void PrintConf(void)
{
   ofstream Conf;

   Conf.open("output.conf.0",ios::app);
   Conf << x << endl; 
   Conf.close(); 

   return;
}	
