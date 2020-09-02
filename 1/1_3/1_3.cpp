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
int Nthr,Nblk, Lblk;
double d,L;
//Main
double teta,y_c;
//Averages
double blk_av,blk_norm,err;
double glob_av,glob_av2,stima;
double Nh;  //Number of hit

//Functions
void Input(void);
void Sample(void);
void Averages(int,int);
double Error(double,double,int);

int main()
{  
   Input();
   for(int iblk=1;iblk<=Nblk;iblk++)
   {
      Nh=0;
      for(int j=0;j<Lblk;j++)
      {
         Sample();    
         if(y_c < double(L) / 2 * sin(teta) ) Nh +=1; //Accumulate if hit the line
      }
      Averages(iblk,Lblk);
   }
   return 0;
}

void Input(void)
{
  //Lettura input
  ifstream InFile;
  InFile.open("./3.in");
  InFile >> Nthr; //Numero di lanci
  InFile >> Nblk; //Numero di blocchi
  InFile >> d; //Distanza tra linee orizzontali
  InFile >> L; //Lunghezza ago
  InFile.close();

  Lblk=int(Nthr/Nblk);

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

void Sample(void)  //campiona coord y e angolo teta dell'ago
{ 
  double x_cam, y_cam;
  double r_cam;

  //Campionamento cord y del centro dell'ago
   y_c = double(d) / 2 * rnd.Rannyu();
   teta = 0;
   //Campionamento angolo dell'ago rispetto dir orizzontale
   while(teta == 0)
   {
      x_cam = -1. + 2 * rnd.Rannyu();
      y_cam = -1. + 2 * rnd.Rannyu();
      r_cam= sqrt( x_cam *x_cam + y_cam * y_cam);
      if(r_cam < 1 && y_cam>0 && x_cam>0)
      {
         teta = acos( x_cam / r_cam ) ;
      }
   }

   return;
}

void Averages(int iblk, int Lblk) //Print results for current block
{

   ofstream Pi;
   const int wd=12;

    cout << "Block number " << iblk << endl;

    //<r>
    Pi.open("3Pi.out",ios::app);
    stima = (2 * L * Lblk ) / ( d * Nh );
    glob_av  += stima;
    glob_av2 += stima*stima;
    err=Error(glob_av,glob_av2,iblk);
    Pi << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av/(double)iblk << setw(wd) << err << endl;
    Pi.close();


    cout << "----------------------------" << endl << endl;

    return;
}


double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

