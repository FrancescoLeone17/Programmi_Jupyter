#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include "random.h"
using namespace std;

//Input
const int M=1000000;
//Random numbers
int seed[4];
Random rnd;
//main
const double pi= acos(double(-1));
double x_0[M],x_1[M],x_2[M]; //Rispettivamente num distr casuale,exp e lorentzi
double r;
//Functions
void Input(void);
void Sample(void);
void S(int);

int main(){

   Input();
   Sample();
   S(1);
   S(2);
   S(10);
   S(100);

  return 0;
}

void Input()
{
   //Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   return;
   input.close();
}

void Sample()
{
  for(int i=0;i<M;i++)
  {                
     r = rnd.Rannyu();
     x_0[i] = r;  //sampling uniform distribution
     x_1[i] = -log( 1 - r);   //sampling exponential distribution
     x_2[i] = tan( pi * ( r - 0.5) );   //sampling Cauchy-Lorentz distribution
  }
   return;
}
void S(int N)
{
   double sum=0;
   const int wd=12;
   int k;
   ofstream  OutFile;

   string str="S_";
   string a=to_string(N);
   str.append(a);
   OutFile.open(str);

   //S_N for uniform distribution
   for(int i=0;i<10000;i++)
   {
      for(int j=0;j<N;j++)
      {
         k= j + i * N;
         sum += x_0[k];
      }
      OutFile << setw(wd) << sum / N << endl;
      sum = 0;
   }


   //S_N for exponential distrib
   for(int i=0;i<10000;i++)
   {
      for(int j=0;j<N;j++)
      {
         k= j + i * N;
         sum += x_1[k];
      }
      OutFile << setw(wd) << sum / N << endl;
      sum = 0;
   }

   //S_N for Cauchy distrib
   for(int i=0;i<10000;i++)
   {
      for(int j=0;j<N;j++)
      {
         k= j + i * N;
         sum += x_2[k];
      }
      sum = 0;
      OutFile << setw(wd) << sum / N << endl;
   OutFile.close();
   }
}
