#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // Exp
#include <iomanip>
#include <vector>
#include <algorithm>    // random_shuffle, rotate
#include "TSP.h"        // global variable

using namespace std;


int main()
{ 
   Input();               
   CityPosition(key);
   InitPerc();
   for(int ipr=0;ipr<Npr;ipr++) l1[ipr]=Fitness(pop[ipr]);
   Order();
    
   ifstream TempStep;
   TempStep.open("TempStep.dat"); 

   for(int j=0;j<N_temp;j++) //ciclo sulle temperature
   {  
      TempStep >> Temp >> nstep; 
      for(int istep=0; istep < nstep; ++istep)  //ciclo sugli step del SA
      {  
         pop_new=pop;     
         for(int ipr=0; ipr<Npr; ++ipr)  //ciclo sui percorsi
         {   
            v_new=pop[ipr];
            Mutation(v_new); 
      	    Metropolis(j,ipr);
         }
	 pop=pop_new;
         for(int k=0;k<Npr;k++) l1[k]=Fitness(pop[k]);  
         Order();              
         Averages(istep);     
      }    
   }
   BestPath();
  return 0;
}

//Funzioni per inizializzare
void Input(void)
{

  ifstream ReadInput;
  cout << "Traveling Salesman Problem" << endl ;
  cout << "Numero di città:   " << Nct << endl;
  cout << "Numero di percorsi:   "  << Npr << endl << endl;
  
  //Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();

  cout << "Lettura input" << endl;
  ReadInput.open("input.dat");
  
  ReadInput >> N_temp;
  ReadInput >> pb_mut;
  ReadInput >> key;

  ReadInput.close();
  
  return; 
}

void CityPosition(int key)
{
  double r;
  double teta,cx,cy;

  //Città su un cerchio
  if( key == 0 )
    {
       for(int i=0; i<Nct; i++)
       {
          r = rnd.Rannyu(); //Modified
	  teta = 2*pi * r;

	  c_pos[i][0] = cos(teta);
	  c_pos[i][1] = sin(teta);  
       }
    }
 
  //Città in un quadrato
  else
    {
     for(int i=0; i<Nct; i++)
       {  

          r = rnd.Rannyu(); //Modified
          cx = 2*r - 1.;

	  r = rnd.Rannyu(); //Modified
	  cy = 2*r - 1.; 

          c_pos[i][0] = cx;
          c_pos[i][1] = cy;
       }
    }  

  return;
}


void InitPerc(void)
{ 
   vector <int>  tmp;
    cout << "Inizializza i percorsi" << endl;

    //Matrix inizialization
    for(int i = 0; i < Npr; i ++) pop[i].resize(Nct,0);
    for(int i = 0; i < Npr; i ++) pop_new[i].resize(Nct,0);  

   //Generazione dei percorsi
   for(int j=0; j<Nct; j++) pop[0][j] = j;
    
   for(int j=0; j<Npr; j++)
     {
        pop[j] = pop[0];
        Perm(1, Nct, pop[j]);
     }
   
   //Check percorsi generati corretti
   for(int j=0; j<Npr; j++)   Check(pop[j]); 
   
   cout << "-------------------------------" << endl; 
   return;
}

//Controlla che il primo elemento sia zero
//Controlla che la città sia chiamata una e una sola volta

void Check(vector<int> & vett)
{
   int index;
   vector <int> tmp1,tmp2;  

   tmp1.resize(Nct,0);
   tmp2.resize(Nct,0);
   for(int i=0; i<Nct; i++) tmp1[i]=i;
   
   //Primo check
   if( vett[0] != 0 )
   {
      cout << "Errore: Generato percorso non lecito: Check 1" << endl;
   }

   //Secondo check
   for( int i=1; i<Nct; i++)
   {
      index = vett[i];
      tmp2[index] = index;
   }
    
   if(tmp1 != tmp2)
   {
      cout << "Errore: Generato percorso non lecito: Check 2" << endl;
   }
   
   return;
}


//Fitness Function
double Fitness(vector <int> & vect)
{
  double d;
  double fit=0;

  //Calcolo funzioni costo
  for(int j=0; j<Nct; ++j )  //Ciclo fino a penultimo step
  {
     d = Distance(j,vect);
     fit += d;
  }
  
  return fit;
}

double Distance(int j,vector <int> & vect)
{
   double dist=0;
   int ind1 = vect[j], ind2 = vect[j+1]; //Indice della città

   if( j == Nct-1)
   {
      for(int k=0; k<2; ++k )
      {
         dist += pow( c_pos[ind1][k] - c_pos[0][k], 2);
      }

      return sqrt(dist);
   }

   else
   {
      for(int k=0; k<2; ++k )
      {
         dist += pow( c_pos[ind1][k] - c_pos[ind2][k], 2);
      }

      return sqrt(dist);
   }
}

void Order(void)
{  
   vector< pair <double, vector <int> > > zip;
   
   //zipvector inizialization
   for(int i=0; i<Npr; i++) zip.push_back( make_pair( l1[i], pop[i]));

   //Ordering l2
   sort(zip.begin(),zip.end());

   for(int i=0; i<Npr; i++)
   {
      l1[i] = zip[i].first;
      pop[i] = zip[i].second;
   }
}

//Mutation Operator
void Mutation(std::vector<int> & vect)
{
  double r;
  int a,b,c;
  vector <int> tmp;
  
  //Pair Perm
  r = rnd.Rannyu(); //Modified
  if(r<pb_mut)
  {
     a = rnd.Rannyu(1,31); //Modified
     b = rnd.Rannyu(1,31); //Modified
     iter_swap(vect.begin()+a, vect.begin()+b);
  }

  //Shift di 4 città per 5 posizioni
  r = rnd.Rannyu(); //Modified
  if(r<pb_mut)
  {
     a = rnd.Rannyu(1,23); //Modified
     b = a+4;
     c = b+5;
     rotate(vect.begin()+a, vect.begin()+b, vect.begin()+c);
  }

  //Permutation of 4 cities with other 4 cities
  r = rnd.Rannyu(); //Modified
  if(r<pb_mut)
  {
     r = rnd.Rannyu(); //Modified
     a = rnd.Rannyu(1,19); //Modified //19
  
     tmp.resize(8,0);    //Conversion //8
     for(int i=0;i<4;i++) tmp[i]=vect[a+i]; //4
     for(int i=7;i<11;i++) tmp[i-3]=vect[a+i]; // 7   11   3
     Perm(0,8,tmp);      //Perm
  
     for(int i=0;i<4;i++) vect[a+i]=tmp[i];  //Conversion inversa
     for(int i=7;i<11;i++) vect[a+i]=tmp[i-3];
  }
  //Inversion
  r = rnd.Rannyu(); //Modified
  if(r<pb_mut)
  {
     a = rnd.Rannyu(1,27); //Modified //era 27
     b = a+5;              //era 5
     reverse(vect.begin()+a,vect.begin()+b);
  }

  Check(vect);
  return;
}

void Perm(int a, int b, vector<int> & vvv)
{
   random_shuffle ( vvv.begin()+a, vvv.begin()+b, myrandom);  // [a, b)
   return;
}

int myrandom (int i) { return rnd.Rannyu(0,i);} //Modify

void Averages(int istep)
{
   ofstream BsL1;  //L1 valued over the best path
   const int wd=12;

    BsL1.open("BsL1.out",ios::app);

    BsL1 << setw(wd) << istep <<  setw(wd) << l1[0] << endl;
   
     return;
}

void BestPath(void)
{
   ofstream BsPath;
   const int wd=12;
   int ind1; 
   
   BsPath.open("BsPath.out");

   for(int i=0;i<Nct;i++)
   {  
      ind1=pop[0][i];	   
      BsPath << setw(wd) << c_pos[ind1][0] <<  setw(wd) << c_pos[ind1][1]  << endl;
   }
   BsPath << setw(wd) << c_pos[0][0] <<  setw(wd) << c_pos[0][1]  << endl;

   BsPath.close();

   cout << "Accettazione   " << acc/att << endl;
   return;
} 

//Attenzione alle variabili
void Metropolis(int j,int ipr)
{
   vector <int> v_old(Nct);
   double pb;
   
   l1_new=Fitness(v_new);
   pb = exp(-l1_new/Temp) / exp(-l1[ipr]/Temp);

   if(pb > rnd.Rannyu())
   {
      for(int k=0;k<Nct;k++) pop_new[ipr][k] = v_new[k];
      acc = acc + 1;
   }
   att = att + 1;
   return;
}

