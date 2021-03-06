#include "mpi.h"
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // Exp
#include <iomanip>
#include <vector>
#include <algorithm>    // random_shuffle, rotate
#include <string>       // to append in Averages
#include "TSP.h"        // global variable

using namespace std;


int main(int argc, char* argv[])
{ 
  //Parallelize code
  int size, rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Input(rank);
  CityPosition(key);
  InitPerc();
  Fitness();

for(int imigr=0;imigr<nmigr;imigr++)
  {
     Order();
     for(int istep=0;istep<nstep;++istep)
     {    
        for(int ipr=0;ipr<Npr;ipr=ipr+2)
        {   
           Selection();

           if(rnd.Rannyu()<0.7) Crossover(int(rnd.Rannyu(1,15)),v_c1,v_c2); 
           Mutation(v_c1);
           Mutation(v_c2);

           v_new[ipr] = v_c1;
           v_new[ipr+1] = v_c2;	 
        }
        v_old = v_new;
        Fitness();
        Order();
        Averages(istep,rank);
     }
     Migration(rank);
  }    
BestPath(rank); 
MPI_Finalize();
  return 0;
 
}

//Funzioni per inizializzare
void Input(int rank)
{

  ifstream ReadInput;
  
  //Read seed for random numbers
   int p1, p2;
   int ind_rank=2+rank;
   ifstream Primes("Primes");
   for(int i=0;i<ind_rank;i++) Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();

  ReadInput.open("input.dat");
  
  ReadInput >> nstep;
  ReadInput >> key;
  ReadInput >> nmigr;

  ReadInput.close();
  
  srand(37167316);
  return; 
}

void CityPosition(int key)
{
  double r;
  double cx,cy;
  ifstream CityPosition;


  //Città su un cerchio
  if( key == 0 )
    {    
       CityPosition.open("CityPosition_Circle.dat");
       for(int i=0;i<Nct;i++) 
       {
	  CityPosition >> c_pos[i][0] >> c_pos[i][1];
	  r = rnd.Rannyu();
       }
       CityPosition.close();
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

    //Matrix inizialization
    for(int i = 0; i < Npr; i ++) v_old[i].resize(Nct+1,0);
    for(int i = 0; i < Npr; i ++) v_new[i].resize(Nct+1,0);  

   //Generazione dei percorsi
   for(int j=0; j<Nct; j++) v_old[0][j] = j;
   v_old[0][Nct] = 0; 

   for(int j=0; j<Npr; j++)
     {
        v_old[j] = v_old[0];
        Perm(1, Nct, v_old[j]);
     }
   
   //Check percorsi generati corretti
   for(int j=0; j<Npr; j++)   Check(v_old[j]); 
   
   return;
}

//Controlla che il primo elemento sia zero
//Controlla che la città sia chiamata una e una sola volta
//Controlla che ultimo elem sia la prima città
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

   //Terzo check
   if( vett[Nct] != 0 )
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
void Fitness(void)
{
  double d;

  for(int i=0;i<Npr;i++)  l2[i]=0;
  for(int i=0;i<Npr;i++)  l1[i]=0;

  //Calcolo funzioni costo
  for(int i=0; i<Npr; ++i )
  {
     for(int j=0; j<Nct; ++j )  //Ciclo fino a penultimo step
     {
        d = Distance(i,j);
        l1[i] += d;
        l2[i] += d*d;
     }
  }
  
  return;
}

double Distance(int i, int j)
{
   double dist=0;
   int ind1 = v_old[i][j], ind2 = v_old[i][j+1]; //Indice della città


   for(int k=0; k<2; ++k )
   {
      dist += pow( c_pos[ind1][k] - c_pos[ind2][k], 2);
   }

   return sqrt(dist);
}

void Order(void)
{  
   vector< pair <double,vector<int>> > zip;
   
   //zipvector inizialization
   for(int i=0; i<Npr; i++) zip.push_back( make_pair( l2[i], v_old[i] ));

   //Ordering l2
   sort(zip.begin(),zip.end());

   for(int i=0; i<Npr; i++)
   {
      l2[i] = zip[i].first;
      v_old[i] = zip[i].second;
   }
   //Ordering l1
   sort(l1.begin(),l1.end());
}

//Selection Operators
void Selection(void)
{  
   int ic,jc;
   double r;  //Numero casuale tra 0 e 1
   double p=2; //Esponente
   

   //Scelta del primo percorso
   r = rnd.Rannyu();  //Modified
   ic = int(Npr * pow(r,p));

   //Scelta del secondo percorso
   r = rnd.Rannyu(); //Modified
   jc = int( Npr * pow(r,p) );
   
   //Coppia di percorsi scelta --> (vc1,vc2)
   for(int i=0;i<Nct;i++)
   {
     v_c1[i] = v_old[ic][i];
     v_c2[i] = v_old[jc][i];
   }
 
   Check(v_c1);
   Check(v_c2);
   return;
}


//Mutation Operator
void Mutation(std::vector<int> & vect)
{
  double r;
  double pb;
  int a,b,c;
  vector <int> tmp;
  
  //Pair Perm
  pb = 0.08;
  r = rnd.Rannyu(); //Modified
  if(r<pb)
  {
     a = int( rnd.Rannyu(1,31) ); //Modified
     b = int( rnd.Rannyu(1,31) ); //Modified
     iter_swap(vect.begin()+a, vect.begin()+b);
  }

  //Shift di 4 città per 5 posizioni
  pb = 0.08;
  r = rnd.Rannyu(); //Modified
  if(r<pb)
  {
     a = int( rnd.Rannyu(1,23) ); //Modified
     b = a+4;
     c = b+5;
     rotate(vect.begin()+a, vect.begin()+b, vect.begin()+c);
  }

  //Permutation of 4 cities with other 4 cities
  pb = 0.08;
  r = rnd.Rannyu(); //Modified
  if(r<pb)
  {
     r = rnd.Rannyu(); //Modified
     a = int( rnd.Rannyu(1,19) ); //Modified //19
  
     tmp.resize(8,0);    //Conversion //8
     for(int i=0;i<4;i++) tmp[i]=vect[a+i]; //4
     for(int i=7;i<11;i++) tmp[i-3]=vect[a+i]; // 7   11   3
     Perm(0,8,tmp);      //Perm
  
     for(int i=0;i<4;i++) vect[a+i]=tmp[i];  //Conversion inversa
     for(int i=7;i<11;i++) vect[a+i]=tmp[i-3];
  }
  //Inversion
  pb = 0.08;
  r = rnd.Rannyu(); //Modified
  if(r<pb)
  {
     a = int( rnd.Rannyu(1,27) ); //Modified //era 27
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

int myrandom (int i) { return int( rnd.Rannyu(0,i) );} //Modify

void Crossover(int m,vector <int> & v1, vector <int> & v2)
{
   int counter=0;
   int tmp_1[v1.size()], tmp_2[v2.size()];

   //Taglia primi m pezzi di v2 e aggiunge da v1
   for(unsigned int i=0; i<v1.size(); i++)
   {
      for(unsigned int j=m; j<v2.size(); j++)
      {
         if(v1[i] == v2[j])
         {
            tmp_2[m+counter] = v1[i];
            counter += 1;
            break;
         }
      }
   }

   counter=0;

   //Taglia primi m pezzi di v1 e aggiunge da v2
   for(unsigned int i=0; i<v2.size(); i++)
   {
      for(unsigned int j=m; j<v1.size(); j++)
      {
         if(v2[i] == v1[j])
         {
            tmp_1[m+counter] = v2[i];
            counter += 1;
            break;
         }
      }
   }

   for( unsigned int j=m; j<v1.size(); j++)
   {
      v1[j] = tmp_1[j];
      v2[j] = tmp_2[j];
   }

}

void Averages(int istep, int rank)
{
   ofstream AvL1,AvL2;  //L1 (or L2) averaged over the best path
   ofstream BsL1,BsL2;  //L1 (or L2) valued over the best path
   const int wd=12;
   double sum_cum1=0,sum_cum2=0; //Somma cumulativa
   int hp = Npr/2;
   string str1="AvL1_";
   string str2="AvL2_";
   string str3="BsL1_";
   string str4="BsL2_";
   string a= to_string(rank);

   str1.append(a);
   str2.append(a);
   str3.append(a);
   str4.append(a);

    AvL1.open(str1,ios::app);
    AvL2.open(str2,ios::app);
    BsL1.open(str3,ios::app);
    BsL2.open(str4,ios::app);
    
    for(int i=0;i< hp;i++)
    {
       sum_cum1 += l1[i];
       sum_cum2 += l2[i];
    }

    AvL1 << setw(wd) << istep <<  setw(wd) << sum_cum1 / hp  << endl;
    AvL2 << setw(wd) << istep <<  setw(wd) << sum_cum2 / hp  << endl;
    BsL1 << setw(wd) << istep <<  setw(wd) << l1[0] << endl;
    BsL2 << setw(wd) << istep <<  setw(wd) << l2[0] << endl;
   
     return;
}

void BestPath(int rank)
{
   ofstream BsPath;
   const int wd=12;
   int ind1; 
   string str="BsPath_";
   string a=to_string(rank);
   str.append(a);

   BsPath.open(str);

   for(int i=0;i<Nct;i++)
   {  
      ind1=v_old[0][i];	   
      BsPath << setw(wd) << c_pos[ind1][0] <<  setw(wd) << c_pos[ind1][1]  << endl;
   }
   BsPath << setw(wd) << c_pos[0][0] <<  setw(wd) << c_pos[0][1]  << endl;

   BsPath.close();
   return;
}


void Migration(int rank)
{  

   //Scelta degli accoppiamenti per le migrazioni	
   int couple[4]={};
   int counter=2;
  
   couple[0] = 0;
   couple[1] = rand()%3+1;
   for(int i=1; i<=3; i++)
   {
      if( i != couple[1] )
      {
         couple[counter] = i;
         counter += 1;
      } 
   }
    
   //Inizializzazione array da mandare e ricevere
   int tmp_0[32],tmp_1[32],tmp_2[32],tmp_3[32];

   //Inizio migrazioni
   MPI_Status stat1, stat2, stat3, stat4;
   int itag=1;
   int itag2=2;
   int itag3=3;
   int itag4=4;

   for(int j=0; j<40; j++)   //Migrazione dei primi 40 cromosomi della popolazione
   {  
      //Inizializza array da mandare e ricevere 
      for(int i=0;i<32; i++)
      {
         tmp_0[i] = v_old[j][i] ;
         tmp_1[i] = v_old[j][i] ;
         tmp_2[i] = v_old[j][i] ;
         tmp_3[i] = v_old[j][i] ;
      }

      //Migrazioni
      if(rank==couple[0])
      {
          MPI_Send(&tmp_1[0],Nct,MPI_INTEGER,couple[1],itag,MPI_COMM_WORLD);
          MPI_Recv(&tmp_0[0],Nct,MPI_INTEGER,couple[1],itag2,MPI_COMM_WORLD,&stat2);
      }

      else if(rank==couple[1])
      {
          MPI_Recv(&tmp_1[0],Nct,MPI_INTEGER,couple[0],itag,MPI_COMM_WORLD,&stat1);
          MPI_Send(&tmp_0[0],Nct,MPI_INTEGER,couple[0],itag2,MPI_COMM_WORLD);
      }

     // if(rank==couple[0]) for(int k=0;k<Nct;k++)  v_old[j][k] = tmp_0[k];
     // if(rank==couple[1]) for(int k=0;k<Nct;k++)  v_old[j][k] = tmp_1[k];
   
      if(rank==couple[2])
      {  
         MPI_Send(&tmp_3[0],Nct,MPI_INTEGER,couple[3],itag3,MPI_COMM_WORLD);
         MPI_Recv(&tmp_2[0],Nct,MPI_INTEGER,couple[3],itag4,MPI_COMM_WORLD,&stat4);
      }

      else if(rank==couple[3])
      {  
         MPI_Recv(&tmp_3[0],Nct,MPI_INTEGER,couple[2],itag3,MPI_COMM_WORLD,&stat3);
         MPI_Send(&tmp_2[0],Nct,MPI_INTEGER,couple[2],itag4,MPI_COMM_WORLD);
      }

   //   if(rank==couple[2]) for(int k=0;k<Nct;k++)  v_old[j][k] = tmp_2[k];
   //   if(rank==couple[3]) for(int k=0;k<Nct;k++)  v_old[j][k] = tmp_3[k];

   }
   cout << "Fine migrazione   " << endl;
   return;
}

