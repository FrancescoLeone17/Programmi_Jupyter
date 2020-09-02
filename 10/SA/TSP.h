//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//Input 
const int Nct=32; 		      //Numero di città
const int Npr=500;                    //Please pari
int N_temp; 		      //Numero di temperature
int N_term;
int key;                              //Se 0 città su cerchio, se 1 città su un quadrato
int nstep;         	              //Numero di step che algoritmo viene ripetuto
double Temp;                         //Temperatura
double pb_mut;

//Inizialization
std::vector <std::vector<int>> pop(Npr);    // vettore della popolazione
double c_pos[32][2];                        // c[i][k] posizione della i-esima città con dir k

//Fitness
std::vector <double> l1(Npr);    // l1

//Metropolis
double acc=0,att=0;
std::vector <int> v_new(Nct);                     // Nuovo percorso
std::vector <std::vector<int>> pop_new(Npr);      // vettore della nuova popolazione
double l1_new;

//pigreco
const double pi=3.1415927;

//Functions
void Input(void);
void CityPosition(int);
void Perm(int,int, std::vector<int> &);
int myrandom(int);
void InitPerc(void);
void Check(std::vector<int> &);
double Fitness(std::vector <int> & );
double Distance(int,std::vector <int> &);
void Order(void);
void Mutation(std::vector<int> &);
void Averages(int);
void BestPath(void);
void Metropolis(int,int);
