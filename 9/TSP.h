//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//Input 
const int Nct=32; //Numero di città
const int Npr=200;  //Please pari
int key;          //Se 0 città su cerchio, se 1 città su un quadrato
int nstep;        //Numero di step che algoritmo viene ripetuto

//Inizialization
std::vector <std::vector<int>> v_old(Npr);  // vettore della vecchia popolazione
std::vector <std::vector<int>> v_new(Npr);      // vettore della nuova popolazione
double c_pos[32][2];                        // c[i][k] posizione della i-esima città con la k-sima dir
std::vector <int> v_c1(Nct),v_c2(Nct);      // v_c1--> percorso della prima coppia selezionata, v_c2 della seconda

//Fitness
std::vector <double> l1(Npr);              // l1
std::vector <double> l2(Npr);              //l2

//pigreco
const double pi=3.1415927;

//Functions
void Input(void);
void CityPosition(int);
void Perm(int,int, std::vector<int> &);
int myrandom(int);
void InitPerc(void);
void Check(std::vector<int> &);
void Fitness(void);
double Distance(int,int);
void Order(void);
void Selection(void);
void Mutation(std::vector<int> &);
void Crossover(int,std::vector <int> &,std::vector <int> &);
void Averages(int);
void BestPath(void);

