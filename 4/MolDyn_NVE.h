/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=200;
int n_props;
int iv,ik,it,ie;
int igofr; //Aggiunta
double bin_size,nbins; //Aggiunta
double stima_pot, stima_kin, stima_etot, stima_temp;
double walker[m_props]; //Aggiunta

//Aggiunta
// averages
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double err_pot, err_kin, err_etot, err_temp, err_gdir;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;
//Aggiunta
int nblk;

//pigreco
const double pi=3.1415927; //Aggiunta

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
//Aggiunta
void Reset(int);
void Accumulate(void);
void Averages(int);
double Error(double,double,int);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
