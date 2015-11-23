#ifndef CLASS_H
#define CLASS_H

#include <iostream>
#include <ctime>
#include <fstream>
#include <string>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
//multivariate normal libraries
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

using namespace std;
int poisson(double lambda);
double expTrunc(double, double);

class Parameters { //store all parameters in one class to stop them being everywhere. Will also make it easier to change parameters when needed.
  gsl_rng * rando;
  const gsl_rng_type * T;
  double lbda_original,v_to_h_original,sig_original; //original parameters stored so when changed from intervention still have in memory.
  public:
    Parameters(); //imports parameters from text file.
    Parameters(int); //imports parameters from text file and use index to seed random number generator.
    Parameters(double,double); //bite exposure gamma a, bite exposure gamma b
    Parameters(double,double,double); //bite exposure gamma a, bite exposure gamma b, v_to_h ratio.
    ~Parameters(void);
    //functions
    void setParameters(void);
    void importParameters(void);
    void setBR(bool intervention);
    void setVH(bool intervention);
    void setMu(bool intervention);
    double L3Uptake(double);
    double readParamLine(ifstream& myfile);
    double gamma_dist(double, double);
    double normal_dist(double, double);
    double uniform_dist(void);
    double expTrunc(double, double);
    int poisson_dist(double);
    int rmvnorm(const gsl_rng *, const int, const gsl_vector *, const gsl_matrix *, gsl_vector *);
    int setUB(double *, double *, double *);
    double immuneRK4Step(double,double);
    //parameters
    double a;
    double b;
    double shapeRisk; //shape of probability distribution for risk.
    double riskMu1; //mean of risk for age group less than 16
    double riskMu2; //mean of risk for age group  17-29
    double riskMu3; //mean of risk for age group  30+
    double mu; //death rate of worms
    double theta; //immune system response parameter. 0.112
    double gamma; //mf death rate
    double alpha; //mf birth rate
    double v_to_h; //vector to human ratio.
    double xi; //constant bite rate (5760.0 / 30.0)
    double kappas1; //vector uptake and development
    double r1; //vector uptake and development
    double tau; 
    double z; //waning immunity
    double nu; //poly-monogamy parameter   
    //shared state parameters    
    double L3; //larvae density.
    double g; //Proportion of mosquitoes which pick up infection when biting an infected host
    double sig; //death rate of mosquitos
    double psi1; //Proportion of L3 leaving mosquito per bite
    double psi2; //Proportion of L3 leaving mosquito that enter host
    double lbda; //bite rate per mosquito per month
    double dt;
    double lbdaR, v_to_hR, sigR; //reduction in bite-rate and v_to_h ratio.
    int nMDA; //number of rounds of MDA.
    int mdaFreq;
    double covMDA;
    double s2;
    double mfPropMDA; //proportion of mf killed due to MDA
    double wPropMDA; //proportion of worms killed due to MDA
    double sysComp,u0Comp, sigComp; //parameters for systematic non-compliance
    double sysCompN, u0CompN, sigCompN; //parameters for systematic non-compliance of bed-nets.
    double covN; //coverage of bednets.
    double rhoCN; //correlation of non-compliance between bed-nets and MDA (Chemotherapy).
    int mosquitoSpecies; //species of mosquito 0 - Anopheles, 1 - Culex
    double rhoBU; //correlation between non-compliance and bite risk.
    int aWol; //using doxycycline as intervention 0- not used. 1- is used.
    double sN,dN; //probability of mosquito successfully biting given use of bed nets. prob of death due to bed net.
    int IDAControl; //switch to IDA after five rounds of normal MDA treatment. 0- not used, 1 - is used.
} ;

class Host {
    
    
    //double dt; //time-step for tau.
    void initialise(void); //initialise state variables.
    
    
  public:
    
    Host(double);
    Host(void);
    void react(void);
    void evolve(double);
    double mfConc(void);
    double biteRate(void);
    double repRate(void);
    void updateRisk(bool);
    //state variables
    int W,WM,WF; //number of worms, number of male worms, number of female worms
    double I; //immune response (assumed to be deterministic)
    double M; //mf produced. //int
    double t; //time of host
    double a; //age of host.
    int aWol; //host is treated with Doxycycline.
    int bedNet; //host used bednet.
    //host parameters
    Parameters *param; //pointer to params.
    double b; //risk of infection
    double pRisk; //p-value for risk of infection (updates b)
    double uComp; //random parameter used to calculate compliance with MDA
    double uCompN; //random parameter used to calculate compliance with bed nets.
    //const double mu = 0.1/(30.0); //death rate of worms
    //double beta = 0.112; //immune system response parameter.
    //double gamma = 0.1/(30.0);
    //double alpha = 0.2/(30.0);
    int births;
    int deaths;
    //double xi = (5760.0 / 30.0) * 0.414 * 0.32 * 0.2;
    //double tau = 0.04 / (365.0); //0.02
    //double L3 = 1.0; //keeps track of L3 in total population.
} ;


class Model {

    Host *host_pop;
    int n;  
    bool aWolInt; //undergoing aWol intervention.
    int bedNetInt; //undergoing bednet intervention.
  public:
    Model(int,double,double);
    Model(int, Parameters*);
    Parameters *param;
    void react(void);
    void evolve(double);
    void save(string);
    void outputs(void);
    void saveOngoing(string); //appends WM,WF,M to the bottom of filename
    void evolveAndSaves(double tot_t, string filename); //functionally same as evolve, but saves outputs at regular intervals.
    double L3(void);
    void MDAEvent(void);
    void aWolEvent(bool);
    void bedNetEvent(void);
    
    
} ;

#endif
