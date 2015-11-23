#include <iostream>
#include <ctime>
#include <fstream>
#include <string>
#include <cmath>
#include "modelclasses.h"
//------------------------------------------------------------//
//------------------random number generation-----------------//
//----------------------------------------------------------//
//std::default_random_engine generator;
//std::uniform_real_distribution<double> ud(0.0,1.0);

//---------------------------------------------------//
//------------------functions misc.-----------------//
//-------------------------------------------------//

//--------------------------------------------------------//
//------------------Host class functions-----------------//
//------------------------------------------------------//
Host::Host (double x){
  b = x;
  t = 0;
 
  initialise();
}

Host::Host(void){
  t = 0;
  
  initialise();
}

void Host::react( void ) {
  //generator.seed(time(0));
  //std::poisson_distribution<int> birth_dist(5.0);
  //std::poisson_distribution<int> death_dist(5.0);
  double bNReduction = 1 - (1-param->sN) * (double) bedNet;
  //immune state update
  
  //I +=  (param->dt) *( (double) W - param->z * I);  
  I = param->immuneRK4Step((double) W , I);
  //male worm update
  births = param->poisson_dist(0.5 * bNReduction * param->xi  * biteRate() * param->L3 * exp(-1 * param->theta * I) * b *  param->dt); //exp(-1 * beta * I)
  //births = param->poisson_dist(0.5 * param->xi  * biteRate() * param->L3 * exp(-1 * param->theta * I) * b *  param->dt); //exp(-1 * beta * I)
  deaths = param->poisson_dist( param->mu * (1.0+ (double) aWol) * (double) WM * param->dt );
  WM += births - deaths;

  //female worm update
  births = param->poisson_dist(0.5 * bNReduction * param->xi * biteRate() * param->L3 * exp(-1 * param->theta * I) * b  * param->dt); //* exp(-1 * beta * I)
  //births = param->poisson_dist(0.5  * param->xi  * biteRate() * param->L3 * exp(-1 * param->theta * I) * b *  param->dt); //exp(-1 * beta * I)
  deaths = param->poisson_dist( param->mu * (1.0+ (double) aWol) * (double) WF * param->dt );
  WF += births - deaths;

  //Mf update
  //births = poisson(param->alpha * WF * WM);
  //deaths = poisson(param->gamma * M);
  //M += births - deaths;
  M += param->dt * (repRate() * (1.0 - (double) aWol) - param->gamma * M);
  //M += param->dt * (repRate() - param->gamma * M);
  //total worm count
  W = WM+WF;
  //time-step
  t += param->dt;
  a += param->dt;
  //ensure all positive state variables remain positive
  if (W < 0){ W=0; }
  if (WM < 0){ WM = 0; }
  if (WF < 0){ WF = 0; }
  if (I < 0){ I = 0.0; }
  if (M < 0){ M = 0.0; }
  //simulate event where host dies and is replaced by a new host.
  if (param->uniform_dist() < (1 - exp(-1 * param->tau * param->dt) ) || a>1200.0){ //if over age 100 
    initialise();
    a = 0; //birth event so age is 0.
    
  }
}



void Host::evolve(double tot_t ){
  while(t<tot_t){
    react();  
  }
}

void Host::initialise(void){
  
  W = 0;
  WM = 0;
  WF = 0;
  I = 0.0;
  M = 0.0; //0
  bedNet = 0;
  aWol=0;
}

double Host::mfConc(void){ //return concentration of mf as opposed to absolute number
return (double) M;// * 0.005;

}

double Host::biteRate(void){
  
  if (a < 108.0){ //less than 9 * 12 = 108.0
    return a/108.0; 
  }  else {
    return 1.0;  
  }
    
}

double Host::repRate(void){
  if (param->nu == 0){

    if (WM>0){
      return (double) WF;
    } else {
      return 0.0;
    }
   
  } else {
    return param->alpha * std::min((double) WF ,(1/param->nu) * (double) WM);
  }
}

void Host::updateRisk(bool init=false){
  if(init){
    if(a<240.0) { //age under 20
      b = gsl_cdf_gamma_Pinv(pRisk,param->shapeRisk,param->riskMu1/param->shapeRisk);  
    } else if(a < 348.0) { //age between 20 and 29
      b = gsl_cdf_gamma_Pinv(pRisk,param->shapeRisk,param->riskMu2/param->shapeRisk);
    } else { //age 30+
      b = gsl_cdf_gamma_Pinv(pRisk,param->shapeRisk,param->riskMu3/param->shapeRisk);
    }
  } else {
    if(a<12.0) { //age at 0
      b = gsl_cdf_gamma_Pinv(pRisk,param->shapeRisk,param->riskMu1/param->shapeRisk);  
    } else if(a >= 240.0 and a < 252.0) { //age 25
      b = gsl_cdf_gamma_Pinv(pRisk,param->shapeRisk,param->riskMu2/param->shapeRisk);
    } else if(a>= 360.0 and a < 372.0){ //age 30+
      b = gsl_cdf_gamma_Pinv(pRisk,param->shapeRisk,param->riskMu3/param->shapeRisk);
    }
  }
}
//---------------------------------------------------------//
//------------------Model class functions-----------------//
//-------------------------------------------------------//


Model::Model(int size, double a, double b){ //deprecated, need to link to parameters class
  n = size;
  //std::gamma_distribution<double> b_dist(a,b);
  //std::exponential_distribution<double> e_dist(0.112); //0.112
  //std::exponential_distribution<double> age_dist( 0.04);
  //std::gamma_distribution<double>::param_type p{a,b};
  //cool_dist.param(p);
  host_pop = new Host[n];
  for(int i =0; i < n; i++){
    host_pop[i].b = param->gamma_dist(a,b);
    //host_pop[i].beta = 0.112;//e_dist(generator);
    host_pop[i].a = param->expTrunc(param->tau,100.0) * 12.0;//age_dist(generator) * 365.0;
  }
}

Model::Model(int size, Parameters *newparam){
  double a,b;  
  double sU=0;
  double sB=0;
  double sN = 0;
  
  n = size;
  param = newparam;
  a = param-> a;
  b = param-> b;  
  //std::gamma_distribution<double> b_dist(a,b);
  //std::exponential_distribution<double> e_dist(0.112); //0.112
  //std::exponential_distribution<double> age_dist( param->tau * 12.0); //*365.0 as tau is in (days)^-1 and need it in (yrs)^-1
  //std::gamma_distribution<double>::param_type p{a,b};
  //cool_dist.param(p);
  host_pop = new Host[n];
  for(int i =0; i < n; i++){
    host_pop[i].b = param->gamma_dist(a,b);//this is unneccessary now that b is updated via update risk method. Keeping for legacy.
    
    host_pop[i].a = param->expTrunc(param->tau,100.0) * 12.0;//age_dist(generator) * 365.0;
    host_pop[i].param = param;
    //host_pop[i].pRisk = param->uniform_dist(); //create prisk
    //host_pop[i].updateRisk(true); //update risk with initialisation true.
    //host_pop[i].uComp = param->normal_dist(param->u0Comp,param->sigComp); //set individual's compliance parameter u_i
    param->setUB(&sU, &sB, &sN);
    host_pop[i].pRisk = sB;
    host_pop[i].updateRisk(true); //update risk with initialisation true.
    host_pop[i].uComp = sU;
    host_pop[i].uCompN = sN;
    host_pop[i].bedNet = 0;
  }
}

void Model::react(void){ //deprecated. does not include vector dynamics!!!!
  for(int i =0; i < n; i++){
    host_pop[i].react();
  }
}

void Model::evolve(double tot_t){ //deprecated. does not include vector dynamics!!!!

  for(int i =0; i < n; i++){
    host_pop[i].evolve(tot_t * 12.0);
  }
}

void Model::evolveAndSaves(double tot_t, string filename){
double t = 0;
int icount = 0;
double maxMDAt = 1200.0;
double maxoldMDAt; //used in triple drug treatment. 
aWolInt = false;
bedNetInt = 0;
if (param->aWol==0){
  maxMDAt = (1200.0+ (double) param->nMDA * (double) param->mdaFreq);
  if(param->IDAControl==1){ //if switching to IDA after five treatment rounds.
    maxoldMDAt = (1200.0+ 5.0 * (double) param->mdaFreq);
  } else {
    maxoldMDAt = 2*maxMDAt; //this just makes maxoldMDAt larger than total treatment time so there is never a switch.
  }
} else {
  maxMDAt = (1200.0+ 36.0); //if aWol intervention then effects only last two years.
}
//double currentL3 = 0.5;
cout << "mosquito species: " << param->mosquitoSpecies << "\n";
param->L3 = 5.0;
//open file and save a blank line to make sure file is empty before appending to it
ofstream mf;
mf.open(filename.c_str());
mf << " ";
mf.close();
cout << "0----------100\n-";
  while(t< tot_t * 12.0){ //for 100 years update annually, then update monthly when recording and intervention is occuring.
    if (t<960.0){ //1200.0
      param->dt = 12.0;
    }else{
      param->dt = 1.0;
    }
    for(int i =0; i < n; i++){
      host_pop[i].react();
    }
    //update
    t = host_pop[0].t; 

    param->L3=L3();    
    if(((int) t % 2 == 0) && (t < floor(t) +param->dt)){
          //cout << "t = " << (double) t/12.0 << "\n";
          saveOngoing(filename);
    }
    if ( ((int) t % (int) (tot_t * 12.0/10.0) == 0) && (t < floor(t) + param->dt)){ //every 10% of time run.
      cout << "-";    
    }
    if (t>=1200.0 && t < 1200.0 +param->dt){ //events that occur at start of treatment after 100 years.
      cout << "bednet event at " << t;
      bedNetEvent();
      bedNetInt = 1;
    }

    if (((int) t % param->mdaFreq == 0) && (t < floor(t) + param->dt)){ //things that need to occur annually
      if(t>maxoldMDAt){
        param->mfPropMDA = 0.0;
        param->wPropMDA = 0.0;
      }
      if( (t>1200.0) && (t<= maxMDAt) ){ //if after one hundred years and less than 125 years.
        if(param->aWol==0){
          MDAEvent();  
        } else {
          aWolEvent(true);
          aWolInt = true; //switch so know that aWol intervention is occurring and to not update hosts that may not have been treated intially.
        }
        param->setBR(true); //intervention true.
        param->setVH(true);
        param->setMu(true);
      } else {
        param->setBR(false); //intervention false.
        param->setVH(false);
        param->setMu(false);
        if (param->aWol==1){
          aWolEvent(false);
          aWolInt = false;
        }
      }

      for(int i =0; i < n; i++){
        host_pop[i].updateRisk(); //update risk annually as don't need to check this as frequently as reactions.
      }  
    }
    icount++;
  }
cout << "\n";
}

void Model::save(string filename){
  ofstream mf;
  mf.open (filename.c_str());
  for(int i = 0; i < n; i++){
    mf << host_pop[i].W << " " << host_pop[i].WM << " " << host_pop[i].WF << " " << host_pop[i].M << "\n";
  }
  mf.close();
}

void Model::saveOngoing(string filename){
  ofstream mf;
  mf.open (filename.c_str(), std::ios_base::app);
  for(int i = 0; i < n; i++){
    mf << host_pop[i].WM << " ";
  }
  mf << "\n";
  for(int i = 0; i < n; i++){
    mf << host_pop[i].WF << " ";
  }
  mf << "\n";
  for(int i = 0; i < n; i++){
    mf << host_pop[i].M << " ";
  }
  mf << "\n";
  for(int i = 0; i < n; i++){
    mf << param->L3 << " ";
  }
  mf << "\n";
  for(int i = 0; i < n; i++){
    mf << host_pop[i].a << " ";
  }
  mf << "\n";
  for(int i = 0; i < n; i++){
    mf << host_pop[i].I << " ";
  }
  mf << "\n";
  mf.close();
}

void Model::outputs(void){
  for(int i=0; i < n; i++){
    cout << host_pop[i].b << " ";  
  }
  cout << "\n";
}

double Model::L3(void){
  double mf = 0.0;
  double bTot = 0.0;
  for(int i=0; i < n; i++){
    //mf += param->kappas1 * pow(1 - exp(-param->r1 *( host_pop[i].mfConc() * host_pop[i].b)/param->kappas1), 2.0);  
    mf += host_pop[i].b * param->L3Uptake(host_pop[i].mfConc());
    bTot += host_pop[i].b;
  }
  mf = mf / bTot; //(double) n;
  return mf * (1 + (double) bedNetInt * param->covN * (param->sN - 1)) * param->lbda * param-> g /(param->sig + param->lbda * param->psi1);
  
}

void Model::MDAEvent(void){
    for(int i = 0; i<n; i++){
      if (param->normal_dist(host_pop[i].uComp,1.0)<0){    //param->uniform_dist()<param->covMDA
        host_pop[i].M = param->mfPropMDA * host_pop[i].M;
        host_pop[i].WM = (int) floor(param->wPropMDA * (double) host_pop[i].WM ); 
        host_pop[i].WF = (int) floor(param->wPropMDA * (double) host_pop[i].WF ); 
      }
    }
}

void Model::bedNetEvent(void){
    param->sig = param->sig + param->lbda * param->dN * param->covN;
    for(int i = 0; i<n; i++){
      if (param->normal_dist(host_pop[i].uCompN,1.0)<0){    //param->uniform_dist()<param->covMDA
        host_pop[i].bedNet = 1; //using bed-net
      } else {
        host_pop[i].bedNet = 0; //not using bed-net
      }
    }
}

void Model::aWolEvent(bool intervention){
  for(int i = 0; i<n; i++){
    if(intervention & !aWolInt){
      if (param->normal_dist(host_pop[i].uComp,1.0)<0){    //param->uniform_dist()<param->covMDA
        host_pop[i].M = 0.0;
        host_pop[i].WM = (int) floor(0.1 * (double) host_pop[i].WM ); 
        host_pop[i].WF = (int) floor(0.1 * (double) host_pop[i].WF ); 
        host_pop[i].aWol = 1; 
      }
    } else if (!intervention){
      host_pop[i].aWol = 0;
    }
  }
}
//--------------------------------------------------------------//
//------------------Parameters class functions-----------------//
//------------------------------------------------------------//

Parameters::Parameters(double ai, double bi) {
   //gsl random stuff
   gsl_rng_env_setup();
   T = gsl_rng_default;
   rando = gsl_rng_alloc(T);
   gsl_rng_set(rando,time(0));

   setParameters();
   //param input
   a = ai;
   b = bi;


}

Parameters::Parameters(double ai, double bi, double v_to_h) {
   //gsl random stuff
   gsl_rng_env_setup();
   T = gsl_rng_default;
   rando = gsl_rng_alloc(T);
   gsl_rng_set(rando,time(0));
   setParameters();
   //param input
   a = ai;
   b = bi;
   xi = lbda * v_to_h * 0.414 * 0.32 * 0.2; //constant bite rate (5760.0 / 30.0)
  
}

Parameters::Parameters(){ //constructor where all parameters are imported from parameters.txt
  //param input
  importParameters();

  //gsl random stuff
   gsl_rng_env_setup();
   T = gsl_rng_default;
   rando = gsl_rng_alloc(T);
   gsl_rng_set(rando,time(0));
   

}

Parameters::Parameters(int index){ //constructor where all parameters are imported from parameters.txt and use index to seed random number generator
  importParameters();

  //gsl random stuff
   gsl_rng_env_setup();
   T = gsl_rng_default;
   rando = gsl_rng_alloc(T);
   gsl_rng_set(rando,time(0)+index);
   

}

Parameters::~Parameters(void){
  cout << "destructing random number generator \n";  
  gsl_rng_free(rando);
}

void Parameters::setParameters(void){
  riskMu1 = 0.53;
  riskMu2 = 0.29;
  riskMu3 = 0.29;
  shapeRisk = 1.2;
  //param no input
   mu = 0.1; //death rate of worms
   theta = 0.0;//0.001 //immune system response parameter. 0.112
   gamma = 0.1; //mf death rate
   alpha = 0.58; //mf birth rate per fertile worm per 20 uL of blood.
   lbda = 10; //number of bites per mosquito per month.
   xi = lbda * (500.0) * 0.414 * 0.32 * 0.2; //constant bite rate (5760.0 / 30.0)
   kappas1 = 4.395; //vector uptake and development anophelene
   r1 = 0.055; //vector uptake and development anophelene
   tau = 0.02 / (12.0); 
   z = 0.0; //waning immunity
   nu = 0.0; //poly-monogamy parameter   
    //shared state parameters    
   L3 = 0.0; //larvae density.
   g = 0.37; //Proportion of mosquitoes which pick up infection when biting an infected host
   sig = 5.0; //death rate of mosquitos
   psi1 = 0.414; //Proportion of L3 leaving mosquito per bite
   psi2 = 0.32; //Proportion of L3 leaving mosquito that enter host
   s2 = 0.2; //Prop of L3 that enter human host developing into adult.
   dt = 1.0;
}

double Parameters::readParamLine(ifstream& myfile){
    string line;
    getline (myfile,line);
    cout << "reading line :" << line << '\n';
    string token = line.substr(0, line.find(':'));
    cout << "token :" << token << '\n';
    return stod(token);
}

void Parameters::importParameters(void){
  ifstream myfile ("parameters.txt");
  if (myfile.is_open())
  {
    
    riskMu1 = readParamLine(myfile);
    riskMu2 = readParamLine(myfile);
    riskMu3 = readParamLine(myfile);
    shapeRisk = readParamLine(myfile);
    mu = readParamLine(myfile);
    theta = readParamLine(myfile);
    gamma = readParamLine(myfile);
    alpha = readParamLine(myfile);
    lbda = readParamLine(myfile);
    v_to_h = readParamLine(myfile);
    kappas1 = readParamLine(myfile);
    r1 = readParamLine(myfile);
    tau = readParamLine(myfile);
    z = readParamLine(myfile);
    nu = readParamLine(myfile);
    L3 = readParamLine(myfile);
    g = readParamLine(myfile);
    sig = readParamLine(myfile);
    psi1 = readParamLine(myfile);
    psi2 = readParamLine(myfile);
    dt = readParamLine(myfile);
    lbdaR = readParamLine(myfile);
    v_to_hR = readParamLine(myfile);
    nMDA = (int) readParamLine(myfile);
    mdaFreq = (int) readParamLine(myfile);
    covMDA = readParamLine(myfile);
    s2 = readParamLine(myfile);
    mfPropMDA = 1 - readParamLine(myfile);
    wPropMDA = 1 - readParamLine(myfile);
    sysComp = readParamLine(myfile);
    mosquitoSpecies = (int) readParamLine(myfile);
    rhoBU = readParamLine(myfile);
    aWol = (int) readParamLine(myfile);
    sigR = readParamLine(myfile);
    covN = readParamLine(myfile);
    sysCompN = readParamLine(myfile);
    rhoCN = readParamLine(myfile);
    IDAControl = readParamLine(myfile);
    myfile.close();
    //calculate other parameters
    lbda_original = lbda;
    v_to_h_original = v_to_h;
    sig_original = sig;
    xi = lbda * v_to_h * psi1 * psi2 * s2; //constant bite rate (5760.0 / 30.0)
    a = shapeRisk; //shape parameter (can vary)
    b = 1/a; //scale parameter determined so mean is 1.
    //bed net parameters
    sN = 0.03;
    dN=0.41;
    //non-compliance parameters for MDA
    sigComp = sqrt(sysComp/(1+sysComp));
    u0Comp = -1*gsl_cdf_gaussian_Pinv(covMDA, 1.0) * sqrt(1 + sigComp * sigComp);

    //non-compliance parameters for bed-nets
    sigCompN = sqrt(sysCompN/(1+sysCompN));
    u0CompN = -1*gsl_cdf_gaussian_Pinv(covN, 1.0) * sqrt(1 + sigCompN * sigCompN);


  }
}

void Parameters::setBR(bool intervention){
  if (intervention){
    lbda = lbdaR * lbda_original;
    xi = lbda * v_to_h * psi1 * psi2 * s2;
  } else {
    lbda = lbda_original;
    xi = lbda * v_to_h * psi1 * psi2 * s2;
  }
}

void Parameters::setVH(bool intervention){
  if (intervention){
    v_to_h = v_to_hR * v_to_h_original;
    xi = lbda * v_to_h * psi1 * psi2 * s2;
  } else {
    v_to_h = v_to_h_original;
    xi = lbda * v_to_h * psi1 * psi2 * s2;
  }
}

void Parameters::setMu(bool intervention){
  if (intervention){
    sig = sigR; //increase mortality due to bed nets. dN = 0.41 death rate 
  } else {
    sig = sig_original;
  }  
}

double Parameters::L3Uptake(double mf){
  if(mosquitoSpecies==0){
    return kappas1 * pow(1 - exp(-r1 *( mf )/kappas1), 2.0);
  }
  else
  {
    return kappas1 * ( 1 - exp(-r1 *( mf )/kappas1) );
  }
}

double Parameters::gamma_dist(double a, double b){
  
  return gsl_ran_gamma(rando,a,b);
}

int Parameters::poisson_dist(double rate){
  return gsl_ran_poisson(rando,rate);
}

double Parameters::uniform_dist(void){
  return gsl_ran_flat(rando,0,1);
}

double Parameters::expTrunc(double lambda, double trunc){
  return (-1/ lambda)* log(1- uniform_dist() *(1 - exp(- lambda * trunc)) );
}

double Parameters::normal_dist(double mu, double sigma){
  return mu + gsl_ran_gaussian (rando, sigma);
}

int Parameters::rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result){
/* multivariate normal distribution random number generator */
/*
* n dimension of the random vetor
* mean  vector of means of size n
* var variance matrix of dimension n x n
* result  output variable with a sigle random vector normal distribution generation
*/
int k;
gsl_matrix *work = gsl_matrix_alloc(n,n);

gsl_matrix_memcpy(work,var);
gsl_linalg_cholesky_decomp(work);

for(k=0; k<n; k++)
  gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
gsl_vector_add(result,mean);

gsl_matrix_free(work);

return 0;
}

int Parameters::setUB(double *u, double *b, double *n){
  //create correlation matrix and mean for multivariate normal distribution.
  if (sigComp==0){ //if sigComp zero then no correlation and u can only be one value.
    *u = u0Comp;
    *b = uniform_dist();
    *n = normal_dist(u0CompN,sigCompN);
  }else{
    double b0; //this is the normally-distributed result that needs converting into the uniform bite-risk.
    gsl_matrix * m = gsl_matrix_alloc (3, 3);
    gsl_matrix_set (m, 0, 0, 1.0);
    gsl_matrix_set (m, 0, 1, sigComp * rhoBU);
    gsl_matrix_set (m, 1, 0, sigComp * rhoBU);
    gsl_matrix_set (m, 1, 1, sigComp * sigComp);
    gsl_matrix_set (m, 2, 0, sigComp * rhoCN * sigCompN);
    gsl_matrix_set (m, 0, 2, sigComp * rhoCN * sigCompN);
    gsl_matrix_set (m, 2, 1, 0.0);
    gsl_matrix_set (m, 1, 2, 0.0);
    gsl_matrix_set (m, 2, 2, sigCompN * sigCompN);

    //for (int i = 0; i < 3; i++){  
    //  for (int j = 0; j < 3; j++){
    //    printf ("m(%d,%d) = %g\n", i, j, gsl_matrix_get (m, i, j));
    //  }
    //}

    //set mean vector
    gsl_vector * v = gsl_vector_alloc (3);
    gsl_vector_set (v, 0, 0.0);
    gsl_vector_set(v,1, u0Comp);
    gsl_vector_set(v,2, u0CompN);

    //set result vector
    gsl_vector * result = gsl_vector_alloc (3);
    rmvnorm(rando, 3, v, m, result);

    b0 = (double) gsl_vector_get (result, 0);
    *u = (double) gsl_vector_get (result, 1);
    *n = (double) gsl_vector_get (result, 2);
    
    *b = gsl_cdf_gaussian_P (b0, 1.0);
    gsl_matrix_free (m);
    gsl_vector_free (v);
  }
  return 0;
}

double Parameters::immuneRK4Step(double W,double I){
  double k1,k2,k3,k4;
  k1 = dt * (W - z*I);
  k2 = dt * (W - z *(I + 0.5 * k1));
  k3 = dt * (W - z *(I + 0.5 * k2));
  k4 = dt * (W - z *(I + k3));
  return I + 0.1666667 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

