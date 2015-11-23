// pointer to classes example

#include <iostream>
#include <ctime>
#include <fstream>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
//multivariate normal libraries
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <cmath>
#include "modelclasses.h"
using namespace std;

//extern std::default_random_engine generator;



int main(int argc, char **argv) {
  

  //generator.seed(time(0));

  
  double duration;
  double tot_y;
  int n;
  int index;
  string filename = "test3.txt";
  Parameters *param;  
  if (argc==5){
    cout<< "four arguments imported. \n";
    tot_y=atof(argv[1]);
    n = atoi(argv[2]);
    filename = argv[3];
    index = atoi(argv[4]);
    param = new Parameters(index);

  } else if(argc==4){
    cout << "three arguments imported. \n";
    tot_y=atof(argv[1]);
    n = atoi(argv[2]);
    filename = argv[3];
    param = new Parameters();

  } else if(argc==3){
    cout << "two arguments imported. \n";
    tot_y=atof(argv[1]);
    n = atoi(argv[2]);
    param = new Parameters();
  } else {
    tot_y = 120.0;
    n = 1000;
    param = new Parameters();
  }
  double u=0;
  double b=0;
  double pnet = 0;
  param->setUB(&u, &b, &pnet);
  cout << "random u , b and n" << u << " " << b << " " << pnet <<"\n";
  cout << "random gamma example: " << param->gamma_dist(1.0,1.0) << "\n";
  cout << "random poisson example: " << param->poisson_dist(100.0) << "\n";
  cout << "random gaussian example: " << param->normal_dist(1.0,0.1) << "\n";
  cout << "random gamma example inverse method" << gsl_cdf_gamma_Pinv(param->uniform_dist(),1.0,1.0) << "\n";
  std::clock_t start;
  
  
  //define new model with params pop, a,b for gamma-distribution.
  start = std::clock();
  cout << "Initialising array... \n";
  //Model m1(n,2.0,1.0);
  Model m1(n,param);
  m1.evolveAndSaves(tot_y,filename);
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout << "Finished simulation in " << duration << "s" << "\n";
  
  cout << "write to file... \n";
  m1.save("test2.txt");

  //Host * pop;
  //pop = new Host[n];
  //one host time-series example.
  //std::gamma_distribution<double> distribution(2.0,3.0); //2.0, 1.0
  //Host h1(distribution(generator));
  //ofstream tsfile;
  //tsfile.open("tsexample.txt");
  //double t = 0.0;
  //while(t < (tot_y * 365.0)){
   // h1.react();
   // t +=1;
   // tsfile << h1.W << " " << h1.WM << " " << h1.WF << " "  << h1.M << " " << h1.births << " " << h1.deaths <<"\n";
  //} 
  //tsfile.close();

 
  
  delete param;
  
  return 0;
}
