/* This code uses FEM to upscale SMR Reduced Order Lab Scale Model to an Industrial Plug Flow Reactor. The model being used is defined by the text files 'plot_config.txt. The inflow values and other configuration detais such as temperature are grabbed from 'config_Model_Inputs.txt' file.
 */

//Dependency Librariers: Boost, OpenBLAS, Armadillo
//compile command:
//Debugging: g++ -g -std=c++11 -I /home/nirjhar/Libraries/boost_1_71_0/ PFR.cpp -o PFR -larmadillo
//Optimized: g++ -O2 -std=c++11 -I /home/nirjhar/Libraries/boost_1_71_0/ PFR.cpp -o PFR -larmadillo





#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
#include "func_bssanova_2.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <csignal>
#include "debug.hpp" // Contains functions used for GDB debuggig.
#include <chrono>
#include <iostream>
#include <fstream>
#include <utility>
#include <boost/numeric/odeint.hpp>
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/io.hpp"
#include <boost/phoenix/core.hpp>
#include "func_bssanova_2.hpp"
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include "plotter.hpp"
#include <random>
#include <ctime>




using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;

typedef boost::numeric::ublas::vector<double> uvec_t;




int main()

{



uvec_t par(8); //passes theta paramters(physical parameters) to model


//if theta parameters are in logscale or not
bool logscale = false;
//DiskOn variable define in plotter.hpp
int no_parameters; //total number of parameters(both beta and theta)  

  std::string time_stamp;
  int no_plot,upper_bound,lower_bound,break_step;

  //Grabs relevant model and plot details plot_config.txt located in source directory.
   
  std::string confile = "plot_config.txt";//this error number is used to give warning for problem reading from file
  std::ifstream confilestream;
  confilestream.open(confile.c_str(), std::ios::in);
  std::string line;
  while (!confilestream.eof())  
    {
      getline(confilestream, line, '\n');
      std::vector<std::string> str;
      if (line[0] != '#' && line[0] != ' ' && line.length() > 0) 
	{
	  boost::split(str, line, boost::is_any_of(" "));

	  if (str.size() == 9)  //weird
	    {
	      no_parameters = boost::lexical_cast<int>(str[0]);
	      time_stamp = str[1];
	      upper_bound = boost::lexical_cast<double>(str[2]);
	      lower_bound = boost::lexical_cast<int>(str[3]);
	      no_plot  = boost::lexical_cast<int>(str[4]);
	      break_step = boost::lexical_cast<int>(str[5]);
	      logscale = boost::lexical_cast<int>(str[6]);
	      DiskOn = boost::lexical_cast<int>(str[7]);
	    }  
	       
	}
    }
  



  rows = no_plot; //'rows'will specify no of lines in output text files



  

  boost::numeric::ublas::matrix<double> parameters(100100,8);
   
  if (DiskOn)
    parameters.resize(100100,no_parameters); //Change this 



  double N_CH4,N_H2O;


  //Choosing the parameter distribution file specified in plot_config.txt

  confile = "BayesParameterResults_";//"BayesParameterResults_'time_stamp_here'.txt";BayesParameterResults_1516899681
  confile.append(time_stamp);
  confile.append(".txt");      
                 
  std::ifstream constream;
  constream.open(confile.c_str(), std::ios::in);

  if (!constream.is_open())
    {std::cout << "Error opening parameter files" << std::endl;  //Error specifying failure to open parameter containing file
      throw 20;}



  int no_parameter_sets_loaded = 0; // Stores no of parameter sets loaded from text file.

  while (!constream.eof()) {
    getline(constream, line, '\n');
    std::vector<std::string> str;
    if (line[0] != '#' && line[0] != ' ' && line.length() > 0) {
      boost::split(str, line, boost::is_any_of(" "));
           
	   
      if (str.size() == no_parameters + 1) {
                  
		  
	for (int j = 0; j < no_parameters; j++)
	  {parameters(no_parameter_sets_loaded,j) = boost::lexical_cast<double>(str[j]);}         
      }  
                
	     
    }

    no_parameter_sets_loaded++;}

  std::cout << "Number of Parameter Sets Loaded = " << no_parameter_sets_loaded << std::endl; // Declares no of parameters sets loaded from text file.
   
  parameters.resize(no_parameter_sets_loaded,parameters.size2()); // Matrix 'parameters' resized accordingly

  boost::numeric::ublas::vector<double> betas(4);//arbitrary size

  std::default_random_engine RNG; //Setting up random number generator.
    
  //Read config and Data files

  //Reading Config Files
  confile = "config_Model_Inputs.txt";
  std::ifstream constream1; 
  constream1.open(confile.c_str(), std::ios::in);

  if (!constream1.is_open())
    std::cout << "Error Reading Config File" << std::endl; //Error specifying failure to open config file

  int j = 0;
  while (!constream1.eof()) {
    getline(constream1, line, '\n');
    std::vector<std::string> str;
    if (line[0] != '#' && line[0] != ' ' && line.length() > 0) {
      boost::split(str, line, boost::is_any_of(" "));
      if (j == 0) {
	if (str.size() == 2) {
	  no_steady_samples = boost::lexical_cast<int>(str[0]);
	  solver_k2 = boost::lexical_cast<double>(str[1]); // Reaction II assumed infinitely fast hence large k2 value
                   
	}                 
      }
    }
  }


  //Reading Data File
    
  input.resize(no_steady_samples,2);
  confile = "SMR_800_0.005.txt"; // Remove no of entries at top
    
  std::ifstream constream2; 
  constream2.open(confile.c_str(), std::ios::in);
  j = 0;
    
  while (!constream2.eof()) {
    getline(constream2, line, '\n');
    std::vector<std::string> str;
    if (line[0] != '#' && line[0] != ' ' && line.length() > 0) {
      boost::split(str, line, boost::is_any_of(" "));
            
	   
      input(j,0) = boost::lexical_cast<double>(str[0]); input(j,1) = boost::lexical_cast<double>(str[1]); // used to grab temp and CH4 inflow
		
      j++;
	    
    }
  }

  //Globally defined matrices to store PFR steady state values - resized according to configuration
  successful_runs_CH4.resize(rows,no_steady_samples);successful_runs_H2O.resize(rows,no_steady_samples); successful_runs_CO2.resize(rows,no_steady_samples);successful_runs_H2.resize(rows,no_steady_samples);successful_runs_CO.resize(rows,no_steady_samples);

  int index = 0;
  for ( int q = 0; q < no_plot; q++)                             
    {

      //Random number generation for sampling parametrs from posterior distribution

      //std::uniform_real_distribution<> dis(0.0, 1.0);
      //std::uniform_int_distribution<int> distribution(1,100);  
      
      //std::cout << dis(RNG) << std::endl;
      //double number = dis(RNG);
      //int number2 = floor(number*(upper_bound - lower_bound));
      srand((unsigned)time(0)+q);
      int rand_number = (rand()%1000);
      index =  lower_bound + rand_number;
      //index =  lower_bound + number2;


      //make the par and beta vectors to be passed into the plotter function
      for ( int i = 0; i < 8; i++)
	{par[i] = parameters(index,i);}	  
              
      //store beta parameters values if model contains discrepancy terms 
      if (DiskOn)
	{ 
	  betas.resize(no_parameters - 8);
	  for (int i = 0; i < (no_parameters - 8); i++)
	    betas[i] = parameters(index,8 + i);
	}



      std::cout << "Index is " << index << std::endl;  //Display parameter index during run time


      //apply conversion is parameters are in logscale
      if (logscale)
	{
	  for ( int i = 0; i < 8; i++)
	    {par[i] = pow(10,par[i]);}	  
	}

      // function defined in 'plotter.hpp' to solve ODE system and PFR for one parameter se  and store output values in steady_state_all matrix. steady_state_all will incrementally fill up successful_runs_{species_name} matrices.
      plotter(par,0.00037 ,0.000123,8.314,0.00001,800,10,betas);      
    }

  //Resize output matrix before model runs 
  successful_runs_CH4.resize(success_index,no_steady_samples);successful_runs_H2O.resize(success_index,no_steady_samples);
  successful_runs_CO2.resize(success_index,no_steady_samples);successful_runs_H2.resize(success_index,no_steady_samples);successful_runs_CO.resize(success_index,no_steady_samples); 


  //Write results to text files named after each individual chemical species in reactor
  FILE * temp1 = fopen("CH4.txt", "w");
 
  for (int a = 0; a < rows; a++) //was successful_runs_index
    {
      fprintf(temp1,"\n");
      for (int i = 0; i < successful_runs_CH4.size2(); i++)
	{
	  fprintf(temp1, " %lf", successful_runs_CH4(a,i)); //Write the data to a temporary file
	}
    }
 
  fclose(temp1);



  FILE * temp2 = fopen("H2O.txt", "w");
 
  for (int a = 0; a < rows; a++)
    {
      fprintf(temp2,"\n");
      for (int i = 0; i < successful_runs_H2O.size2(); i++)
	{
	  fprintf(temp2, " %lf", successful_runs_H2O(a,i)); //Write the data to a temporary file
	}
    }
 
  fclose(temp2);


  FILE * temp3 = fopen("CO2.txt", "w");
 
  for (int a = 0; a < rows; a++)
    {
      fprintf(temp3,"\n");
      for (int i = 0; i < successful_runs_CO2.size2(); i++)
	{
	  fprintf(temp3, " %lf", successful_runs_CO2(a,i)); //Write the data to a temporary file
	}
    }
 
  fclose(temp3);



  FILE * temp4 = fopen("H2.txt", "w");
 
  for (int a = 0; a <  rows; a++)
    {
      fprintf(temp4,"\n");
      for (int i = 0; i < successful_runs_H2.size2(); i++)
	{
	  fprintf(temp4, " %lf", successful_runs_H2(a,i)); //Write the data to a temporary file
	}
    }
 
  fclose(temp4);


  FILE * temp5 = fopen("CO.txt", "w");
 
  for (int a = 0; a <  rows; a++)
    {
      fprintf(temp5,"\n");
      for (int i = 0; i < successful_runs_CO.size2(); i++)
	{
	  fprintf(temp5, " %lf", successful_runs_CO(a,i)); //Write the data to a temporary file
	}
    }
 
  fclose(temp5);
}

