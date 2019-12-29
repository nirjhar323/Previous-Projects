/* This code uses FEM to upscale SMR Reduced Order Lab Scale Model to an Industrial Plug Flow Reactor. The model being used is defined by the text files 'plot_config.txt', 
 */
//compile command: g++ -g -std=c++11 -I /home/nirjhar/Libraries/boost_1_71_0/ test.cpp -o test -larmadillo



//#include "funeval_base.hpp"
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/vector.hpp>
 #include "boost/algorithm/string.hpp"
 #include "boost/lexical_cast.hpp"
//#include "psopoint_ublas_2.hpp"
 #include "func_bssanova_2.hpp"
#include <iostream>
 #include <fstream>
#include <string>
#include <sstream>
//#include "model.hpp"
#include <csignal>
 #include "debug.hpp"

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




using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;

typedef boost::numeric::ublas::vector<double> uvec_t;

uvec_t par(8); //passes theta paramters(physical parameters) to model


//if theta parameters are in logscale or not
bool logscale = false;
//DiskOn variable define in plotter.hpp
int no_parameters; //total number of parameters(both beta and theta)  


int main()

{

   std::string time_stamp;
   int no_plot,upper_bound,lower_bound,break_step;

   //Grabs relevant model and plot details plot_config.txt located in source directory.
   
   std::string confile = "plot_config.txt";//this error number is used to give warning for problem reading from file
   std::ifstream confilestream;
    confilestream.open(confile.c_str(), std::ios::in);
    //int i = 0;
    std::string line;
    while (!confilestream.eof())  
   {
        getline(confilestream, line, '\n');
        std::vector<std::string> str;
        if (line[0] != '#' && line[0] != ' ' && line.length() > 0) 
{
            boost::split(str, line, boost::is_any_of(" "));
//             if (i == 0) 
// {
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
		//	    }
	}
    }
  



  rows = no_plot; //Sets no of rows in output text files



  //successful_runs_CH4_transient.resize(rows,150); //transient run = 150s so 1s sampling used  

  boost::numeric::ublas::matrix<double> parameters(100100,8);
   
  if (DiskOn)
    parameters.resize(100100,no_parameters); //Change this 



  double N_CH4,N_H2O;

  confile = "BayesParameterResults_";//"BayesParameterResults_'time_stamp_here'.txt";BayesParameterResults_1516899681
  confile.append(time_stamp);
  confile.append(".txt");    
  
                 
    std::ifstream constream;
    constream.open(confile.c_str(), std::ios::in);

    if (!constream.is_open())
      {std::cout << "Error opening parameter files" << std::endl;
	throw 20;}


    int i = 0;
    //std::string line;
    while (!constream.eof()) {
        getline(constream, line, '\n');
        std::vector<std::string> str;
        if (line[0] != '#' && line[0] != ' ' && line.length() > 0) {
            boost::split(str, line, boost::is_any_of(" "));
           
	   
                if (str.size() == no_parameters + 1) {
                  
		  
		  for (int j = 0; j < no_parameters; j++)
		    {parameters(i,j) = boost::lexical_cast<double>(str[j]);}         
}  
                
	     
	}

	i++;}

   

    std::cout << "Number of Parameter Sets Loaded = " << i << std::endl;
   
     parameters.resize(i,parameters.size2());



     boost::numeric::ublas::vector<double> betas(4);//arbitrary size

   

	std::default_random_engine RNG; //I think this is always initialized to same value
    
	//Read config and Data files

	//Reading Config Files
    confile = "config_Model_Inputs.txt";
//error = 0;                    //this error number is used to give warning for problem reading from file
    




    std::ifstream constream1; 
    constream1.open(confile.c_str(), std::ios::in);

    if (!constream1.is_open())
      std::cout << "Error Reading Config File" << std::endl;


    int j = 0;
    //std::string line;
    while (!constream1.eof()) {
        getline(constream1, line, '\n');
        std::vector<std::string> str;
        if (line[0] != '#' && line[0] != ' ' && line.length() > 0) {
            boost::split(str, line, boost::is_any_of(" "));
            if (j == 0) {
                if (str.size() == 2) {
                    no_steady_samples = boost::lexical_cast<int>(str[0]);
                    solver_k2 = boost::lexical_cast<double>(str[1]);
                   
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
            
	   
		  input(j,0) = boost::lexical_cast<double>(str[0]); input(j,1) = boost::lexical_cast<double>(str[1]);
		
                j++;
	    
	}
    }




successful_runs_CH4.resize(rows,no_steady_samples);successful_runs_H2O.resize(rows,no_steady_samples);
 successful_runs_CO2.resize(rows,no_steady_samples);successful_runs_H2.resize(rows,no_steady_samples);successful_runs_CO.resize(rows,no_steady_samples);




   int index = 0;
   for ( int q = 0; q < no_plot; q++)                             
      {


        std::uniform_real_distribution<> dis(0.0, 1.0);
        std::uniform_int_distribution<int> distribution(1,100);  
      

	std::cout << dis(RNG) << std::endl;
        double number = dis(RNG);
       

        int number2 = floor(number*(upper_bound - lower_bound));
                
	  index =  lower_bound + number2;


	//make the par and beta vectors to be passed into the plotter function
	  for ( int i = 0; i < 8; i++)
	    {par[i] = parameters(index,i);}	  



              

              if (DiskOn)
	    { 
              betas.resize(no_parameters - 8);
	      for (int i = 0; i < (no_parameters - 8); i++)
		betas[i] = parameters(index,8 + i);
            }



	std::cout << "Index is " << index << std::endl;



        if (logscale)
	  {
	  for ( int i = 0; i < 8; i++)
	    {par[i] = pow(10,par[i]);}	  
	  }

	
         plotter(par,0.00037 ,0.000123,8.314,0.00001,800,10,betas,betas,betas,betas,betas);


      
 }








uvec_t xvals(rows);

for ( int a = 0; a < rows; a++)
      {
    	       xvals[a] = 0.01*a;
               
      }


//Resize successful runs matrix
 successful_runs_CH4.resize(success_index,no_steady_samples);successful_runs_H2O.resize(success_index,no_steady_samples);
 successful_runs_CO2.resize(success_index,no_steady_samples);successful_runs_H2.resize(success_index,no_steady_samples);successful_runs_CO.resize(success_index,no_steady_samples); //Resize this before model is run



 //Write results to text file data.txt
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


// FILE * temp6 = fopen("CH4_transient.txt", "w");
 
//  for (int a = 0; a < rows; a++)
//    {
//      fprintf(temp6,"\n");
//      for (int i = 0; i < successful_runs_CH4_transient.size2(); i++)
//     {
//       fprintf(temp6, " %lf", successful_runs_CH4_transient(a,i)); //Write the data to a temporary file
//     }
//    }
 
//     fclose(temp6);



}

