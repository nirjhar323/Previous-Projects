#pragma once 
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
//#include "psopoint_ublas_2.hpp"
#include "func_bssanova_2.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
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

#include <stdlib.h>
#include <stdio.h>

//Steady state testing variables
int steady_state_test_counter = 1;

//type definitions used in ODEINT library
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;
typedef boost::numeric::ublas::vector<double> uvec_t;


//Model Update
double total_pressure = 0;
double original_kappa1, original_kappa2, original_kappa3, original_k1, original_k3; 
double r1_temp,r2_temp,r3_temp; 
double solver_kappa1,solver_kappa2,solver_kappa3,solver_k1,solver_k2,solver_k3,solver_N_CH4,solver_N_H2O,solver_N_CO2,solver_N_H2,solver_N_CO,R,V,T,t,N_TOTAL;

//New theta variables
double delH1,delH2,delH3,delS1,delS2,delS3,A1,A3;

//Index variables and Solver Error Handling Variables - to store successful runs only
int success_index = 0;
int col_no = 0;
bool solver_success = false;


//Discrepanyc(Disk) variables  
bool DiskOn = false; // Variable to control if discrepancy is used
boost::numeric::ublas::vector<double> solver_betas(4);
double disk_k1_res,disk_k3_res,disk_kappa1_res,disk_kappa2_res,disk_kappa3_res; //Stores value of discrepancy function
func_bssanova  disk_k1,disk_k3,disk_kappa1,disk_kappa2,disk_kappa3;//func_bssanova - discrepanyc object, one for each of the parameters k1,k2, kappa1,2,3
boost::numeric::ublas::vector<double>  disk_k1_input_model(5),disk_k3_input_model(5),disk_kappa1_input_model(5),disk_kappa2_input_model(5),disk_kappa3_input_model(5); 


//Observer and Result Processing Variables
double sampling_interval = 0.01;//in seconds
double observing_time = 10;//in seconds
int rows;// = 1 + (observing_time/sampling_interval); // No of rows in CH4,H20 . . . containers

//Globally defined matrices to store PFR steady state values 
matrix_type successful_runs_CH4(rows,1),successful_runs_H2O(rows,1),successful_runs_CO2(rows,1),successful_runs_H2(rows,1),successful_runs_CO(rows,1);



//Variables to assess transient/steady state
matrix_type successful_runs_CH4_transient(rows,1);
bool transient_data;
double transient_output;
boost::numeric::ublas::matrix<double> time_sample(15000,5);//Holds transient state values, assumes steady state in 15s max
bool steady_state_condition = false;
matrix_type  steady_state_all(6,5); //Stores the output values for all gases inside tank for all lab reactors - initialized arbitrarily - resized later on  
matrix_type steady_state_all_temp(1,5); //Temp container to hold values for 1 solver run
matrix_type input; //Stores the methane inflow values. F_H2O = F_CH4 * 3;
int no_steady_samples;




void plotter(uvec_t passed_par, double passed_N_CH4,double passed_N_H2O,double passed_R,double passed_V,double passed_T,double passed_t,uvec_t passed_betas)
{

  //Model Update
  delH1 = passed_par[0]; delH2 = passed_par[1];delH3 = passed_par[2];delS1 = passed_par[3];
  delS2 = passed_par[4];delS3 = passed_par[5]; A1 = passed_par[6]; A3 = passed_par[7];

  solver_success = true;  
  steady_state_all.resize(no_steady_samples,5);



  transient_data = false;
 
  //Re-initialize time_sample to zero
  for (int i = 0; i < time_sample.size1();i++)
    {time_sample(i,0) = time_sample(i,1) = time_sample(i,2) = time_sample(i,3) = time_sample(i,4) = 0;}


  //Input starting values of the states - Make this text grabbed
  for (int i = 0; i < time_sample.size2();i++)
    {time_sample(0,i) = 2;}


  total_pressure = 0;

  for (int i = 0; i < time_sample.size2();i++)
    {total_pressure += time_sample(0,i);} 

  R = passed_R;
  V = passed_V;
  T = passed_T;
  t = passed_t;
  //These are all the same vectors!
  solver_betas = passed_betas;

  //initialize discrepancy objects if model has discrepancy terms. 
  if (DiskOn)
    {
      disk_k1.SplineConfig("config/config_Disc_smrEquilibrium.txt"); 
      disk_k3.SplineConfig("config/config_Disc_smrEquilibrium.txt");
      disk_kappa1.SplineConfig("config/config_Disc_smrEquilibrium.txt");
      disk_kappa2.SplineConfig("config/config_Disc_smrEquilibrium.txt");
      disk_kappa3.SplineConfig("config/config_Disc_smrEquilibrium.txt");
    }


  //Structure used by ODEINT library
  struct stiff_system

  {
                                                                                        
    double atm_pres = 101325;  //Atmospheric pressure in pascals

    void operator()( const vector_type &x , vector_type &dxdt , double t  ) // vector 'x' stores transient values of the states at different times

    {
      int index = 0;

      //If any state exceed 10 bar(i.e there is solver error) displays 'Stop'
      for (int i = 0; i < x.size(); i++)
	{
	  if (x[i] > 10)
	    std::cout << "Stop" << std::endl;
        }

      
      if (!transient_data)
	{
	  steady_state_condition = false; //If True,tank is in steady state


	  //Checks for steady state condition after 1s (i.e t>1)
	  if (t > 1)
	    {
           
	      index = int((t - 1)/0.01); 
	      steady_state_condition =  abs((x[0] - time_sample(index,0))) < 0.001 ? true:false;
	      steady_state_condition =  abs((x[1] - time_sample(index,1))) < 0.001 ? true:false;
	      steady_state_condition =  abs((x[2] - time_sample(index,2))) < 0.001 ? true:false;
	      steady_state_condition =  abs((x[3] - time_sample(index,3))) < 0.001 ? true:false;       
	      steady_state_condition =  abs((x[4] - time_sample(index,4))) < 0.001 ? true:false;
	    }

          //New code for steady state testing
          if(steady_state_condition && steady_state_test_counter == 1)
	    {
	      time_sample.resize(t/0.01,time_sample.size2());
	      FILE * temp6 = fopen("steady_state_testing.txt", "w");
 
              for (int a = 0; a <  time_sample.size1(); a++)
		{
		  fprintf(temp6,"\n");
		  for (int i = 0; i < time_sample.size2(); i++)
		    {
		      fprintf(temp6, " %lf", time_sample(a,i)); //Write the data to a temporary file
		    }
		}
 
	      fclose(temp6);
	      
	      steady_state_test_counter++;

	    }



	  
	  //New Code for steady state testing
	  if(steady_state_condition && steady_state_test_counter > 1) 
	    {
	  
	      for (int i = 0; i < steady_state_all.size2(); i++)
		{
		  //if ( i == 0)
		  steady_state_all_temp(0,i) = x[i]; //steady_state_all_temp = stores the the steady state partial pressure values for the specific tank for all species          
		  
		}

	      //Calculate outflow to adjancent reactor
	      double total_outflow = (N_TOTAL + 2*(V*r1_temp) + 2*(V*r3_temp));
	      solver_N_CH4 = x[0]*total_outflow;
	      solver_N_CH4 = (x[0]/total_pressure)*total_outflow;
	      solver_N_H2O =(x[1]/total_pressure)*total_outflow;
	      solver_N_CO2 = (x[2]/total_pressure)*total_outflow;
	      solver_N_H2 = (x[3]/total_pressure)*total_outflow;
	      solver_N_CO = (x[4]/total_pressure)*total_outflow;



	      throw 20; // if steady state condition reached move to next finite element(i.e lab scale reactor) 
        

	    }

	  
	}


    
      //Stores transient state values in matrix time_sample    
      for ( int i = 1; i < time_sample.size1(); i++)//change 1 to 0

        {

  	  if (t > (i-1)*0.01 && t < i*0.01)
	    { time_sample(i,0) = x[0]; time_sample(i,1) = x[1]; time_sample(i,2) = x[2];time_sample(i,3) = x[3];time_sample(i,4) = x[4]; }
    
        }


      if (DiskOn)
	{

	  //Specifying inputs to discrepancy function
	  disk_k1_input_model[0] = x[0];
	  disk_k1_input_model[1] = x[1];           
	  disk_k1_input_model[2] = x[4]; 
	  disk_k1_input_model[3] = x[3];
	  disk_k1_input_model[4] = 1/T;

	  disk_k3_input_model[0] = x[0];
	  disk_k3_input_model[1] = x[1];           
	  disk_k3_input_model[2] = x[2];
	  disk_k3_input_model[3] = x[3];
	  disk_k3_input_model[4] = 1/T;          
  
	  disk_kappa1_input_model[0] = x[0];
	  disk_kappa1_input_model[1] = x[1];           
	  disk_kappa1_input_model[2] = x[4];
	  disk_kappa1_input_model[3] = x[3];
	  disk_kappa1_input_model[4] = 1/T;

	  disk_kappa2_input_model[0] = x[4];;
	  disk_kappa2_input_model[1] = x[1];           
	  disk_kappa2_input_model[2] = x[2]; 
	  disk_kappa2_input_model[3] = x[3];
	  disk_kappa2_input_model[4] = 1/T;

	  disk_kappa3_input_model[0] = x[0];
	  disk_kappa3_input_model[1] = x[1];           
	  disk_kappa3_input_model[2] = x[2];
	  disk_kappa3_input_model[3] = x[3];
	  disk_kappa3_input_model[4] = 1/T;


	    
	  //SplienEval is a func_bssanova type object that updates the parmaeters with discrepancy
          disk_k1.SplineEval(disk_k1_input_model, solver_betas, disk_k1_res);  
	  disk_k3.SplineEval(disk_k3_input_model, solver_betas, disk_k3_res);
	  disk_kappa1.SplineEval(disk_kappa1_input_model, solver_betas, disk_kappa1_res);
	  disk_kappa2.SplineEval(disk_kappa2_input_model, solver_betas, disk_kappa2_res);
          disk_kappa3.SplineEval(disk_kappa3_input_model, solver_betas, disk_kappa3_res);    
	  

	 
	  // Update physical parameters with discrepancy
	  solver_k1 =  original_k1*exp(disk_k1_res);
	  solver_k3 =  original_k3*exp(disk_k3_res);
	  solver_kappa1 =  original_kappa1*exp(disk_kappa1_res);
	  solver_kappa2 =  original_kappa2*exp(disk_kappa2_res);
	  solver_kappa3 =  original_kappa3*exp(disk_kappa3_res);
	}

      else 
	{
	  solver_k1 =  original_k1;
	  solver_k3 =  original_k3;
	  solver_kappa1 = original_kappa1;
	  solver_kappa2 = original_kappa2;
	  solver_kappa3 = original_kappa3;
	}
	
      //Calculate reaction rates 
      double r1 =  r1_temp = solver_k1*((x[0]*x[1]) - ((x[4]*x[3]*x[3]*x[3])/solver_kappa1));
      double r2 =  r2_temp = solver_k2*((x[4]*x[1]) - ((x[2]*x[3])/solver_kappa2));
      double r3 =  r3_temp = (solver_k3*((x[0]*x[1]*x[1]) - ((x[2]*x[3]*x[3]*x[3]*x[3])/solver_kappa3)));
 
     

      //ODEs defining system
      dxdt[0] = (1/atm_pres)*((R*T)/V)*(solver_N_CH4 - (V*r1) - (V*r3) - (x[0]/(total_pressure))*(N_TOTAL + 2*(V*r1) + 2*(V*r3)));
      dxdt[1] = (1/atm_pres)*((R*T)/V)*(solver_N_H2O - (V*r1) - V*r2 - 2*(V*r3) - (x[1]/(total_pressure))*(N_TOTAL + 2*(V*r1) + 2*(V*r3)));  
      dxdt[2] = (1/atm_pres)*((R*T)/V)*(solver_N_CO + (V*r2) + (V*r3) - (x[2]/(total_pressure))*(N_TOTAL + 2*(V*r1) + 2*(V*r3)));
      dxdt[3] = (1/atm_pres)*((R*T)/V)*(solver_N_H2 + 3*(V*r1) + V*r2 + 4*(V*r3) - (x[3]/(total_pressure))*(N_TOTAL + 2*(V*r1) + 2*(V*r3)));
      dxdt[4] = (1/atm_pres)*((R*T)/V)*(solver_N_CO + V*r1 - V*r2 - (x[4]/(total_pressure))*(N_TOTAL + 2*(V*r1) + 2*(V*r3)));


    }

  };


  //Structure used by ODEINT library
 
  struct stiff_system_jacobi
  {

    double atm_pres = 101325;  //Atmospheric pressure in Pa


    void operator()( const vector_type &x  , matrix_type &J , const double & /* t */ , vector_type &dfdt )
    {

      //Calculate reaction rates    
      double r1 = solver_k1*((x[0]*x[1]) - ((x[4]*x[3]*x[3]*x[3])/solver_kappa1));
      double r2 = solver_k2*((x[4]*x[1]) - ((x[2]*x[3])/solver_kappa2));
      double r3 = solver_k3*((x[0]*x[1]*x[1]) - ((x[2]*x[3]*x[3]*x[3]*x[3])/solver_kappa3));  


      // Partial derivatives of reaction rates
      double r1d0 = solver_k1*(x[1]);                   
      double r1d1 = solver_k1*(x[0]);                             
      double r1d2 = 0;
      double r1d3 = -solver_k1*((3*x[4]*x[3]*x[3])/solver_kappa1);
      double r1d4 = -solver_k1*((x[3]*x[3]*x[3])/solver_kappa1);

      double r2d0 = 0;                            
      double r2d1 = solver_k2*(x[4]);                            
      double r2d2 = -solver_k2*(x[3]/solver_kappa2);
      double r2d3 = -solver_k2*(x[2]/solver_kappa2);
      double r2d4 = solver_k2*x[1];
 
      double r3d0 = solver_k3*(x[1]*x[1]);
      double r3d1 = solver_k3*2*x[0]*x[1];                                   
      double r3d2 = -solver_k3*((x[3]*x[3]*x[3]*x[3])/solver_kappa3); 
      double r3d3 = -solver_k3*((4*x[2]*x[3]*x[3]*x[3])/solver_kappa3);
      double r3d4 = 0;


 
      J(0,0) = (1/atm_pres)*((R*T)/V)*(-(V*r1d0) - (V*r3d0) - (1/total_pressure)*(N_TOTAL + 2*(V*r1d0*x[0] + V*r1) + 2*((V*r3d0*x[0])+(V*r3))));       
      J(0,1) = (1/atm_pres)*((R*T)/V)*(-(V*r1d1) - (V*r3d1) - (x[0]/(total_pressure))*(2*(V*r1d1) + 2*(V*r3d1)));                                   
      J(0,2) = (1/atm_pres)*((R*T)/V)*(-(V*r1d2) - (V*r3d2) - (x[0]/(total_pressure))*(2*(V*r1d2) + 2*(V*r3d2)));
      J(0,3) = (1/atm_pres)*((R*T)/V)*(-(V*r1d3) - (V*r3d3) - (x[0]/(total_pressure))*(2*(V*r1d3) + 2*(V*r3d3)));
      J(0,4) = (1/atm_pres)*((R*T)/V)*(-(V*r1d4) - (V*r3d4) - (x[0]/(total_pressure))*(2*(V*r1d4) + 2*(V*r3d4)));
  

      J(1,0) = (1/atm_pres)*((R*T)/V)*(-2*(V*r1d0) - V*r2d0 - 2*(V*r3d0) - (x[1]/(total_pressure))*(2*(V*r1d0) + 2*(V*r3d0))); 
      J(1,1) = (1/atm_pres)*((R*T)/V)*(-2*(V*r1d1) - V*r2d1 - 2*(V*r3d1) - (1/(total_pressure))*(N_TOTAL + 2*(V*r1d1*x[1]+V*r1) + 2*(V*r3d1*x[1]+ V*r3)));
      J(1,2) = (1/atm_pres)*((R*T)/V)*(-2*(V*r1d2) - V*r2d2 - 2*(V*r3d3) - (x[1]/(total_pressure))*(2*(V*r1d3) + 2*(V*r3d3)));
      J(1,3) = (1/atm_pres)*((R*T)/V)*(-2*(V*r1d3) - V*r2d3 - 2*(V*r3d3) - (x[1]/(total_pressure))*(2*(V*r1d3) + 2*(V*r3d3)));
      J(1,4) = (1/atm_pres)*((R*T)/V)*(-2*(V*r1d4) - V*r2d4 - 2*(V*r3d4) - (x[1]/(total_pressure))*(2*(V*r1d4) + 2*(V*r3d4)));
  
      J(2,0) = (1/atm_pres)*((R*T)/V)*((V*r2d0) + (V*r3d0) - (x[2]/(total_pressure))*(2*(V*r1d0) + 2*(V*r3d0)));
      J(2,1) = (1/atm_pres)*((R*T)/V)*((V*r2d1) + (V*r3d1) - (x[2]/(total_pressure))*(2*(V*r1d1) + 2*(V*r3d1)));
      J(2,2) = (1/atm_pres)*((R*T)/V)*((V*r2d2) + (V*r3d2) - (1/(total_pressure))*(N_TOTAL + 2*(V*r1d2*x[2] + V*r1) + 2*(V*r3d2*x[2] + V*r3)));
      J(2,3) = (1/atm_pres)*((R*T)/V)*((V*r2d3) + (V*r3d3) - (x[2]/(total_pressure))*(2*(V*r1d3) + 2*(V*r3d3)));
      J(2,4) = (1/atm_pres)*((R*T)/V)*((V*r2d4) + (V*r3d4) - (x[2]/(total_pressure))*(2*(V*r1d4) + 2*(V*r3d4)));
  
      J(3,0) = (1/atm_pres)*((R*T)/V)*( 3*(V*r1d0) +  V*r2d0 + 4*(V*r3d0) - (x[3]/(total_pressure))*(2*(V*r1d0) + 2*(V*r3d0)));
      J(3,1) = (1/atm_pres)*((R*T)/V)*( 3*(V*r1d1) +  V*r2d1 + 4*(V*r3d1) - (x[3]/(total_pressure))*(2*(V*r1d1) + 2*(V*r3d1)));
      J(3,2) = (1/atm_pres)*((R*T)/V)*( 3*(V*r1d2) +  V*r2d2 + 4*(V*r3d2) - (x[3]/(total_pressure))*(2*(V*r1d2) + 2*(V*r3d2)));
      J(3,3) = (1/atm_pres)*((R*T)/V)*( 3*(V*r1d3) +  V*r2d3 + 4*(V*r3d3) - ((1/(total_pressure)))*(N_TOTAL + 2*(V*r1d3*x[3]+ V*r1) + 2*(V*r3 + V*r3d3*x[3])));
      J(3,4) = (1/atm_pres)*((R*T)/V)*( 3*(V*r1d4) +  V*r2d4 + 4*(V*r3d4) - (x[3]/(total_pressure))*(2*(V*r1d4) + 2*(V*r3d4)));

      J(4,0) = (1/atm_pres)*((R*T)/V)*((V*r1d0) - (V*r2d0)  - (x[4]/(total_pressure))*(2*(V*r1d0) + 2*(V*r3d0)));
      J(4,1) = (1/atm_pres)*((R*T)/V)*((V*r1d1) - (V*r2d1)  - (x[4]/(total_pressure))*(2*(V*r1d1) + 2*(V*r3d1)));
      J(4,2) = (1/atm_pres)*((R*T)/V)*((V*r1d2) - (V*r2d2)  - (x[4]/(total_pressure))*(2*(V*r1d2) + 2*(V*r3d2)));
      J(4,3) = (1/atm_pres)*((R*T)/V)*((V*r1d3) - (V*r2d3)  - (x[4]/(total_pressure))*(2*(V*r1d3) + 2*(V*r3d3)));
      J(4,4) = (1/atm_pres)*((R*T)/V)*((V*r1d4) -  V*r2d4 -  4*(V*r3d4) - (1/(total_pressure))*(N_TOTAL + 2*(((V*r1d4)*x[4]) + V*r1)+ 2*(V*r3d4*x[4] + V*r3))); 
   
         
     
      dfdt[0] = 0.0;
      dfdt[1] = 0.0;
      dfdt[2] = 0.0;
      dfdt[3] = 0.0;
      dfdt[4] = 0.0;


    }
  };



  

  for ( int k = 0; k < steady_state_all.size1(); k++)

    {

      //setting inflow values
      T = input(k,1);
      if ( k == 0)                  //i.e first reactor - k = reactor no?
	{
	  solver_N_CH4 = input(k,0);
	  solver_N_H2O = solver_N_CH4*3;
	  solver_N_H2 = 0;
	}

      N_TOTAL = solver_N_CH4 + solver_N_H2O + solver_N_H2 + solver_N_CO2 + solver_N_CO; 

      original_k1 =  A1*exp(-delH1/(R*T));  
      original_k3 =  A3*exp(-delH3/(R*T));  
      original_kappa1 = exp(delS1/R)*exp(-delH1/(R*T)); 
      original_kappa2 = exp(delS2/R)*exp(-delH2/(R*T));  
      original_kappa3 = exp(delS3/R)*exp(-delH3/(R*T));  


      //Display parameter and configuration values for each model run
      std::cout << original_kappa1 << " "  << original_kappa2 << " " << original_kappa3 << " " << original_k1 << " " << original_k3<< " "  <<  solver_betas[0] <<  " " <<  solver_betas[1] << " " << solver_betas[2] <<" " << solver_betas[3] << " " << "CH4=" << solver_N_CH4 <<  " " << "H2O="<< solver_N_H2O <<  " " <<  R  <<  " " <<  V  <<  " " <<  T << " " << t <<  " " << std::endl;
 
      vector_type x( 5, 2);

      /* For debugging 
	 if ( k > 22)
	 std::cout << "Stop" << std::endl;
      */
 
      try
	{
	  size_t num_of_steps = integrate_const( make_dense_output< rosenbrock4< double > >( 1.0e-3 , 1.0e-3,1.0e-3) , make_pair( stiff_system() , stiff_system_jacobi() ) ,x , 0.0 , observing_time , 1.0 ,cout << phoenix::arg_names::arg2 << " " << phoenix::arg_names::arg1[0] << " " << phoenix::arg_names::arg1[1] << " " << phoenix::arg_names::arg1[2] << " " << phoenix::arg_names::arg1[3] << " " << phoenix::arg_names::arg1[4]  << "\n" );    

	}

      catch(...)
	{
 

	  //display output matrix
	  for (int a = 0; a < steady_state_all.size2(); a++)
	    {steady_state_all(k,a) = steady_state_all_temp(0,a);}

	  get_matrix(steady_state_all);

	}
  

    }

  if (steady_state_condition) 
    {
      for(int q = 0; q < steady_state_all.size1(); q++)
	{
	  successful_runs_CH4(success_index,q) = steady_state_all(q,0);
	  successful_runs_H2O(success_index,q) = steady_state_all(q,1);
	  successful_runs_CO2(success_index,q) = steady_state_all(q,2);
	  successful_runs_H2(success_index,q) = steady_state_all(q,3);
	  successful_runs_CO(success_index,q) = steady_state_all(q,4);

        }
      
      // success_index++;        incremented after running model for transient data
    }
     


  //Run model for getting transient data
  if ( no_steady_samples == input.size1())
    {success_index++;
      return; }// Terminate function if no transient data
}









