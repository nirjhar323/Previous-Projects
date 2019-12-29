#pragma once 
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

#include <stdlib.h>
#include <stdio.h>


#define NUM_COMMANDS 2 //no of commands passed to gnuplot





 typedef boost::numeric::ublas::vector< double > vector_type;
 typedef boost::numeric::ublas::matrix< double > matrix_type;

using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;
typedef boost::numeric::ublas::vector<double> uvec_t;


//steady state  and PFR variables

bool steady_state = false;
vector_type steady_state_output(5);//returns output  with 1st 4 elements partial pressures and last one outflow




int col_no = 0;
int success_index = 0;

//Set number of plot
int no_of_plots = 3;

//Will store result of successful runs only









double solver_kappa1,solver_kappa2,solver_kappa3,solver_k1,solver_k3,solver_N_CH4,solver_N_H2O,solver_N_CO2,solver_N_H2,R,V,T,t,N_TOTAL;



  
boost::numeric::ublas::vector<double> solver_betas_kappa1(4), solver_betas_kappa2(4), solver_betas_kappa3(4), solver_betas_k1(4), solver_betas_k3(4);
    
double disk_k1_res,disk_k3_res,disk_kappa1_res,disk_kappa2_res,disk_kappa3_res; //Stores value of discrepancy function



double sampling_interval = 0.01;//in seconds
double observing_time = 15;//in seconds

int rows = 1 + (observing_time/sampling_interval);

matrix_type yvals(rows,no_of_plots);//stores state variables for final result







boost::numeric::ublas::matrix<double> time_sample(rows,4);

func_bssanova  disk_k1,disk_k3,disk_kappa1,disk_kappa2,disk_kappa3;//func_bssanova  disk_k1,disk_k3,disk_kap1,disk_kap2,disk_kap3; //discrepanyc object ,there should be one for each of the parameters k1,k2, kappa1,2,3
boost::numeric::ublas::vector<double>  disk_k1_input_model(4),disk_k3_input_model(4),disk_kappa1_input_model(4),disk_kappa2_input_model(4),disk_kappa3_input_model(4); 










vector_type plotter(double passed_kappa1,double passed_kappa2,double passed_kappa3,double passed_k1,double passed_k3, double passed_R,double passed_V,double passed_T,double passed_t,uvec_t passed_betas_kappa1,uvec_t passed_betas_kappa2,uvec_t passed_betas_kappa3,uvec_t passed_betas_k1,uvec_t passed_betas_k3 )
{

  bool success = true;
  
typedef boost::numeric::ublas::vector <double> vector_type; 
  typedef boost::numeric::ublas::matrix <double> matrix_type;




//Setting Values
solver_kappa1 = passed_kappa1;
solver_kappa2 = passed_kappa2;
solver_kappa3 = passed_kappa3; 
solver_k1 = passed_k1;
solver_k3 = passed_k3;

// solver_N_CH4 = passed_N_CH4;
// solver_N_H2O = passed_N_H2O;
// solver_N_CO2 = passed_N_CH4;
// solver_N_H2 = passed_N_H2O;

 solver_N_CH4 = steady_state_output[4]*steady_state_output[0];
solver_N_H2O = steady_state_output[4]*steady_state_output[1];
solver_N_CO2 = steady_state_output[5]*steady_state_output[2];
solver_N_H2 =  steady_state_output[5]*steady_state_output[3];




R = passed_R;
V = passed_V;
T = passed_T;
t = passed_t;
solver_betas_k1 = passed_betas_k1;
solver_betas_k1 = passed_betas_k3;
solver_betas_kappa1 = passed_betas_kappa1;
solver_betas_kappa2 = passed_betas_kappa2;
solver_betas_kappa3 = passed_betas_kappa3;

N_TOTAL = steady_state_output[4];




 // set_values(passed_kappa1,passed_kappa2,passed_kappa3,passed_k1,passed_k3,passed_N_CH4,passed_N_H2O,passed_R,passed_V,passed_T,passed_t,passed_betas_k1,passed_betas_k3,passed_betas_kappa1,passed_betas_kappa2,passed_betas_kappa3,global_kappa1,global_kappa2,global_kappa3,global_k1,global_k3,model_N_CH4,model_N_H2O,R,V,T,t,global_betas_k1,global_betas_k3,global_betas_kappa1,global_betas_kappa2,global_betas_kappa3);



time_sample(0,0) = time_sample(0,1) = time_sample(0,2) = time_sample(0,3) = 0.25;  //Initial state value


//initialize discrepancy objects


      disk_k1.SplineConfig("config/config_Disc_smrEquilibrium.txt"); 
      disk_k3.SplineConfig("config/config_Disc_smrEquilibrium.txt");
      disk_kappa1.SplineConfig("config/config_Disc_smrEquilibrium.txt");
      disk_kappa2.SplineConfig("config/config_Disc_smrEquilibrium.txt");
      disk_kappa3.SplineConfig("config/config_Disc_smrEquilibrium.txt");













struct stiff_system

{
                                                                                        
  
double atm_pres = 101325;  //Atmospheric pressure in pascals

    void operator()( const vector_type &x , vector_type &dxdt , double t  )

    {
       
     
       //This needs fixing. Time intervals are all wrong


      //Condition for applying discrepancy
      //if (unniversal_counter%600 == 0)
      
      for ( int i = 1; i< time_sample.size1(); i++)   
	{
	  if (t > (i-1)*0.01 & t < i*0.01)
	    { time_sample(i,0) = x[0]; time_sample(i,1) = x[1]; time_sample(i,2) = x[2]; time_sample(i,3) = x[3];}

        }

    





            disk_k1_input_model[0] = x[0];
            disk_k1_input_model[1] = x[1];           
            disk_k1_input_model[2] = x[2]; //*time_sample(15,3))/(kappa2*time_sample(15,1));;//This is not CO
            disk_k1_input_model[3] = x[3];


            disk_k3_input_model[0] = x[0];
            disk_k3_input_model[1] = x[1];           
            disk_k3_input_model[2] = x[2]; //*time_sample(15,3))/(kappa2*time_sample(15,1));;//This is not CO
            disk_k3_input_model[3] = x[3];

            disk_kappa1_input_model[0] = x[0];
            disk_kappa1_input_model[1] = x[1];           
            disk_kappa1_input_model[2] = (x[2]*x[3])/(solver_kappa2*x[1]); //*time_sample(15,3))/(kappa2*time_sample(15,1));;//This is not CO
            disk_kappa1_input_model[3] = x[3];


            disk_kappa2_input_model[0] = (x[2]*x[3])/(solver_kappa2*x[1]);
            disk_kappa2_input_model[1] = x[1];           
            disk_kappa2_input_model[2] = x[2]; //*time_sample(15,3))/(kappa2*time_sample(15,1));;//This is not CO
            disk_kappa2_input_model[3] = x[3];


            disk_kappa3_input_model[0] = x[0];
            disk_kappa3_input_model[1] = x[1];           
            disk_kappa3_input_model[2] = x[2]; //*time_sample(15,3))/(kappa2*time_sample(15,1));;//This is not CO
            disk_kappa3_input_model[3] = x[3];





          disk_k1.SplineEval(disk_k1_input_model, solver_betas_k1, disk_k1_res);  
	  disk_k3.SplineEval(disk_k3_input_model, solver_betas_k3, disk_k3_res);
	  disk_kappa1.SplineEval(disk_kappa1_input_model, solver_betas_kappa1, disk_kappa1_res);
	  disk_kappa2.SplineEval(disk_kappa2_input_model, solver_betas_kappa2, disk_kappa2_res);
          disk_kappa3.SplineEval(disk_kappa3_input_model, solver_betas_kappa3, disk_kappa3_res);    
	  

   solver_k1 =  solver_k1*exp(disk_k1_res);
   solver_k3 =  solver_k3*exp(disk_k3_res);
   solver_kappa1 =  solver_kappa1*exp(disk_kappa1_res);
   solver_kappa2 =  solver_kappa2*exp(disk_kappa2_res);
   solver_kappa3 =  solver_kappa3*exp(disk_kappa3_res);
	

  double p_CO = (x[2]*x[3])/(x[1]*solver_kappa2);

  double r1 =  solver_k1*((x[0]*x[1]) - ((p_CO*x[3]*x[3]*x[3])/solver_kappa1));
  double r3 = (solver_k3*((x[0]*x[1]*x[1]) - ((x[2]*x[3]*x[3]*x[3]*x[3])/solver_kappa2)));
     
  
  dxdt[0] = (1/atm_pres)*((R*T)/V)*(solver_N_CH4 - (V*r1) - (V*r3) - x[0]*(N_TOTAL + (V*r1) + 2*(V*r3)));
  dxdt[1] = (1/atm_pres)*((R*T)/V)*(solver_N_H2O - (V*r1) - 2*(V*r3) - x[1]*(N_TOTAL + (V*r1) + 2*(V*r3)));  
  dxdt[2] = (1/atm_pres)*((R*T)/V)*(solver_N_CO2 + (V*r1) + (V*r3) - x[2]*(N_TOTAL + (V*r1) + 2*(V*r3)));
  dxdt[3] = (1/atm_pres)*((R*T)/V)*(solver_N_H2 + (V*r1) + 3*(V*r1) + 4*(V*r3) - x[3]*(N_TOTAL + (V*r1) + 2*(V*r3)));



      //Check for steady-state
      steady_state = false;
      if (t > 1)
        {
          int index = 0; 
	  index = (t - 1)/0.01; 
	  steady_state =  abs((x[0] - time_sample(index,0))) < 0.001 ? true:false;
          steady_state =  abs((x[1] - time_sample(index,1))) < 0.001 ? true:false;
          steady_state =  abs((x[2] - time_sample(index,2))) < 0.001 ? true:false;
          steady_state =  abs((x[3] - time_sample(index,3))) < 0.001 ? true:false;       



}
   
      if (steady_state)
	{std::cout << "ha" << std::endl;
	            steady_state_output[0] = x[0]; steady_state_output[1] = x[1];steady_state_output[2] = x[2];                     steady_state_output[3] = x[3]; steady_state_output[4] = (N_TOTAL + (V*r1) + 2*(V*r3));
                    N_TOTAL = steady_state_output[4];
        }




    }

};

struct stiff_system_jacobi
{



double atm_pres = 101325;  //Atmospheric pressure in pascals


    void operator()( const vector_type &x  , matrix_type &J , const double & /* t */ , vector_type &dfdt )
    {
 
     
  double p_CO = (x[2]*x[3])/(x[1]*solver_kappa2);

  double r1 = solver_k1*((x[0]*x[1]) - ((p_CO*x[3]*x[3]*x[3])/solver_kappa1));
  double r3 = solver_k3*((x[0]*x[1]*x[1]) - ((x[2]*x[3]*x[3]*x[3]*x[3])/solver_kappa2));  //maybe kappa3?? check




  double r1d1 = solver_k1*(x[1]);                            
  double r1d2 = solver_k1*(x[0]);                            
  double r1d3 = 0;
  double r1d4 = -solver_k1*((3*p_CO*x[3]*x[3])/solver_kappa1);
 
  double r3d1 = solver_k3*(x[1]*x[1]);
  double r3d2 = solver_k3*2*x[0]*x[1];                                   
  double r3d3 = -solver_k3*((x[3]*x[3]*x[3]*x[3])/solver_kappa2); 
  double r3d4 = -solver_k3-((4*x[2]*x[3]*x[3]*x[3])/solver_kappa2);


 
  J(0,0) = (1/atm_pres)*((R*T)/V)*(-(V*r1d1) - (V*r3d1) - (N_TOTAL + (V*r1d1*x[0] + V*r1) + 2*((V*r3d1*x[0])+(V*r3))));       
  J(0,1) = (1/atm_pres)*((R*T)/V)*(-(V*r1d2) - (V*r3d2) - x[0]*((V*r1d2) + 2*(V*r3d2)));                                   //
  J(0,2) = (1/atm_pres)*((R*T)/V)*(-(V*r1d3) - (V*r3d3) - x[0]*((V*r1d3) + 2*(V*r3d3)));
  J(0,3) = (1/atm_pres)*((R*T)/V)*(-(V*r1d4) - (V*r3d4) - x[0]*((V*r1d4) + 2*(V*r3d4)));
  
  J(1,0) = (1/atm_pres)*((R*T)/V)*(-(V*r1d1) - 2*(V*r3d1) - x[1]*((V*r1d1) + 2*(V*r3d1))); 
  J(1,1) = (1/atm_pres)*((R*T)/V)*(-(V*r1d2) - 2*(V*r3d2) - (N_TOTAL + (V*r1d2*x[1]+V*r1) + 2*(V*r3d2*x[1]+ V*r3)));
  J(1,2) = (1/atm_pres)*((R*T)/V)*(-(V*r1d3) - 2*(V*r3d3) - x[1]*((V*r1d3) + 2*(V*r3d3)));
  J(1,3) = (1/atm_pres)*((R*T)/V)*(-(V*r1d4) - 2*(V*r3d4) - x[1]*((V*r1d4) + 2*(V*r3d4)));
  
  J(2,0) = (1/atm_pres)*((R*T)/V)*((V*r1d1) + (V*r3d1) - x[2]*((V*r1d1) + 2*(V*r3d1)));
  J(2,1) = (1/atm_pres)*((R*T)/V)*((V*r1d2) + (V*r3d2) - x[2]*((V*r1d2) + 2*(V*r3d2)));
  J(2,2) = (1/atm_pres)*((R*T)/V)*((V*r1d3) + (V*r3d3) - (N_TOTAL + (V*r1d3*x[2] + V*r1) + 2*(V*r3d3*x[2] + V*r3)));
  J(2,3) = (1/atm_pres)*((R*T)/V)*((V*r1d4) + (V*r3d4) - x[2]*((V*r1d4) + 2*(V*r3d4)));
  
  J(3,0) = (1/atm_pres)*((R*T)/V)*((V*r1d1) + 3*(V*r1d1) + 4*(V*r3d1) - x[3]*((V*r1d1) + 2*(V*r3d1)));
  J(3,1) = (1/atm_pres)*((R*T)/V)*((V*r1d2) + 3*(V*r1d2) + 4*(V*r3d2) - x[3]*((V*r1d2) + 2*(V*r3d2)));
  J(3,2) = (1/atm_pres)*((R*T)/V)*((V*r1d3) + 3*(V*r1d3) + 4*(V*r3d3) - x[3]*((V*r1d3) + 2*(V*r3d3)));
  J(3,3) = (1/atm_pres)*((R*T)/V)*((V*r1d4) + 3*(V*r1d4) + 4*(V*r3d4) - (N_TOTAL + (V*r1d4*x[3]+ V*r1) + 2*(V*r3 + V*r3d4*x[3])));   
         
        //WHAT ARE THESE? // Like sundials
        dfdt[0] = 0.0;
        dfdt[1] = 0.0;
        dfdt[2] = 0.0;
        dfdt[3] = 0.0;
    }
};

  
 std::cout << solver_kappa1 << " "  << solver_kappa2 << " " << solver_kappa3<< " " << solver_k1 << " " << solver_k3<< " "  <<  solver_betas_kappa1[0] <<  " " <<  solver_betas_kappa1[1] << " " << solver_betas_kappa1[2] <<" " << solver_betas_kappa1[3] << " " << solver_N_CH4<<  " " << solver_N_H2O <<  " " <<  R  <<  " " <<  V  <<  " " <<  T << " " << t <<  " " << std::endl;




vector_type x( 4, 0.25);

           






  

try
  {
  size_t num_of_steps = integrate_const( make_dense_output< rosenbrock4< double > >( 1.0e-4 , 1.0e-4, 1.0e-4) , make_pair( stiff_system() , stiff_system_jacobi() ) ,x , 0.0 , observing_time , 0.1,cout << phoenix::arg_names::arg2 << " " << phoenix::arg_names::arg1[0] << " " << phoenix::arg_names::arg1[1] << " " << phoenix::arg_names::arg1[2] << " " << phoenix::arg_names::arg1[3] << "\n" );  
  }

 catch(...)
   { 
  std::cout << "Solver Failure"<< std::endl;
  
     
  // for (int i =0; i < time_sample.size1(); i++)
  //      {time_sample(i,0) = time_sample(i,1)= time_sample(i,2)= time_sample(i,3) = 0;} 
 

  success = false;

 }






 
 
//filling up xvals and yvals with time_sample values
     

    if (steady_state)    
  for ( int a = 0; a < time_sample.size1(); a++)
      {
	yvals(a,col_no) = time_sample(a,0);//store state values in yvals
      }
    
    else
      std::cout << "Still Unsteady" << std::endl;

    
    col_no++;








}









