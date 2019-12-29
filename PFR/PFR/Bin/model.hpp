
/*
Takes in the theta parameters(i.e k1,k3,kap1,2,3) and beta parameters and solves the ODE system, applying discrepancy at each 0.1s interval(i.e theta parameter values are updated each 0.01 second)
*/




#pragma once
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


typedef boost::numeric::ublas::vector<double> uvec_t;


//Global Variables
double  global_kappa1,global_kappa2,global_kappa3,global_k1,global_k3,N_CH4,N_H2O,R,V,T,t; 
uvec_t global_betas_k1,global_betas_k3,global_betas_kap1,global_betas_kap2,global_betas_kap3;

double  N_TOTAL;
time_t start_time;
int unniversal_counter;

/*
//Discrepancy Variables
double k1_mod,k3_mod,kappa1_mod,kappa2_mod,kappa3_mod;

//func_bssanova  disk_k1,disk_k3,disk_kap1,disk_kap2,disk_kap3; //discrepanyc object ,there should be one for each of the parameters k1,k2, kappa1,2,3
double disk_k1_res,disk_k3_res,disk_kappa1_res,disk_kappa2_res,disk_kappa3_res; //Stores value of discrepancy function

func_bssanova  disk_k1,disk_k3,disk_kap1,disk_kap2,disk_kap3;

boost::numeric::ublas::vector<double>  disk_k1_input_model(4),disk_k3_input_model(4),disk_kap1_input_model(4),disk_kap2_input_model(4),disk_kap3_input_model(4);         
bool diskconsuccess = false;
*/

      // disk_k1.SplineConfig("config/config_Disc_smrEquilibrium.txt"); 
      // disk_k3.SplineConfig("config/config_Disc_smrEquilibrium.txt");
      // disk_kap1.SplineConfig("config/config_Disc_smrEquilibrium.txt");
      // disk_kap2.SplineConfig("config/config_Disc_smrEquilibrium.txt");
      // disk_kap3.SplineConfig("config/config_Disc_smrEquilibrium.txt");

 

boost::numeric::ublas::matrix<double> time_sample(16,4);





void set_values(double passed_kappa1,double passed_kappa2,double passed_kappa3,double passed_k1,double passed_k3,double passed_N_CH4,double passed_N_H2O,double passed_R,double passed_V,double passed_T,double passed_t,uvec_t passed_betas_k1,uvec_t passed_betas_k3,uvec_t passed_betas_kap1,uvec_t passed_betas_kap2,uvec_t passed_betas_kap3,double &ref_kaapa1,double &ref_kaapa2,double &ref_kaapa3,double &ref_k1,double &ref_k3,double &ref_N_CH4,double &ref_N_H2O,double &ref_R,double &ref_V,double &ref_T,double &ref_t,uvec_t &ref_betas_k1,uvec_t &ref_betas_k3,uvec_t &ref_betas_kap1,uvec_t &ref_betas_kap2,uvec_t &ref_betas_kap3)
{
ref_kaapa1 = passed_kappa1;
ref_kaapa2 = passed_kappa2;
ref_kaapa3 = passed_kappa3; 
ref_k1 = passed_k1;
ref_k3 = passed_k3;
ref_N_CH4 = passed_N_CH4;
ref_N_H2O = passed_N_H2O;
ref_R = passed_R;
ref_V = passed_V;
ref_T = passed_T;
ref_t = passed_t;
ref_betas_k1 = passed_betas_k1;
ref_betas_k1 = passed_betas_k1;
ref_betas_k1 = passed_betas_k1;
ref_betas_kap2 = passed_betas_kap2;
ref_betas_kap3 = passed_betas_kap3;

N_TOTAL = N_CH4 + N_H2O;
}







boost::numeric::ublas::matrix <double>  model(double passed_kappa1,double passed_kappa2,double passed_kappa3,double passed_k1,double passed_k3,double passed_N_CH4,double passed_N_H2O,double passed_R,double passed_V,double passed_T,double passed_t,uvec_t passed_betas_k1,uvec_t passed_betas_k3,uvec_t passed_betas_kap1,uvec_t passed_betas_kap2,uvec_t passed_betas_kap3 )
{


 


  
using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;




typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;


 set_values(passed_kappa1,passed_kappa2,passed_kappa3,passed_k1,passed_k3,passed_N_CH4,passed_N_H2O,passed_R,passed_V,passed_T,passed_t,passed_betas_k1,passed_betas_k3,passed_betas_kap1,passed_betas_kap2,passed_betas_kap3,global_kappa1,global_kappa2,global_kappa3,global_k1,global_k3,N_CH4,N_H2O,R,V,T,t,global_betas_k1,global_betas_k3,global_betas_kap1,global_betas_kap2,global_betas_kap3);

time_sample(0,0) = time_sample(0,1) = time_sample(0,2) = time_sample(0,3) = 0.25;  //Initial state value


//initialize discrepancy objects


      // disk_k1.SplineConfig("config/config_Disc_smrEquilibrium.txt"); 
      // disk_k3.SplineConfig("config/config_Disc_smrEquilibrium.txt");
      // disk_kap1.SplineConfig("config/config_Disc_smrEquilibrium.txt");
      // disk_kap2.SplineConfig("config/config_Disc_smrEquilibrium.txt");
      // disk_kap3.SplineConfig("config/config_Disc_smrEquilibrium.txt");













struct stiff_system

{
                                                                                        
 
double atm_pres = 101325;  //Atmospheric pressure in pascals

    void operator()( const vector_type &x , vector_type &dxdt , double t  )

    {

       unniversal_counter++;

       //This needs fixing. Time intervals are all wrong
      if ( t > 0.00 && t < 0.010001 ){  
	time_sample(1,0) = x[0];time_sample(1,1) = x[1];time_sample(1,2) = x[2];time_sample(1,3) = x[3];}
      else if ( t > 0.010001 && t < 0.020001){
  	time_sample(2,0) = x[0];time_sample(2,1) = x[1];time_sample(2,2) = x[2];time_sample(2,3) = x[3];}
      else if (t > 0.020001 && t < 0.030001){
  	time_sample(3,0) = x[0];time_sample(3,1) = x[1];time_sample(3,2) = x[2];time_sample(3,3) = x[3];}
      else if ( t > 0.030001 && t < 0.040001){
  	time_sample(4,0) = x[0];time_sample(4,1) = x[1];time_sample(4,2) = x[2];time_sample(4,3) = x[3];}
      else if ( t > 0.040001 && t < 0.050001){
  	time_sample(5,0) = x[0];time_sample(5,1) = x[1];time_sample(5,2) = x[2];time_sample(5,3) = x[3];}
      else if ( t > 0.050001 && t < 0.060001){
  	time_sample(6,0) = x[0];time_sample(6,1) = x[1];time_sample(6,2) = x[2];time_sample(6,3) = x[3];}
      else if ( t > 0.060001 && t < 0.070001){
  	time_sample(7,0) = x[0];time_sample(7,1) = x[1];time_sample(7,2) = x[2];time_sample(7,3) = x[3];}
      else if ( t > 0.070001 && t < 0.080001){
  	time_sample(8,0) = x[0];time_sample(8,1) = x[1];time_sample(8,2) = x[2];time_sample(8,3) = x[3];}
      else if ( t > 0.080001 && t < 0.090001){
  	time_sample(9,0) = x[0];time_sample(9,1) = x[1];time_sample(9,2) = x[2];time_sample(9,3) = x[3];}
      else if ( t > 0.090001 && t < 0.100001){
  	time_sample(10,0) = x[0];time_sample(10,1) = x[1];time_sample(10,2) = x[2];time_sample(10,3) = x[3];}
      else if ( t > 0.100001 && t < 0.110001){
  	time_sample(11,0) = x[0];time_sample(11,1) = x[1];time_sample(11,2) = x[2];time_sample(11,3) = x[3];}
      else if ( t > 0.110001 && t < 0.120001){
  	time_sample(12,0) = x[0];time_sample(12,1) = x[1];time_sample(12,2) = x[2];time_sample(12,3) = x[3];}
      else if ( t > 0.120001 && t < 0.130001){
  	time_sample(13,0) = x[0];time_sample(13,1) = x[1];time_sample(13,2) = x[2];time_sample(13,3) = x[3];}
      else if ( t > 0.130001 && t < 0.140001){
  	time_sample(14,0) = x[0];time_sample(14,1) = x[1];time_sample(14,2) = x[2];time_sample(14,3) = x[3];}
      else if (t > 0.140001 && t < 0.150001){
  	time_sample(15,0) = x[0];time_sample(15,1) = x[1];time_sample(15,2) = x[2];time_sample(15,3) = x[3];}

      //Condition for applying discrepancy
      //if (unniversal_counter%600 == 0)
      

    






            disk_k1_input_model[0] = x[0];
            disk_k1_input_model[1] = x[1];           
            disk_k1_input_model[2] = x[2]; //*time_sample(15,3))/(kappa2*time_sample(15,1));;//This is not CO
            disk_k1_input_model[3] = x[3];


            disk_k3_input_model[0] = x[0];
            disk_k3_input_model[1] = x[1];           
            disk_k3_input_model[2] = x[2]; //*time_sample(15,3))/(kappa2*time_sample(15,1));;//This is not CO
            disk_k3_input_model[3] = x[3];

            disk_kap1_input_model[0] = x[0];
            disk_kap1_input_model[1] = x[1];           
            disk_kap1_input_model[2] = (x[2]*x[3])/(global_kappa2*x[1]); //*time_sample(15,3))/(kappa2*time_sample(15,1));;//This is not CO
            disk_kap1_input_model[3] = x[3];


            disk_kap2_input_model[0] = (x[2]*x[3])/(global_kappa2*x[1]);
            disk_kap2_input_model[1] = x[1];           
            disk_kap2_input_model[2] = x[2]; //*time_sample(15,3))/(kappa2*time_sample(15,1));;//This is not CO
            disk_kap2_input_model[3] = x[3];


            disk_kap3_input_model[0] = x[0];
            disk_kap3_input_model[1] = x[1];           
            disk_kap3_input_model[2] = x[2]; //*time_sample(15,3))/(kappa2*time_sample(15,1));;//This is not CO
            disk_kap3_input_model[3] = x[3];





          disk_k1.SplineEval(disk_k1_input_model, global_betas_k1, disk_k1_res);  
	  disk_k3.SplineEval(disk_k3_input_model, global_betas_k3, disk_k3_res);
	  disk_kap1.SplineEval(disk_kap1_input_model, global_betas_kap1, disk_kappa1_res);
	  disk_kap2.SplineEval(disk_kap2_input_model, global_betas_kap2, disk_kappa2_res);
          disk_kap3.SplineEval(disk_kap3_input_model, global_betas_kap3, disk_kappa3_res);    
	  

  global_k1 = global_k1*exp(disk_k1_res);
  global_k3 = global_k3*exp(disk_k3_res);
  global_kappa1 = global_kappa1*exp(disk_kappa1_res);
  global_kappa2 = global_kappa2*exp(disk_kappa2_res);
  global_kappa3 = global_kappa3*exp(disk_kappa3_res);
	

  double p_CO = (x[2]*x[3])/(x[1]*global_kappa2);

  double r1 =  global_k1*((x[0]*x[1]) - ((p_CO*x[3]*x[3]*x[3])/global_kappa1));
  double r3 = (global_k3*((x[0]*x[1]*x[1]) - ((x[2]*x[3]*x[3]*x[3]*x[3])/global_kappa2)));
     
  
  dxdt[0] = (1/atm_pres)*((R*T)/V)*(N_CH4 - (V*r1) - (V*r3) - x[0]*(N_TOTAL + (V*r1) + 2*(V*r3)));
  dxdt[1] = (1/atm_pres)*((R*T)/V)*(N_H2O - (V*r1) - 2*(V*r3) - x[1]*(N_TOTAL + (V*r1) + 2*(V*r3)));  
  dxdt[2] = (1/atm_pres)*((R*T)/V)*((V*r1) + (V*r3) - x[2]*(N_TOTAL + (V*r1) + 2*(V*r3)));
  dxdt[3] = (1/atm_pres)*((R*T)/V)*((V*r1) + 3*(V*r1) + 4*(V*r3) - x[3]*(N_TOTAL + (V*r1) + 2*(V*r3)));


    }

};

struct stiff_system_jacobi
{

double atm_pres = 101325;  //Atmospheric pressure in pascals


    void operator()( const vector_type &x  , matrix_type &J , const double & /* t */ , vector_type &dfdt )
    {
 
      
  double p_CO = (x[2]*x[3])/(x[1]*global_kappa2);

  double r1 = global_k1*((x[0]*x[1]) - ((p_CO*x[3]*x[3]*x[3])/global_kappa1));
  double r3 = global_k3*((x[0]*x[1]*x[1]) - ((x[2]*x[3]*x[3]*x[3]*x[3])/global_kappa2));  //maybe kappa3?? check




  double r1d1 = global_k1*(x[1]);                            
  double r1d2 = global_k1*(x[0]);                            
  double r1d3 = 0;
  double r1d4 = -global_k1*((3*p_CO*x[3]*x[3])/global_kappa1);
 
  double r3d1 = global_k3*(x[1]*x[1]);
  double r3d2 = global_k3*2*x[0]*x[1];                                   
  double r3d3 = -global_k3*((x[3]*x[3]*x[3]*x[3])/global_kappa2); 
  double r3d4 = -global_k3-((4*x[2]*x[3]*x[3]*x[3])/global_kappa2);


 
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

  
cout << global_kappa1 << " "  << global_kappa2 << " " << global_kappa3<< " " << global_k1 << " " << global_k3 <<  " " << N_CH4<<  " " << N_H2O<<  " " <<  R  <<  " " <<  V  <<  " " <<  T << " " << t <<  " " << endl;


  vector_type x( 4, 0.25);

try
  {
  size_t num_of_steps = integrate_const( make_dense_output< rosenbrock4< double > >( 1.0e-4 , 1.0e-4, 1.0e-4) , make_pair( stiff_system() , stiff_system_jacobi() ) ,x , 0.0 , t , 0.01,cout << phoenix::arg_names::arg2 << " " << phoenix::arg_names::arg1[0] << " " << phoenix::arg_names::arg1[1] << " " << phoenix::arg_names::arg1[2] << " " << phoenix::arg_names::arg1[3] << "\n" );  
  }

 catch(...)
   {std::cout << "Caught it !"<< std::endl;}

//std::cout << num_of_steps << std::endl;
    
    
  return time_sample;
}
