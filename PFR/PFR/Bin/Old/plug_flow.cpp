//compile command: g++ -g -std=c++11 -I /home/nirjhar/boost_source_code/boost_1_62_0/ plug_flow.cpp -o test -larmadillo


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

// #define NUM_POINTS 5  // no of data points
// #define NUM_COMMANDS 3 //no of commands passed to gnuplot

using namespace std;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;

typedef boost::numeric::ublas::vector<double> uvec_t;



  
int main()

{


  boost::numeric::ublas::matrix<double> parameters(50,9); // size of rows and columns in file here



  double N_CH4,N_H2O;

  std::string confile = "Sample_Parameters.txt";//"BayesParameterResults_'time_stamp_here'.txt";BayesParameterResults_1516899681
                     
    std::ifstream constream;
    constream.open(confile.c_str(), std::ios::in);
    int i = 0;
    std::string line;
    while (!constream.eof()) {
        getline(constream, line, '\n');
        std::vector<std::string> str;
        if (line[0] != '#' && line[0] != ' ' && line.length() > 0) {
            boost::split(str, line, boost::is_any_of(" "));
           
                if (str.size() == 10) {
                  
		  
		  parameters(i,0) = boost::lexical_cast<double>(str[0]); parameters(i,1) = boost::lexical_cast<double>(str[1]); parameters(i,2) = boost::lexical_cast<double>(str[2]);parameters(i,3) = boost::lexical_cast<double>(str[3]);parameters(i,4) = boost::lexical_cast<double>(str[4]);parameters(i,5) = boost::lexical_cast<double>(str[5]);parameters(i,6) = boost::lexical_cast<double>(str[6]);parameters(i,7) = boost::lexical_cast<double>(str[7]);parameters(i,8) = boost::lexical_cast<double>(str[8]);
                }
	    

	}

	i++;}

    std::cout << i << std::endl;
   



   	boost::numeric::ublas::vector<double> betas(4);

    //solve and plot in gnuplot for selected paramter values.

	steady_state_output[0] = 0.25;steady_state_output[1] = 0.25;steady_state_output[2] = 0.25;steady_state_output[3] = 0.25;steady_state_output[4] = 0.25;
    
   int index = 0;
    for ( int q = 0; q < no_of_plots; q++)                             
      {
	index =  q;//enter desired index here;


	//make the beta vectors to be passed into the plotter function

	betas[0] = parameters(index,5); betas[1] = parameters(index,6); betas[2] = parameters(index,7); betas[3] = parameters(index,8);





	//Fills up yvals with values
	plotter(parameters(index,0),parameters(index,1),parameters(index,2),parameters(index,3),parameters(index,4),8.314,0.005,800,0.15,betas,betas,betas,betas,betas);
      }//see above comment to see what each passed number stands for. 


    if(steady_state)
      std::cout << "ha" << std::endl;




uvec_t xvals(rows);

for ( int a = 0; a < rows; a++)
      {
    	       xvals[a] = 0.01*a;
               
      }





 //Write results to text file data.txt
FILE * temp = fopen("data.txt", "w");
 


 // for (int a = 0; a < rows; a++)
 //   {
 //     fprintf(temp,"\n");
 //   for (int i = 0; i < no_of_plots; i++)
 //    {
 //      fprintf(temp, " %lf",yvals(a,i)); //Write the data to a temporary file
 //    }
 //   }





 // for (int a = 0; a < rows; a++)
 //   {
 //     fprintf(temp,"\n");
 //     for (int i = 0; i < successful_runs.size2(); i++)
 //    {
 //      fprintf(temp, " %lf",successful_runs(a,i)); //Write the data to a temporary file
 //    }
 //   }


//plot
 // FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");


 // char * commandsForGnuplot[] = {"set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5   # --- blue","set title \"TITLEEEEE\"", "plot 'data.temp' using 1:2 with linespoints ls 1", "replot 'data.temp' using 1:2 with linespoints ls 1"};//This acts like a stri


 // for (int i=0; i < 2; i++) //implements first 2 commands from the command string
 //    {
 //    fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
 //    }
   






 //     for (int q = 0; q < no_of_plots; q++)
 //   {


 //  FILE * temp = fopen("data.temp", "w");
 //    /*Opens an interface that one can use to send commands as if they were typing into the
 //     *     gnuplot command line.  "The -persistent" keeps the plot open even after your
 //     *     C program terminates.
 //     */
   
    
 //    for (int i = 0; i < NUM_POINTS; i++)
 //    {
 //      fprintf(temp, "%lf %lf \n",xvals[i],yvals(i,q)); //Write the data to a temporary file
 //    }


   
   





//PLOT

  
    
    // if ( q == 0)
    // fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[2]); //Send commands to gnuplot one by one.
    // else
    // fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[3]);

     fclose(temp);


     //}
    // printer(parameters(index,0),parameters(index,1),parameters(index,2),parameters(index,3),parameters(index,4),0.00278 ,0.00556,8.314,0.005,800,0.15,betas,betas,betas,betas,betas);

    //printer(betas);
 



 }
