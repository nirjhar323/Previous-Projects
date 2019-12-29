
// Copyright 2012 David S. Mebane.  Function evaluator for
// BSS-ANOVA-decomposed stochastic functions.

#ifndef _FUNC_BSSANOVA_
#define _FUNC_BSSANOVA_
#include "debug.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include "boost/algorithm/string.hpp"
#include <boost/lexical_cast.hpp>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "BasisSplineLoader.hpp" // Spline Coefficients for the Basis Functions


//Delete!
int screwup = 0; 


//Global Variables
double normxpar_issues = 0;





class func_bssanova {
    
public:                          // Public members can be freely accessed anywhere
    
   
private:
    
};

bool func_bssanova::configset(std::string const &confile) {
    
    
}

bool func_bssanova::evaluate(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<double> const &xpass, double &ypass) {
    
    
    
}

bool func_bssanova::function(boost::numeric::ublas::vector<double> const &par, boost::numeric::ublas::vector<double> const &xpass, double const &delt, double const &init, int const &ind, double &ypass) {
    
    
   

bool func_bssanova::SplineConfig(std::string const &SplineConf) {
    /* Method to parse through a config file that has comments in it. File also has flags that enable inversion of parameters. */
   
    //}
    
    SplineLoader BasLoad(numfunc, numdisc);
    BasLoad.Loader(BasisCoeffs);
    
  
   
    }
    return 1;
}
#endif
