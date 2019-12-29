#pragma once 
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
//#include "psopoint_ublas_2.hpp"
//#include "func_bssanova_2.hpp"
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


#include <stdlib.h>
#include <stdio.h>

typedef boost::numeric::ublas::vector< double > vector_type;
using namespace std;

int main()

{
  vector_type k(2)[0.01,0];
  // k[0] = 1; k[1] = 2;
  // cout << k << endl;
  //k[0] = 0.01;
  cout << k << endl;



}
