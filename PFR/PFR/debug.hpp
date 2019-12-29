#pragma once
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <csignal>

//Function sets used for debugging

void get_vector(std::vector<std::string> &m)
{
  for (int i = 0; i < m.size(); i++)
    std::cout << m[i] << std::endl;
}

void get_vector(boost::numeric::ublas::vector<size_t> &m)
{
  for (int i = 0; i < m.size(); i++)
    std::cout << m[i] << std::endl;
}

float get_vector_element_(boost::numeric::ublas::vector<double> &m,size_t a)
{
  std::cout << m(a)<< std::endl;
}

// float get_vector_element_(boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double>> &m,size_t a)
// {
//   std::cout << m(a)<< std::endl;
// }

void get_vector(boost::numeric::ublas::vector<double> &m)

{
  std::cout<< m << std::endl;
  
}

void get_vector(boost::numeric::ublas::vector<int> &m)

{
  std::cout<< m << std::endl;
  
}

void get_all_vector_elements(boost::numeric::ublas::vector<float> &m)

{
  std::cout<< m << std::endl;
  
}

void get_vector(boost::numeric::ublas::vector<boost::numeric::ublas::matrix<double> >  &m)

{
  std::cout<< m << std::endl;
  
}

void get_vector(boost::numeric::ublas::vector<bool> &m)

{
  std::cout<< m << std::endl;
  
}

void get_vector(boost::numeric::ublas::vector<std::string> &m)

{
  std::cout<< m << std::endl;
  
}

float get_matrix_element(boost::numeric::ublas::matrix<double> &m,size_t a,size_t b)

{
  std::cout<< m(a,b) << std::endl;
  return m(a,b);
}

void get_matrix(boost::numeric::ublas::matrix<double> &m)

{
  std::cout<< m << std::endl;
  
}

void get_matrix(boost::numeric::ublas::matrix<int> &m)

{
  std::cout<< m << std::endl;
  
}

void get_matrix(boost::numeric::ublas::matrix<float> &m)

{
  std::cout<< m << std::endl;
  
}

void get_matrix(boost::numeric::ublas::matrix<std::string> &m)

{
  std::cout<< m << std::endl;
  
}

