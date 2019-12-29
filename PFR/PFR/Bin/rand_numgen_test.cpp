#include <iostream>
#include <random>
//#include <cstdlib>
#include <ctime>

int main()

{

  
  
  for (int i = 0; i < 20; i++){

    int index = 0;

   ;

   //these generate the same random number always, gottas use random seed(time)   
    
    
    
    

    int upper_bound = 47000;
    int lower_bound = 46000;
    
    //int number2 = floor(number*(upper_bound - lower_bound));
    
    srand((unsigned)time(0)+i);
    int rand_number = (rand()%1000);
    index =  lower_bound + rand_number;


    //std::cout << rand_number << std::endl;
    


     std::cout << index << "\n" << std::endl;

  }

  

}
