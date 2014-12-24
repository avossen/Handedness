#include "FitResults.h"
#include <iostream>


using namespace std;



void FitResults::print()
{

  cout <<"=--------------FitResults Contents------"<<endl;
  cout <<"A1: " << A1 <<" +- " <<eA1 <<" A2: " << A2 <<" +- " <<eA2 <<" A3= " <<A3 << " +- " << eA3 <<endl;
  cout <<"-----------------------------------------------"<<endl;


};
ClassImp(FitResults);
