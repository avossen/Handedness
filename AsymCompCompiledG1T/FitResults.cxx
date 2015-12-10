#include "FitResults.h"
#include <iostream>


using namespace std;



void FitResults::print()
{

  cout <<"=--------------FitResults Contents------"<<endl;
  cout <<"A1: " << A1 <<" +- " <<eA1<<"A1PP: " << A1PP <<" +- " <<eA1PP <<" A2: " << A2 <<" +- " <<eA2 <<"A2PP: " << A2PP <<" +- " <<eA2PP<<" A3= " <<A3 << " +- " << eA3 << " A3PP: "<< A3PP<<  " +- " << eA3PP<<endl;
  cout <<"-----------------------------------------------"<<endl;


};
ClassImp(FitResults);
