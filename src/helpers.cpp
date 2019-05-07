#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export()]]
LogicalVector determineGroup(CharacterVector names, CharacterVector compnames)
{
  LogicalVector results(names.size());
  for(int i = 0; i < names.size(); i++)
  {
//Rcpp::Rcout << i << std::endl;
    if(names[i] == "(_End_)"|| names[i] == "(_Sink_)")
    {
      results[i] = LogicalVector::get_na();
    }
    else
    {
      for(int j = 0; j < compnames.size(); j++)
      {
  //      Rcpp::Rcout << j << std::endl;
        if(names[i] == compnames[j])
        {
         results[i] = true;
         break;
         }
      }
    }
  }
  return results;
}


