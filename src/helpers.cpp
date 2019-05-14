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

// [[Rcpp::export()]]
NumericVector handle_new_cs_cpp(NumericMatrix problem, NumericVector prices, //need to make sure prices are in the same order as they are in the problem matrix
                                double reso, NumericVector new_c_index)
{
  NumericVector suggestions(new_c_index.size());
  for (int i = 0; i < new_c_index.size(); ++i)
  {
    NumericVector distvec = problem(_, new_c_index[i]);
    int val_to_use;

    distvec = round((reso * distvec), 0);
    NumericVector comp_vals(distvec.size());
    for (int j = 0; j < distvec.size(); ++j)
    {
      if(Rcpp::traits::is_infinite<REALSXP>(distvec[j]))
      {
        comp_vals[j] = R_NegInf;

      }
      else if (Rcpp::internal::Rcpp_IsNA(prices[j]))
      {
        comp_vals[j] = NumericVector::get_na();
      }
      else
      {
        comp_vals[j] = prices[j] - distvec[j];
      }
      int maxA = max(comp_vals);
      val_to_use = maxA - 1;
    }
    suggestions[i] = val_to_use;
  }
  return(suggestions);
}

// [[Rcpp::export()]]
NumericVector handle_new_ts_cpp(NumericMatrix problem, NumericVector prices, //need to make sure prices are in the same order as they are in the problem matrix
                                double reso, NumericVector new_t_index)
{
  NumericVector suggestions(new_t_index.size());
  for (int i = 0; i < new_t_index.size(); ++i)
  {
    NumericVector distvec = problem(new_t_index[i], _);
    int val_to_use;

    distvec = round((reso * distvec), 0);
    NumericVector comp_vals(distvec.size());
    for (int j = 0; j < distvec.size(); ++j)
    {
      if(Rcpp::traits::is_infinite<REALSXP>(distvec[j]))
      {
        comp_vals[j] = R_PosInf;

      }
      else if (Rcpp::internal::Rcpp_IsNA(prices[j]))
      {
        comp_vals[j] = NumericVector::get_na();
      }
      else
      {
        comp_vals[j] = distvec[j] + prices[j];
      }
      int minA = min(comp_vals);
      val_to_use = minA + 1;
    }
    suggestions[i] = val_to_use;
  }
  return(suggestions);
}


