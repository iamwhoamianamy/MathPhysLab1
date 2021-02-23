#pragma once
#include <vector>
using namespace std;

class Test
{
public:
   double f(const double& X, const double& Y)
   {
      return X + Y;
   }

   double lambda()
   {
      return 1 ;
   }

   double gamma()
   {
      return 1 ;
   }

   double u(const double& X, const double& Y)
   {
      return X + Y;
   }
};