#pragma once
#include <vector>
using namespace std;

class Test
{
public:
   double f(const double& X, const double& Y)
   {
      return (0) * lambda() + u(X, Y) * gamma();
   }

   double lambda()
   {
      return 1;
   }

   double gamma()
   {
      return 1000;
   }

   double u(const double& X, const double& Y)
   {
      return X + Y;
   }
};