#pragma once
#include <vector>
using namespace std;

class Test
{
public:

   // Функция - полином какого порядка
   int N;

   Test(const int& t_N) : N(t_N) {};

   Test() : N(0) {};

   double f(const double& X, const double& Y)
   {
      switch(N)
      {
      case(0): return (0)* lambda() + u(X, Y) * gamma();
      case(1): return (0)* lambda() + u(X, Y) * gamma();
      case(2): return (-4)* lambda() + u(X, Y) * gamma();
      case(3): return (-6 * X - 6 * Y)* lambda() + u(X, Y) * gamma();
      };
   }

   double lambda()
   {
      return 1;
   }

   double gamma()
   {
      return 1;
   }

   double theta()
   {
      return 1;
   }

   double u(const double& X, const double& Y)
   {
      switch(N)
      {
      case(0): return 2.0;
      case(1): return X + Y;
      case(2): return X * X + Y * Y;
      case(3): return X * X * X + Y * Y * Y;
      };
   }
};