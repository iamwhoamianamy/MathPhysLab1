#pragma once
#include <vector>
using namespace std;

class Test
{
public:

   int N;

   Test(const int& t_N) : N(t_N) {};

   Test() : N(0) {};

   double f(const double& x, const double& y)
   {
      switch(N)
      {
      case(0): return (0)* lambda() + u(x, y) * gamma();
      case(1): return (0)* lambda() + u(x, y) * gamma();
      case(2): return (-4)* lambda() + u(x, y) * gamma();
      case(3): return (-6 * x - 6 * y) * lambda() + u(x, y) * gamma();
      case(4): return (-12 * x * x - 12 * y * y) * lambda() + u(x, y) * gamma();
      case(5): return 2 * sin(x + y) * lambda() + u(x, y) * gamma();
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

   vector<double> theta(const double& x, const double& y)
   {
      // Нормали вниз, вправо, вверх, влево
      switch(N)
      {
      case(0): return vector<double>(4, 0);
      case(1): return { -1, 1, 1, -1 };
      case(2): return { -2*x, 2*y, 2*x, -2*y };
      };
   }

   double u(const double& x, const double& y)
   {
      switch(N)
      {
      case(0): return 2.0;
      case(1): return x + y;
      case(2): return x * x + y * y;
      case(3): return x * x * x + y * y * y;
      case(4): return x * x * x * x + y * y * y * y;
      case(5): return sin(x + y);
      };
   }
};