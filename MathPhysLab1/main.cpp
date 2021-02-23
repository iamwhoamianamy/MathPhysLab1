#include <iostream>
#include "EllipticalProblem.h"

using namespace std;

void main()
{
   EllipticalProblem ep = EllipticalProblem();
   
   ep.FormMatrix();

   ep.slae->GaussSeidel(1000, 10e-10, 0.1);

   cout << "x           y           calc        prec       dif" << endl;
   
   for(int j = 0; j < ep.Yn.size(); j++)
   {
      for(int i = 0; i < ep.Xn.size(); i++)
      {
         cout << fixed << ep.Yn[j];
         cout << setw(12) << ep.Xn[i];
         double t = ep.slae->xk[j * ep.Xn.size() + i];
         cout << setw(12) << t;
         double tt = ep.test.u(ep.Xn[i], ep.Yn[j]);
         cout << setw(12) << tt;
         cout << setw(15) << scientific <<
            abs(t - tt) << endl;
      }
   }
   
}