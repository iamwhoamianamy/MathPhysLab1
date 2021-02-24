#include <iostream>
#include "EllipticalProblem.h"

using namespace std;

void main()
{
   EllipticalProblem ep = EllipticalProblem();
   
   ep.FormMatrix();

   ep.slae->GaussSeidel(10000, 10e-14, 0.65);

   /*vector<double> prec(ep.slae->N);

   for (int j = 0; j < ep.Yn.size(); j++)
   {
      for (int i = 0; i < ep.Xn.size(); i++)
      {
         
         prec[j * ep.Xn.size() + i] = ep.test.u(ep.Xn[i], ep.Yn[j]);
      }
   }

   vector<double> calc(ep.slae->N);

   ep.slae->Multiplication(prec, calc);*/

   //for (int i = 0; i < ep.slae->N; i++)
   //{
   //   cout << calc[i] << endl;
   //}

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