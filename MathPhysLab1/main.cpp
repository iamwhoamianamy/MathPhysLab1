#include <iostream>
#include "EllipticalProblem.h"

using namespace std;

void main()
{
   EllipticalProblem ep = EllipticalProblem();
   
   ep.FormMatrix();

   ep.slae->GaussSeidel(10000, 1e-14, 0.65);

   ep.PrintSolution();
   
}