#include <iostream>
#include "EllipticalProblem.h"

using namespace std;

void main()
{
   EllipticalProblem ep = EllipticalProblem();
   


   ep.slae->GaussSeidel(10e-10, 1000, 1);

}