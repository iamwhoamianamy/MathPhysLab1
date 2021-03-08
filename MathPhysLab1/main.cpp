#include <iostream>
#include "EllipticalProblem.h"

using namespace std;

void main()
{
   EllipticalProblem ep = EllipticalProblem();

   ep.FormGridEven("coords.txt");

   ep.ReadBordConditions("borders.txt");
   ep.FormBordConditions();
   
   // Инициализация СЛАУ
   ep.slae = new SLAE(ep.N_X * ep.N_Y, ep.N_X);

   // Инициализация тестовых данных
   ep.test = Test(3);
   
   ep.FormMatrix();

   ep.slae->GaussSeidel(10000, 1e-14, 0.65);

   ep.PrintSolution("res.txt");
   
}