#include <iostream>
#include "EllipticalProblem.h"

using namespace std;

void main()
{
   EllipticalProblem ep = EllipticalProblem();

   ep.ReadCoordLinesEven("coords.txt");
   ep.FormGridEven();

   ep.ReadBordConditions("borders.txt");
   ep.FormBordConditions();
   
   // ������������� ����
   ep.slae = new SLAE(ep.N_X * ep.N_Y, ep.N_X);

   // ������������� �������� ������
   ep.test = Test(2);
   
   ep.FormMatrix();

   ep.slae->GaussSeidel(10000, 1e-14, 0.65);

   ep.PrintSolution();
   
}