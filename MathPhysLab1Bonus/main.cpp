#include <iostream>
#include "EllipticalProblem.h"

using namespace std;

void main()
{
   EllipticalProblem ep = EllipticalProblem();

   ep.ReadFormGrid("regions.txt");

   // ������������� ����
   ep.slae = new SLAE(ep.n_nodes);

   for(int reg_i = 0; reg_i < ep.n_regions; reg_i++)
   {
      Region* r = &ep.regions[reg_i];

      for(int node_i = 0; node_i < r->n_nodes; node_i++)
      {
         int global_i = node_i + r->first;

         ep.slae->index[0][global_i] = -r->n_x;
         ep.slae->index[1][global_i] = -1;
         ep.slae->index[2][global_i] = 0;
         ep.slae->index[3][global_i] = 1;
         ep.slae->index[4][global_i] = r->n_x;
      }
   }

   // ������������� �������� ������
   ep.test = Test(3);
   
   // ������������ ������� �������
   ep.FormMatrix();

   ep.slae->GaussSeidel(10000, 1e-14, 0.65);

   ep.PrintSolution("res.txt");
}