#pragma once
#include <vector>
#include <fstream>
#include <string>
#include "SLAE.h"
#include "Test.h"

using namespace std;

class EllipticalProblem
{
public:
   vector<double> x_reg;               // ������� ������� �� X
   vector<double> y_reg;               // ������� ������� �� Y

   vector<double> x_node;              // ���������� ����� �� X
   vector<double> y_node;              // ���������� ����� �� Y

   int N_X;                            // ���������� ����� �� X
   int N_Y;                            // ���������� ����� �� Y

   int x_bord;                         // ������ ���������� �������
                                       // L-������� �� X
   int y_bord;                         // ������ ���������� �������
                                       // L-������� �� Y

   const int N_BORD = 6;               // ���������� �����

   vector<vector<int>> borders;        // ���������� � ��������� ��������

   SLAE* slae;                         // �������
   Test test;                          // �������� ����������

   EllipticalProblem()
   {
      x_reg = vector<double>(3);
      y_reg = vector<double>(3);

      ReadCoordLines("coords.txt");

      double h;
      // ��������� ��������� ����� �� X
      
      h = (x_reg[1] - x_reg[0]) / x_bord;

      for(int i = 0; i <= x_bord; i++)
         x_node[i] = x_reg[0] + h * i;

      h = (x_reg[2] - x_reg[1]) / (N_X - x_bord - 1);

      for(int i = 0; i <= N_X - x_bord - 1; i++)
         x_node[i + x_bord] = x_reg[1] + h * i;

      // ��������� ��������� ����� �� Y
      h = (y_reg[1] - y_reg[0]) / y_bord;

      for(int i = 0; i <= y_bord; i++)
         y_node[i] = y_reg[0] + h * i;

      h = (y_reg[2] - y_reg[1]) / (N_Y - y_bord - 1);

      for(int i = 0; i <= (N_Y - y_bord - 1); i++)
         y_node[i + y_bord] = y_reg[1] + h * i;

      ReadBordConditions("borders.txt");

      // ������������ �������� ������ ����� � ���������������
      // �������� ���������
      for(int i = 0; i < N_BORD; i++)
      {
         borders[i][1] = CorrespondX(borders[i][1]);
         borders[i][2] = CorrespondX(borders[i][2]);
         borders[i][3] = CorrespondY(borders[i][3]);
         borders[i][4] = CorrespondY(borders[i][4]);
      }

      // ������������� ����
      slae = new SLAE(N_X * N_Y, N_X);

      // ������������� �������� ������
      test = Test(2);
   }

   ~EllipticalProblem()
   {
      delete slae;
   }

   // ������� ���������� ������ ������� �� ����� FILE_NAME
   void ReadCoordLines(const string& FILE_NAME)
   {
      ifstream fin(FILE_NAME);

      // ��������� ������� �������
      for(int i = 0; i < 3; i++)
         fin >> x_reg[i];

      for(int i = 0; i < 3; i++)
         fin >> y_reg[i];

      // ���������� ������������ ����� �� X
      int count = 0;    // ����� ����� �� X
      int t1, t2;

      fin >> t1;
      count += t1;

      fin >> t2;
      count += t2;

      N_X = count + 1;
      x_bord = t1;

      x_node.resize(N_X);

      // ���������� ������������ ����� �� Y
      count = 0;        // ����� ����� �� Y

      fin >> t1;
      count += t1;

      fin >> t2;
      count += t2;

      N_Y = count + 1;
      y_bord = t1;

      y_node.resize(N_Y);

      fin.close();
   }

   int CorrespondX(const int& I)
   {
      switch(I)
      {
      case(0): return 0;
      case(1): return x_bord;
      case(2): return N_X - 1;
      }
   }
   int CorrespondY(const int& I)
   {
      switch(I)
      {
      case(0): return 0;
      case(1): return y_bord;
      case(2): return N_Y - 1;
      }
   }

   // ������� ���������� �������� ��� �������� ������� �������
   // �� ����� FILE_NAME
   void ReadBordConditions(const string& FILE_NAME)
   {
      ifstream fin(FILE_NAME);

      borders.resize(N_BORD);

      for(int i = 0; i < N_BORD; i++)
      {
         borders[i].resize(5);

         for(int j = 0; j < 5; j++)
            fin >> borders[i][j];
      }

      fin.close();
   }

   // ������������ ������� �������
   void FormMatrix()
   {
      for(int n = 0; n < slae->N; n++)
      {
         // ������� ������������ ����
         int x_cent = n % N_X;
         int y_cent = floor(n / N_X);

         // ��������� ��������� ����� ������ L-�����
         if(x_cent < N_X - 1 && x_cent > 0 &&
            y_cent < y_bord && y_cent > 0 || 
            x_cent < x_bord && x_cent > 0 &&
            y_cent < N_Y - 1 && y_cent > 0)
         {
            // �������� �� X
            double hi = x_node[x_cent + 1] - x_node[x_cent + 0];
            double hi1 = x_node[x_cent - 0] - x_node[x_cent - 1];

            // �������� �� Y
            double hj = y_node[y_cent + 1] - y_node[y_cent + 0];
            double hj1 = y_node[y_cent - 0] - y_node[y_cent - 1];
            
            // ������ ����
            slae->matrix[0][n] = -test.lambda() *
               (2.0 / (hj1 * (hj + hj1)));

            // ����� ����
            slae->matrix[1][n] = -test.lambda() *
               (2.0 / (hi1 * (hi + hi1)));

            // ����������� ����
            slae->matrix[2][n] = +test.lambda() *
               (2.0 / (hi1 * hi) + 2.0 / (hj1 * hj)) + test.gamma();

            // ������ ����
            slae->matrix[3][n] = -test.lambda() *
               (2.0 / (hi * (hi + hi1)));

            // ������� ����
            slae->matrix[4][n] = -test.lambda() *
               (2.0 / (hj * (hj + hj1)));

            // ������ ������ �����
            slae->f[n] = test.f(x_node[x_cent], y_node[y_cent]);
         }
         // ��������� �������� ����
         else if(x_cent <= x_bord || y_cent <= y_bord)
         {
            // ����� �� ���� ������
            for(int b = 0; b < N_BORD; b++)
            {
               // ������� ���� ���������� �� ����� b
               if(x_cent >= borders[b][1] && x_cent <= borders[b][2] &&
                  y_cent >= borders[b][3] && y_cent <= borders[b][4])
               {
                  // ������ ������� �������
                  if(borders[b][0] == 0)
                  {
                     slae->matrix[2][n] = 1.0;
                     slae->f[n] = test.u(x_node[x_cent], y_node[y_cent]);
                  }
                  // ������ ������� �������
                  else if(borders[b][0] == 1)
                  {
                     // ���� ����� ����������� ��� X
                     if(borders[b][3] == borders[b][4])
                     {
                        double h = x_node[x_cent + 1] - x_node[x_cent - 1];

                        slae->matrix[1][n] = -test.lambda() / h;
                        slae->matrix[3][n] = +test.lambda() / h;

                        // ���� ������� ���������� ����
                        if(y_cent == 0)
                        {
                           slae->matrix[1][n] *= -1;
                           slae->matrix[3][n] *= -1;
                        }
                     }
                     // ���� ����� ����������� ��� Y
                     else if(borders[b][1] == borders[b][2])
                     {
                        double h = y_node[y_cent + 1] - y_node[y_cent - 1];

                        slae->matrix[0][n] = -test.lambda() / h;
                        slae->matrix[4][n] = +test.lambda() / h;

                        // ���� ������� ���������� �����
                        if(x_cent == 0)
                        {
                           slae->matrix[0][n] *= -1;
                           slae->matrix[4][n] *= -1;
                        }
                     }

                     slae->matrix[2][n] = test.gamma();
                     slae->f[n] = test.f(x_node[x_cent], y_node[y_cent]);
                  }

                  break;
               }
            }
         }
         // ��������� �������� �� ��������� L-�����
         else
            slae->matrix[2][n] = 1.0;
      }
   }

   // ����� �������
   void PrintSolution()
   {
      int w = ceil(log10(N_X * N_Y)) + 2;

      cout << " x          y              calc           prec      dif         ";
      
      for(int i = 0; i < w - 1; i++)
         cout << " ";

      cout << "N  location" << endl << fixed;
      for(int j = 0; j < N_Y; j++)
      {
         for(int i = 0; i < N_X; i++)
         {
            int n = j * N_X + i;
            //if (i % 2 == 0 && j % 2 == 0)
            {
               cout << setw(9) << y_node[j];
               cout << setw(11) << x_node[i];
               double t = slae->xk[n];
               cout << setw(15) << t;
               double tt = 0;
               if (i <= x_bord || j <= y_bord) tt = test.u(x_node[i], y_node[j]);
               cout << setw(15) << tt;
               cout << setw(14) << scientific <<
                  abs(t - tt);
               cout << fixed << setw(w) << n;

               if (i < N_X - 1 && i > 0 &&
                  j < y_bord && j > 0 ||
                  i < x_bord && i > 0 &&
                  j < N_Y - 1 && j > 0)
                  cout << "  inner";
               else if (i <= x_bord || j <= y_bord)
                  cout << "  border";
               else
                  cout << "  outer";

               cout << endl;
            }
         }
      }
   }
};