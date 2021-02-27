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
   vector<double> Xr;                  // ������� ������� �� X
   vector<double> Yr;                  // ������� ������� �� Y

   vector<double> Xn;                  // ���������� ����� �� X
   vector<double> Yn;                  // ���������� ����� �� Y

   int N_X;                            // ���������� ����� �� X
   int N_Y;                            // ���������� ����� �� Y

   int x_bord;                         // ������ ���������� �������
                                       // L ������� �� X
   int y_bord;                         // ������ ���������� �������
                                       // L ������� �� Y

   SLAE* slae;                         // �������
   Test test;                          // �������� ����������

   EllipticalProblem()
   {
      Xr = vector<double>(3);
      Yr = vector<double>(3);

      ReadCoordLines("coords.txt");

      slae = new SLAE(N_X * N_Y, N_X);
   }

   ~EllipticalProblem()
   {
      delete slae;
   }

   // ������� ���������� ������ ������� �� ����� FILE_NAME,
   // ��������� ��������� �����
   void ReadCoordLines(const string& FILE_NAME)
   {
      std::ifstream fin(FILE_NAME);

      // ��������� ������� �������
      for(int i = 0; i < 3; i++)
         fin >> Xr[i];

      for(int i = 0; i < 3; i++)
         fin >> Yr[i];

      // ��������� ��������� ����� �� X
      int count = 0;        // ����� ����� �� X
      int t1, t2;

      fin >> t1;
      count += t1;

      fin >> t2;
      count += t2;

      N_X = count + 1;
      x_bord = t1;

      Xn.resize(N_X);

      for(int i = 0; i <= t1; i++)
         Xn[i] = Xr[0] + (Xr[1] - Xr[0]) / t1 * i;

      for(int i = 0; i <= t2; i++)
         Xn[i + t1] = Xr[1] + (Xr[2] - Xr[1]) / t2 * i;

      // ��������� ��������� ����� �� Y
      count = 0;        // ����� ����� �� Y

      fin >> t1;
      count += t1;

      fin >> t2;
      count += t2;

      N_Y = count + 1;
      y_bord = t1;

      Yn.resize(N_Y);

      for(int i = 0; i <= t1; i++)
         Yn[i] = Yr[0] + (Yr[1] - Yr[0]) / t1 * i;

      for(int i = 0; i <= t2; i++)
         Yn[i + t1] = Yr[1] + (Yr[2] - Yr[1]) / t2 * i;
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
            y_cent < y_bord - 1 && y_cent > 0 || 
            x_cent < x_bord - 1 && x_cent > 0 &&
            y_cent < N_Y - 1 && y_cent > 0)
         {
            // �������� �� X
            double hi = Xn[x_cent + 1] - Xn[x_cent + 0];
            double hi1 = Xn[x_cent - 0] - Xn[x_cent - 1];

            // �������� �� Y
            double hj = Yn[y_cent + 1] - Yn[y_cent + 0];
            double hj1 = Yn[y_cent - 0] - Yn[y_cent - 1];
            
            // ������ ����
            slae->matrix[0][n] = -test.lambda() *
               (2.0 / (hj1 * (hj + hj1)));

            // ����� ����
            slae->matrix[1][n] = -test.lambda() *
               (2.0 / (hi1 * (hi + hi1)));

            // ����������� ����
            slae->matrix[2][n] = test.lambda() *
               (2.0 / (hi1 * hi) + 2.0 / (hj1 * hj)) + test.gamma();

            // ������ ����
            slae->matrix[3][n] = -test.lambda() *
               (2.0 / (hi * (hi + hi1)));

            // ������� ����
            slae->matrix[4][n] = -test.lambda() *
               (2.0 / (hj * (hj + hj1)));

            // ������ ������ �����
            slae->f[n] = test.f(Xn[x_cent], Yn[y_cent]);
         }
         // ��������� �������� ����
         else if(x_cent <= x_bord || y_cent <= y_bord)
         {
            slae->matrix[2][n] = 1.0;
            slae->f[n] = test.u(Xn[x_cent], Yn[y_cent]);
         }
         else
            slae->matrix[2][n] = 1.0;
      }
   }

   // ����� �������
   void PrintSolution()
   {
      cout << "x           y           calc        prec       dif        N" << endl;
      cout << fixed;
      for(int j = 0; j < N_Y; j++)
      {
         for(int i = 0; i < N_X; i++)
         {
            int n = j * N_Y + i;
            cout << Yn[j];
            cout << setw(12) << Xn[i];
            double t = slae->xk[n];
            cout << setw(12) << t;
            double tt = test.u(Xn[i], Yn[j]);
            cout << setw(12) << tt;
            cout << setw(15) << scientific <<
               abs(t - tt);
            cout << fixed << setw(5) << n;

            if(j > y_bord && i > x_bord)
               cout << "  +";

            cout << endl;

         }
      }
   }
};