#pragma once
#include <vector>
#include "Vector.h"

using namespace std;

class SLAE
{
public:
   vector<vector<double>> matrix; // ������� �������
   vector<int> index;             // ������� ��������
   vector<double> f;              // ������ ������ �����

   const int D = 5;               // ���������� ���������� �������
   int N = 0;                     // ����������� �������
   int GAP = 0;                   // ���������� �� ������� ����������

   vector<double> xk, xk1;  // ��������������� �������

   //const double EPS = 10e-10;     // �������� �������
   //const double MAX_ITER = 1000;  // ������������ ����� ��������

   SLAE(const int& N, const int& N_X)
   {
      GAP = N_X - 1;
      this->N = N;

      index.resize(D);
      index[0] = -(2 + GAP);
      index[1] = -1;
      index[2] = 0;
      index[3] = 1;
      index[4] = 2 + GAP;

      for(int i = 0; i < D; i++)
         matrix.resize(N);

      xk.resize(N);
      xk1.resize(N);
   }

   // ��������� ������� ������� �� ������ vec,
   // ��������� � res
   void Multiplication(vector<double>& vec, vector<double>& res)
   {
      int n = vec.size(), k = 0;
      for(int i = 0; i < D; i++)
      {
         k = index[i];
         if(k < 0)
            for(int j = abs(k); j < n; j++)
               res[j] += matrix[i][j] * vec[k + j];
         else
            for(int j = 0; j < n - k; j++)
               res[j] += matrix[i][j] * vec[k + j];
      }
   }

   // ��������� ������������� ������� �������
   double RelativeResidual(vector<double>& vec)
   {
      vector<double> mult(N);

      Multiplication(vec, mult);
       mult = f - mult;

      return Norm(mult) / Norm(f);
   }

   // ������������ ������� ������ ������-�������
   void IterativeProcess(const int& j, double& sum)
   {
      int k = 0, n = xk.size();
      for(int i = 0; i < D; i++)
      {
         k = index[i];
         if(k + j >= 0 && k + j < n)
         {
            if(i < 3) // ������ �����������
               sum += matrix[i][j] * xk1[k + j];
            else // ������� �����������
               sum += matrix[i][j] * xk[k + j];
         }
      }
   }

   // ������� ������� ������� ������-�������
   void GaussSeidel(const int& MAX_ITER, const double& EPS,
      const double& RELAX)
   {
      double residual = 0.0, sum = 0.0;
      residual = RelativeResidual(xk);
      for(int k = 0; k < MAX_ITER && residual > EPS; k++)
      {
         for(int j = 0; j < N; j++)
         {
            IterativeProcess(j, sum);
            xk1[j] = xk[j] + (RELAX / matrix[3][j]) * (f[j] - sum);
            sum = 0.;
         }
         xk.swap(xk1);
         residual = RelativeResidual(xk);
         cout << k << " " << residual << endl;
      }
   }
};