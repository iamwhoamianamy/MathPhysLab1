#pragma once
#include <vector>

using namespace std;

class SLAE
{
public:
   vector<vector<double>> matrix; // ������� �������
   vector<vector<int>> index;     // ������� ��������
   vector<double> f;              // ������ ������ �����

   const int D = 5;               // ���������� ���������� �������
   int N = 0;                     // ����������� �������

   vector<double> xk, xk1;  // ��������������� �������

   SLAE(const int& t_n)
   {
      N = t_n;

      index.resize(D);

      for(int d = 0; d < D; d++)
         index[d].resize(N);

      matrix.resize(D);

      for(int d = 0; d < D; d++)
         matrix[d].resize(N);

      xk.resize(N);
      xk1.resize(N);
      f.resize(N);
   }

   // ��������� ������� ������� �� ������ vec,
   // ��������� � res
   void Multiplication(vector<double>& vec, vector<double>& res)
   {
      int n = vec.size(), k = 0;
      for(int j = 0; j < n; j++)
      {
         for(int d = 0; d < D; d++)
         {
            k = index[d][j];
            if(k < 0 && j + k > 0 ||
               k > 0 && j + k < n)
               res[j] += matrix[d][j] * vec[k + j];
         }
      }
   }

   // ����� �������
   double Norm(const vector<double>& vec)
   {
      double res = 0;
      for(int i = 0; i < N; i++)
         res += vec[i] * vec[i];

      return sqrt(res);
   }

   // ��������� ������������� ������� �������
   double RelativeResidual(vector<double>& vec)
   {
      vector<double> mult(N);

      Multiplication(vec, mult);
      for(size_t i = 0; i < N; i++)
         mult[i] = f[i] - mult[i];

      return Norm(mult) / Norm(f);
   }

   // ������������ ������� ������ ������-�������
   void IterativeProcess(const int& j, double& sum)
   {
      int k = 0, n = xk.size();
      for(int i = 0; i < D; i++)
      {
         k = index[i][j];
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
   void GaussSeidel(const int& max_iter, const double& eps,
      const double& relax)
   {
      double residual = 0.0, sum = 0.0;
      residual = RelativeResidual(xk);
      for(int k = 0; k < max_iter && residual > eps; k++)
      {
         for(int j = 0; j < N; j++)
         {
            IterativeProcess(j, sum);
            xk1[j] = xk[j] + (relax / matrix[2][j]) * (f[j] - sum);
            sum = 0.;
         }
         xk.swap(xk1);
         residual = RelativeResidual(xk);
      }
   }
};