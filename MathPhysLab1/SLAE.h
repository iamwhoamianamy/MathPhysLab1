#pragma once
#include <vector>
#include "Vector.h"

using namespace std;

class SLAE
{
public:
   vector<vector<double>> matrix; // Матрица системы
   vector<int> index;             // Индексы столбцов
   vector<double> f;              // Вектор правой части

   const int D = 5;               // Количество диагоналей матрицы
   int N = 0;                     // Размерность матрицы
   int M = 0;                     // Расстояние до крайних диагоналей

   vector<double> xk, xk1;  // Вспомогательные векторы

   SLAE(const int& N, const int& N_X)
   {
      M = N_X - 2;
      this->N = N;

      index.resize(D);
      index[0] = -(2 + M);
      index[1] = -1;
      index[2] = 0;
      index[3] = 1;
      index[4] = 2 + M;

      matrix.resize(D);

      for(int i = 0; i < D; i++)
         matrix[i].resize(N);

      xk.resize(N);
      xk1.resize(N);
      f.resize(N);
   }

   // Умножение матрицы системы на вектор vec,
   // результат в res
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

   // Норма вектора
   double Norm(const vector<double>& vec)
   {
      double res = 0;
      for (int i = 0; i < N; i++)
         res += vec[i] * vec[i];

      return sqrt(res);
   }

   // Получение относительной невязки системы
   double RelativeResidual(vector<double>& vec)
   {
      vector<double> mult(N);

      Multiplication(vec, mult);
      for (size_t i = 0; i < N; i++)
         mult[i] = f[i] - mult[i];

      return Norm(mult) / Norm(f);
   }

   // Итерационный процесс метода Гаусса-Зейделя
   void IterativeProcess(const int& j, double& sum)
   {
      int k = 0, n = xk.size();
      for(int i = 0; i < D; i++)
      {
         k = index[i];
         if(k + j >= 0 && k + j < n)
         {
            if(i < 3) // нижний треугольник
               sum += matrix[i][j] * xk1[k + j];
            else // верхний треугольник
               sum += matrix[i][j] * xk[k + j];
         }
      }
   }

   // Решение системы методом Гаусса-Зейделя
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
            xk1[j] = xk[j] + (RELAX / matrix[2][j]) * (f[j] - sum);
            sum = 0.;
         }
         xk.swap(xk1);
         residual = RelativeResidual(xk);
         //cout << k << " " << residual << endl;
      }
   }
};