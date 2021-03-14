#pragma once
#include <vector>

using namespace std;

class SLAE
{
public:
   vector<vector<double>> matrix; // Матрица системы
   vector<vector<int>> index;     // Индексы столбцов
   vector<double> f;              // Вектор правой части

   const int D = 5;               // Количество диагоналей матрицы
   int N = 0;                     // Размерность матрицы

   vector<double> xk, xk1;  // Вспомогательные векторы

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

   // Умножение матрицы системы на вектор vec,
   // результат в res
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

   // Норма вектора
   double Norm(const vector<double>& vec)
   {
      double res = 0;
      for(int i = 0; i < N; i++)
         res += vec[i] * vec[i];

      return sqrt(res);
   }

   // Получение относительной невязки системы
   double RelativeResidual(vector<double>& vec)
   {
      vector<double> mult(N);

      Multiplication(vec, mult);
      for(size_t i = 0; i < N; i++)
         mult[i] = f[i] - mult[i];

      return Norm(mult) / Norm(f);
   }

   // Итерационный процесс метода Гаусса-Зейделя
   void IterativeProcess(const int& j, double& sum)
   {
      int k = 0, n = xk.size();
      for(int i = 0; i < D; i++)
      {
         k = index[i][j];
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