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
         //cout << k << " " << residual << endl;
      }
   }
};


//
//   // Вывод решения в файл FILE_NAME
//   void PrintSolution(const string& FILE_NAME)
//   {
//      ofstream fout(FILE_NAME);
//      int w = ceil(log10(N_X * N_Y)) + 2;
//      double norm = 0., norm_u = 0.;
//
//      fout << " y          x              calc           prec      dif         ";
//      
//      for(int i = 0; i < w - 1; i++)
//         fout << " ";
//
//      fout << "N  location" << endl << fixed;
//      for(int j = 0; j < N_Y; j++)
//      {
//         for(int i = 0; i < N_X; i++)
//         {
//            int n = j * N_X + i;
//            //if (i % 8 == 0 && j % 8 == 0)
//            {
//               fout << setw(9) << y_node[j];
//               fout << setw(11) << x_node[i];
//               double t = slae->xk[n];
//               fout << setw(15) << t;
//               double tt = 0;
//               if (i <= x_bord || j <= y_bord) tt = test.u(x_node[i], y_node[j]);
//               fout << setw(15) << tt;
//               fout << setw(14) << scientific << abs(t - tt);
//               fout << fixed << setw(w) << n;
//
//               if (i < N_X - 1 && i > 0 &&
//                  j < y_bord && j > 0 ||
//                  i < x_bord && i > 0 &&
//                  j < N_Y - 1 && j > 0)
//                  fout << "  inner";
//               else if (i <= x_bord || j <= y_bord)
//                  fout << "  border";
//               else
//                  fout << "  outer";
//               fout << endl;
//
//               norm_u += tt * tt;
//               norm += abs(t - tt) * abs(t - tt);
//            }
//         }
//      }
//      fout << "||u-u*||/||u*|| = " << scientific << sqrt(norm) / sqrt(norm_u) << endl;
//      fout << "||u-u*|| = " << scientific << sqrt(norm);
//      fout.close();
//   }
//;