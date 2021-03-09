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
   vector<double> x_reg;               // Границы области по X
   vector<double> y_reg;               // Границы области по Y

   vector<double> x_node;              // Координаты узлов по X
   vector<double> y_node;              // Координаты узлов по Y

   int N_X;                            // Количество узлов по X
   int N_Y;                            // Количество узлов по Y

   int x_bord;                         // Индекс внутренней границы
                                       // L-области по X
   int y_bord;                         // Индекс внутренней границы
                                       // L-области по Y

   const int N_BORD = 6;               // Количество ребер

   vector<vector<int>> borders;        // Информация о граничных условиях

   SLAE* slae;                         // Система
   Test test;                          // Тестовая информация

   EllipticalProblem()
   {
      x_reg = vector<double>(3);
      y_reg = vector<double>(3);
   }

   ~EllipticalProblem()
   {
      delete slae;
   }

   // Функция считывания границ области из файла FILE_NAME
   // и формирования сетки
   void FormGridEven(const string& FILE_NAME)
   {
      ifstream fin(FILE_NAME);

      // Считываем границы области
      for(int i = 0; i < 3; i++)
         fin >> x_reg[i];

      for(int i = 0; i < 3; i++)
         fin >> y_reg[i];

      // Генерация координат узлов по X
      int n;
      double h, q;

      fin >> q >> n;
      //q = sqrt(q);
      //q = sqrt(q);
      N_X = n + 1;
      x_bord = n;
      x_node.resize(N_X);

      h = (x_reg[1] - x_reg[0]);
      if (q != 1)
         h *= (1 - q) / (1 - pow(q, n));
      else
         h /= n;

      x_node[0] = x_reg[0];
      for (int i = 0; i < n; i++)
         x_node[i + 1] = x_node[i] + h * pow(q, i);

      fin >> q >> n;
      //q = sqrt(q);
      //q = sqrt(q);
      N_X += n;
      x_node.resize(N_X);

      h = (x_reg[2] - x_reg[1]);
      if (q != 1)
         h *= (1 - q) / (1 - pow(q, n));
      else
         h /= n;

      for (int i = 1; i <= n; i++)
         x_node[i + x_bord] = x_node[i + x_bord - 1] + h * pow(q, i - 1);

      // Генерация координат узлов по Y
      fin >> q >> n;
      //q = sqrt(q);
      //q = sqrt(q);
      N_Y = n + 1;
      y_bord = n;
      y_node.resize(N_Y);

      h = (y_reg[1] - y_reg[0]);
      if (q != 1)
         h *= (1 - q) / (1 - pow(q, n));
      else
         h /= n;

      y_node[0] = y_reg[0];
      for (int i = 0; i < n; i++)
         y_node[i+1] = y_node[i] + h * pow(q, i);

      fin >> q >> n;
      //q = sqrt(q);
      //q = sqrt(q);
      N_Y += n;
      y_node.resize(N_Y);

      h = (y_reg[2] - y_reg[1]);
      if (q != 1)
         h *= (1 - q) / (1 - pow(q, n));
      else
         h /= n;

      for (int i = 1; i <= n; i++)
         y_node[i + y_bord] = y_node[i + y_bord - 1] + h * pow(q, i - 1);

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

   // Функция считывания индексов для описания краевых условий
   // из файла FILE_NAME
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

   // Функция формирования индексов границ ребер с соответсвующими
   // краевыми условиями
   void FormBordConditions()
   {
      for(int i = 0; i < N_BORD; i++)
      {
         borders[i][1] = CorrespondX(borders[i][1]);
         borders[i][2] = CorrespondX(borders[i][2]);
         borders[i][3] = CorrespondY(borders[i][3]);
         borders[i][4] = CorrespondY(borders[i][4]);
      }
   }

   // Формирование матрицы системы
   void FormMatrix()
   {
      for(int n = 0; n < slae->N; n++)
      {
         // Индексы центрального узла
         int x_cent = n % N_X;
         int y_cent = floor(n / N_X);

         // Обработка некраевых узлов внутри L-формы
         if(x_cent < N_X - 1 && x_cent > 0 &&
            y_cent < y_bord && y_cent > 0 || 
            x_cent < x_bord && x_cent > 0 &&
            y_cent < N_Y - 1 && y_cent > 0)
         {
            // Приросты по X
            double hi = x_node[x_cent + 1] - x_node[x_cent + 0];
            double hi1 = x_node[x_cent - 0] - x_node[x_cent - 1];

            // Приросты по Y
            double hj = y_node[y_cent + 1] - y_node[y_cent + 0];
            double hj1 = y_node[y_cent - 0] - y_node[y_cent - 1];
            
            // Нижний узел
            slae->matrix[0][n] = -test.lambda() *
               (2.0 / (hj1 * (hj + hj1)));

            // Левый узел
            slae->matrix[1][n] = -test.lambda() *
               (2.0 / (hi1 * (hi + hi1)));

            // Центральный узел
            slae->matrix[2][n] = +test.lambda() *
               (2.0 / (hi1 * hi) + 2.0 / (hj1 * hj)) + test.gamma();

            // Правый узел
            slae->matrix[3][n] = -test.lambda() *
               (2.0 / (hi * (hi + hi1)));

            // Верхний узел
            slae->matrix[4][n] = -test.lambda() *
               (2.0 / (hj * (hj + hj1)));

            // Вектор правой части
            slae->f[n] = test.f(x_node[x_cent], y_node[y_cent]);
         }
         // Обработка краевого узла
         else if(x_cent <= x_bord || y_cent <= y_bord)
         {
            // Обход по всем ребрам
            for(int b = 0; b < N_BORD; b++)
            {
               // Условие узла нахождения на ребре b
               if(x_cent >= borders[b][1] && x_cent <= borders[b][2] &&
                  y_cent >= borders[b][3] && y_cent <= borders[b][4])
               {
                  // Первое краевое условие
                  if(borders[b][0] == 0)
                  {
                     slae->matrix[2][n] = 1.0;
                     slae->f[n] = test.u(x_node[x_cent], y_node[y_cent]);
                  }
                  // Второе краевое условие
                  else if(borders[b][0] == 1)
                  {
                     // Если ребро параллельно оси X
                     if(borders[b][3] == borders[b][4])
                     {
                        double h = x_node[x_cent + 1] - x_node[x_cent - 1];

                        slae->matrix[1][n] = -test.lambda() / h;
                        slae->matrix[3][n] = +test.lambda() / h;

                        // Если нормаль направлена вниз
                        if(y_cent == 0)
                        {
                           slae->matrix[1][n] *= -1;
                           slae->matrix[3][n] *= -1;
                        }
                        // Если нормаль направлена вверх
                        else
                        {

                        }
                     }
                     // Если ребро параллельно оси Y
                     else if(borders[b][1] == borders[b][2])
                     {
                        double h = y_node[y_cent + 1] - y_node[y_cent - 1];

                        slae->matrix[0][n] = -test.lambda() / h;
                        slae->matrix[4][n] = +test.lambda() / h;

                        // Если нормаль направлена влево
                        if(x_cent == 0)
                        {
                           slae->matrix[0][n] *= -1;
                           slae->matrix[4][n] *= -1;
                        }
                        // Если нормаль направлена вправо
                        else
                        {

                        }
                     }

                     slae->matrix[2][n] = test.gamma();
                     slae->f[n] = test.f(x_node[x_cent], y_node[y_cent]);
                  }

                  break;
               }
            }
         }
         // Обработка квадрата за пределами L-формы
         else
            slae->matrix[2][n] = 1.0;
      }
   }

   // Вывод решения в файл FILE_NAME
   void PrintSolution(const string& FILE_NAME)
   {
      ofstream fout(FILE_NAME);
      int w = ceil(log10(N_X * N_Y)) + 2;
      double norm = 0., norm_u = 0.;

      fout << " x          y              calc           prec      dif         ";
      
      for(int i = 0; i < w - 1; i++)
         fout << " ";

      fout << "N  location" << endl << fixed;
      for(int j = 0; j < N_Y; j++)
      {
         for(int i = 0; i < N_X; i++)
         {
            int n = j * N_X + i;
            //if (i % 8 == 0 && j % 8 == 0)
            {
               fout << setw(9) << x_node[j];
               fout << setw(11) << y_node[i];
               double t = slae->xk[n];
               fout << setw(15) << t;
               double tt = 0;
               if (i <= x_bord || j <= y_bord) tt = test.u(x_node[i], y_node[j]);
               fout << setw(15) << tt;
               fout << setw(14) << scientific << abs(t - tt);
               fout << fixed << setw(w) << n;

               if (i < N_X - 1 && i > 0 &&
                  j < y_bord && j > 0 ||
                  i < x_bord && i > 0 &&
                  j < N_Y - 1 && j > 0)
                  fout << "  inner";
               else if (i <= x_bord || j <= y_bord)
                  fout << "  border";
               else
                  fout << "  outer";
               fout << endl;

               norm_u += tt * tt;
               norm += abs(t - tt) * abs(t - tt);
            }
         }
      }
      fout << "||u-u*||/||u*|| = " << scientific << sqrt(norm) / sqrt(norm_u) << endl;
      fout << "||u-u*|| = " << scientific << sqrt(norm);
      fout.close();
   }
};