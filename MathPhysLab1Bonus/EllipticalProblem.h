#pragma once
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include "SLAE.h"
#include "Test.h"
#include "Region.h"

using namespace std;

class EllipticalProblem
{
public:

   vector<Region> regions;             // Регионы расчетной области

   int n_regions = 0;                  // Количество регионов
   int n_nodes = 0;                    // Общее количество узлов

   vector<vector<int>> borders;        // Информация о граничных условиях

   SLAE* slae;                         // Система
   Test test;                          // Тестовая информация

   EllipticalProblem()
   {

   }

   ~EllipticalProblem()
   {
      delete slae;
   }

   // Функция считывания областей из файла FILE_NAME
   // и формирования сетки
   void ReadFormGrid(const string& FILE_NAME)
   {
      ifstream fin(FILE_NAME);

      fin >> n_regions;
      string s;

      regions.resize(n_regions);

      for(int reg_i = 0; reg_i < n_regions; reg_i++)
      {
         fin >> s;
         Region* r = &regions[reg_i];

         // Считывание границы области
         fin >> r->left;
         fin >> r->right;
         fin >> r->bot;
         fin >> r->top;

         // Генерация координат узлов по X
         int n;
         double h, q;

         fin >> q >> n;

         r->n_x = n + 1;
         r->x_node.resize(r->n_x);

         h = r->right - r->left;

         if(q != 1)
            h *= (1 - q) / (1 - pow(q, n));
         else
            h /= n;

         r->x_node[0] = r->left;

         for(int i = 0; i < n; i++)
            r->x_node[i + 1] = r->x_node[i] + h * pow(q, i);

         // Генерация координат узлов по Y

         fin >> q >> n;

         r->n_y = n + 1;
         r->y_node.resize(r->n_y);

         h = r->top - r->bot;

         if(q != 1)
            h *= (1 - q) / (1 - pow(q, n));
         else
            h /= n;

         r->y_node[0] = r->bot;

         for(int i = 0; i < n; i++)
            r->y_node[i + 1] = r->y_node[i] + h * pow(q, i);

         if(reg_i != 0)
            r->first = regions[reg_i - 1].last + 1;
         else
            r->first = 0;

         r->n_nodes = r->n_x * r->n_y;

         r->last = r->first + r->n_nodes - 1;

         n_nodes += r->n_nodes;

         r->borders.resize(4);
         // Считывание информации о краевых условиях
         for(int bord_i = 0; bord_i < 4; bord_i++)
            fin >> r->borders[bord_i];
      }

      fin.close();
   }

   //   int CorrespondX(const int& I)
   //   {
   //      switch(I)
   //      {
   //      case(0): return 0;
   //      case(1): return x_bord;
   //      case(2): return N_X - 1;
   //      }
   //   }
   //   int CorrespondY(const int& I)
   //   {
   //      switch(I)
   //      {
   //      case(0): return 0;
   //      case(1): return y_bord;
   //      case(2): return N_Y - 1;
   //      }
   //   }
   //
   //   // Функция считывания индексов для описания краевых условий
   //   // из файла FILE_NAME
   //   void ReadBordConditions(const string& FILE_NAME)
   //   {
   //      ifstream fin(FILE_NAME);
   //
   //      borders.resize(N_BORD);
   //
   //      for(int i = 0; i < N_BORD; i++)
   //      {
   //         borders[i].resize(5);
   //
   //         for(int j = 0; j < 5; j++)
   //            fin >> borders[i][j];
   //      }
   //
   //      fin.close();
   //   }
   //
   //   // Функция формирования индексов границ ребер с соответсвующими
   //   // краевыми условиями
   //   void FormBordConditions()
   //   {
   //      for(int i = 0; i < N_BORD; i++)
   //      {
   //         borders[i][1] = CorrespondX(borders[i][1]);
   //         borders[i][2] = CorrespondX(borders[i][2]);
   //         borders[i][3] = CorrespondY(borders[i][3]);
   //         borders[i][4] = CorrespondY(borders[i][4]);
   //      }
   //   }
   //

   // Формирование матрицы системы
   void FormMatrix()
   {
      // Проход по всем регионам
      for(int reg_i = 0; reg_i < n_regions; reg_i++)
      {
         Region* r = &regions[reg_i];

         // Проход по всем узлам региона
         for(int node_i = 0; node_i < r->n_nodes; node_i++)
         {
            // Индекс узла в глобальной нумерации
            int global_i = node_i + r->first;

            // Индексы центрального узла
            int x_cent = node_i % r->n_x;
            int y_cent = floor(node_i / r->n_x);

            // Обработка некраевых узлов
            if(0 < x_cent && x_cent < r->n_x - 1 &&
               0 < y_cent && y_cent < r->n_y - 1)
            {
               // Приросты по X
               double hi = r->x_node[x_cent + 1] - r->x_node[x_cent + 0];
               double hi1 = r->x_node[x_cent - 0] - r->x_node[x_cent - 1];

               // Приросты по Y
               double hj = r->y_node[y_cent + 1] - r->y_node[y_cent + 0];
               double hj1 = r->y_node[y_cent - 0] - r->y_node[y_cent - 1];

               // Нижний узел
               slae->matrix[0][global_i] = -test.lambda() *
                  (2.0 / (hj1 * (hj + hj1)));

               // Левый узел
               slae->matrix[1][global_i] = -test.lambda() *
                  (2.0 / (hi1 * (hi + hi1)));

               // Центральный узел
               slae->matrix[2][global_i] = +test.lambda() *
                  (2.0 / (hi1 * hi) + 2.0 / (hj1 * hj)) + test.gamma();

               // Правый узел
               slae->matrix[3][global_i] = -test.lambda() *
                  (2.0 / (hi * (hi + hi1)));

               // Верхний узел
               slae->matrix[4][global_i] = -test.lambda() *
                  (2.0 / (hj * (hj + hj1)));

               // Вектор правой части
               slae->f[global_i] = test.f(r->x_node[x_cent], r->y_node[y_cent]);
            }
            // Обработка краевых узлов
            else
            {
               int border_i;

               if(x_cent == 0)
                  border_i = 3;
               else if(x_cent == r->n_x - 1)
                  border_i = 1;

               if(y_cent == 0)
                  border_i = 0;
               else if(y_cent == r->n_y - 1)
                  border_i = 2;

               // Первое краевое
               switch(r->borders[border_i])
               {
                  // Первое краевое
                  case 0:
                  {
                     slae->matrix[0][global_i] = 0;
                     slae->matrix[1][global_i] = 0;
                     slae->matrix[2][global_i] = 1.0;
                     slae->matrix[3][global_i] = 0;
                     slae->matrix[4][global_i] = 0;
                     slae->f[global_i] = test.u(r->x_node[x_cent], r->y_node[y_cent]);
                     break;
                  }
                  // Второе краевое
                  case 1:
                  {
                     //switch(border_i)
                     //{
                     //// Нормаль направлена вниз
                     //case 0:
                     //   double dy = r->y_node[y_cent + 1] - r->y_node[y_cent];
                     //   slae->matrix[2][global_i] = test.lambda() / dy;
                     //   slae->matrix[4][global_i] = -test.lambda() / dy;
                     //   slae->f[global_i] = test.theta(r->x_node[x_cent], r->y_node[y_cent])[0];
                     //   break;
                     //// Если нормаль направлена вправо
                     //case 1:
                     //   double dx = x_node[x_cent] - x_node[x_cent - 1];
                     //   slae->matrix[1][n] = -test.lambda() / dx;
                     //   slae->matrix[2][n] = test.lambda() / dx;
                     //   slae->f[n] = test.theta(x_node[x_cent], y_node[y_cent])[1];
                     //   break;
                     //}
                     break;
                  }
                  // Нет краевого условия
                  case -1:
                  {
                     double hi = 0, hi1 = 0, hj = 0, hj1 = 0;
                     
                     if(x_cent == 0)
                        hi1 = hi = r->x_node[x_cent + 1] - r->x_node[x_cent];
                     else if(x_cent == r->n_x - 1)
                        hi = hi1 = r->x_node[x_cent] - r->x_node[x_cent - 1];
                     else
                     {
                        hi = r->x_node[x_cent + 1] - r->x_node[x_cent];
                        hi1 = r->x_node[x_cent] - r->x_node[x_cent - 1];
                     }

                     if(y_cent == 0)
                        hj1 = hj = r->y_node[y_cent + 1] - r->y_node[y_cent];
                     else if(y_cent == r->n_y - 1)
                        hj = hj1 = r->y_node[y_cent] - r->y_node[y_cent - 1];
                     else
                     {
                        hj = r->y_node[y_cent + 1] - r->y_node[y_cent];
                        hj1 = r->y_node[y_cent] - r->y_node[y_cent - 1];
                     }

                     //// Нижний узел
                     //if(y_cent != 0)
                     //   slae->matrix[0][global_i] = -test.lambda() *
                     //      (2.0 / (hj1 * (hj + hj1)));

                     //// Левый узел
                     //if(x_cent != 0)
                     //   slae->matrix[1][global_i] = -test.lambda() *
                     //      (2.0 / (hi1 * (hi + hi1)));

                     //// Центральный узел
                     //slae->matrix[2][global_i] = +test.lambda() *
                     //   (2.0 / (hi1 * hi) + 2.0 / (hj1 * hj)) + test.gamma();

                     //// Правый узел
                     //if(x_cent != r->n_x - 1)
                     //   slae->matrix[3][global_i] = -test.lambda() *
                     //      (2.0 / (hi * (hi + hi1)));

                     //// Верхний узел
                     //if(y_cent != r->n_y - 1)
                     //   slae->matrix[4][global_i] = -test.lambda() *
                     //      (2.0 / (hj * (hj + hj1)));

                     // Нижний узел
                     slae->matrix[0][global_i] = -test.lambda() *
                        (2.0 / (hj1 * (hj + hj1)));

                     // Левый узел
                     slae->matrix[1][global_i] = -test.lambda() *
                        (2.0 / (hi1 * (hi + hi1)));

                     // Центральный узел
                     slae->matrix[2][global_i] = +test.lambda() *
                        (2.0 / (hi1 * hi) + 2.0 / (hj1 * hj)) + test.gamma();

                     // Правый узел
                     slae->matrix[3][global_i] = -test.lambda() *
                        (2.0 / (hi * (hi + hi1)));

                     // Верхний узел
                     slae->matrix[4][global_i] = -test.lambda() *
                        (2.0 / (hj * (hj + hj1)));

                     // Вектор правой части
                     slae->f[global_i] = test.f(r->x_node[x_cent], r->y_node[y_cent]);
                     break;
                  }
               }
            }
         }
      }
   }


   void PrintSolution(const string& file_name)
   {
      ofstream fout(file_name);
      double norm = 0.0, norm_u = 0.0;

      fout << " y          x              calc           prec";
      fout << "      dif           N   reg  location" << endl << fixed;

      // Проход по всем регионам
      for(int reg_i = 0; reg_i < n_regions; reg_i++)
      {
         Region* r = &regions[reg_i];

         // Проход по всем узлам региона
         for(int node_i = 0; node_i < r->n_nodes; node_i++)
         {
            // Индекс узла в глобальной нумерации
            int global_i = node_i + r->first;

            // Индексы центрального узла
            int x_cent = node_i % r->n_x;
            int y_cent = floor(node_i / r->n_x);

            fout << setw(9) << r->y_node[y_cent];
            fout << setw(11) << r->x_node[x_cent];

            double calc = slae->xk[global_i];
            fout << setw(15) << calc;
            double prec = test.u(r->x_node[x_cent], r->y_node[y_cent]);
            fout << setw(15) << prec;

            fout << setw(14) << scientific << abs(prec - calc);

            fout << fixed << setw(5) << global_i << setw(4) << reg_i;

            // Обработка некраевых узлов
            if(0 < x_cent && x_cent < r->n_x - 1 &&
               0 < y_cent && y_cent < r->n_y - 1)
            {
               fout << "  inner";
            }
            else
            {
               fout << "  border";
            }

            fout << endl;
         }
      }
   }
};