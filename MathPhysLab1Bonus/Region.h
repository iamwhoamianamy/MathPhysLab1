#pragma once
#include <vector>

using namespace std;

struct Region
{
   double left, right, top, bot;       // Границы областей

   vector<double> x_node;              // Координаты узлов по X
   vector<double> y_node;              // Координаты узлов по Y

   int n_nodes;                         // Количество узлов

   int n_x;                            // Количество узлов по X
   int n_y;                            // Количество узлов по Y

   int first, last;                    // Индексы первого и последнего
                                       // узлов в глобальной нумерации
   // Массив c информацией о краевых условиях региона
   // 0 - нижнее
   // 1 - правое
   // 2 - верхнее
   // 3 - левое
   vector<int> borders;
};