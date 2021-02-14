#pragma once
#include <vector>
#include <fstream>
#include <string>
#include "matrix.h"

using namespace std;

class EllipticalProblem
{
public:
   vector<double> Xr;            // Границы области по X
   vector<double> Yr;            // Границы области по Y

   vector<double> Xn;            // Координаты узлов по X
   vector<double> Yn;            // Координаты узлов по Y

   matrix* MainMatrix;           // Матрица системы

   EllipticalProblem()
   {
      Xr = vector<double>(3);
      Yr = vector<double>(3);

      read_coord_lines("coords.txt");

      MainMatrix = new matrix(pow(Xn.size() + Yn.size(), 2));
   }

   ~EllipticalProblem()
   {
      delete MainMatrix;
   }

   // Функция считывания границ области из файла FILE_NAME,
   // генерация координат узлов
   void read_coord_lines(const string& FILE_NAME)
   {
      std::ifstream fin(FILE_NAME);

      // Считываем границы области
      for(int i = 0; i < 3; i++)
         fin >> Xr[i];

      for(int i = 0; i < 3; i++)
         fin >> Yr[i];

      // Генерация координат уззлов по X
      int count = 0;        // Число узлов по X
      int t1, t2;

      fin >> t1;
      count += t1;

      fin >> t2;
      count += t2;

      Xn.resize(count + 1);

      for(int i = 0; i <= t1; i++)
         Xn[i] = Xr[0] + (Xr[1] - Xr[0]) / t1 * i;

      for(int i = 0; i <= t2; i++)
         Xn[i + t1] = Xr[1] + (Xr[2] - Xr[1]) / t2 * i;

      // Генерация координат уззлов по Y
      count = 0;        // Число узлов по Y

      fin >> t1;
      count += t1;

      fin >> t2;
      count += t2;

      Yn.resize(count + 1);

      for(int i = 0; i <= t1; i++)
         Yn[i] = Yr[0] + (Yr[1] - Yr[0]) / t1 * i;

      for(int i = 0; i <= t2; i++)
         Yn[i + t1] = Yr[1] + (Yr[2] - Yr[1]) / t2 * i;
   }
};