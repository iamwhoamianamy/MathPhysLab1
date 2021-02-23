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
   vector<double> Xr;                  // Границы области по X
   vector<double> Yr;                  // Границы области по Y

   vector<double> Xn;                  // Координаты узлов по X
   vector<double> Yn;                  // Координаты узлов по Y

   SLAE* slae;                         // Система

   EllipticalProblem()
   {
      Xr = vector<double>(3);
      Yr = vector<double>(3);

      ReadCoordLines("coords.txt");

      slae = new SLAE(Xn.size() * Yn.size(), Xn.size());
   }

   ~EllipticalProblem()
   {
      delete slae;
   }

   // Функция считывания границ области из файла FILE_NAME,
   // генерация координат узлов
   void ReadCoordLines(const string& FILE_NAME)
   {
      std::ifstream fin(FILE_NAME);

      // Считываем границы области
      for(int i = 0; i < 3; i++)
         fin >> Xr[i];

      for(int i = 0; i < 3; i++)
         fin >> Yr[i];

      // Генерация координат узлов по X
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

      // Генерация координат узлов по Y
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

   // Формирование матрицы системы
   void FormMatrix()
   {
      for(int n = 0; n < slae->N; n++)
      {
         
      }
   }
};