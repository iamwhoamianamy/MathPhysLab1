#pragma once
#include <vector>
#include <fstream>
#include <string>
#include "matrix.h"

using namespace std;

class EllipticalProblem
{
public:
   vector<double> Xr;            // ������� ������� �� X
   vector<double> Yr;            // ������� ������� �� Y

   vector<double> Xn;            // ���������� ����� �� X
   vector<double> Yn;            // ���������� ����� �� Y

   matrix* MainMatrix;           // ������� �������

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

   // ������� ���������� ������ ������� �� ����� FILE_NAME,
   // ��������� ��������� �����
   void read_coord_lines(const string& FILE_NAME)
   {
      std::ifstream fin(FILE_NAME);

      // ��������� ������� �������
      for(int i = 0; i < 3; i++)
         fin >> Xr[i];

      for(int i = 0; i < 3; i++)
         fin >> Yr[i];

      // ��������� ��������� ������ �� X
      int count = 0;        // ����� ����� �� X
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

      // ��������� ��������� ������ �� Y
      count = 0;        // ����� ����� �� Y

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