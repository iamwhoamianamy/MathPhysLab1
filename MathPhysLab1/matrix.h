#pragma once
#include <vector>

using namespace std;

class matrix
{
public:
   int N;               // ������ �������
   int M;               // ���������� ��������� � ������������

   vector<int> ig;      // ��������� ������ �����
   vector<int> jg;      // ������ �������� ��������������� ���������

   vector<double> ggl;    // ������� �����������
   vector<double> ggu;    // ������ �����������

   vector<double> di;     // ���������

   // ������������� ������� � ������� �����������
   matrix(int _N)
   {
      N = _N;
      M = N * (N - 1) / 2;

      ggl.resize(M);
      ggu.resize(M);
      jg.resize(M);
      di.resize(N);
      ig.resize(N + 1);

      ig[0] = ig[1] = 0;

      di[0] = 1.0;

      for(int i = 1; i < N; i++)
      {
         int i0 = ig[i + 0];
         int i1 = ig[i + 1] = ig[i + 0] + i;

         for(int j = 0, k = i0; j < i; j++, k++)
            jg[k] = j;
      }
   }
};