#pragma once
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;

typedef double real;

// ��������� ������� �� �����
vector<real> operator * (real val, const vector<real>& vec)
{
   vector<real> res(vec.size());

   for (size_t i = 0; i < vec.size(); ++i)
      res[i] = val * vec[i];
   return res;
}

// �������� ��������
vector<real> operator + (const vector<real>& vec1, const vector<real>& vec2)
{
   if (vec1.size() != vec2.size())
      throw("a.size() != b.size()");

   vector<real> res(vec1.size());

   for (size_t i = 0; i < vec1.size(); ++i)
      res[i] = vec1[i] + vec2[i];
   return res;
}

// ��������� ��������
vector<real> operator - (const vector<real>& vec1, const vector<real>& vec2)
{
   if (vec1.size() != vec2.size())
     throw("a.size() != b.size()");

   vector<real> res(vec1.size());

   for (size_t i = 0; i < vec1.size(); ++i)
      res[i] = vec1[i] - vec2[i];
   return res;
}

// ��������� ������������ ��������
real operator *(const vector<real>& vec1, const vector<real>& vec2)
{
   if (vec1.size() != vec2.size())
      throw("vec1.size() != vec2.size()");

   int n = vec1.size();
   real res = 0;

   for (int i = 0; i < n; i++)
      res += vec1[i] * vec2[i];

   return res;
}

// ����� �������
real Norm(const vector<real>& vec)
{
   return sqrt(vec * vec);
}