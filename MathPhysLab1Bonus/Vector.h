#pragma once
#include <vector>
#include <iomanip>
#include <fstream>
using namespace std;


// Умножение вектора на число
vector<double> operator * (double val, const vector<double>& vec)
{
   size_t n = vec.size();
   vector<double> res(n);
   
   for (size_t i = 0; i < n; ++i)
      res[i] = val * vec[i];
   return res;
}

// Сложение векторов
vector<double> operator + (const vector<double>& vec1, const vector<double>& vec2)
{
   size_t n = vec1.size();
   if (n != vec2.size())
      throw("a.size() != b.size()");

   vector<double> res(n);

   for (size_t i = 0; i < n; ++i)
      res[i] = vec1[i] + vec2[i];
   return res;
}

// Вычитание векторов
vector<double> operator - (const vector<double>& vec1, const vector<double>& vec2)
{
   size_t n = vec1.size();
   if (n != vec2.size())
     throw("a.size() != b.size()");

   vector<double> res(n);

   for (size_t i = 0; i < n; ++i)
      res[i] = vec1[i] - vec2[i];
   return res;
}

// Скалярное произведение векторов
double operator *(const vector<double>& vec1, const vector<double>& vec2)
{
   size_t n = vec1.size();
   if (n != vec2.size())
      throw("vec1.size() != vec2.size()");

   double res = 0;
   for (int i = 0; i < n; i++)
      res += vec1[i] * vec2[i];
   return res;
}

// Норма вектора
double Norm(const vector<double>& vec)
{
   return sqrt(vec * vec);
}