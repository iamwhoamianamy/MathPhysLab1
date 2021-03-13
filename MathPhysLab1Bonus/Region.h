#pragma once
#include <vector>

using namespace std;

struct Region
{
   double left, right, top, bot;       // ������� ��������

   vector<double> x_node;              // ���������� ����� �� X
   vector<double> y_node;              // ���������� ����� �� Y

   int n_nodes;                         // ���������� �����

   int n_x;                            // ���������� ����� �� X
   int n_y;                            // ���������� ����� �� Y

   int first, last;                    // ������� ������� � ����������
                                       // ����� � ���������� ���������
   // ������ c ����������� � ������� �������� �������
   // 0 - ������
   // 1 - ������
   // 2 - �������
   // 3 - �����
   vector<int> borders;
};