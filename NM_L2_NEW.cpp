#include <iostream>
#include <vector>
#include <stdio.h>
#include <string>
#include <io.h>
#include <math.h>
#include <time.h>

using namespace std;
FILE* in;
const int d = 5;

void InputVec(vector<double>& vec, int const &beg, int const &end)
{
	for (int i = beg; !feof(in) && i < end; i++)
		fscanf_s(in, "%lf", &vec[i]);
}

int InputMatrix(vector<int>& index, vector<vector<double>>& matrix, int& n, int& m)
{
	if (fopen_s(&in, "matrix.txt", "r") == 0)
	{
		fscanf_s(in, "%d%d", &n, &m);
		index[0] = -(3 + m);
		index[1] = -(2 + m);
		index[2] = -1;
		index[3] = 0;
		index[4] = 1;
		index[5] = 2 + m;
		index[6] = 3 + m;

		for (int i = 0; !feof(in) && i < d; i++)
		{
			matrix[i].resize(n);
			if (index[i] < 0)
				InputVec(matrix[i], -index[i], n);
			else
				InputVec(matrix[i], 0, n - index[i]);
		}
	}
	else
	{
		cout << "file matrix.txt cannot open";
		return 1;
	}
	fclose(in);
	return 0;
}

int InputLimits(double& relax, double& eps, int& maxiter)
{
	if (fopen_s(&in, "limits.txt", "r") == 0)
		fscanf_s(in, "%lf%lf%d", &relax, &eps, &maxiter);
	else
	{
		cout << "file limits.txt cannot open";
		return 1;
	}
	fclose(in);
	return 0;
}

int InputData(vector<double> &X0, vector<double>&F)
{
	int n = F.size();
	if (fopen_s(&in, "vectors.txt", "r") == 0)
	{
		InputVec(X0, 0, n);
		InputVec(F, 0, n);
	}
	else
	{
		cout << "file vectors.txt cannot open";
		return 1;
	}
	fclose(in);
	return 0;
}

void Multiplication( vector<vector<double>> &matrix, 
	 vector<int> &index,  vector<double> &V, vector<double> &res)
{
	int n = V.size(), k = 0;
	for (int i = 0; i < d; i++)
	{
		k = index[i];
		if (k < 0)
			for (int j = abs(k); j < n; j++)
				res[j] += matrix[i][j] * V[k + j];
		else
			for (int j = 0; j < n - k; j++)
				res[j] += matrix[i][j] * V[k + j];
	}
}

double VecNorm( vector<double> &vec)
{
	double sum = 0.;
	for (int i = 0; i < vec.size(); i++)
		sum += vec[i] * vec[i];
	return sqrt(sum);
}

double RelativeResidual( vector<double> &F,  vector<vector<double>> &matrix, 
	 vector<int> &index,  vector<double> &X)
{
	int n = F.size();
	vector<double> mult(n);

	Multiplication(matrix, index, X, mult);
	for (int i = 0; i < n; i++)
		mult[i] = F[i] - mult[i];

	return VecNorm(mult) / VecNorm(F);
}

void IterativeProcess(vector<double>& Xk, vector<double>& Xk1, vector<int>& index,
	vector<vector<double>> &matrix, int &j, double &sum)
{
	int k = 0, n = Xk.size();
	for (int i = 0; i < d; i++)
	{
		k = index[i];
		if (k + j >= 0 && k + j < n)
		{
			if (i < 3) // нижний треугольник
				sum += matrix[i][j] * Xk1[k + j];
			else // верхний треугольник
				sum += matrix[i][j] * Xk[k + j];
		}
	}
}

void Jacobi(vector<double> &Xk, vector<double> &Xk1, vector<int> &index, vector<vector<double>> &matrix,
	vector<double> &F, int &maxiter, double &eps, double &relax)
{
	int n = F.size();
	double r = 0., sum = 0.;
	r = RelativeResidual(F, matrix, index, Xk);
	for (int k = 0; k < maxiter && r > eps; k++)
	{
		for (int j = 0; j < n; j++)
		{
			IterativeProcess(Xk, Xk, index, matrix, j, sum);
			Xk1[j] = Xk[j] + (relax / matrix[3][j]) * (F[j] - sum);
			sum = 0.;
		}
		Xk.swap(Xk1);
		r = RelativeResidual(F, matrix, index, Xk);
		cout << k << " " << r << endl;
	}
}

void GaussZeidel(vector<double>& Xk, vector<double>& Xk1, vector<int>& index, vector<vector<double>>& matrix,
	vector<double>& F, int& maxiter, double& eps, double& relax)
{
	int n = F.size();
	double r = 0., sum = 0.;
	r = RelativeResidual(F, matrix, index, Xk);
	for (int k = 0; k < maxiter && r > eps; k++)
	{
		for (int j = 0; j < n; j++)
		{
			IterativeProcess(Xk, Xk1, index, matrix, j, sum);
			Xk1[j] = Xk[j] + (relax / matrix[3][j]) * (F[j] - sum);
			sum = 0.;
		}
		Xk.swap(Xk1);
		r = RelativeResidual(F, matrix, index, Xk);
		cout << k << " " << r << endl;
	}
}

int Output(vector<double> res)
{
	if (fopen_s(&in, "out.txt", "w") == 0)
		for (int i = 0; i < res.size(); i++)
			fprintf_s(in, "%.14e \n", res[i]);
	else
	{
		cout << "file out.txt cannot open";
		return 1;
	}
	return 0;
}

int main()
{
	int n = 0, maxiter = 0, m = 0;
	double relax = 0., eps = 0.;

	vector<vector<double>> matrix(d);
	vector<int> index(d);
	InputMatrix(index, matrix, n, m);
	InputLimits(relax, eps, maxiter);
	vector<double> Xk(n), Xk1(n), F(n);
	InputData(Xk, F);
	//Jacobi(Xk, Xk1, index, matrix, F, maxiter, eps, relax);
	clock_t start = clock();
	GaussZeidel(Xk, Xk1, index, matrix, F, maxiter, eps, relax);
	clock_t end = clock();
	cout << "Time: " << (double)(end - start) * 1000 / CLOCKS_PER_SEC << " msec\n";

	Output(Xk);
	return 0;
}