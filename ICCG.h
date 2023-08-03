
#include "pch.h"

int ICCG(const vector<vector<double>> &press_coeff, const vector<double> &data_all, vector<double> &output, int n, int &max_iter, double &eps);
int IncompleteCholeskyDecomp2(const vector<vector<double>> &A, vector<vector<double>> &L, vector<double> &d, int n);
void ICRes(const vector<vector<double>> &L, const vector<double> &d, const vector<double> &r, vector<double> &u, int n);
double dot(vector<double> A, vector<double> B, int n);

int ICCG(const vector<vector<double>> &press_coeff, const vector<double> &data_all, vector<double> &output, int n, int &max_iter, double &eps)
{
	if (n <= 0)
		return 0;

	vector<double> r(n), p(n), y(n), r2(n);
	output.assign(n, 0.0);
	vector<double> d(n);
	vector<vector<double>> L(n, vector<double>(n, 0.0));
	IncompleteCholeskyDecomp2(press_coeff, L, d, n);

	// ��0�ߎ����ɑ΂���c���̌v�Z
	for (int i = 0; i < n; ++i)
	{
		double ax = 0.0;
		for (int j = 0; j < n; ++j)
		{
			ax += press_coeff[i][j] * output[j];
		}
		r[i] = data_all[i] - ax;
	}

	// p_0 = (LDL^T)^-1 r_0 �̌v�Z
	ICRes(L, d, r, p, n);

	double rr0 = dot(r, p, n);
	double rr1;
	double alpha, beta;

	double e = 0.0;
	int k;
	for (k = 0; k < max_iter; ++k)
	{
		// y = AP �̌v�Z
		for (int i = 0; i < n; ++i)
		{
			y[i] = dot(press_coeff[i], p, n);
		}

		// alpha = r*r/(P*AP)�̌v�Z
		alpha = rr0 / dot(p, y, n);

		// ��x�A�c��r�̍X�V
		for (int i = 0; i < n; ++i)
		{
			output[i] += alpha * p[i];
			r[i] -= alpha * y[i];
		}

		// (r*r)_(k+1)�̌v�Z
		ICRes(L, d, r, r2, n);
		rr1 = dot(r, r2, n);

		// �������� (||r||<=eps)
		e = sqrt(rr1);
		if (e < eps)
		{
			k++;
			break;
		}

		// ���̌v�Z��P�̍X�V
		beta = rr1 / rr0;
		for (int i = 0; i < n; ++i)
		{
			p[i] = r2[i] + beta * p[i];
		}

		// (r*r)_(k+1)�����̃X�e�b�v�̂��߂Ɋm�ۂ��Ă���
		rr0 = rr1;
	}

	max_iter = k;
	eps = e;

	return 1;
}

void ICRes(const vector<vector<double>> &L, const vector<double> &d, const vector<double> &r, vector<double> &u, int n)
{
	vector<double> y(n);
	for (int i = 0; i < n; ++i)
	{
		double rly = r[i];
		for (int j = 0; j < i; ++j)
		{
			rly -= L[i][j] * y[j];
		}
		y[i] = rly / L[i][i];
	}

	for (int i = n - 1; i >= 0; --i)
	{
		double lu = 0.0;
		for (int j = i + 1; j < n; ++j)
		{
			lu += L[j][i] * u[j];
		}
		u[i] = y[i] - d[i] * lu;
	}
}

int IncompleteCholeskyDecomp2(const vector<vector<double>> &A, vector<vector<double>> &L, vector<double> &d, int n)
{
	if (n <= 0)
		return 0;

	L[0][0] = A[0][0];
	d[0] = 1.0 / L[0][0];

	for (int i = 1; i < n; ++i)
	{
		for (int j = 0; j <= i; ++j)
		{
			if (fabs(A[i][j]) < 1.0e-10)
				continue;

			double lld = A[i][j];
			for (int k = 0; k < j; ++k)
			{
				lld -= L[i][k] * L[j][k] * d[k];
			}
			L[i][j] = lld;
		}

		d[i] = 1.0 / L[i][i];
	}

	return 1;
}

double dot(vector<double> A, vector<double> B, int n)
{
	int i;
	double sum = 0.0;

	for (i = 0; i < n; ++i)
	{
		sum += A[i] * B[i];
	}

	return sum;
}