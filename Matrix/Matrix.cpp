#include "Matrix.h"

template<typename Element>
void Swap(Element& left, Element& right)
{
	Element temp = left;
	left = right;
	right = temp;
}

void swapRow(std::vector<std::vector<double>>& A, int row1, int row2)
{
	size_t size = A.size();

	for (size_t i(0); i < size; i++)
		Swap(A[row1][i], A[row2][i]);
}

void swapColumn(std::vector<std::vector<double>>& A, int column1, int column2)
{
	size_t size = A.size();

	for (size_t i(0); i < size; i++)
		Swap(A[i][column1], A[i][column2]);
}

bool isZeroDiagonal(std::vector<std::vector<double>> A)
{
	size_t size = A.size();

	for (size_t i(0); i < size; i++)
		if (A[i][i] == 0)
			return true;
	return false;
}

size_t Matrix::size()
{
	return matrix.size();
}

Matrix::Matrix(size_t N, size_t M)
{
	std::vector<double> row;
	for (int i(0); i < N; i++)
	{
		for (int j(0); j < M; j++)
			row.push_back(0);

		matrix.push_back(row);
		row.clear();
	}
}

Matrix::Matrix(size_t N)
{
	std::vector<double> row;
	for (int j(0); j < N; j++)
		row.push_back(0);

	matrix.push_back(row);
	row.clear();
}

Matrix::Matrix(double** A, size_t N, size_t M)
{
	for (int i(0); i < N; i++)
	{
		for (int j(0); j < M; j++)
		{
			matrix[i][j] = A[i][j];
		}
	}
}

Matrix::Matrix(double* A, size_t N)
{
	for (int i(0); i < N; i++)
		matrix[i][0] = A[i];
}

Matrix::Matrix(std::vector<std::vector<double>> A)
{
	for (auto row : A)
		matrix.push_back(row);
}

Matrix::Matrix(std::vector<double> A)
{
	matrix.push_back(A);
}

void Matrix::Output(bool Absolute)
{
	for (auto i : matrix)
	{
		for (auto j : i)
			std::cout << std::setw(3) << (Absolute ? abs(j) : j) << " ";

		std::cout << "\n";
	}
}

double Matrix::Norm1()
{
	double S = 0.;

	for (auto i : matrix)
	{
		double St = 0.;
		for (auto j : i)
			St += abs(j);

		if (S < St)
			S = St;
	}

	return S;
}

double& Matrix::operator() (size_t n, size_t m)
{
	return matrix[n][m];
}

double& Matrix::operator() (size_t n)
{
	return matrix[0][n];
}

Matrix operator - (Matrix& A, Matrix& B)
{
	size_t row = A.matrix.size();
	size_t column = A.matrix.at(0).size();

	Matrix ans(row, column);

	for (size_t i(0); i < row; i++)
		for (size_t j(0); j < column; j++)
			ans(i, j) = A(i, j) - B(i, j);

	return ans;
}

Matrix operator + (Matrix& A, Matrix& B)
{
	Matrix ans(A.size(), B.size());

	for (size_t i(0); i < A.matrix.size(); i++)ö
		for (size_t j(0); j < B.matrix.size(); j++)
			ans(i, j) = A(i, j) + B(i, j);

	return ans;
}

Matrix operator * (Matrix& A, Matrix& B)
{
	Matrix ans(A.size(), B.size());

	size_t size = A.size();

	for (size_t i(0); i < size; i++)
		for (size_t j(0); j < size; j++)
			for (size_t l(0); l < size; l++)
				ans(i, j) += A(i, l) * B(l, j);

	return ans;
}


double Matrix::Determinant()
{
	int sign = 1;
	double ans = 1.; // determinant matrix

	size_t N = size(); // size matrix NxN

	std::vector<std::vector<double>> matrix_copy(matrix);

	bool DiagonalZero = isZeroDiagonal(matrix);

	int iter = 0;
	while (DiagonalZero)
	{
		for (int i(0); i < N - 1; i++)
			for (int j(i + 1); j < N; j++)
			{
				swapRow(matrix_copy, i, j);
				sign *= -1;

				DiagonalZero = isZeroDiagonal(matrix_copy);
				if (DiagonalZero)
					break;
			}
	
		if (++iter > 1000)
			return 0.;
	}
	
	for (int i(0); i < N - 1; i++)
	{	
		for (int j(i); j < N - 1; j++)
		{
			double d = -matrix_copy[j + 1][i] / matrix_copy[i][i];

			for(int l(0); l < N; l++)
				matrix_copy[j + 1][l] += matrix_copy[i][l] * d;	
		}
	}

	for (int i(0); i < N; i++)
		ans *= matrix_copy[i][i];

	return ans * sign;
}

void Matrix::Copy(Matrix A)
{
	size_t n = A.size(), m = A.size();

	for (int i(0); i < n; i++)
	{
		for (int j(0); j < n; j++)
			matrix[i][j] = A.matrix[i][j];
	}
}

double Matrix::Max()
{
	double max = matrix[0][0];

	for (auto i : matrix)
	{
		for (auto j : i)
			if (j > max)
				max = j;
	}

	return max;
}

double& Matrix::at(size_t index)
{
	return matrix[0][index];
}