#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

class Matrix
{
private:
	std::vector<std::vector<double>> matrix;

public:
	Matrix(double** A, size_t N, size_t M);
	Matrix(size_t N, size_t M);
	Matrix(std::vector<std::vector<double>> A);

	Matrix(double* A, size_t N);
	Matrix(size_t N);
	Matrix(std::vector<double> A);

	size_t size();

	void Output(bool Absolute = false);

	void Copy(Matrix A);

	double Norm1();

	double Determinant();

	double Max();

	///override operator	
	friend Matrix operator - (Matrix& A, Matrix& B);
	friend Matrix operator + (Matrix& A, Matrix& B);
	friend Matrix operator * (Matrix& A, Matrix& B);

	double& operator() (size_t n, size_t m);
	double& operator() (size_t n);

	double& at(size_t index);

	///end override operator

};
	/*friend Matrix operator += (const Matrix& A, Matrix& B);*/