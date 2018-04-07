#include <vector>
#include "mtrx_functions.h"

int compareLd(const ld x, const ld y) {
	if (x > EPS + y)
		return 1;
	if (fabs(x - y) < EPS)
		return 0;
	if (x < y - EPS)
		return -1;
	return 0;
}
void resizeMatrix(vvld &matrix, const int n, const int m) {
	matrix.resize(n);
	for (int i = 0; i < n; ++i)
		matrix[i].resize(m);
}
void transposeMatrix(const vvld &matrix, vvld &transposedMatrix) {
	for (int i = 0; i < (int)transposedMatrix.size(); ++i) {
		for (int j = 0; j < (int)transposedMatrix[i].size(); ++j)
			transposedMatrix[i][j] = matrix[j][i];
	}
}
void multSquareMatrixs(const vvld &matrix1, const vvld &matrix2, vvld &multMatrix, const int size) {
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			for (int q = 0; q < size; ++q)
				multMatrix[i][j] += matrix1[i][q] * matrix2[q][j];
		}
	}
}
void inverseMatrix(const vvld &matrix, vvld &invMatrix) {
	//Gauss
	//Needs to be fixed
	//Has problems with zeros matrices
	//There will not be zeros matrices
	vvld temp;
	resizeMatrix(temp, (int)matrix.size(), (int)matrix.size() * 2);
	for (int i = 0; i < (int)matrix.size(); ++i) {
		for (int j = 0; j < (int)matrix.size(); ++j) {
			temp[i][j] = matrix[i][j];
		}
		for (int j = (int)matrix.size(); j < (int)matrix.size() * 2; ++j) {
			if (i + (int)matrix.size() == j)
				temp[i][j] = 1.0;
			else
				temp[i][j] = 0.0;
		}
	}

	for (int i = 0; i < (int)matrix.size(); ++i) {
		ld x = temp[i][i];
		if (x == 0) {
			for (int j = 0; j < (int)matrix.size(); ++j) {
				if (j == i)
					continue;
				if (temp[j][i] != 0) {
					for (int q = 0; q < (int)matrix.size() * 2; ++q)
						swap(temp[i][q], temp[j][q]);
					break;
				}
			}
			x = temp[i][i];
			if (x == 0)
				return;
		}
		for (int j = 0; j < (int)matrix.size() * 2; ++j)
			temp[i][j] /= x;
		for (int k = i + 1; k < (int)matrix.size(); ++k) {
			ld y = temp[k][i];
			for (int j = 0; j < (int)matrix.size() * 2; ++j)
				temp[k][j] -= temp[i][j] * y;
		}
	}
	for (int i = (int)matrix.size() - 1; i >= 0; --i) {
		for (int k = i - 1; k >= 0; --k) {
			ld x = temp[k][i];
			for (int j = 0; j < (int)matrix.size() * 2; ++j)
				temp[k][j] -= temp[i][j] * x;
		}
	}
	for (int i = 0; i < (int)matrix.size(); ++i) {
		int t = 0;
		for (int j = (int)matrix.size(); j < (int)matrix.size() * 2; ++j)
			invMatrix[i][t++] = temp[i][j];
	}
}