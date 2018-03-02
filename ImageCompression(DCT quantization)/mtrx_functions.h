#pragma once
#include <vector>

using namespace std;
typedef long double ld;
typedef vector < vector <ld> > vvld;

const ld EPS = (ld)1E-18;

int compareLd(const ld x, const ld y);
void resizeMatrix(vvld &matrix, const int n, const int m);
void transposeMatrix(const vvld &matrix, vvld &transposedMatrix);
void multSquareMatrixs(const vvld &matrix1, const vvld &matrix2, vvld &multMatrix, const int size);
void inverseMatrix(const vvld &matrix, vvld &invMatrix);