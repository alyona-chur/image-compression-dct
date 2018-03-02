#include <vector>
#include "quantization.h"
#include "mtrx_functions.h"

void form_matrix_quantos(vvld &matrix_quantos, ld q) {
	for (int i = 0; i < (int)matrix_quantos.size(); ++i) {
		for(int j = 0; j < (int)matrix_quantos[0].size(); ++j)
			matrix_quantos[i][j] = 1.0 + ((1.0 + (ld)i + (ld)j) * q);
	}
}
void quant(const vvld &matrixPostDCT, vvld &quanted_matrix, ld q) {
	vvld matrix_quantos;
	resizeMatrix(matrix_quantos, (int)matrixPostDCT.size(), (int)matrixPostDCT[0].size());
	form_matrix_quantos(matrix_quantos, q);

	for (int i = 0; i < (int)matrixPostDCT.size(); ++i) {
		for (int j = 0; j < (int)matrixPostDCT[i].size(); ++j) {
			quanted_matrix[i][j] = round(matrixPostDCT[i][j] / matrix_quantos[i][j]);
			if (quanted_matrix[i][j] == 0.0) //looks better without -0
				quanted_matrix[i][j] = 0.0;
		}
	}
}
void quant_back(const vvld& quanted_matrix, vvld& matrixPostDCT, ld q) {
	vvld matrix_quantos;
	resizeMatrix(matrix_quantos, (int)matrixPostDCT.size(), (int)matrixPostDCT[0].size());
	form_matrix_quantos(matrix_quantos, q);

	for (int i = 0; i < (int)matrixPostDCT.size(); ++i) {
		for (int j = 0; j < (int)matrixPostDCT[i].size(); ++j) {
			matrixPostDCT[i][j] = matrix_quantos[i][j] * quanted_matrix[i][j];
			if (matrixPostDCT[i][j] == 0.0) //looks better without -0
				matrixPostDCT[i][j] = 0.0;
		}
	}
}
