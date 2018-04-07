#include <vector>
#include "mtrx_functions.h"
#include "dct.h"

void formDCTmatrix(vvld &DCTmatrix) {
	for (int i = 0; i < DCT_BLOCK_SIZE; ++i) {
		for (int j = 0; j < DCT_BLOCK_SIZE; ++j) {
			if (i == 0)
				DCTmatrix[i][j] = 1.0 / sqrt((ld)DCT_BLOCK_SIZE);
			else
				DCTmatrix[i][j] = sqrt(2.0 / (ld)DCT_BLOCK_SIZE) * cos((2.0 * (ld)j + 1.0) * (ld)i * PI / 2.0 / (ld)DCT_BLOCK_SIZE);
		}
	}
}
void doBlockDCT(const vvld &block, vvld &DCTblock, const vvld &DCTmatrix, const vvld &transposedDCT) {
	vvld temp;
	resizeMatrix(temp, DCT_BLOCK_SIZE, DCT_BLOCK_SIZE);
	multSquareMatrixs(DCTmatrix, block, temp, DCT_BLOCK_SIZE);
	multSquareMatrixs(temp, transposedDCT, DCTblock, DCT_BLOCK_SIZE);
}
void doBlockDCTback(vvld &block, const vvld &DCTblock, const vvld &DCTmatrix, const vvld &transposedDCT) {
	vvld invDCTmatrix;
	resizeMatrix(invDCTmatrix, DCT_BLOCK_SIZE, DCT_BLOCK_SIZE);
	inverseMatrix(DCTmatrix, invDCTmatrix);
	vvld invTransposedDCT;
	resizeMatrix(invTransposedDCT, DCT_BLOCK_SIZE, DCT_BLOCK_SIZE);
	inverseMatrix(transposedDCT, invTransposedDCT);

	vvld temp;
	resizeMatrix(temp, DCT_BLOCK_SIZE, DCT_BLOCK_SIZE);
	multSquareMatrixs(invDCTmatrix, DCTblock, temp, DCT_BLOCK_SIZE);
	multSquareMatrixs(temp, invTransposedDCT, block, DCT_BLOCK_SIZE);
}

void DCT(const vvld &matrix, vvld &postDCT) {
	//postDCT = matrix;

	vvld DCTmatrix;
	resizeMatrix(DCTmatrix, DCT_BLOCK_SIZE, DCT_BLOCK_SIZE);
	formDCTmatrix(DCTmatrix);

	vvld transposedDCT;
	resizeMatrix(transposedDCT, DCT_BLOCK_SIZE, DCT_BLOCK_SIZE);
	transposeMatrix(DCTmatrix, transposedDCT);

	for (int st_i = 0; st_i + DCT_BLOCK_SIZE < (int)matrix.size(); st_i += DCT_BLOCK_SIZE) {
		for (int st_j = 0; st_j + DCT_BLOCK_SIZE < (int)matrix[st_i].size(); st_j += DCT_BLOCK_SIZE) {
			vvld block;
			resizeMatrix(block, DCT_BLOCK_SIZE, DCT_BLOCK_SIZE);
			int q = 0, w = 0;

			for (int i = st_i; i < st_i + DCT_BLOCK_SIZE; ++i) {
				for (int j = st_j; j < st_j + DCT_BLOCK_SIZE; ++j)
					block[q][w++] = matrix[i][j] - 128;
				q++; w = 0;
			}

			vvld DCTblock;
			resizeMatrix(DCTblock, DCT_BLOCK_SIZE, DCT_BLOCK_SIZE);
			doBlockDCT(block, DCTblock, DCTmatrix, transposedDCT);
			q = 0; w = 0;

			for (int i = st_i; i < st_i + DCT_BLOCK_SIZE; ++i) {
				for (int j = st_j; j < st_j + DCT_BLOCK_SIZE; ++j)
					postDCT[i][j] = DCTblock[q][w++];
				q++; w = 0;
			}
		}
	}
}
void DCT_back(const vvld &postDCT, vvld &matrix) {
	vvld DCTmatrix;
	resizeMatrix(DCTmatrix, DCT_BLOCK_SIZE, DCT_BLOCK_SIZE);
	formDCTmatrix(DCTmatrix);

	vvld transposedDCT;
	resizeMatrix(transposedDCT, DCT_BLOCK_SIZE, DCT_BLOCK_SIZE);
	transposeMatrix(DCTmatrix, transposedDCT);

	for (int st_i = 0; st_i + DCT_BLOCK_SIZE < (int)postDCT.size(); st_i += DCT_BLOCK_SIZE) {
		for (int st_j = 0; st_j + DCT_BLOCK_SIZE < (int)postDCT[st_i].size(); st_j += DCT_BLOCK_SIZE) {
			vvld DCTblock;
			resizeMatrix(DCTblock, DCT_BLOCK_SIZE, DCT_BLOCK_SIZE);
			int q = 0, w = 0;

			for (int i = st_i; i < st_i + DCT_BLOCK_SIZE; ++i) {
				for (int j = st_j; j < st_j + DCT_BLOCK_SIZE; ++j)
					DCTblock[q][w++] = postDCT[i][j];
				q++; w = 0;
			}

			vvld block;
			resizeMatrix(block, DCT_BLOCK_SIZE, DCT_BLOCK_SIZE);
			doBlockDCTback(block, DCTblock, DCTmatrix, transposedDCT);
			q = 0; w = 0;

			for (int i = st_i; i < st_i + DCT_BLOCK_SIZE; ++i) {
				for (int j = st_j; j < st_j + DCT_BLOCK_SIZE; ++j)
					matrix[i][j] = block[q][w++] + 128;
				q++; w = 0;
			}
		}
	}
}

int tr(int x) {
	return x + x;
}