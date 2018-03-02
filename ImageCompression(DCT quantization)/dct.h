#pragma once
#include <vector>

using namespace std;
typedef long double ld;
typedef vector < vector <ld> > vvld;
const ld PI = (ld)3.141592653589793238462643;
const int DCT_BLOCK_SIZE = 8;

void formDCTmatrix(vvld& DCTmatrix);
void doBlockDCT(const vvld& block, vvld& DCTblock, const vvld& DCTmatrix, const vvld& transposedDCT);
void doBlockDCTback(vvld& block, const vvld& DCTblock, const vvld& DCTmatrix, const vvld& transposedDCT);

void DCT(const vvld &matrix, vvld &postDCT);
void DCT_back(const vvld &postDCT, vvld &matrix);