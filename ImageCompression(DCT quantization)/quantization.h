#pragma once
#include <vector>

using namespace std;
typedef long double ld;
typedef vector < vector <ld> > vvld;

void quant(const vvld &matrixPostDCT, vvld &quanted_matrix, ld q);
void quant_back(const vvld& quanted_matrix, vvld& matrixPostDCT, ld q);