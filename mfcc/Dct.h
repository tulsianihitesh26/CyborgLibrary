#ifndef DCT_H
#define DCT_H
#include <iostream>
#include "Dsp.h"
using namespace std;
class Dct {
public:
	Dct() {
	}
	;
	void dctInit(
			int numberFilters,
			int numberCoefficients,
			int dctType);
	~Dct() {
	}
	;
	NumberArray applyDct(
			NumberArrayRef data);
	NumberArray applyDct_legacy(
			NumberArrayRef data);

private:

	Number2DArray m_dctBank;
	double m_sqrt_inv_n, m_sqrt_inv_2n;
	int m_dctType;
};

#endif // DCT_H
