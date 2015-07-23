#ifndef MFCC_H
#define MFCC_H

#include "Dsp.h"
#include "MelFilterBank.h"
#include "Dct.h"
#include <cmath>
#include <math.h>
#include <cstddef>
#include <vector>

using namespace std;

class Mfcc {
public:
	Mfcc() {
	}
	;
	void mfccInit(
			int inputSize,
			int Fs = 8000,
			Number fLow = 300,
			Number fHigh = 3500,
			int numFilt = 31,
			int numFeatures = 13,
			int nfft = 512,
			int dctType = 0);
	~Mfcc();
	NumberArray calculate(
			NumberArrayRef source);

private:
	int m_inputSize;
	int m_pow2Size;
	int m_samplingFreq;
	int m_fftSize;
	int m_numFeatures;
	MelFilterBank m_bank;
	int m_dctType;
	Dct m_dct;

};

#endif // MFCC_H
