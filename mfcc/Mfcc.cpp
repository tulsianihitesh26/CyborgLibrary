/**
 * Mfcc.cpp
 *
 * Calculation of MFCC signal features.
 *
 */

#include "Mfcc.h"

Mfcc::~Mfcc() {
}

void Mfcc::mfccInit(
		int inputSize,
		int Fs,
		Number fLow,
		Number fHigh,
		int numFilt,
		int numFeatures,
		int nfft,
		int dctType) {
	m_inputSize = inputSize;
	m_pow2Size = Dsp::nextPow2(nfft);
	m_samplingFreq = Fs;
	m_fftSize = nfft;
	m_numFeatures = numFeatures;
	m_bank.melFilterBankInit(m_samplingFreq, fLow, fHigh, numFilt, Dsp::round(((Number) m_fftSize / (Number)2.0) + 1));
	m_dct.dctInit(numFilt, m_numFeatures, dctType);

	int nfftReqd = (1 << Dsp::nextPow2(m_inputSize));
	if (nfftReqd > nfft) {
		cout << "ERROR: nfft set is less than the required nfft of " << nfftReqd << endl;
		exit(0);
	}
}

/**
 * Calculates MFCC given audio frame
 */
NumberArray Mfcc::calculate(
		NumberArrayRef audioFrame) {
	//------------calculate FFT-spectrum------------------//
	NumberArray realPart = audioFrame;
	realPart.resize(m_fftSize, 0.0); // zero padding

	NumberArray imagPart(m_fftSize, 0.0);
	Dsp::fft(Dsp::FORWARD_FFT, m_pow2Size, m_fftSize, realPart, imagPart);

	//---------Calculate LOG of( Mel-spectrum)------------------//
	NumberArray filterOutput = m_bank.applyFilter(realPart, imagPart);

	//---------Calculate Mel-cepstrum------------------//
	return m_dct.applyDct(filterOutput);
}

