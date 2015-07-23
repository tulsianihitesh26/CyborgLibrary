/*
 * Dsp.h
 *
 *  Created on: 20-Mar-2013
 *      Author: sujeet
 */

#ifndef DSP_H_
#define DSP_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <iterator>
#include "../cyborg/commonUtils.h"

using namespace std;
class Dsp {
public:

	const static Number INVALID_CENT;
	const static int SAMPLING_RATE = 1024;
	const static Number PI;
	Dsp();
	virtual ~Dsp();

	enum FFT_DIRECTION {
		FORWARD_FFT = 1, INVERSE_FFT = -1
	};

	static void zeroMeanUnitVarianceNormalization(
			Number2DArrayRef oldMFCCMatrix);
	static void linspace(
			Number a,
			Number b,
			NumberArrayRef array);
	static bool levinson(
			NumberArrayRef t,
			NumberArrayRef b,
			NumberArrayRef x);
	static void changePitchCentOctave(
			NumberArrayRef userPcopy,
			int octaveShift);
	static void autoCorrCoeffs(
			NumberArrayConstRef inputSamples,
			NumberArrayRef acfCoeffs);
	static void autoCorrMatrix(
			NumberArrayConstRef adaptationRegionRefSamples,
			Number2DArrayRef acfCoeffs2D);
	static void crossCorrCoeffs(
			NumberArrayConstRef referenceInput,
			NumberArrayConstRef userInput,
			NumberArrayRef ccfCoeffs);
	static void crossCorrMatrix(
			NumberArrayConstRef adaptationRegionRefSamples,
			NumberArrayConstRef adaptationRegionUserSamples,
			Number2DArrayRef ccfCoeffs);
	static void dtwPath(
			PitchesConstRef refMIDINotesInterpolated,
			PitchesConstRef userP,
			unsigned int padVectorSizeForDTW,
			vector<int>& path,
			int constraintLengthFrames);
	static void dtwPath(
			PitchesConstRef refMIDINotesInterpolated,
			PitchesConstRef userP,
			unsigned int padVectorSizeForDTW,
			vector<int>& path);
	static Number computeDCOfVibrato(
			NumberArrayConstRef pdV);
	static Number sum(
			NumberArrayConstRef esd,
			int startIdx,
			int endIdx);
	static int sum(
			IntArrayRef esd,
			int startIdx,
			int endIdx);
	static Number energy(
			NumberArrayRef userFile,
			int startSample,
			int endSample);
	static Number max(
			NumberArrayConstRef v,
			int startIdx,
			int endIdx,
			unsigned int * maxIdx);
	static Number median(
			NumberArrayConstRef pdV,
			Number dMedianFactor);
	static Number sumSquareSub(
			NumberArrayConstRef data,
			Number value);
	static Number sumSquareSubm(
			NumberArrayConstRef data,
			NumberArrayConstRef data2);
	static void hanning(
			NumberArrayRef hannWin);
	static void hamming(
			NumberArrayRef hammingWin);
	static Number sinc(
			Number x);
	static int gcd(
			int a,
			int b);
	static void sincResample(
			NumberArrayRef userFile,
			Number shift,
			int newFreq);
	static void esd(
			NumberArrayRef vEsd,
			NumberArrayConstRef vXr,
			NumberArrayConstRef vXi);
	static void magnitude(
			NumberArrayRef vEsd,
			NumberArrayConstRef vXr,
			NumberArrayConstRef vXi);
	static void esdFFTS(
			NumberArrayRef vEsd,
			NumberArrayConstRef vX);
	static int nextPow2(
			int length);
	static void fft(
			FFT_DIRECTION dir,
			unsigned const int m,
			unsigned const int n,
			NumberArrayRef xr,
			NumberArrayRef xi);
	static void scale(
			NumberArrayRef vEsd,
			Number scale);
	static void divideAByB_thresholded(
			NumberArrayRef vResult,
			vector<Number> const & vA,
			vector<Number> const & vB,
			Number threshold);
	static void autoCorr(
			NumberArrayRef currWin,
			int order,
			NumberArrayRef aCorrPts);
	static void filter(
			int ord,
			NumberArrayRef a,
			NumberArrayRef b,
			NumberArrayRef x,
			NumberArrayRef y);
	static void prEmpFilter(
			NumberArrayRef x,
			NumberArrayRef y);
	static void conv(
			NumberArrayRef x,
			int L,
			NumberArrayRef h);
	static void convolveAdaptiveFilter(
			NumberArrayRef x,
			NumberArrayRef h);
	static void firResample(
			NumberArrayConstRef coeffs,
			NumberArrayRef input,
			NumberArrayRef output,
			int downFactor);
	static void fftConv(
			NumberArrayRef aud,
			NumberArrayConstRef rir,
			NumberArrayRef acousticVec);
	static void fftConvFFTvec(
			NumberArrayRef aud,
			NumberArrayConstRef rirFFT,
			NumberArrayConstRef rirFFTI,
			NumberArrayRef acousticVec);
	static void getLpcCoeffs(
			NumberArrayRef acfCoeffs,
			NumberArrayRef lpCoeffs,
			int lpOrder);
	static void warpBpole(
			Number alpha,
			NumberArrayRef Bcoeff);
	static void warpApole(
			NumberArrayRef a,
			Number alpha,
			NumberArrayRef Acoeff);

	static void addEcho(
			NumberArrayRef awesomeVec,
			NumberArrayRef filteredVec);
	static void addRingMod(
			NumberArrayRef userAudVec,
			NumberArrayRef awesomeVec);
	static void addFlanger(
			NumberArrayRef userAudVec,
			NumberArrayRef awesomeVec);
	static void changeFormant(
			NumberArrayRef userAudio,
			Number alpha,
			NumberArrayRef audioOut);
	static int round(Number val);

};

#endif /* DSP_H_ */
