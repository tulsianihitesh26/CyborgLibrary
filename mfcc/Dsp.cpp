/*
 * Dsp.cpp
 *
 *  Created on: 20-Mar-2013
 *      Author: sujeet
 */

#include "Dsp.h"
#include "../cyborg/NeuralNetwork.h"
#include <math.h>
#include <algorithm>
#include <numeric>

const Number Dsp::INVALID_CENT = -9999.0;
const Number Dsp::PI = 3.14159265359;


Dsp::Dsp() {
}

void Dsp::zeroMeanUnitVarianceNormalization(
		Number2DArrayRef oldMFCCMatrix) {
	NumberArray meanVec, stdVec;

	NeuralNetwork::matrixMean(oldMFCCMatrix, 1, meanVec);
	NeuralNetwork::matrixStd(oldMFCCMatrix, 1, stdVec);
	for (int i = 0; i < (int) oldMFCCMatrix.size(); i++) {
		for (int j = 0; j < (int) oldMFCCMatrix[i].size(); j++) {
			oldMFCCMatrix[i][j] = (oldMFCCMatrix[i][j] - meanVec[j]) / stdVec[j];
		}
	}

}
void Dsp::firResample(
		NumberArrayConstRef coeffs,
		NumberArrayRef input,
		NumberArrayRef output,
		int downFactor) {

	NumberArray insamp(input.size() + coeffs.size() - 1, 0.0);
	int offset = coeffs.size() - 1;
	// put the new samples at the high end of the buffer
	for (int i = 0; i < (int) input.size(); i++) {
		insamp[offset + i] = input[i];
	}

	// apply the filter to each input sample
	for (int n = ((int) coeffs.size()) >> 1, m = 0; n < (int) input.size(); n = n + downFactor, m++) // filterLength>>1 does the filter delay compensation
			{
		Number acc = 0;
		for (int k = 0; k < (int) coeffs.size(); k++) {
			acc += coeffs[k] * insamp[offset + n - k];
		}
		output.push_back(acc);
	}
}

void Dsp::linspace(
		Number a,
		Number b,
		NumberArrayRef array) {
	int n = array.size();
	Number step = (b - a) / (n - 1);
	int i;
	for (i = 0; i < n; i++) {
		array[i] = a;
		a += step;
	}
}
//--------------------------------------------------------------------------//
// Implementation for Levinson-Durbin algorithm. Zhang Ming, 2010-11, Xi'an Jiaotong University. http://m.oschina.net/blog/8515
//--------------------------------------------------------------------------//
// To solve tx=b equation where t is a Toeplitz matrix
// input :
// 		t    : t(0), t(1), ..., t(n-1) of Toeplitz coefficient matrix
// 		b    : constant vector
// output :
// 		x    : unknown coeffs
//		boolean indicating Solution was reached
bool Dsp::levinson(
		NumberArrayRef t,
		NumberArrayRef b,
		NumberArrayRef x) {
	int n = t.size();
	Number alpha, beta, q, c, omega;
	NumberArray y(n, 0.0);
	NumberArray yy(n, 0.0);
	//NumberArray x(n,0.0);

	alpha = t[0];
	if (alpha == 0) {
		cout << "Levinson: Alpha Value Zero. The matrix is ill-conditioned!" << endl;
		return false;
	}
	y[0] = 1;
	x[0] = b[0] / alpha;

	for (int k = 1; k < n; ++k) {
		q = 0;
		beta = 0;
		for (int j = 0; j < k; ++j) {
			q += x[j] * t[k - j];
			beta += y[j] * t[j + 1];
		}
		c = -beta / alpha;

		yy[0] = c * y[k - 1];
		y[k] = y[k - 1];
		for (int i = 1; i < k; ++i) {
			yy[i] = y[i - 1] + c * y[k - i - 1];
		}
		yy[k] = y[k - 1];

		alpha += c * beta;
		if (alpha == 0) {
			cout << "Levinson: Alpha Value Zero. The matrix is ill-conditioned!" << endl;
			return false;
		}

		omega = (b[k] - q) / alpha;
		for (int i = 0; i < k; ++i) {
			x[i] += omega * yy[i];
			y[i] = yy[i];
		}
		x[k] = omega * y[k];
	}
	return true;
}

//--------------------------------------------------------------------------//
// Function to shift pitch in Cents by different octave
//--------------------------------------------------------------------------//
// input :
// 	userPcopy: userPitch Vector in Cents
// output :
// 	octaveShift - Octave Shift Factor
//--------------------------------------------------------------------------//
void Dsp::changePitchCentOctave(
		NumberArrayRef userPcopy,
		int octaveShift) {
	for (unsigned int i = 0; i < userPcopy.size(); i++) {
		userPcopy[i] = (userPcopy[i] != INVALID_CENT) ? userPcopy[i] + (1200 * octaveShift) : INVALID_CENT;
	}
}

//--------------------------------------------------------------------------//
// Function to calculate autocorrelation coefficients of a signal
//--------------------------------------------------------------------------//
// input :
// 	inputSamples: input sequence for which Autocorrelation is to be computed
// output :
// 	acfCoeffs - Autocorrelation Coefficient
//--------------------------------------------------------------------------//

void Dsp::autoCorrCoeffs(
		NumberArrayConstRef inputSamples,
		NumberArrayRef acfCoeffs) {

	unsigned int lengthOfAdaptationData = inputSamples.size(); // Length of Signal
	unsigned int lengthOfFilter = acfCoeffs.size(); // Autocorrelation Matrix Size
	Number sum = 0;

	// Create the autocorrelation vector first column
	for (unsigned int i = 0; i < lengthOfFilter; i++) {
		for (unsigned int j = i; j < lengthOfAdaptationData; j++) {
			sum += (inputSamples[j] * inputSamples[j - i]);
		}
		acfCoeffs[i] = sum / (Number) lengthOfAdaptationData; //normalization

		sum = 0;
	}
}

//--------------------------------------------------------------------------//
// Function to calculate autocorrelation matrix of a signal
//--------------------------------------------------------------------------//
// input :
// 	inputSamples: input sequence for which Autocorrelation is to be computed
// output :
// 	acfCoeffs - Autocorrelation Coefficients Matrix
//--------------------------------------------------------------------------//

void Dsp::autoCorrMatrix(
		NumberArrayConstRef inputSamples,
		Number2DArrayRef acfCoeffs) {

	unsigned int lengthOfAdaptationData = inputSamples.size(); // Length of Signal
	unsigned int lengthOfFilter = acfCoeffs.size(); // Autocorrelation Matrix Size
	Number sum = 0;

	// Create the autocorrelation vector first column
	for (unsigned int i = 0; i < lengthOfFilter; i++) {
		for (unsigned int j = i; j < lengthOfAdaptationData; j++) {
			sum += (inputSamples[j] * inputSamples[j - i]);
		}
		acfCoeffs[0][i] = sum / (Number) lengthOfAdaptationData; //normalization

		sum = 0;
	}

	// Create the toeplitz matrix from the vector
	for (unsigned int i = 1; i < lengthOfFilter; i++) //First column is already created so start from second
			{
		for (unsigned int j = 0; j < lengthOfFilter; j++) {
			acfCoeffs[i][j] = acfCoeffs[0][abs((int) i - (int) j)];
		}
	}

}

//--------------------------------------------------------------------------//
// Function to calculate cross-correlation coefficients of two signals reference and user input
//--------------------------------------------------------------------------//
// input :
// 	referenceInput: reference input sequence for which cross-correlation is to be computed
// 	userInput: user input sequence for which cross-correlation is to be computed
// output :
// 	ccfCoeffs - Cross-correlation Coefficients Coeffs
//--------------------------------------------------------------------------//

void Dsp::crossCorrCoeffs(
		NumberArrayConstRef referenceInput,
		NumberArrayConstRef userInput,
		NumberArrayRef ccfCoeffs) {
	unsigned int lengthOfAdaptationData = referenceInput.size(); // Length of Signals
	unsigned int lengthOfFilter = ccfCoeffs.size(); // Correlation Matrix Size
	Number sum = 0;
	for (unsigned int i = 0; i < lengthOfFilter; i++) {
		for (unsigned int j = i; j < lengthOfAdaptationData; j++) {
			sum += (userInput[j] * referenceInput[j - i]);
		}
		ccfCoeffs[i] = sum / (Number) lengthOfAdaptationData;
		sum = 0;
	}

}

//--------------------------------------------------------------------------//
// Function to calculate cross-correlation matrix of two signals reference and user input
//--------------------------------------------------------------------------//
// input :
// 	referenceInput: reference input sequence for which cross-correlation is to be computed
// 	userInput: user input sequence for which cross-correlation is to be computed
// output :
// 	ccfCoeffs - Cross-correlation Coefficients Matrix
//--------------------------------------------------------------------------//

void Dsp::crossCorrMatrix(
		NumberArrayConstRef referenceInput,
		NumberArrayConstRef userInput,
		Number2DArrayRef ccfCoeffs) {
	unsigned int lengthOfAdaptationData = referenceInput.size(); // Length of Signals
	unsigned int lengthOfFilter = ccfCoeffs.size(); // Correlation Matrix Size
	Number sum = 0;
	for (unsigned int i = 0; i < lengthOfFilter; i++) {
		for (unsigned int j = i; j < lengthOfAdaptationData; j++) {
			sum += (userInput[j] * referenceInput[j - i]);
		}
		ccfCoeffs[i][0] = sum / (Number) lengthOfAdaptationData;
		sum = 0;
	}

}

//--------------------------------------------------------------------------//
// Function to calculate DTW path matrix having Sakoe-Chiba constraints between two same length signals: reference and user input
//--------------------------------------------------------------------------//
// input :
// 	refMIDINotesInterpolated: reference input sequence for which dtw path is to be computed
// 	userP: user input sequence for which dtw path is to be computed
// output :
// 	path - Vector containing column number of warped index
//--------------------------------------------------------------------------//

void Dsp::dtwPath(
		PitchesConstRef refMIDINotesInterpolated,
		PitchesConstRef userP,
		unsigned int padVectorSizeForDTW,
		vector<int> & path,
		int oneSideConstraintLengthFrames) {

	// For Adding CONSTRAINT TO DTW, 2-D Similarity Matrix Initialized with maximum Cost 9999
	// TODO [Critical] - Maximum Cost should be greater than 9999???
	NumberArray similarityMatrix1D(refMIDINotesInterpolated.size(), 9999.0);
	Number2DArray similarityMatrix(refMIDINotesInterpolated.size(), similarityMatrix1D);

	// 2D Matrix containing cost of path upto a particular point
	NumberArray costMatrix1D(refMIDINotesInterpolated.size(), 0.0);
	Number2DArray costMatrix(refMIDINotesInterpolated.size(), costMatrix1D);

	// 2D Matrix Containing direction to be moved to
	IntArray moveMatrix1D(refMIDINotesInterpolated.size(), 0);
	Int2DArray moveMatrix(refMIDINotesInterpolated.size(), moveMatrix1D);

	// 2-D DTW Path Matrix Internal
	IntArray path1D(refMIDINotesInterpolated.size(), 0);
	Int2DArray path2D(refMIDINotesInterpolated.size(), path1D);

	unsigned int bandSakoeChiba = min(oneSideConstraintLengthFrames, int(refMIDINotesInterpolated.size() / 2));

	// Using Euclidean Distance to compute Similarity Matrix within Sakoe Chiba Band
	for (unsigned int i = 0; i < bandSakoeChiba; i++) {
		for (unsigned int j = 0; j < bandSakoeChiba + i; j++) {
			similarityMatrix[i][j] = fabs(refMIDINotesInterpolated[i] - userP[j]);
		}
	}
	for (unsigned int i = bandSakoeChiba; i < refMIDINotesInterpolated.size() - bandSakoeChiba; i++) {
		for (unsigned int j = i - bandSakoeChiba + 1; j < i - bandSakoeChiba + 1 + (2 * bandSakoeChiba); j++) {
			similarityMatrix[i][j] = fabs(refMIDINotesInterpolated[i] - userP[j]);
		}
	}
	for (unsigned int i = refMIDINotesInterpolated.size() - bandSakoeChiba; i < refMIDINotesInterpolated.size(); i++) {
		for (unsigned int j = i - bandSakoeChiba + 1; j < refMIDINotesInterpolated.size(); j++) {
			similarityMatrix[i][j] = fabs(refMIDINotesInterpolated[i] - userP[j]);
		}
	}

	/*
	 for (unsigned int i=0;i<refMIDINotesInterpolated.size();i++)
	 {
	 for (unsigned int j=0;j<refMIDINotesInterpolated.size();j++)
	 {
	 SimilarityMatrix[i][j]=fabs(refMIDINotesInterpolated[i]-userP[j]);
	 //SimilarityMatrix[i][j]=sqrt((Number)abs((NoteDiscreteValuesForDTW[i]*NoteDiscreteValuesForDTW[i])-(UserPitchValuesForDTW[j]*UserPitchValuesForDTW[j])));
	 }

	 }
	 */

	// Compute cost across different paths
	costMatrix[0][0] = similarityMatrix[0][0];
	moveMatrix[0][0] = 2;

	// Compute Cost and Path along the edges of the matrix
	for (unsigned int i = 1; i < refMIDINotesInterpolated.size(); i++) {
		costMatrix[0][i] = costMatrix[0][i - 1] + similarityMatrix[0][i];
		moveMatrix[0][i] = 1;
	}
	for (unsigned int i = 1; i < refMIDINotesInterpolated.size(); i++) {
		costMatrix[i][0] = costMatrix[i - 1][0] + similarityMatrix[i][0];
		moveMatrix[i][0] = 3;
	}

	// Compute Cost and path within the matrix
	Number MinValueofTrio = 0;
	for (unsigned int i = 1; i < refMIDINotesInterpolated.size(); i++) {
		for (unsigned int j = 1; j < refMIDINotesInterpolated.size(); j++) {
			if (costMatrix[i][j - 1] + 0.3 <= costMatrix[i - 1][j - 1]) {
				moveMatrix[i][j] = 1;
				MinValueofTrio = costMatrix[i][j - 1] + 0.3;
			} else {
				moveMatrix[i][j] = 2;
				MinValueofTrio = costMatrix[i - 1][j - 1];
			}
			if (MinValueofTrio > costMatrix[i - 1][j] + 0.3) {
				moveMatrix[i][j] = 3;
				MinValueofTrio = costMatrix[i - 1][j] + 0.3;
			}
			costMatrix[i][j] = similarityMatrix[i][j] + MinValueofTrio;
		}
	}

	// Compute 2D DTW Path from cost matrix
	// TODO [Optimization] can the 1D Path Vector be filled in this loop instead of a separate loop
	int XCord = refMIDINotesInterpolated.size() - 1;
	int YCord = refMIDINotesInterpolated.size() - 1;
	while (1) {
		path2D[YCord][XCord] = 1;
		if ((XCord == 0) && (YCord == 0)) {
			break;
		}

		if (moveMatrix[YCord][XCord] == 1) {
			XCord--;
		} else if (moveMatrix[YCord][XCord] == 2) {
			XCord--;
			YCord--;
		} else if (moveMatrix[YCord][XCord] == 3) {
			YCord--;
		}
	}

	// Convert 2D DTW Path matrix to 1D vector containing warped index
	for (unsigned int index = padVectorSizeForDTW; index < path2D[index].size() - padVectorSizeForDTW; index++) {

		for (unsigned int i = path2D[index].size() - padVectorSizeForDTW; i >= padVectorSizeForDTW; i--) {
			// TODO [Critical] Change to include +-50 Cent on either side of each note
			if (path2D[index][i] == 1) {
				path[index - padVectorSizeForDTW] = i - padVectorSizeForDTW;
				break;
			}
		}
	}

}

//--------------------------------------------------------------------------//
// Function to calculate constrained DTW path matrix (constraint is hard-coded to 8 frames) between two same length signals: reference and user input
//--------------------------------------------------------------------------//
// input :
// 	refMIDINotesInterpolated: reference input sequence for which dtw path is to be computed
// 	userP: user input sequence for which dtw path is to be computed
// output :
// 	path - Vector containing column number of warped index
//--------------------------------------------------------------------------//
void Dsp::dtwPath(
		PitchesConstRef refMIDINotesInterpolated,
		PitchesConstRef userP,
		unsigned int padVectorSizeForDTW,
		vector<int> & path) {

	// For Adding CONSTRAINT TO DTW Similarity Matrix Initialized with maximum Cost 9999
	NumberArray SimilarityMatrix1D(refMIDINotesInterpolated.size(), 9999.0);
	Number2DArray SimilarityMatrix(refMIDINotesInterpolated.size(), SimilarityMatrix1D);
	NumberArray DTW1D(refMIDINotesInterpolated.size(), 0.0);
	Number2DArray DTW(refMIDINotesInterpolated.size(), DTW1D);
	IntArray MoveMatrix1D(refMIDINotesInterpolated.size(), 0);
	Int2DArray MoveMatrix(refMIDINotesInterpolated.size(), MoveMatrix1D);

	IntArray path1D(refMIDINotesInterpolated.size(), 0);
	Int2DArray path2D(refMIDINotesInterpolated.size(), path1D);

	unsigned int bandSakoeChiba = min(8, int(refMIDINotesInterpolated.size() / 2)); // 150ms delay or advance
	// Using Euclidean Distance to compute Similarity Matrix
	for (unsigned int i = 0; i < bandSakoeChiba; i++) {
		for (unsigned int j = 0; j < bandSakoeChiba + i; j++) {
			SimilarityMatrix[i][j] = fabs(refMIDINotesInterpolated[i] - userP[j]);
		}
	}
	for (unsigned int i = bandSakoeChiba; i < refMIDINotesInterpolated.size() - bandSakoeChiba; i++) {
		for (unsigned int j = i - bandSakoeChiba + 1; j < i - bandSakoeChiba + 1 + (2 * bandSakoeChiba); j++) {
			SimilarityMatrix[i][j] = fabs(refMIDINotesInterpolated[i] - userP[j]);
		}
	}
	for (unsigned int i = refMIDINotesInterpolated.size() - bandSakoeChiba; i < refMIDINotesInterpolated.size(); i++) {
		for (unsigned int j = i - bandSakoeChiba + 1; j < refMIDINotesInterpolated.size(); j++) {
			SimilarityMatrix[i][j] = fabs(refMIDINotesInterpolated[i] - userP[j]);
		}
	}

	DTW[0][0] = SimilarityMatrix[0][0];
	MoveMatrix[0][0] = 2;
	for (unsigned int i = 1; i < refMIDINotesInterpolated.size(); i++) {
		DTW[0][i] = DTW[0][i - 1] + SimilarityMatrix[0][i];
		MoveMatrix[0][i] = 1;
	}
	for (unsigned int i = 1; i < refMIDINotesInterpolated.size(); i++) {
		DTW[i][0] = DTW[i - 1][0] + SimilarityMatrix[i][0];
		MoveMatrix[i][0] = 3;
	}
	Number MinValueofTrio = 0;
	for (unsigned int i = 1; i < refMIDINotesInterpolated.size(); i++) {
		for (unsigned int j = 1; j < refMIDINotesInterpolated.size(); j++) {
			if (DTW[i][j - 1] + 0.3 <= DTW[i - 1][j - 1]) {
				MoveMatrix[i][j] = 1;
				MinValueofTrio = DTW[i][j - 1] + 0.3;
			} else {
				MoveMatrix[i][j] = 2;
				MinValueofTrio = DTW[i - 1][j - 1];
			}
			if (MinValueofTrio > DTW[i - 1][j] + 0.3) {
				MoveMatrix[i][j] = 3;
				MinValueofTrio = DTW[i - 1][j] + 0.3;
			}
			DTW[i][j] = SimilarityMatrix[i][j] + MinValueofTrio;
		}
	}
	int XCord = refMIDINotesInterpolated.size() - 1;
	int YCord = refMIDINotesInterpolated.size() - 1;
	while (1) {
		path2D[YCord][XCord] = 1;
		if ((XCord == 0) && (YCord == 0)) {
			break;
		}

		if (MoveMatrix[YCord][XCord] == 1) {
			XCord--;
		} else if (MoveMatrix[YCord][XCord] == 2) {
			XCord--;
			YCord--;
		} else if (MoveMatrix[YCord][XCord] == 3) {
			YCord--;
		}
	}

	for (unsigned int index = padVectorSizeForDTW; index < path2D[index].size() - padVectorSizeForDTW; index++) {

		for (unsigned int i = path2D[index].size() - padVectorSizeForDTW; i >= padVectorSizeForDTW; i--) {
			// TODO Change to include +-50 Cent on either side of each note
			if (path2D[index][i] == 1) {
				path[index - padVectorSizeForDTW] = i - padVectorSizeForDTW;
				break;
			}
		}
	}

}

//--------------------------------------------------------------------------//
// calculates energy from complex vector
//--------------------------------------------------------------------------//
// input : real and imaginary vector
// output : energy vector
//--------------------------------------------------------------------------//
void Dsp::esd(
		NumberArrayRef vEsd,
		NumberArrayConstRef vXr,
		NumberArrayConstRef vXi) {
	unsigned const int size = vEsd.size();
	for (unsigned int i = 0; i < size; i++) {
		vEsd[i] = (vXr[i] * vXr[i]) + (vXi[i] * vXi[i]);
	}
}

//--------------------------------------------------------------------------//
// calculates magnitude of complex vector
//--------------------------------------------------------------------------//
// input : real and imaginary vector
// output : magnitude vector
//--------------------------------------------------------------------------//
void Dsp::magnitude(
		NumberArrayRef vEsd,
		NumberArrayConstRef vXr,
		NumberArrayConstRef vXi) {
	unsigned const int size = vEsd.size();
	for (unsigned int i = 0; i < size; i++) {
		vEsd[i] = sqrt((vXr[i] * vXr[i]) + (vXi[i] * vXi[i]));
	}
}

//--------------------------------------------------------------------------//
// calculates energy from output of ffts function
//--------------------------------------------------------------------------//
// input : output vector from ffts
// output : energy vector
//--------------------------------------------------------------------------//
void Dsp::esdFFTS(
		NumberArrayRef vEsd,
		NumberArrayConstRef vX) {
	unsigned const int size = vEsd.size();
	for (unsigned int i = 0, j = 0; i < size; i = i + 2, j++) {
		vEsd[j] = (vX[i] * vX[i]) + (vX[i + 1] * vX[i + 1]);
	}
}

//--------------------------------------------------------------------------//
// Multiply a vector with scaling factor
//--------------------------------------------------------------------------//
// input : vector and scaling factor
// output : scaled vector
//--------------------------------------------------------------------------//
void Dsp::scale(
		vector<Number> & vect,
		Number scale) {
	std::transform(vect.begin(), vect.end(), vect.begin(), std::bind1st(std::multiplies<Number>(), scale));
}

//--------------------------------------------------------------------------//
// Divides a vector A by B only if value of B is greater than a threshold
//--------------------------------------------------------------------------//
// input : vector A and B, threshold value
// output : thresholded vector A/B
//--------------------------------------------------------------------------//
void Dsp::divideAByB_thresholded(
		vector<Number> & vResult,
		vector<Number> const & vA,
		vector<Number> const & vB,
		Number threshold) {
	for (unsigned int i = 0; i < vResult.size(); i++) {
		if (threshold < vB[i]) {
			vResult[i] = vA[i] / vB[i];
		} else {
			vResult[i] = 0;
		}
	}
}

//--------------------------------------------------------------------------//
// calculates FFT and IFFT according to numerical recipes
//--------------------------------------------------------------------------//
// input :
//	dir: FFT/IFFT flag,
//	m: nfft,
//  n:size of input vector,
//	xr= input vector
// output :
//	xr- contains the real values of output,
//	xi= contains the imaginary values of output
//--------------------------------------------------------------------------//
void Dsp::fft(
		FFT_DIRECTION dir,
		unsigned const int m,
		unsigned const int n,
		NumberArrayRef xr,
		NumberArrayRef xi) {

	unsigned int i, i1, j, k, i2, l, l1, l2;
	Number c1, c2, tx, ty, t1, t2, u1, u2, z;

	// number of points n=2^m
	/*
	 n = 1;
	 for (i = 0; i < m; i++)
	 n *= 2;
	 */

	// bit reversal
	i2 = n >> 1;
	j = 0;
	for (i = 0; i < n - 1; i++) {
		if (i < j) {
			tx = xr[i];
			ty = xi[i];
			xr[i] = xr[j];
			xi[i] = xi[j];
			xr[j] = tx;
			xi[j] = ty;
		}
		k = i2;
		while (k <= j) {
			j -= k;
			k >>= 1;

		}
		j += k;
	} //end of bit reversal loop

	//Compute FFT

	c1 = -1.0;
	c2 = 0.0;
	l2 = 1;

	for (l = 0; l < m; l++) {
		l1 = l2;
		l2 <<= 1;
		u1 = 1.0;
		u2 = 0.0;

		for (j = 0; j < l1; j++) {
			for (i = j; i < n; i += l2) {
				i1 = i + l1;
				t1 = u1 * xr[i1] - u2 * xi[i1];
				t2 = u1 * xi[i1] + u2 * xr[i1];

				xr[i1] = xr[i] - t1;
				xi[i1] = xi[i] - t2;

				xr[i] += t1;
				xi[i] += t2;

			}

			z = u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = z;
		}

		c2 = sqrt((1.0 - c1) / 2.0);
		if (dir == FORWARD_FFT)
			c2 = -c2;
		c1 = sqrt((1.0 + c1) / 2.0);

	}

	if (dir == INVERSE_FFT) {
		for (i = 0; i < n; i++) {
			xr[i] /= n;
			xi[i] /= n;
		}
	}
}
//--------------------------------------------------------------------------//
// gives the next power of 2 for a given input
//--------------------------------------------------------------------------//
// input : input length
// output : nearest power of 2 (if length =10 then output is 4 because 2^4 =16)
//--------------------------------------------------------------------------//

int Dsp::nextPow2(
		int length) {
	return ceil(log10((Number)length) / log10((Number)2.0));
}

//--------------------------------------------------------------------------//
// Compute DC component of vibrato
//--------------------------------------------------------------------------//
// input :
// output :
//--------------------------------------------------------------------------//
Number Dsp::computeDCOfVibrato(
		NumberArrayConstRef pdV) {
	Number dc = 0;
	int totalNonZeroElem = 0;
	unsigned int i;

	for (i = 0; i < pdV.size(); i++) {
		// Checking for INVALID_CENT PITCHES
		if (pdV[i] != INVALID_CENT) {
			dc += pdV[i];
			totalNonZeroElem++;
		}
	}
	dc /= totalNonZeroElem;
	return dc;
}

//--------------------------------------------------------------------------//
// Getting value of median for vibrato
//--------------------------------------------------------------------------//
// input = vector pdV
// output= returns median value
//--------------------------------------------------------------------------//
Number Dsp::median(
		NumberArrayConstRef pdV,
		Number dMedianFactor) {
	Number median;
	NumberArray NonZeropdV;
	// Find INVALID PITCHES IN THIS SO AS TO IGNORE IN MEDIAN COMPUTATION
	for (unsigned int i = 0; i < pdV.size(); i++) {
		if (pdV[i] != INVALID_CENT) {
			NonZeropdV.push_back(pdV[i]);
		}
	}
	sort(NonZeropdV.begin(), NonZeropdV.end());
	median = NonZeropdV.size() > 0 ? NonZeropdV[int((NonZeropdV.size() * dMedianFactor) + 0.5)] : INVALID_CENT;
	return median;
}

//--------------------------------------------------------------------------//
// sums vector values between 2 indices
//--------------------------------------------------------------------------//
// input :
//	esd: input vector
//	startIdx: start index,
//	endIdx: end index
// output : returns sum value
//--------------------------------------------------------------------------//
Number Dsp::sum(
		NumberArrayConstRef esd,
		int startIdx,
		int endIdx) {
	Number mySum = 0;
	for (int i = startIdx; i < endIdx; i++)
		mySum += esd[i];
	return mySum;
//	return std::accumulate(esd.begin()+startIdx,esd.begin()+1+endIdx,0);//#include <numeric>
}

int Dsp::sum(
		IntArrayRef esd,
		int startIdx,
		int endIdx) {
	int mySum = 0;
	for (int i = startIdx; i < endIdx; i++)
		mySum += esd[i];
	return mySum;
//	return std::accumulate(esd.begin()+startIdx,esd.begin()+1+endIdx,0);//#include <numeric>
}

//--------------------------------------------------------------------------//
// returns the max value and the index of it from a vector
//--------------------------------------------------------------------------//
// input :
//	v: input vector
//	startIdx: start index,
//	endIdx: end index
// output : returns the maximum value, and outputs the max index through maxIdx function variable.
//--------------------------------------------------------------------------//
Number Dsp::max(
		NumberArrayConstRef v,
		int startIdx,
		int endIdx,
		unsigned int * maxIdx) {
	Number maxVal = -9e99;
	for (int i = startIdx; i <= endIdx; i++) {
		if (v[i] > maxVal) {
			maxVal = v[i];
			*maxIdx = i;
		}
	}
	return maxVal;
}

//--------------------------------------------------------------------------//
// returns the sum square of a vector subtracted by a constant value
//--------------------------------------------------------------------------//
// input : vector data, constant value to be subtracted
// output : returns sum  square of a vector subtracted by a constant value
//--------------------------------------------------------------------------//
Number Dsp::sumSquareSub(
		NumberArrayConstRef data,
		Number value) {
	int iNumElements = data.size();
	Number iSum = 0;
	for (int i = 0; i < iNumElements; i++) {
		iSum += (data[i] - value) * (data[i] - value);

	}
	return iSum;
}

//--------------------------------------------------------------------------//
// returns the sum square of a vector subtracted by another vector
//--------------------------------------------------------------------------//
// input : vector data and vector data2
// output : returns sum square of vector subtracted by another vector
//--------------------------------------------------------------------------//
Number Dsp::sumSquareSubm(
		NumberArrayConstRef data,
		NumberArrayConstRef data2) {
	int iNumElements = data.size();
	Number iSum = 0;
	for (int i = 0; i < iNumElements; i++) {
		iSum += (data[i] - data2[i]) * (data[i] - data2[i]);
	}
	return iSum;
}

Dsp::~Dsp() {
}

//--------------------------------------------------------------------------//
// hanning window
//--------------------------------------------------------------------------//
// input : hanning window size N
// output : hanning window hannWin
//--------------------------------------------------------------------------//
void Dsp::hanning(
		NumberArrayRef hannWin) {
	int N = hannWin.size();
	for (int i = 0; i < N; i++)
		hannWin[i] = .5 * (1 - cos(2 * PI * i / (N - 1)));
}

//--------------------------------------------------------------------------//
// hamming window
//--------------------------------------------------------------------------//
// input : hamming window size N
// output : hamming window hammingWin
//--------------------------------------------------------------------------//
void Dsp::hamming(
		NumberArrayRef hammingWin) {
	int N = hammingWin.size();
	for (int i = 0; i < N; i++)
		hammingWin[i] = .54 - 0.46 * cos(2 * PI * i / (N - 1));
}

//--------------------------------------------------------------------------//
// sinc value for index
//--------------------------------------------------------------------------//
// input : index x
// output : sinc value for index x
//--------------------------------------------------------------------------//
Number Dsp::sinc(
		Number x) {
	return (0 == x) ? 1 : (sin(PI * x) / (PI * x));
}

// Greatest Common Divisor
// input : values a and b
// output : returns gcd of a and b
int Dsp::gcd(
		int a,
		int b) {
	int c;
	while (a != 0) {
		c = a;
		a = b % a;
		b = c;
	}
	return b;
}

//--------------------------------------------------------------------------//
// carries out autocorrelation in limited area, specified by order (used by LPC)
//--------------------------------------------------------------------------//
// input :
//	currWin: input vector
//	order : search area order
// output : vector of autocorrelated points aCorrPts
//--------------------------------------------------------------------------//
void Dsp::autoCorr(
		NumberArrayRef currWin,
		int order,
		NumberArrayRef aCorrPts) {
	int winLen = currWin.size();
	// compute autocorrelations
	for (int i = 0; i <= order; i++) {
		Number sum = 0;
		for (int k = 0; k < winLen - i; k++)
			sum += currWin[k] * currWin[k + i];
		aCorrPts[i] = sum;
	}
}

//--------------------------------------------------------------------------//
// Fast resampling of input vector using sinc function, with subsample shift if required
//--------------------------------------------------------------------------//
// input :
//	userFile: input vector ,
//	shift: resampled value at a new subsample shift,
//  newFreq: new resampled frequency
// output :
//	userFile: returns resampled audio in the same input vector
//--------------------------------------------------------------------------//
void Dsp::sincResample(
		NumberArrayRef userFile,
		Number shift,
		int newFreq) {
	int sincLen = 0.002 * SAMPLING_RATE; // cut off for sinc ripples, set to 2ms.
	int userAudioLen = userFile.size();
	int gcdValue = gcd(SAMPLING_RATE, newFreq);
	int numFilt = newFreq / gcdValue;
	Number multirateFrac = (Number) SAMPLING_RATE / newFreq;

	NumberArray audioIn(userAudioLen + sincLen, 0.0); // on resampling, the resampled audio is advanced by 'sincLen' number of samples,
	//to compensate this, the input audio is padded with 'sincLen' samples of audio in the beginning
	for (int i = 0; i < userAudioLen; i++) {
		audioIn[i + sincLen] = userFile[i];
	}
	fill(userFile.begin(), userFile.end(), 0);

	//------------- INITIALIZE-----------------//
	NumberArray hann(sincLen, 0);
	hanning(hann);

	Number d;
	NumberArray nInDelayFract(numFilt, 0.0);
	for (int i = 0; i < numFilt; i++) {
		d = i * multirateFrac;
		nInDelayFract[i] = d - floor(d) - shift;
	}

	// calculate Sinc coeffecients at the delay factors calulated above
	vector<vector<Number> > h(numFilt, NumberArray(sincLen, 0.0));

	Number nH;
	for (int i = 0; i < numFilt; i++) {
		for (int j = 0; j < sincLen; j++) {
			nH = (sincLen / 2) - 1.0 - j + nInDelayFract[i];
			h[i][j] = sinc(nH) * hann[j];
		}
	}

	int orgFact = SAMPLING_RATE / gcdValue;
	for (int i = 0; i < userAudioLen; i++) {
		int iQ = (int) floor(i / (Number) numFilt);
		int iR = i % numFilt;
		int delta = (iQ * orgFact) + (int) floor(iR * ((Number) orgFact / (Number) numFilt));
		int filterIdx = i % numFilt;

		for (int j = 0; j < sincLen; j++) {
			if ((j + delta + (sincLen / 2) + 1) < userAudioLen + sincLen) {
				userFile[i] += (audioIn[j + delta + (sincLen / 2) + 1] * h[filterIdx][j]);
			}
		}
	}

	if (newFreq != SAMPLING_RATE) {
		userFile.resize(userAudioLen / multirateFrac, 0.0);
	}
	//int xUpLen= (userAudioLen+sincLen)/multirateFrac;
}

//--------------------------------------------------------------------------//
// Computes LPC coeff
//--------------------------------------------------------------------------//
// Input : autocorrelation coeffs
// Output: LPC coeff, flag to indicate unstable filter
//--------------------------------------------------------------------------//
void Dsp::getLpcCoeffs(
		NumberArrayRef acfCoeffs,
		NumberArrayRef lpCoeffs,
		int lpOrder) {
	NumberArray E(lpOrder + 1, 0.0), K(lpOrder + 1, 0.0);
	vector<vector<Number> > a(lpOrder + 1, NumberArray(lpOrder + 1, 0.0));
	Number sum;
	int i, j;

	// Step S1
	E[0] = acfCoeffs[0];// Energy in the sequence is the first value of the autocorrelation of the iwndowed s[n] i.e r[0]
	a[0][0] = 1;

	// Step S5
	for (i = 1; i <= lpOrder; i++) {
		// Step S1
		sum = 0;
		a[i][0] = 1;

		// Step S2
		for (j = 0; j < i; j++) {
			sum += acfCoeffs[i - j] * a[i - 1][j];
		}
		K[i] = (-1) * sum / E[i - 1];

		if ((K[i] > 1.0) || (K[i] < -1.0)) {
			cout << "Unstable filter for frame using prev frame's LPC coeff & gain\n" << endl;
			// TODO[Critical]: Change this to work for failed case
			// return a_prev;
		}

		// Step S3
		a[i][i] = K[i];
		for (j = 1; j < i; j++) {
			a[i][j] = a[i - 1][j] + K[i] * a[i - 1][i - j];
		}

		// Step S4
		E[i] = (1 - K[i] * K[i]) * E[i - 1];
	}

	a[lpOrder][0] = 1.0;

	// Copy results in output array
	for (i = 0; i <= lpOrder; i++) {
		lpCoeffs[i] = a[lpOrder][i];
	}

}
//--------------------------------------------------------------------------//
// filtering code
//--------------------------------------------------------------------------//
// input:
//	order: filter order is always length(a)-1,
//	b & a: coefficients vector, length of vector (a)==vector(b) fill necessary zeros if required
// output: filtered vector y
//--------------------------------------------------------------------------//
void Dsp::filter(
		int ord,
		NumberArrayRef a,
		NumberArrayRef b,
		NumberArrayRef x,
		NumberArrayRef y) { // ord = is always length(a)-1

	int i, j;
	y[0] = b[0] * x[0];
	for (i = 1; i < ord + 1; i++) {
		y[i] = 0.0;
		for (j = 0; j < i + 1; j++)
			y[i] = y[i] + b[j] * x[i - j];
		for (j = 0; j < i; j++)
			y[i] = y[i] - a[j + 1] * y[i - j - 1];

		if (a[0] != 0 && a[0] != 1)
			y[i] = y[i] / a[0];
	}
	/* end of initial part */
	for (i = ord + 1; i < (int) x.size(); i++) {
		y[i] = 0.0;
		for (j = 0; j < ord + 1; j++)
			y[i] = y[i] + b[j] * x[i - j];
		for (j = 0; j < ord; j++)
			y[i] = y[i] - a[j + 1] * y[i - j - 1];

		if (a[0] != 0 && a[0] != 1)
			y[i] = y[i] / a[0];
	}

} /* end of filter */

void Dsp::prEmpFilter(
		NumberArrayRef x,
		NumberArrayRef y) {
	y[0] = x[0];
	for (int i = 1; i < (int) x.size(); i++) {
		y[i] = x[i] - 0.97 * x[i - 1];
	}
}
//--------------------------------------------------------------------------//
// time domain convolution of two vectors
//--------------------------------------------------------------------------//
//inputs:
//	x  : input vector
//	L  : length of x
//	h  : vector h
//outputs: convolved vector is written in x
//--------------------------------------------------------------------------//

void Dsp::conv(
		NumberArrayRef x,
		int L,
		NumberArrayRef h) {
	int M = h.size();
	//int L=x.size();

	//M is assumed odd
	int n, i, k;
	Number sum = 0;
	int halfspan;

	if (M % 2 == 0)
		halfspan = M / 2;
	else
		halfspan = (M - 1) / 2;

	NumberArray temp((L + M), 0.0);
	for (i = 0; i < L + M; i++) //zero pad the signal at the start and end.
			{
		if (i < halfspan)
			temp[i] = 0;
		else if (i >= halfspan && i < (L + halfspan))
			temp[i] = x[i - halfspan]; //temp[halfspan] = x[0] achieved
		else
			temp[i] = (x[L - 1] + x[L - 2]) / 2;  //to avoid the biphasic giving strong -ve peaks.

	}
	fill(x.begin(), x.end(), 0.0);
	for (n = 0; n < L; n++) {
		sum = 0;
		for (k = 0; k < M; k++)
			sum += temp[n - k + (M - 1)] * h[k];
		x[n] = sum;

		//if(n<L) printf("\n%f\t%f",temp[n], y[n]);
	}
}

/********************************************************************/
/* Input signal : x (length L)
 Filter kernel   : h (length M (presently expected to be an odd number))
 Output signal   : y  (length L+M-1)

 Implementation: eg:
 x = 1 1 1
 h = 1 2 3 (gather kernel)
 y = 1 3 6 5 1
 (same as matlab conv.m)
 */
/********************************************************************/
void Dsp::convolveAdaptiveFilter(
		NumberArrayRef x,
		NumberArrayRef h) {
	int n, i, k;
	Number sum = 0;
	int L = x.size();
	int M = h.size();

	NumberArray temp((L + 2 * M - 2), 0.0);

	for (i = 0; i < L + 2 * M - 2; i++) //zero pad the signal at the start and end.
			{
		if (i < M - 1)
			temp[i] = 0;
		else if (i >= (M - 1) && i < (L + M))
			temp[i] = x[i - (M - 1)];
		else
			temp[i] = 0;

	}
	NumberArray filtereddata(L + M - 1, 0.0);
	for (n = M - 1; n < L + 2 * M - 2; n++) {
		sum = 0;
		for (k = 0; k < M; k++)
			sum += temp[n - k] * h[k];
		filtereddata[n - (M - 1)] = sum;
	}

	for (i = 0; i < L; i++) {
		x[i] = filtereddata[i];
	}

}

//--------------------------------------------------------------------------//
// frequency domain convolution of two vectors
//inputs:
//	aud  : input vector
//	rir  : vector h
//outputs:
//	acousticVec : convolved vector
//--------------------------------------------------------------------------//
void Dsp::fftConv(
		NumberArrayRef aud,
		NumberArrayConstRef rir,
		NumberArrayRef acousticVec) {

	int Len = aud.size() + rir.size() - 1;
	int m = Dsp::nextPow2(Len);
	int outLen = pow((Number) 2, (Number) m);
	NumberArray audR(outLen, 0.0), audI(outLen, 0.0), rirR(outLen, 0.0), rirI(outLen, 0.0), outR(outLen, 0.0), outI(
			outLen, 0.0);

	for (int i = 0; i < (int) aud.size(); i++)
		audR[i] = aud[i];
	for (int i = 0; i < (int) rir.size(); i++)
		rirR[i] = rir[i];

	Dsp::fft(Dsp::FORWARD_FFT, m, outLen, audR, audI);
	Dsp::fft(Dsp::FORWARD_FFT, m, outLen, rirR, rirI);
	for (int i = 0; i < outLen; i++) {
		outR[i] = audR[i] * rirR[i] - audI[i] * rirI[i];
		outI[i] = audR[i] * rirI[i] + audI[i] * rirR[i];
	}
	Dsp::fft(Dsp::INVERSE_FFT, m, outLen, outR, outI);

	for (int i = 0; i < (int) aud.size(); i++)
		acousticVec[i] = outR[i];
}

//--------------------------------------------------------------------------//
// frequency domain convolution of two vectors
//--------------------------------------------------------------------------//
//inputs:
//	aud  : input vector
//	rir  : vector h
//outputs:
//	acousticVec : convolved vector
//--------------------------------------------------------------------------//
//void Dsp::fftConvFFTvec(
//		NumberArrayRef aud,
//		NumberArrayConstRef rirFFT,
//		NumberArrayConstRef rirFFTI,
//		NumberArrayRef acousticVec) {
//
//	int N = rirFFT.size();
//	int M = N / 2;
//	int nfft = log((Number) N) / M_LN2;
//	if (N > (int) aud.size())
//		aud.resize(N, 0.0);
//	for (int frame = 0; frame < (int) aud.size() - M + 1; frame = frame + M) {
//		NumberArray audR(N, 0.0), audI(N, 0.0), outR(N, 0.0), outI(N, 0.0);
//
//		for (int i = 0; i < M; i++)
//			audR[i] = aud[frame + i];
//		Dsp::fft(Dsp::FORWARD_FFT, nfft, N, audR, audI);
//
//		for (int i = 0; i < N; i++) {
//			outR[i] = audR[i] * rirFFT[i] - audI[i] * rirFFTI[i];
//			outI[i] = audR[i] * rirFFTI[i] + audI[i] * rirFFT[i];
//		}
//		Dsp::fft(Dsp::INVERSE_FFT, nfft, N, outR, outI);
//	}
//}

//--------------------------------------------------------------------------//
// warp lpc coefficients for voice transformation. This function gives the numerator coefficients of the new transformed FIR filter
//--------------------------------------------------------------------------//
//inputs:
//	alpha: lpc coefficients
//outputs:
//	Bcoeff: transformed coefficients
//--------------------------------------------------------------------------//
void Dsp::warpBpole(
		Number alpha,
		NumberArrayRef Bcoeff) {
	int order = Bcoeff.size();
	NumberArray c(2, 0);
	c[0] = 1;
	c[1] = -alpha;

	Bcoeff[0] = 1;
	for (int n = 1; n < order; n++) {
		conv(Bcoeff, order, c);
	}

}

//--------------------------------------------------------------------------//
// warp lpc coefficients for voice transformation. This function gives the denominator coefficients of the new transformed FIR filter
//--------------------------------------------------------------------------//
//inputs:
//  a: transformation coefficient - Negative alpha shifts poles up in frequency
//	alpha: lpc coefficients
//outputs:
//	Acoeff: transformed coefficients
//--------------------------------------------------------------------------//
void Dsp::warpApole(
		NumberArrayRef a,
		Number alpha,
		NumberArrayRef Acoeff) {

	int order = a.size();
	NumberArray d(2, 0), dd(order + 1, 0);
	d[0] = -alpha;
	d[1] = 1;
	dd[0] = -alpha;
	dd[1] = 1;
	NumberArray c(2, 0), cc(order + 1, 0);
	c[0] = 1;
	c[1] = -alpha;
	cc[0] = 1;
	cc[1] = -alpha;

	Acoeff[0] = a[0];

	for (int n = 1; n < order; n++) {
		conv(Acoeff, order, c);
		for (int j = 0; j < order; j++)
			Acoeff[j] = Acoeff[j] + a[n] * dd[j];

		conv(dd, order, d);
		conv(cc, order, c);
	}
}

//--------------------------------------------------------------------------//
// energy of vector
//--------------------------------------------------------------------------//
//inputs:
//  userFile: input vector
//	startSample: start calculating energy from this point
//	endSample: end calculating energy at this point
//outputs:
//	returns energy value
//--------------------------------------------------------------------------//
Number Dsp::energy(
		NumberArrayRef userFile,
		int startSample,
		int endSample) {
	// Function to calculate the energy of a segment in a file (given by Start and End Samples)
	Number energy = 0;
	for (int i = startSample; i < endSample; i++)
		energy += userFile[i] * userFile[i];
	energy = energy / (endSample - startSample);
	return energy;
}

//--------------------------------------------------------------------------//
// change formant of user audio using LPC
//--------------------------------------------------------------------------//
//inputs:
//  userAudio: audio of user
//	alpha: move formants by this factor, Negative alpha shifts poles up in frequency
//outputs:
//	audioOut: formant shifted audio
//--------------------------------------------------------------------------//
void Dsp::changeFormant(
		NumberArrayRef userAudio,
		Number alpha,
		NumberArrayRef audioOut) {

	//-------------------------- INITIALIZE-----------------------------//
	int lpcOrder, hopLen, winLen;
	lpcOrder = 18;
	hopLen = 256;
	winLen = 2 * hopLen;

	//------------------- READ AUDIO-------------------------//
	int userAudioLen = userAudio.size();
	int numOfFrames = (userAudioLen - winLen) / hopLen;

	NumberArray audioBuff(userAudioLen, 0.0);
	NumberArray audioIn(userAudioLen, 0.0);

	//------------------ PRE-EMPHASIS FILTER -----------------------//
	NumberArray bPreEmp(2, 0.0);
	bPreEmp[0] = 1;
	bPreEmp[1] = -0.9;
	NumberArray aPreEmp(2, 0.0);
	aPreEmp[0] = 1;

	Dsp::filter(1, aPreEmp, bPreEmp, userAudio, audioIn);

	NumberArray analWin(winLen, 0.0);
	NumberArray hanWin(winLen, 0.0);
	Dsp::hanning(hanWin);

	NumberArray Bcoeff(lpcOrder + 1, 0.0);
	Dsp::warpBpole(alpha, Bcoeff);

	NumberArray acfCoeffs(lpcOrder + 1, 0.0);
	NumberArray dummy(lpcOrder + 1, 0.0);
	dummy[0] = 1;
	NumberArray lpCoeffs(lpcOrder + 1, 0.0);

	for (int frame = 0; frame < numOfFrames; frame++) {
		for (int i = 0, j = frame * hopLen; i < winLen; i++, j++) {
			analWin[i] = audioIn[j] * hanWin[i];
		}

		//--------------------- LPC Analysis---------------------//

		Dsp::autoCorr(analWin, lpcOrder, acfCoeffs);
		Dsp::getLpcCoeffs(acfCoeffs, lpCoeffs, lpcOrder);
		NumberArray residual(winLen, 0.0);
		Dsp::filter(lpcOrder, dummy, lpCoeffs, analWin, residual);

		//---------------------- Shift Formant--------------------//
		NumberArray Acoeff(lpcOrder + 1, 0.0);
		Dsp::warpApole(lpCoeffs, alpha, Acoeff);

		//---------------------- LPC Synthesis--------------------//
		NumberArray synWin(winLen, 0.0);
		Dsp::filter(lpcOrder, Acoeff, Bcoeff, residual, synWin);

		for (int i = 0, j = frame * hopLen; i < winLen; i++, j++) {
			audioBuff[j] += synWin[i];
		}
	}
	//------------------ DE-EMPHASIS FILTER ----------------------//
	Dsp::filter(1, bPreEmp, aPreEmp, audioBuff, audioOut);

}

//--------------------------------------------------------------------------//
// add flanger (chorus) effect to input vector
//--------------------------------------------------------------------------//
//inputs:
//  userAudVec: audio of user
//outputs:
//	awesomeVec: flanger added audio
//--------------------------------------------------------------------------//
void Dsp::addFlanger(
		NumberArrayRef userAudVec,
		NumberArrayRef awesomeVec) {

	int max_samp_delay = .01 * SAMPLING_RATE;
	Number amp = 0.7; // suggested coefficient from page 71 DAFX

	int numberOfSamples = awesomeVec.size();
	for (int i = 0; i < numberOfSamples; i++) {
		if (i > max_samp_delay) {
			int cur_delay = round(sin(2.0 * PI * (i >> 1) / SAMPLING_RATE)) * max_samp_delay;
			awesomeVec[i] = amp * userAudVec[i] + amp * userAudVec[i - cur_delay];   // add delayed sample
		}
	}

}


int Dsp::round(Number x)
{
   return (int)x >= (Number)0.0 ? floorf(x + (Number)0.5) : ceilf(x - (Number)0.5);
}

//--------------------------------------------------------------------------//
// add ring modulation(alien) effect to input vector
//--------------------------------------------------------------------------//
//inputs:
//  userAudVec: audio of user
//outputs:
//	awesomeVec: ring modulation added audio
//--------------------------------------------------------------------------//
void Dsp::addRingMod(
		NumberArrayRef userAudVec,
		NumberArrayRef awesomeVec) {
	//Ring Modulate with a sine wave frequency Fc
	int Fc = 150; //Both 150-350 is good
	int numberOfSamples = awesomeVec.size();
	for (int i = 0; i < numberOfSamples; i++) {
		awesomeVec[i] = userAudVec[i] * sin(2.0 * PI * i * (Fc / (Number) SAMPLING_RATE));
	}
}

//--------------------------------------------------------------------------//
// add echo effect to input vector
//--------------------------------------------------------------------------//
//inputs:
//  userAudVec: audio of user
//outputs:
//	awesomeVec: echo added audio
//--------------------------------------------------------------------------//
void Dsp::addEcho(
		NumberArrayRef awesomeVec,
		NumberArrayRef filteredVec) {
	Number echoDecay = 0.2; // decay in ms
	int NumberOfFadeInSamples = 0.18 * SAMPLING_RATE;
	int numberOfSamples = awesomeVec.size();
	for (int i = 0; i < numberOfSamples; i++) {
		if ((i >= NumberOfFadeInSamples) & (i <= numberOfSamples - NumberOfFadeInSamples)) {
			filteredVec[i] = awesomeVec[i] + (Number) echoDecay * awesomeVec[i - NumberOfFadeInSamples];
		} else if ((i < NumberOfFadeInSamples)) {
			filteredVec[i] = awesomeVec[i] * ((Number) (i) / NumberOfFadeInSamples);
		} else {
			filteredVec[i] = awesomeVec[i]
					+ ((Number) echoDecay * awesomeVec[i - NumberOfFadeInSamples])
							* ((Number) (numberOfSamples - i) / NumberOfFadeInSamples);
		}
	}
}

