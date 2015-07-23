#include "MelFilterBank.h"

MelFilterBank::~MelFilterBank() {
}
void MelFilterBank::melFilterBankInit(
		Number samplingFreq,
		Number lowFreq,
		Number highFreq,
		int numFilt,
		int nfftBy2) {
	m_FS = samplingFreq;
	m_lowerF = lowFreq;
	m_upperF = highFreq;
	m_numFilt = numFilt;
	m_nfftBy2 = nfftBy2;
	m_nfft = (nfftBy2 - 1) * 2;
	m_specStart.resize(m_numFilt, 0.0);
	m_filtStart.resize(m_numFilt, 0.0);
	m_filtWidth.resize(m_numFilt, 0.0);

	//create return array
	NumberArray centers(m_numFilt + 2, 0.0);
	Number maxFreqMel, minFreqMel, deltaFreqMel, deltaFreq;

	//compute mel min./max. frequency
	maxFreqMel = linearToMel(m_upperF);
	minFreqMel = linearToMel(m_lowerF);
	deltaFreqMel = (maxFreqMel - minFreqMel) / (m_numFilt + 1);
	deltaFreq = m_FS / (Number) m_nfft;

	/* Count and place filter coefficients. */
	int n_coeffs = 0;
	for (int i = 0; i < m_numFilt; ++i) {
		float freqs[3];

		/* Left, center, right frequencies in Hertz */
		for (int j = 0; j < 3; ++j) {
			freqs[j] = melToLinear((i + j) * deltaFreqMel + minFreqMel);
			/* Round them to DFT points if requested */
			freqs[j] = ((int) (freqs[j] / deltaFreq + 0.5)) * deltaFreq;
		}

		/* spec_start is the start of this filter in the power spectrum. */
		m_specStart[i] = -1;
		/* There must be a better way... */
		for (int j = 0; j < m_nfftBy2; ++j) {
			float hz = j * deltaFreq;
			if (hz < freqs[0])
				continue;
			else if (hz > freqs[2] || j == m_nfftBy2 - 1) {
				/* filt_width is the width in DFT points of this filter. */
				m_filtWidth[i] = j - m_specStart[i];
				/* filt_start is the start of this filter in the filt_coeffs array. */
				m_filtStart[i] = n_coeffs;
				n_coeffs += m_filtWidth[i];
				break;
			}
			if (m_specStart[i] == -1)
				m_specStart[i] = j;
		}
	}

	m_filtCoeff.resize(n_coeffs, 0.0);

	/* And now generate the coefficients. */
	n_coeffs = 0;
	for (int i = 0; i < m_numFilt; ++i) {
		float freqs[3];
		/* Left, center, right frequencies in Hertz */
		for (int j = 0; j < 3; ++j) {
			freqs[j] = melToLinear((i + j) * deltaFreqMel + minFreqMel);
			/* Round them to DFT points if requested */
			freqs[j] = ((int) (freqs[j] / deltaFreq + 0.5)) * deltaFreq;
		}

		for (int j = 0; j < m_filtWidth[i]; ++j) {
			float hz, loslope, hislope;

			hz = (m_specStart[i] + j) * deltaFreq;
			if (hz < freqs[0] || hz > freqs[2]) {
				cout << "ERROR: Failed to create filterbank, frequency range does not match. "
						"Sample rate " << m_FS << ", FFT size " << m_nfft << ", lowerf " << freqs[2] << " < freq " << hz
						<< " > upperf " << freqs[0] << "." << endl;
				exit(0);
			}
			loslope = (hz - freqs[0]) / (freqs[1] - freqs[0]);
			hislope = (freqs[2] - hz) / (freqs[2] - freqs[1]);

			loslope *= 2 / (freqs[2] - freqs[0]);
			hislope *= 2 / (freqs[2] - freqs[0]);

			if (loslope < hislope) {

				m_filtCoeff[n_coeffs] = loslope;
			} else {

				m_filtCoeff[n_coeffs] = hislope;
			}
//              cout <<m_filtCoeff[n_coeffs] << "\t";
			++n_coeffs;
		}
	}
//      cout<<endl;
}

NumberArray MelFilterBank::applyFilter(
		NumberArrayRef realPart,
		NumberArrayRef imagPart) {
	NumberArray filteredVec(m_numFilt, 0.0);
	NumberArray absSpectrum(m_nfftBy2, 0.0);
	for (int i = 0; i < m_nfftBy2; i++) {
		absSpectrum[i] = realPart[i] * realPart[i] + imagPart[i] * imagPart[i];
	}

	for (int whichfilt = 0; whichfilt < m_numFilt; whichfilt++) {
		filteredVec[whichfilt] = 0;
		for (int i = 0; i < m_filtWidth[whichfilt]; i++) {
			filteredVec[whichfilt] += absSpectrum[m_specStart[whichfilt] + i] * m_filtCoeff[m_filtStart[whichfilt] + i];
		}

		if (filteredVec[whichfilt] > 0)
			filteredVec[whichfilt] = log(filteredVec[whichfilt]);
		else
			filteredVec[whichfilt] = -10.0;
	}
//	for(int i=0; i<m_numFilt;i++)cout<<filteredVec[i]<<" ";
//	cout<<endl;
	return filteredVec;
}

