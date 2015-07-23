#ifndef MELFILTERBANK_H
#define MELFILTERBANK_H

#include "Dsp.h"
#include <cstddef>
#include <vector>
#include <math.h>

using namespace std;

class MelFilterBank {
public:
	MelFilterBank() {
	}
	;
	void melFilterBankInit(
			Number samplingFreq,
			Number lowFreq,
			Number highFreq,
			int numFilt,
			int nfftBy2);
	~MelFilterBank();
	NumberArray applyFilter(
			NumberArrayRef realPart,
			NumberArrayRef imagPart);
	static Number linearToMel(
			Number linearFrequency) {
		return 1127.01048 * std::log(1.0 + linearFrequency / 700.0);
	}
	static Number melToLinear(
			Number melFrequency) {
		return 700.0 * (std::exp(melFrequency / 1127.01048) - 1.0);
	}

private:
	IntArray m_specStart;
	IntArray m_filtStart;
	IntArray m_filtWidth;
	NumberArray m_filtCoeff;
	Number m_FS;
	Number m_lowerF;
	Number m_upperF;
	int m_numFilt;
	int m_nfftBy2;
	int m_nfft;
};

#endif // MFCC_H
