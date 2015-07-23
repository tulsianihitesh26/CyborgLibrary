#include "Dct.h"
void Dct::dctInit(
		int numberFilters,
		int numberCoefficients,
		int dctType) {
	m_dctType = dctType;
	//compute constants
	Number k = M_PI / numberFilters;
	m_sqrt_inv_n = 1.0 / (sqrt((Number)numberFilters));
	m_sqrt_inv_2n = sqrt(2.0 / numberFilters);

	//generate dct matrix
	for (int i = 0; i < numberCoefficients; i++) {
		NumberArray dctVect(numberFilters, 0.0);
		for (int j = 0; j < numberFilters; j++) {
			dctVect[j] = cos(k * i * (j + 0.5));
//			cout <<dctVect[j]<<"\t";
		}
		m_dctBank.push_back(dctVect);
//		cout<<endl;
	}
}

// Apply DCT transform
NumberArray Dct::applyDct(
		NumberArrayRef mflogspec) {
	NumberArray mfcep(m_dctBank.size(), 0.0);
	int i, j, beta;
	int numFilt = m_dctBank[0].size();
	int numCep = m_dctBank.size();

	if (m_dctType == 0) // legacy
			{
		mfcep[0] = mflogspec[0] / 2; /* beta = 0.5 */
		for (j = 1; j < numFilt; j++)
			mfcep[0] += mflogspec[j]; /* beta = 1.0 */
		mfcep[0] /= numFilt;

		for (i = 1; i < numCep; ++i) {
			mfcep[i] = 0;
			for (j = 0; j < numFilt; j++) {
				if (j == 0)
					beta = 1; /* 0.5 */
				else
					beta = 2; /* 1.0 */
				mfcep[i] += (mflogspec[j] * m_dctBank[i][j]) * beta;
			}
			/* Note that this actually normalizes by num_filters, like the
			 * original Sphinx front-end, due to the doubled 'beta' factor
			 * above.  */
			mfcep[i] /= numFilt * 2;
		}
	} else if (m_dctType == 1) // dct
			{
		mfcep[0] = mflogspec[0];
		for (j = 1; j < numFilt; j++)
			mfcep[0] += mflogspec[j];

		mfcep[0] = (mfcep[0] * m_sqrt_inv_n);

		for (i = 1; i < numCep; ++i) {
			mfcep[i] = 0;
			for (j = 0; j < numFilt; j++) {
				mfcep[i] += (mflogspec[j] * m_dctBank[i][j]);
			}
			mfcep[i] *= m_sqrt_inv_2n;
		}
	} else if (m_dctType == 2) //htk
			{
		mfcep[0] = mflogspec[0];
		for (j = 1; j < numFilt; j++)
			mfcep[0] += mflogspec[j];

		mfcep[0] = (mfcep[0] * m_sqrt_inv_2n); // HTK

		for (i = 1; i < numCep; ++i) {
			mfcep[i] = 0;
			for (j = 0; j < numFilt; j++) {
				mfcep[i] += (mflogspec[j] * m_dctBank[i][j]);
			}
			mfcep[i] *= m_sqrt_inv_2n;
		}
	}
//	for(int z=0; z<m_dctBank.size(); z++)cout <<output[z]<<"\t";
//	cout<<endl;
	return mfcep;
}
