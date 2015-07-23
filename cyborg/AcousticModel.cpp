/*
 * AcousticModel.cpp
 *
 *  Created on: Jul 23, 2014
 *      Author: sharath
 */

#include "AcousticModel.h"

Number AcousticModel::m_logBase = 2.71828182846;
Number AcousticModel::m_logBaseDenominator = 1;

void AcousticModel::acousticModelInit(
		string amDir,
		Number varFloor,
		Number mixwFloor,
		Number tmatFloor,
		Number logbase,
		Number beamWidth) {
	m_cepLen = 39;
	m_varFloor = varFloor;
	m_mixwFloor = mixwFloor;
	m_tmatFloor = tmatFloor;
	m_logBase = logbase;
	m_logBaseDenominator = 1.0 / log(logbase);
	m_beamWidth = log(beamWidth) * m_logBaseDenominator;

//	loadPhoneThresh(thresholdPath);
	loadMean(amDir + "//mean_bin");
	loadVar(amDir + "//variance_bin");
	preComputeGaussianLikelihoodConstants();
	loadTmat(amDir + "//tmat_bin");
	loadMixWt(amDir + "//mixWt_bin");
	loadMdef(amDir + "//mdef_tab");
}

void AcousticModel::acousticModelInit_NN(
		string amDir,
		Number varFloor,
		Number mixwFloor,
		Number tmatFloor,
		Number logbase,
		Number beamWidth) {
	m_cepLen = 39;
	m_varFloor = varFloor;
	m_mixwFloor = mixwFloor;
	m_tmatFloor = tmatFloor;
	m_logBase = logbase;
	m_logBaseDenominator = 1.0 / log(logbase);
	m_beamWidth = log(beamWidth) * m_logBaseDenominator;

	loadTmat(amDir + "//tmat_bin");
	loadMdef(amDir + "//mdef_tab");
	nn.loadNeuralNetwork(amDir + "//NN_TIFR_STATES_NORM_sigm.txt");
}

void AcousticModel::setBeamWidth(
		Number beamWidth) {
	m_beamWidth = log(beamWidth) * m_logBaseDenominator;
}
/**
 * This function takes the binary mean file given as an
 * argument and stores it in a static 3-d array 'MEAN[][][]' for fast computation.
 *
 * The format of the mean array is as follows:<br>
 * <li>First dimension - senoneID </li>
 * <li>Second dimension - GAUSSIAN </li>
 * <li>Third dimension - feature number </li> <br>
 * e.g. MEAN[500][4][38] refers to 500th senone, 4th gaussian and
 * 38th feature.<br>
 * SenoneID, gaussian and feature vectors starts from index <b>0</b>
 * and goes on till the number of senones, gaussians and
 * feature vectors <b> minus 1.</b>
 * @param path path of the binary file
 * @throws IOException
 */
int AcousticModel::loadMean(
		string filePath) {
	ifstream meanFile(filePath.c_str(), ios::binary);
	if (meanFile.is_open()) {
		meanFile.read((char*) &SENONES, sizeof(int));
		meanFile.read((char*) &GAUSSIAN, sizeof(int));
//		cout << "SENONES "<< SENONES<<endl;
//		cout << "GAUSSIAN "<< GAUSSIAN<<endl;

		MEAN = Number3DArray(SENONES, Number2DArray(GAUSSIAN, NumberArray(m_cepLen)));

		for (int i = 0; i < SENONES; i++) {
			for (int j = 0; j < GAUSSIAN; j++) {
				meanFile.read((char*) &MEAN[i][j][0], sizeof(Number) * m_cepLen);
			}

		}
		meanFile.close();
	} else {
		cout << "ERROR: Unable to open : " << filePath << endl;
		exit(0);
	}

	return 0;
}

/**
 * * This function takes the binary variance file given as an
 * argument and stores it in a static 3-d array 'VAR[][][]' for
 *  fast computation.
 *
 * The format of the variance array is as follows:<br>
 * <li>First dimension - senoneID </li>
 * <li>Second dimension - gaussian </li>
 * <li>Third dimension - feature number </li> <br>
 * e.g. VAR[500][4][38] refers to 500th state, 4th gaussian and
 * 38th feature.<br>
 * senoneID, gaussian and feature vectors starts from index <b>0</b>
 * and goes on till the number of senones, gaussians and
 * feature vectors <b> minus 1.</b>
 * @param path path of the binary file
 * @throws IOException
 */

int AcousticModel::loadVar(
		string filePath) {
	ifstream varFile(filePath.c_str(), ios::binary);
	if (varFile.is_open()) {

		varFile.read((char*) &SENONES, sizeof(int));
		varFile.read((char*) &GAUSSIAN, sizeof(int));
//		cout << "SENONES "<< SENONES<<endl;
//		cout << "GAUSSIAN "<< GAUSSIAN<<endl;

		VAR = Number3DArray(SENONES, Number2DArray(GAUSSIAN, NumberArray(m_cepLen)));
		for (int i = 0; i < SENONES; i++) {
			for (int j = 0; j < GAUSSIAN; j++) {
				varFile.read((char*) &VAR[i][j][0], sizeof(Number) * m_cepLen);
				for (int k = 0; k < m_cepLen; k++) {
					if (VAR[i][j][k] < m_varFloor) {
						VAR[i][j][k] = m_varFloor;
					}
//					cout << VAR[i][j][k]<<" ";

				}
//				cout << endl;
			}
		}
		varFile.close();
	} else {
		cout << "ERROR: Unable to open : " << filePath << endl;
		exit(0);
	}
	return 0;
}

int AcousticModel::preComputeGaussianLikelihoodConstants() {
	PRECOMP = Number2DArray(SENONES, NumberArray(GAUSSIAN));
	int m, c, i;
	for (m = 0; m < SENONES; m++) {
		for (c = 0; c < GAUSSIAN; c++) {

			Number lrd = 0;
			for (i = 0; i < (int) VAR[m][c].size(); i++) {

				Number val = (Number) ((VAR[m][c][i] == 0.0) ? LOGPROB_ZERO : log(VAR[m][c][i]));
				lrd += val;

				/* Precompute this part of the exponential */
				VAR[m][c][i] = (1 / (VAR[m][c][i] * 2));
			}

			lrd += VAR[m][c].size() * log(2.0 * M_PI); /* (2pi)^velen */
			PRECOMP[m][c] = -0.5 * lrd; /* Reciprocal, sqrt */
//            cout << PRECOMP[m][c] << " ";
		}
//        cout <<endl;
	}

	return 0;
}

/**
 * This function takes the binary transition matrix file given
 * as an argument and stores it in a static 3-d array 'TMAT[][][]'
 * for fast computation.
 *
 * The format of the tmat array is as follows:<br>
 * <li>First dimension - phone ID </li>
 * <li>Second dimension - transition from </li>
 * <li>Third dimension - transition to.</li> <br>
 * e.g. TMAT[60][1][2] refers to 60th CI phone with a transition
 * from first to  second state.<br>
 * The assumption is that the HMM cant skip state. <br>
 * So TMAT[60][1][3], TMAT[60][0][2] will be zero and so on.<br>
 * State, transition from and transition to starts from index <b>0</b>
 * and goes on till the number of CI phones and the number of
 * states per senone/CI phone.
 * @param path path of the binary file
 * @throws IOException
 */
int AcousticModel::loadTmat(
		string filePath) {
	ifstream tMatFile(filePath.c_str(), ios::binary);
	if (tMatFile.is_open()) {
		int states_per_triphone = 0;
		tMatFile.read((char*) &CI_STATES, sizeof(int));
		tMatFile.read((char*) &states_per_triphone, sizeof(int));
//		cout << "CI_STATES "<< CI_STATES<<endl;
//		cout << "states_per_triphone "<< states_per_triphone<<endl;

		TMAT = Int3DArray(CI_STATES, Int2DArray(states_per_triphone, IntArray(states_per_triphone)));
		for (int i = 0; i < CI_STATES; i++) {
			for (int j = 0; j < states_per_triphone - 1; j++) {
				NumberArray tMatTmp(states_per_triphone, 0.0);
				for (int k = j; k <= j + 1; k++) {
					tMatFile.read((char*) &tMatTmp[k], sizeof(Number));
					if (tMatTmp[k] < m_tmatFloor) {
						tMatTmp[k] = m_tmatFloor;
					}
				}
				vectorNormalize(tMatTmp);
				for (int k = j; k <= j + 1; k++) {
					TMAT[i][j][k] = myLog(tMatTmp[k]);
//					cout << TMAT[i][j][k]<<" ";
				}
//				cout << endl;

			}
		}
		tMatFile.close();
	} else {
		cout << "ERROR: Unable to open : " << filePath << endl;
		exit(0);
	}
	return 0;
}

/**
 * This function takes the binary mixture weight file given
 * as an argument and stores it in a static 2-d array 'MIXWT[][]'
 * for fast computation.
 *
 * The format of the MIXWT array is as follows:<br>
 * <li>First dimension - senoneID;</li>
 * <li>Second dimension - Gaussian mixture number; </li><br>
 * e.g. MIXWT[400][8] refers to 400th senone and 8th gaussian.<br>
 * SenoneID and gaussian mix no. starts from index <b>0</b>
 * and goes on till the number of senones and gaussians
 * minus 1.</b>
 * @param path
 * @throws IOException
 */
int AcousticModel::loadMixWt(
		string filePath) {
	ifstream mixWtFile(filePath.c_str(), ios::binary);
	if (mixWtFile.is_open()) {

		mixWtFile.read((char*) &SENONES, sizeof(int));
		mixWtFile.read((char*) &GAUSSIAN, sizeof(int));
//		cout << "SENONES "<< SENONES<<endl;
//		cout << "GAUSSIAN "<< GAUSSIAN<<endl;

		MIXWT = Int2DArray(SENONES, IntArray(GAUSSIAN));
		for (int i = 0; i < SENONES; i++) {
			NumberArray mixWtTmp(GAUSSIAN, 0.0);
			mixWtFile.read((char*) &mixWtTmp[0], sizeof(Number) * GAUSSIAN);
			for (int j = 0; j < GAUSSIAN; j++) {
				if (mixWtTmp[j] < m_mixwFloor) {
					mixWtTmp[j] = m_mixwFloor;
				}
			}

			vectorNormalize(mixWtTmp);
			for (int j = 0; j < GAUSSIAN; j++) {
				MIXWT[i][j] = myLog(mixWtTmp[j]);
//				cout << MIXWT[i][j]<<" ";
			}
//			cout << endl;
		}
		mixWtFile.close();
	} else {
		cout << "ERROR: Unable to open : " << filePath << endl;
		exit(0);
	}
	return 0;
}

int AcousticModel::loadMdef(
		string filePath) {

	hashMapMdef = map<string, string>();
	ifstream mdefFile(filePath.c_str(), ios::binary);
	string line;
	bool startProcessingLines = false;
	if (mdefFile.is_open()) {
		while (getline(mdefFile, line)) {
			if (startProcessingLines) {
				StringArray tokens = Dictionary::string2words(line);
				if (tokens.size() != 10) {
					cout << "ERROR: Number of row elements in mdef should be 10, only " << tokens.size()
							<< " elements found. " << endl;
					return -1;
					break;
				}
				string triphone, states;
				triphone = tokens[0] + "\t" + tokens[1] + "\t" + tokens[2] + "\t" + tokens[3];
				states = tokens[4] + "\t" + tokens[5] + "\t" + tokens[6] + "\t" + tokens[7] + "\t" + tokens[8] + "\t"
						+ tokens[9];
//			    cout << triphone<< " = " << states<<endl;
				hashMapMdef[triphone] = states;
				phoneList.push_back(tokens[0]);
//				cout << phoneList.back()<< " ";
			}

			if (line.find("#base") == 0)
				startProcessingLines = true;

		}
		mdefFile.close();
		phoneCnt = hashMapMdef.size();
//		cout <<endl<<"Number of Phones: " <<phoneCnt <<endl;
	} else {
		cout << "ERROR: " << filePath << " file does not exist." << endl;
		exit(0);
	}
	return 0;
}

/**
 *
 * @param String of triPhones where each phone is separated by tab ('\t') <br> eg : a	SIL	d'	b
 * @return Integer array containing corresponding attrib, tmat and states of triPhone <br>
 * eg: <br>
 * states[0] = 1    // attrib[0: filler and 1: n/a] <br>
 * states[1] = 11   // tmat <br>
 * states[2] = 265  // state 0 <br>
 * states[3] = 296  // state 1 <br>
 * states[4] = 324  // state 2 <br>
 *
 * state 0, 1, 2 are states of tri-phone 'a SIL b'
 * @throws Exception
 */

IntArray AcousticModel::getStates(
		string triPhones) {

	IntArray states(5, 0);
	StringArray statesString;

	if (hashMapMdef.find(triPhones) == hashMapMdef.end()) {
//			cout << "ERROR : "<< triPhones <<" : TriPhone doesn't exist in 'mdef' training data." <<endl;
		return states;
//			exit(0);
	} else {
		statesString = Dictionary::string2words(hashMapMdef[triPhones]);
	}

	states[0] = (statesString[0].compare("filler") == 0) ? 0 : 1;
	states[1] = atoi(statesString[1].c_str());
	states[2] = atoi(statesString[2].c_str());
	states[3] = atoi(statesString[3].c_str());
	states[4] = atoi(statesString[4].c_str());

	return states;
}

void AcousticModel::vectorNormalize(
		NumberArrayRef vec) {
	Number sum = 0, f;
	int i;
	for (i = 0; i < (int) vec.size(); i++)
		sum += vec[i];

	if (sum != 0.0) {
		f = 1 / sum;
		for (i = 0; i < (int) vec.size(); i++)
			vec[i] *= f;
	}
}

int AcousticModel::myLog(
		Number val) {
	return (val == 0.0) ? LOGPROB_ZERO : (int) (log(val) * m_logBaseDenominator);
}

void AcousticModel::loadPhoneThresh(
		string thresholdPath) {
	hashMapThresh = map<string, int>();
	ifstream thresholdFile(thresholdPath.c_str());
	string line;
	if (thresholdFile.is_open()) {
		while (getline(thresholdFile, line)) {
			istringstream iss(line);
			string phone, threshold;
			iss >> phone;
			phone.erase(std::remove(phone.begin(), phone.end(), '\t'), phone.end()); // remove tab if any

			ostringstream oss;
			oss << iss.rdbuf();
			threshold = oss.str();
			hashMapThresh[phone] = atoi(threshold.c_str());
//			cout << phone<< "  " << hashMapThresh[phone]<<endl;

		}
		thresholdFile.close();
	} else {
		cout << "ERROR: " << thresholdPath << " file does not exist." << endl;
		exit(0);
	}
}

int AcousticModel::getThreshold(
		string triPhones) {

	if (hashMapThresh.find(triPhones) == hashMapThresh.end()) {
		cout << "ERROR : " << triPhones << " : TriPhone doesn't exist in threshold data." << endl;
		exit(0);
	} else {
		return hashMapThresh[triPhones];
	}
}
