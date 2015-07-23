/*
 * AcousticModel.h
 *
 *  Created on: Jul 23, 2014
 *      Author: sharath
 */

#ifndef ACOUSTICMODEL_H_
#define ACOUSTICMODEL_H_
#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "commonUtils.h"
#include "Dictionary.h"
#include "NeuralNetwork.h"

class AcousticModel {
public:

	Number3DArray MEAN, VAR, CLUSTERS;
	Int3DArray TMAT;
	Number2DArray PRECOMP, CODEBOOK;
	Int2DArray MIXWT;
	map<string, string> hashMapMdef;
	map<string, int> hashMapThresh;
	StringArray phoneList;
	int phoneCnt;
	int GAUSSIAN;
	int SENONES;
	int CI_STATES;
	int NO_CLUSTERS;
	int m_cepLen;
	Number m_varFloor;
	Number m_mixwFloor;
	Number m_tmatFloor;
	static Number m_logBase;
	static Number m_logBaseDenominator;
	Number m_beamWidth;

	NeuralNetwork nn;

	AcousticModel() {
	}
	;
	void acousticModelInit(
			string amDir,
			Number varFloor,
			Number mixwFloor,
			Number tmatFloor,
			Number logbase,
			Number beamWidth);
	void acousticModelInit_NN(
			string amDir,
			Number varFloor,
			Number mixwFloor,
			Number tmatFloor,
			Number logbase,
			Number beamWidth);
	void setBeamWidth(
			Number beamWidth);
	~AcousticModel() {
	}
	;

	int loadMean(
			string path);
	int loadVar(
			string path);
	int preComputeGaussianLikelihoodConstants();
	int loadTmat(
			string path);
	int loadMixWt(
			string path);
	int loadMdef(
			string path);
	IntArray getStates(
			string triPhones);
	void vectorNormalize(
			NumberArrayRef vec);
	static int myLog(
			Number val);
	void loadPhoneThresh(
			string thresholdPath);
	int getThreshold(
			string triPhones);

};

#endif /* ACOUSTICMODEL_H_ */
