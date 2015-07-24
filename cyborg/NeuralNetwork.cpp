/*
 * NeuralNetwork.cpp
 *
 *  Created on: May 20, 2015
 *      Author: sharath
 */

#include "NeuralNetwork.h"

namespace std {

NeuralNetwork::NeuralNetwork(
		IntArrayRef architecture) {
//	NNSETUP creates a Feedforward Backpropagate Neural Network
//	 nn = nnsetup(architecture) returns an neural network structure with n=numel(architecture)
//	 layers, architecture being a n x 1 vector of layer sizes e.g. [784 100 10]

	size = architecture;
	n = architecture.size();

	activation_function = "tanh_opt"; //  Activation functions of hidden layers: 'sigm' (sigmoid) or 'tanh_opt' (optimal tanh). Use sigm if data is scaled from [0 1] and tanh for data [-1 1]
	learningRate = 2; //  learning rate Note: typically needs to be lower when using 'sigm' activation function and non-normalized inputs.
	momentum = 0.5;          //  Momentum
	scaling_learningRate = 1;            //  Scaling factor for the learning rate (each epoch)
	weightPenaltyL2 = 0;            //  L2 regularization
	nonSparsityPenalty = 0;            //  Non sparsity penalty
	sparsityTarget = 0.05;         //  Sparsity target
	inputZeroMaskedFraction = 0;            //  Used for Denoising AutoEncoders
	dropoutFraction = 0;            //  Dropout level (http://www.cs.toronto.edu/~hinton/absps/dropout.pdf)
	testing = 0;            //  Internal variable. nntest sets this to one.
	output = "sigm";       //  output unit 'sigm' (=logistic), 'softmax' and 'linear'

	for (int i = 1; i < n; i++) {
		// weights and weight momentum
		Number2DArray randomMatrix;
		randomMatrixGenerator(size[i], size[i - 1] + 1, randomMatrix);
		Number constMultiplier = 2.0 * 4 * sqrt(6.0 / (size[i] + size[i - 1]));

		for (int m = 0; m < (int) randomMatrix.size(); m++) {
			for (int n = 0; n < (int) randomMatrix[m].size(); n++) {
				randomMatrix[m][n] = (randomMatrix[m][n] - 0.5) * constMultiplier;
			}
		}
		W.push_back(randomMatrix);
		vW.push_back(Number2DArray(size[i], NumberArray(size[i - 1] + 1, 0)));
		// average activations (for use with sparsity)
		p.push_back(NumberArray(size[i], 0));
	}
}

int NeuralNetwork::loadNeuralNetwork(
		string nnFileName) {

	ifstream nnFile(nnFileName.c_str());
	if (nnFile.is_open()) {
		string line;
		StringArray words;

		//size
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("size") == 0) {
			for (int i = 1; i < (int) words.size(); i++) {
				size.push_back(atoi(words[i].c_str()));
			}
		} else {
			LOGE("Incorrect NN model format");
		}

		//n
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("n") == 0) {
			n = atoi(words[1].c_str());
		} else {
			LOGE("Incorrect NN model format");
		}

		//activation_function
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("activation_function") == 0) {
			activation_function = words[1];
		} else {
			LOGE("Incorrect NN model format");
		}

		//learningRate
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("learningRate") == 0) {
			learningRate = atof(words[1].c_str());
		} else {
			LOGE("Incorrect NN model format");
		}

		//momentum
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("momentum") == 0) {
			momentum = atof(words[1].c_str());
		} else {
			LOGE("Incorrect NN model format");
		}

		//scaling_learningRate
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("scaling_learningRate") == 0) {
			scaling_learningRate = atof(words[1].c_str());
		} else {
			LOGE("Incorrect NN model format");
		}

		//weightPenaltyL2
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("weightPenaltyL2") == 0) {
			weightPenaltyL2 = atof(words[1].c_str());
		} else {
			LOGE("Incorrect NN model format");
		}

		//nonSparsityPenalty
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("nonSparsityPenalty") == 0) {
			nonSparsityPenalty = atof(words[1].c_str());
		} else {
			LOGE("Incorrect NN model format");
		}

		//sparsityTarget
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("sparsityTarget") == 0) {
			sparsityTarget = atof(words[1].c_str());
		} else {
			LOGE("Incorrect NN model format");
		}

		//inputZeroMaskedFraction
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("inputZeroMaskedFraction") == 0) {
			inputZeroMaskedFraction = atof(words[1].c_str());
		} else {
			LOGE("Incorrect NN model format");
		}

		//dropoutFraction
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("dropoutFraction") == 0) {
			dropoutFraction = atof(words[1].c_str());
		} else {
			LOGE("Incorrect NN model format");
		}

		//testing
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("testing") == 0) {
			testing = atoi(words[1].c_str());
		} else {
			LOGE("Incorrect NN model format");
		}

		//output
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("output") == 0) {
			output = words[1];
		} else {
			LOGE("Incorrect NN model format");
		}

		//W
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("W") == 0) {
			for (int i = 1; i < n; i++) {
				// weights and weight momentum
				Number2DArray tmpMatrix(size[i], NumberArray((size[i - 1] + 1), 0.0));
				for (int m = 0; m < size[i]; m++) {
					getline(nnFile, line);
					words = string2words(line);
					for (int n = 0; n < (size[i - 1] + 1); n++) {
						tmpMatrix[m][n] = atof(words[n].c_str());
					}
				}
				W.push_back(tmpMatrix);
			}

		} else {
			LOGE("Incorrect NN model format");
		}

		//vW
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("vW") == 0) {
			for (int i = 1; i < n; i++) {
				// weights and weight momentum
				Number2DArray tmpMatrix(size[i], NumberArray((size[i - 1] + 1), 0.0));
				for (int m = 0; m < size[i]; m++) {
					getline(nnFile, line);
					words = string2words(line);
					for (int n = 0; n < (size[i - 1] + 1); n++) {
						tmpMatrix[m][n] = atof(words[n].c_str());
					}
				}
				vW.push_back(tmpMatrix);
			}

		} else {
			LOGE("Incorrect NN model format");
		}

		//p
		getline(nnFile, line);
		words = string2words(line);
		if (words[0].compare("p") == 0) {
			for (int i = 1; i < n; i++) {
				// weights and weight momentum
				NumberArray tmpMatrix(size[i], 0.0);
				getline(nnFile, line);
				words = string2words(line);
				for (int n = 0; n < size[i]; n++) {
					tmpMatrix[n] = atof(words[n].c_str());
				}
				p.push_back(tmpMatrix);
			}

		} else {
			LOGE("Incorrect NN model format");
		}

		nnFile.close();
	} else {
		LOGE("File not found : %s", nnFileName.c_str());
		return -1;
	}

	return 0;

}

void NeuralNetwork::nnff(
		Number2DArrayRef featMat,
		Number2DArrayRef dummy) {
	//NNFF performs a feedforward pass
	// nn = nnff(nn, x, y) returns an neural network structure with updated
	// layer activations, error and loss (nn.a, nn.e and nn.L)

	for (int i = 0; i < a.size(); i++) {
		for (int j = 0; j < a[i].size(); j++) {
			a[i][j].clear();
		}
		a[i].clear();
	}
	a.clear();

	Number2DArray tmpMatrix(featMat.size(), NumberArray((featMat[0].size() + 1), 1.0));
	for (int i = 0; i < (int) featMat.size(); i++) {
		for (int j = 0; j < (int) featMat[0].size(); j++) {
			tmpMatrix[i][j + 1] = featMat[i][j];
		}
	}
	a.push_back(tmpMatrix);

	//feedforward pass
	for (int i = 1; i < (n - 1); i++) {
		Number2DArray outMatMult, updatedA, transposeMat;
		if (activation_function.compare("sigm") == 0) {
			// Calculate the unit's outputs (including the bias term)
			matrixTranspose(W[i - 1], transposeMat);
			matrixMult(a[i - 1], transposeMat, outMatMult);
			sigm(outMatMult, updatedA);

		} else if (activation_function.compare("tanh_opt") == 0) {
			matrixTranspose(W[i - 1], transposeMat);
			matrixMult(a[i - 1], transposeMat, outMatMult);
			tanh_opt(outMatMult, updatedA);

		} else {
			LOGE("NeuralNetwork::nnff Unknown activation function.");
		}

		//dropout
		if (dropoutFraction > 0) {
			if (testing) {
				matrixConstMult(a[i], (1 - dropoutFraction));
			}
			// TODO implement dropoutFraction for training
		}

		//Add the bias term
		for (int m = 0; m < (int) updatedA.size(); m++) {
			updatedA[m].insert(updatedA[m].begin(), 1.0);
		}
		a.push_back(updatedA);
	}

	Number2DArray outMatMult, updatedA, transposeMat;
	if (output.compare("sigm") == 0) {
		matrixTranspose(W[n - 2], transposeMat);
		matrixMult(a[n - 2], transposeMat, outMatMult);
		sigm(outMatMult, updatedA);
		a.push_back(updatedA);
	} else if (output.compare("linear") == 0) {
		matrixTranspose(W[n - 2], transposeMat);
		matrixMult(a[n - 2], transposeMat, outMatMult);
		a.push_back(outMatMult);

	} else if (output.compare("softmax") == 0) {
		matrixTranspose(W[n - 2], transposeMat);
		matrixMult(a[n - 2], transposeMat, outMatMult);
		NumberArray maxVec, sumVec;
		matrixMax(outMatMult, 0, maxVec);
		matrixVectorMinus(outMatMult, maxVec);
		matrixExp(outMatMult);
		matrixSum(outMatMult, 0, sumVec);
		matrixVectorDivide(outMatMult, sumVec);
		a.push_back(outMatMult);
	} else {
		LOGE("NeuralNetwork::nnff Unknown output layer");
	}

}

void NeuralNetwork::matrixVectorMinus(
		Number2DArrayRef inOutMatMult,
		NumberArrayRef maxVec) {

	if (inOutMatMult.size() != maxVec.size()) {
		LOGE("NeuralNetwork::matrixVectorMinus Size mismatch");
	}

	for (int i = 0; i < (int) inOutMatMult.size(); i++) {
		for (int j = 0; j < (int) inOutMatMult[i].size(); j++) {
			inOutMatMult[i][j] -= maxVec[i];
		}
	}
}

void NeuralNetwork::matrixVectorDivide(
		Number2DArrayRef inOutMatMult,
		NumberArrayRef sumVec) {

	if (inOutMatMult.size() != sumVec.size()) {
		LOGE("NeuralNetwork::matrixVectorMinus Size mismatch");
	}

	for (int i = 0; i < (int) inOutMatMult.size(); i++) {
		for (int j = 0; j < (int) inOutMatMult[i].size(); j++) {
			inOutMatMult[i][j] /= sumVec[i];
		}
	}
}

void NeuralNetwork::matrixMax(
		Number2DArrayRef P,
		int dim,
		NumberArrayRef maxVec) {
	//dim =0 max of rows
	//       =1 max of col
	if (dim == 0) {
		maxVec.clear();
		maxVec.resize(P.size(), 0);
		for (int i = 0; i < (int) P.size(); i++) {
			Number maxVal = P[i][0];
			for (int j = 1; j < (int) P[i].size(); j++) {
				if (P[i][j] > maxVal) {
					maxVal = P[i][j];
				}
			}
			maxVec[i] = maxVal;
		}
	} else {
		maxVec.clear();
		maxVec.resize(P[0].size(), 0);
		for (int i = 0; i < (int) P[0].size(); i++) {
			Number maxVal = P[0][i];
			for (int j = 1; j < (int) P.size(); j++) {
				if (P[j][i] > maxVal) {
					maxVal = P[j][i];
				}
			}
			maxVec[i] = maxVal;
		}
	}
}

void NeuralNetwork::matrixMean(
		Number2DArrayRef P,
		int dim,
		NumberArrayRef meanVec) {
	//dim =0 max of rows
	//       =1 max of col
	if (dim == 0) {
		meanVec.clear();
		meanVec.resize(P.size(), 0);
		for (int i = 0; i < (int) P.size(); i++) {
			double sumVal = 0;
			for (int j = 0; j < (int) P[i].size(); j++) {
				sumVal += P[i][j];
			}
			meanVec[i] = sumVal / P[i].size();
		}
	} else {
		meanVec.clear();
		meanVec.resize(P[0].size(), 0);
		for (int i = 0; i < (int) P[0].size(); i++) {
			double sumVal = 0;
			for (int j = 0; j < (int) P.size(); j++) {
				sumVal += P[j][i];
			}
			meanVec[i] = sumVal / P.size();
		}
	}
}

void NeuralNetwork::matrixStd(
		Number2DArrayRef P,
		int dim,
		NumberArrayRef stdVec) {
	//dim =0 max of rows
	//     =1 max of col
	NumberArray meanVec;
	matrixMean(P, dim, meanVec);
	if (dim == 0) {
		stdVec.clear();
		stdVec.resize(P.size(), 0);
		for (int i = 0; i < (int) P.size(); i++) {
			double sumVal = 0;
			for (int j = 0; j < (int) P[i].size(); j++) {
				sumVal += ((P[i][j] - meanVec[i]) * (P[i][j] - meanVec[i]));
			}
			stdVec[i] = sqrt(sumVal);//TODO the std eqn is sum(values-mean(values))/len(values), but we are not using len(values) because NN has not been trained without it. Which is a mistake, change this again when NN is trained with the scaling factor.
		}
	} else {
		stdVec.clear();
		stdVec.resize(P[0].size(), 0);
		for (int i = 0; i < (int) P[0].size(); i++) {
			double sumVal = 0;
			for (int j = 0; j < (int) P.size(); j++) {
				sumVal += ((P[j][i] - meanVec[i]) * (P[j][i] - meanVec[i]));
			}
			stdVec[i] = sqrt(sumVal);
		}
	}
}

void NeuralNetwork::matrixSum(
		Number2DArrayRef P,
		int dim,
		NumberArrayRef sumVec) {
	//dim =0 max of rows
	//       =1 max of col
	if (dim == 0) {
		sumVec.clear();
		sumVec.resize(P.size(), 0);
		for (int i = 0; i < (int) P.size(); i++) {
			double sumVal = 0;
			for (int j = 0; j < (int) P[i].size(); j++) {
				sumVal += P[i][j];
			}
			sumVec[i] = sumVal;
		}
	} else {
		sumVec.clear();
		sumVec.resize(P[0].size(), 0);
		for (int i = 0; i < (int) P[0].size(); i++) {
			double sumVal = 0;
			for (int j = 0; j < (int) P.size(); j++) {
				sumVal += P[j][i];
			}
			sumVec[i] = sumVal;
		}
	}
}

void NeuralNetwork::classify(
		Number2DArrayRef featMat,
		Number2DArrayRef outPosteriors) {
	Number2DArray dummy;
	testing = 1;
	nnff(featMat, dummy);
	testing = 0;
	outPosteriors = a[n - 1];

}

void NeuralNetwork::classify(
		NumberArrayRef featMat,
		NumberArrayRef outPosteriors) {
	Number2DArray dummy;
	testing = 1;
	Number2DArray locFeatMat;
	locFeatMat.push_back(featMat);
	nnff(locFeatMat, dummy);
	testing = 0;
	outPosteriors = a[n - 1][0];

}

void NeuralNetwork::tanh_opt(
		Number2DArrayRef A,
		Number2DArrayRef f) {
	f.resize(A.size(), NumberArray(A[0].size(), 0));
	for (int i = 0; i < (int) A.size(); i++) {
		for (int j = 0; j < (int) A[i].size(); j++) {
			f[i][j] = 1.7159 * (tanh(A[i][j] * (2 / 3.0)));
		}
	}

}

void NeuralNetwork::sigm(
		Number2DArrayRef P,
		Number2DArrayRef X) {
	X.resize(P.size(), NumberArray(P[0].size(), 0));
	for (int i = 0; i < (int) P.size(); i++) {
		for (int j = 0; j < (int) P[i].size(); j++) {
			X[i][j] = 1.0 / (1 + exp(-P[i][j]));
		}
	}

}

void NeuralNetwork::matrixExp(
		Number2DArrayRef P) {
	for (int i = 0; i < (int) P.size(); i++) {
		for (int j = 0; j < (int) P[i].size(); j++) {
			P[i][j] = exp(P[i][j]);
		}
	}
}

void NeuralNetwork::matrixTranspose(
		Number2DArrayRef P,
		Number2DArrayRef transposeMat) {
	transposeMat.resize(P[0].size(), NumberArray(P.size(), 0.0));
	for (int i = 0; i < (int) transposeMat.size(); i++) {
		for (int j = 0; j < (int) transposeMat[i].size(); j++) {
			transposeMat[i][j] = P[j][i];
		}
	}
}

void NeuralNetwork::matrixConstMult(
		Number2DArrayRef P,
		Number constNumber) {
	for (int i = 0; i < (int) P.size(); i++) {
		for (int j = 0; j < (int) P[i].size(); j++) {
			P[i][j] *= constNumber;
		}
	}
}

void NeuralNetwork::matrixMult(
		Number2DArrayRef first,
		Number2DArrayRef second,
		Number2DArrayRef multiply) {

	int m = first.size();
	int n = first[0].size();
	int q = second.size();
	int p = second[0].size();

	if (n != q) {
		LOGE("NeuralNetwork::matrixMult Matrix size mismatch");
		exit(-1);
	}
	multiply.resize(m, NumberArray(p, 0.0));

	for (int c = 0; c < m; c++) {
		for (int d = 0; d < p; d++) {
			Number sum = 0;
			for (int k = 0; k < q; k++) {
				sum = sum + (first[c][k] * second[k][d]);
			}
			multiply[c][d] = sum;
		}
	}
}

NeuralNetwork::~NeuralNetwork() {
}

void NeuralNetwork::randomMatrixGenerator(
		int row,
		int col,
		Number2DArrayRef outRandomMatrix) {

	srand(static_cast<unsigned>(time(0)));
	/* generate number: */
	for (int i = 0; i < row; i++) {
		NumberArray tmp;
		for (int j = 0; j < col; j++) {
			tmp.push_back((rand() & 255) / 256.0);
		}
		outRandomMatrix.push_back(tmp);
	}
}

StringArray NeuralNetwork::string2words(
		string line) {
	StringArray tokens;
	istringstream iss(line);
	copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<StringArray>(tokens));
	return tokens;
}

void NeuralNetwork::readPosteriors(string posteriorpath, Number2DArrayRef posterior)
 {
	//cout << "In Read Posteriors: "<<endl;
	//cout << posteriors << " ";
	ifstream posteriorfile(posteriorpath.c_str());
		if (!posteriorfile.is_open()) {
			cout << "ERROR: Unable to open : " << posteriorpath << endl;
			exit(0);
		}
		//cout << "Posterior file read";
		string posteriorfileLine;
		//Number2DArray posterior;
		NumberArray temp;
		while (getline(posteriorfile, posteriorfileLine)) { // Iterate through filenames in fileid file
				string delimiter1 = " ";
				int pos = 0;
				string token;
				while ((pos = posteriorfileLine.find(delimiter1)) != std::string::npos) {
				        //cout << pos<<endl;
						token = posteriorfileLine.substr(0, pos);
				        posteriorfileLine.erase(0, pos + delimiter1.length());
				        //cout << "Posterior erased";
				        temp.push_back(stof(token));
				}
				posterior.push_back(temp);
				temp.clear();
				//cout << posterior[0]<<endl;
		}
		//return posterior;
		cout << "Posteriors read!!";

 }
} /* namespace std */
