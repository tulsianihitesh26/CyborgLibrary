/*
 * NeuralNetwork.h
 *
 *  Created on: May 20, 2015
 *      Author: sharath
 */

#ifndef NEURALNETWORK_H_
#define NEURALNETWORK_H_
#include "commonUtils.h"
#include <iostream>
#include <ctime>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cmath>
namespace std {

class NeuralNetwork {

private:
	int n;							//Number of layers
	IntArray size; 					//Architecture

	string activation_function; //Activation functions of hidden layers: 'sigm' (sigmoid) or 'tanh_opt' (optimal tanh). Use sigm if data is scaled from [0 1] and tanh for data [-1 1]
	Number learningRate; //learning rate Note: typically needs to be lower when using 'sigm' activation function and non-normalized inputs.
	Number momentum; //Momentum
	Number scaling_learningRate; //Scaling factor for the learning rate (each epoch)
	Number weightPenaltyL2; //L2 regularization
	Number nonSparsityPenalty; //Non sparsity penalty
	Number sparsityTarget; //Sparsity target
	Number inputZeroMaskedFraction;	//Used for Denoising AutoEncoders
	Number dropoutFraction;	//Dropout level (http://www.cs.toronto.edu/~hinton/absps/dropout.pdf)
	int testing;	//Internal variable. nntest sets this to one.
	string output;	//output unit 'sigm' (=logistic), 'softmax' and 'linear'

	Number3DArray W;
	Number3DArray vW;
	Number2DArray p;
	Number3DArray a;
	Number2DArray e;
	Number L;
	Number3DArray dW;

public:
	NeuralNetwork() {
	}
	;
	NeuralNetwork(
			IntArrayRef architecture);
	int loadNeuralNetwork(
			string nnFileName);
	virtual ~NeuralNetwork();

	void NeuralNetwork::readPosteriors(
			string posteriorpath,
			Number2DArrayRef posterior);

	void classify(
			Number2DArrayRef featMat,
			Number2DArrayRef outPosteriors);

	void classify(
			NumberArrayRef featMat,
			NumberArrayRef outPosteriors);

	void nnff(
			Number2DArrayRef featMat,
			Number2DArrayRef dummy);

	void matrixVectorMinus(
			Number2DArrayRef outMatMult,
			NumberArrayRef maxVec);

	void matrixVectorDivide(
			Number2DArrayRef inOutMatMult,
			NumberArrayRef sumVec);

	void matrixMax(
			Number2DArrayRef P,
			int dim,
			NumberArrayRef maxVec);

	static void matrixMean(
			Number2DArrayRef P,
			int dim,
			NumberArrayRef meanVec);

	static void matrixStd(
			Number2DArrayRef P,
			int dim,
			NumberArrayRef stdVec);

	void matrixSum(
			Number2DArrayRef P,
			int dim,
			NumberArrayRef sumVec);

	void randomMatrixGenerator(
			int row,
			int col,
			Number2DArrayRef outRandomMatrix);
	static StringArray string2words(
			string line);

	void tanh_opt(
			Number2DArrayRef A,
			Number2DArrayRef f);

	void sigm(
			Number2DArrayRef P,
			Number2DArrayRef X);

	void matrixExp(
			Number2DArrayRef P);

	void matrixTranspose(
			Number2DArrayRef P,
			Number2DArrayRef transposeMat);

	void matrixConstMult(
			Number2DArrayRef P,
			Number constNumber);

	void matrixMult(
			Number2DArrayRef first,
			Number2DArrayRef second,
			Number2DArrayRef multiply);

	};

} /* namespace std */
#endif /* NEURALNETWORK_H_ */
