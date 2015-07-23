/*
 * Fst.h
 *
 *  Created on: Jul 24, 2014
 *      Author: sharath
 */

#ifndef FST_H_
#define FST_H_
#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "commonUtils.h"
#include "Dictionary.h"
namespace std {

class Fst {
public:
	//lexiconPos 0->1 are reserved for silence and fillers
	int lexiconPos;
	int symsPos;
	ofstream wordGraph;
	ofstream lexicon;
	ofstream syms;
	map<string, string> hashMapSyms;
	string m_fstFolderPath;

	Fst(
			string fstFolderPath,
			string fileName,
			string transcription,
			Dictionary dict,
			bool backwardFA,
			bool decoderMode);
	~Fst() {
	}
	;
	int getFstTxtFilesMultipleTrans(
			string wfstFile,
			string transcription,
			Dictionary dict);
	int getFstTxtFiles(
			string transcription,
			Dictionary dict);
	void printSilencesAndFillers(
			StringArray fillerWords);
	void printSyms(
			string c);
	void printWordGraph(
			int pos,
			string c);
	void printWordGraphFromTo(
			int pos,
			int nextPos,
			string c);
	void printFillerInWordGraph(
			int& pos,
			StringArray fillerWords);
	void createHash(
			string word,
			string phoneme);
	void writeLexicon(
			string phoneme,
			string word);
	void printLexicon(
			int pos,
			int nextPos,
			string c,
			string d);
	void printSilAndFillersInLexicon(
			StringArray fillerWords);
	bool isFiller(
			string word,
			StringArray fillerWords);
	StringArray removeAllFillers(
			StringArray words,
			StringArray fillerWords);
};

} /* namespace std */
#endif /* FST_H_ */
