/*
 * dictionary.h
 *
 *  Created on: Jul 23, 2014
 *      Author: sharath
 */

#ifndef DICTIONARY_H_
#define DICTIONARY_H_
#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "commonUtils.h"

namespace std {

class Dictionary {
public:
	map<string, string> hashMapDict;
	StringArray fillerWords;

	Dictionary() {
	}
	;
	void dictionaryInit(
			string dictPath,
			string fillerDicPath,
			bool backwardFA);
	void loadDict(
			string dictPath,
			bool backwardFA);
	void loadFillerDict(
			string fillerDicPath);
	StringArray getFillerWords();
	string getPhonemes(
			string word);
	~Dictionary() {
	}
	;
	static StringArray string2words(
			string line);
	static StringArray getFileNameAndTranscript(
			string line);
	static StringArray getFolderFileName(
			string transcriptFileLine);
	static string reverseString(
			string inString);
	static StringArray splitStringUsingDelimiter(
			string s,
			string delimiter);
};

} /* namespace std */
#endif /* DICTIONARY_H_ */
