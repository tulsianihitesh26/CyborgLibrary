/*
 * dictionary.cpp
 *
 *  Created on: Jul 23, 2014
 *      Author: sharath
 */

#include "Dictionary.h"

namespace std {

void Dictionary::dictionaryInit(
		string dictPath,
		string fillerDicPath,
		bool backwardFA) {
	loadDict(dictPath, backwardFA);
	loadFillerDict(fillerDicPath);
}

/**
 * Reads a Dictionary file and stores it into HashMap. Word and phoneme sequence assumed to be separated by tab (\t).<br>
 *
 * <b>Line format</b> : word	phonemes
 * <b>eg</b>: ahamadanagara 	a h a m a d n a g a r
 *
 * <b>Hash Entry example</b>
 * <blockquote>
 * <b>Key:</b> ahamadanagara<br>
 * <b>Value:</b> a h a m a d n a g a r<br>
 * </blockquote>
 */

void Dictionary::loadDict(
		string dictPath,
		bool backwardFA) {
	hashMapDict = map<string, string>();
	ifstream dicFile(dictPath.c_str(), ios::binary);
	string line;
	if (dicFile.is_open()) {
		while (getline(dicFile, line)) {
			istringstream iss(line);
			string phones, word;
			iss >> word;

			ostringstream oss;
			oss << iss.rdbuf();
			phones = oss.str();
			phones.erase(std::remove(phones.begin(), phones.end(), '\t'), phones.end()); // remove tab if any
			phones.erase(std::remove(phones.begin(), phones.end(), '\n'), phones.end()); // remove newline if any
			phones.erase(std::remove(phones.begin(), phones.end(), '\r'), phones.end()); // remove newline if any
			//		cout << word<< "  " << phones<<endl;
			if (backwardFA) {
				hashMapDict[word] = reverseString(phones);
			} else {
				hashMapDict[word] = phones;
			}

		}
		dicFile.close();
	} else {
		cout << "ERROR: dictPath : " << dictPath << " file does not exist." << endl;
		exit(0);
	}
}

void Dictionary::loadFillerDict(
		string fillerDicPath) {

	ifstream fillerDict(fillerDicPath.c_str());
	string line;
	if (fillerDict.is_open()) {
		while (getline(fillerDict, line)) {
			fillerWords.push_back(string2words(line)[0]);
		}
		fillerDict.close();
	} else {
		cout << "ERROR: fillerDicPath: " << fillerDicPath << " file does not exist." << endl;
		exit(0);
	}
}

StringArray Dictionary::getFillerWords() {
	return fillerWords;
}
string Dictionary::getPhonemes(
		string word) {

	if (hashMapDict.find(word) == hashMapDict.end()) {
		return "NULL";
	} else
		return hashMapDict[word];
}

StringArray Dictionary::string2words(
		string line) {
	StringArray tokens;
	istringstream iss(line);
	copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<StringArray>(tokens));
	return tokens;
}

//returns a stringArray of format { 'transcript', 'fileid' }
StringArray Dictionary::getFileNameAndTranscript(
		string transcriptFileLine) {
	StringArray transcriptFileLineBreakUp;

	string delimiter1 = "(";
	string delimiter2 = ")";
	size_t pos = 0;
	string token;
	pos = transcriptFileLine.find(delimiter1);
	transcriptFileLineBreakUp.push_back(transcriptFileLine.substr(0, pos));
	transcriptFileLine.erase(0, pos + delimiter1.length());

	pos = transcriptFileLine.find(delimiter2);
	transcriptFileLineBreakUp.push_back(transcriptFileLine.substr(0, pos));

	return transcriptFileLineBreakUp;
}

//returns a stringArray of format { 'transcript', 'fileid' }
StringArray Dictionary::getFolderFileName(
		string transcriptFileLine) {
	StringArray transcriptFileLineBreakUp;

	string delimiter1 = "/";

	size_t pos = 0;
	string token;
	pos = transcriptFileLine.find(delimiter1);
	if (pos != std::string::npos) {
		transcriptFileLineBreakUp.push_back(transcriptFileLine.substr(0, pos));
		transcriptFileLine.erase(0, pos + delimiter1.length());
		transcriptFileLineBreakUp.push_back(transcriptFileLine);
	} else {
		transcriptFileLineBreakUp.push_back(transcriptFileLine);
	}

	return transcriptFileLineBreakUp;
}

string Dictionary::reverseString(
		string inString) {
	string outString = "";
	StringArray transList = string2words(inString);
	int transListLen = transList.size() - 1;
	for (int i = transListLen; i != -1; i--) {
		if (i == transListLen)
			outString += transList[i];
		else
			outString += (" " + transList[i]);
	}
	return outString;
}
StringArray Dictionary::splitStringUsingDelimiter(
		string s,
		string delimiter) {
	StringArray outVec;
	size_t pos = 0;
	std::string token;
	while ((pos = s.find(delimiter)) != std::string::npos) {
		outVec.push_back(s.substr(0, pos));
		s.erase(0, pos + delimiter.length());
	}
	outVec.push_back(s);
	return outVec;
}

} /* namespace std */
