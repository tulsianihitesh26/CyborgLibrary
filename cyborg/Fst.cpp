/*
 * Fst.cpp
 *
 *  Created on: Jul 24, 2014
 *      Author: sharath
 */

#include "Fst.h"

namespace std {

Fst::Fst(
		string fstFolderPath,
		string fileName,
		string transcription,
		Dictionary dict,
		bool backwardFA,
		bool decoderMode) :
				m_fstFolderPath(fstFolderPath) {

	char fstCompileCommand[2048];
	int ret;

	// CREATE FOLDER IF NECESSARY
	StringArray folderFileName = Dictionary::getFolderFileName(fileName);
//	cout << folderFileName.size() <<endl;
	if (folderFileName.size() > 2) {
		cout << " ERROR: Wrong fileName format " << fileName << endl;
		exit(0);
	} else if (folderFileName.size() == 2) {
		sprintf(fstCompileCommand, "mkdir \"%s/%s\"", m_fstFolderPath.c_str(), folderFileName[0].c_str());
		ret = system(fstCompileCommand);
		if (ret != 0) {
			cout << " ERROR: Folder creation failed " << fstCompileCommand << endl;
//			exit(0);
		}
	}

	if (backwardFA) {
		transcription = Dictionary::reverseString(transcription);
	}

	if (decoderMode) {
		// FIND UNIQUE WORDS TO BE DECODED
		StringArray wordsInTranscription = Dictionary::string2words(transcription);
		StringArray uniqueWords;
		for (int i = 0; i < (int) wordsInTranscription.size(); i++) {
			if (!(std::find(uniqueWords.begin(), uniqueWords.end(), wordsInTranscription[i]) != uniqueWords.end())) {/* uniqueWords does not contain wordsInTranscription[i] */
				uniqueWords.push_back(wordsInTranscription[i]);
			}
		}

		// WRITE JSGF FILE
		ofstream jsgfFile;
		jsgfFile.open((m_fstFolderPath + fileName + ".jsgf").c_str());
		jsgfFile << "#JSGF V1.0;" << endl;
		jsgfFile << "grammar commands;" << endl;
		jsgfFile << "public <decodewords> = (SIL)( ";
		for (int i = 0; i < (int) uniqueWords.size(); i++) {
			jsgfFile << uniqueWords[i];
			if (i != ((int) uniqueWords.size() - 1)) {
				jsgfFile << " | ";
			}
		}
		jsgfFile << ")*[SIL];" << endl;

		//JSGF TO FSG
#if defined(__CYGWIN__) || defined(_WIN32) || defined(_WIN64)
		string sphinx_jsgf2fsg =
				"D:/Hitesh/CYBORG/sphinx_jsgf2fsg_exe/sphinx_jsgf2fsg.exe";

		sprintf(fstCompileCommand, "cmd.exe /c %s -jsgf %s/%s.jsgf -fsm %s/%s.fsm -symtab %s/%s.sym",
				sphinx_jsgf2fsg.c_str(), m_fstFolderPath.c_str(), fileName.c_str(), m_fstFolderPath.c_str(),
				fileName.c_str(), m_fstFolderPath.c_str(), fileName.c_str());
		ret = system(fstCompileCommand);
		if (ret != 0) {
			cout << " JSGF to FSG failed : " << fstCompileCommand << endl;
			exit(0);
		}

		//FSG TO WFST
		sprintf(fstCompileCommand,
				"cmd.exe /c fstcompile --acceptor --isymbols=%s/%s.sym --keep_isymbols %s/%s.fsm | fstdeterminize | fstminimize | fstrmepsilon | fstprint > %s/%s.wfst ",
				m_fstFolderPath.c_str(), fileName.c_str(), m_fstFolderPath.c_str(), fileName.c_str(),
				m_fstFolderPath.c_str(), fileName.c_str());

		ret = system(fstCompileCommand);
		if (ret != 0) {
			cout << " FSG to WFST failed : " << fstCompileCommand << endl;
			exit(0);
		}

		getFstTxtFilesMultipleTrans(m_fstFolderPath + fileName + ".wfst", transcription, dict);

#else
		cout << "ERROR: Decoder mode not supported for LINUX" <<endl;
		exit(0);
#endif
	} else {
		getFstTxtFiles(transcription, dict);
	}

//creating the fst file for wordGraph..
#if defined(__CYGWIN__) || defined(_WIN32) || defined(_WIN64)
	sprintf(fstCompileCommand,
			"cmd.exe /c fstcompile --isymbols=%s/lexicon.syms --osymbols=%s/lexicon.syms %s/wordGraph.stxt | fstrmepsilon | fstarcsort > %s/wordGraph_dir.fst",
			m_fstFolderPath.c_str(), m_fstFolderPath.c_str(), m_fstFolderPath.c_str(), m_fstFolderPath.c_str());
	ret = system(fstCompileCommand);
	if (ret != 0) {
		cout << " FST failed : " << fstCompileCommand << endl;
		exit(0);
	}

#elif defined(__linux) ||defined(__unix)
	sprintf(fstCompileCommand,"fstcompile --isymbols=%s/lexicon.syms --osymbols=%s/lexicon.syms %s/wordGraph.stxt | fstrmepsilon | fstarcsort > %s/wordGraph_dir.fst",
			m_fstFolderPath.c_str() , m_fstFolderPath.c_str() , m_fstFolderPath.c_str() , m_fstFolderPath.c_str() );
	ret=system(fstCompileCommand);
	if (ret!=0)
	{
		cout << " FST failed : " << fstCompileCommand <<endl;
		exit(0);
	}

#else
	cout << " FST not supported for this OS" <<endl;
	exit(0);
#endif

//creating the fst file for lexicon..
#if defined(__CYGWIN__) || defined(_WIN32) || defined(_WIN64)
	sprintf(fstCompileCommand,
			"cmd.exe /c fstcompile --isymbols=%s/lexicon.syms --osymbols=%s/lexicon.syms %s/lexicon.stxt | fstclosure | fstarcsort > %s/lexicon_dir.fst",
			m_fstFolderPath.c_str(), m_fstFolderPath.c_str(), m_fstFolderPath.c_str(), m_fstFolderPath.c_str());
	ret = system(fstCompileCommand);
	if (ret != 0) {
		cout << " FST failed : " << fstCompileCommand << endl;
		exit(0);
	}
#elif defined(__linux) ||defined(__unix)
	sprintf(fstCompileCommand,"fstcompile --isymbols=%s/lexicon.syms --osymbols=%s/lexicon.syms %s/lexicon.stxt | fstclosure | fstarcsort > %s/lexicon_dir.fst",
			m_fstFolderPath.c_str() , m_fstFolderPath.c_str() , m_fstFolderPath.c_str() , m_fstFolderPath.c_str() );
	ret=system(fstCompileCommand);
	if (ret!=0)
	{
		cout << " FST failed : " << fstCompileCommand <<endl;
		exit(0);
	}
#else
	cout << " FST not supported for this OS" <<endl;
	exit(0);
#endif

//print the final  fst file
#if defined(__CYGWIN__) || defined(_WIN32) || defined(_WIN64)
	sprintf(fstCompileCommand,
			"cmd.exe /c fstcompose %s/lexicon_dir.fst %s/wordGraph_dir.fst | fstrmepsilon | fstdeterminize | fstminimize > %s/final_dir.fst",
			m_fstFolderPath.c_str(), m_fstFolderPath.c_str(), m_fstFolderPath.c_str());
	ret = system(fstCompileCommand);
	if (ret != 0) {
		cout << " FST failed : " << fstCompileCommand << endl;
		exit(0);
	}
#elif defined(__linux) ||defined(__unix)
	sprintf(fstCompileCommand,"fstcompose %s/lexicon_dir.fst %s/wordGraph_dir.fst | fstrmepsilon | fstdeterminize | fstminimize > %s/final_dir.fst",
			m_fstFolderPath.c_str() , m_fstFolderPath.c_str() , m_fstFolderPath.c_str());
	ret=system(fstCompileCommand);
	if (ret!=0)
	{
		cout << " FST failed : " << fstCompileCommand <<endl;
		exit(0);
	}
#else
	cout << " FST not supported for this OS" <<endl;
	exit(0);
#endif

//print the final fst txt file
#if defined(__CYGWIN__) || defined(_WIN32) || defined(_WIN64)
	sprintf(fstCompileCommand,
			"cmd.exe /c fstprint --isymbols=%s/lexicon.syms --osymbols=%s/lexicon.syms %s/final_dir.fst %s/%s.fst",
			m_fstFolderPath.c_str(), m_fstFolderPath.c_str(), m_fstFolderPath.c_str(), m_fstFolderPath.c_str(),
			fileName.c_str());
	ret = system(fstCompileCommand);
	if (ret != 0) {
		cout << " FST failed : " << fstCompileCommand << endl;
		exit(0);
	}
#elif defined(__linux) ||defined(__unix)
	sprintf(fstCompileCommand,"fstprint --isymbols=%s/lexicon.syms --osymbols=%s/lexicon.syms %s/final_dir.fst %s/%s.fst",
			m_fstFolderPath.c_str() , m_fstFolderPath.c_str() , m_fstFolderPath.c_str(), m_fstFolderPath.c_str(),fileName.c_str());
	ret=system(fstCompileCommand);
	if (ret!=0)
	{
		cout << " FST failed : " << fstCompileCommand <<endl;
		exit(0);
	}
#else
	cout << " FST not supported for this OS" <<endl;
	exit(0);
#endif

}

/**
 * This function creates three text files in the current directory
 * <b>wordGraph.stxt</b>, <b>lexicon.stxt</b> and <b>lexicon.syms</b> from the transcription and the
 * dictionary.<br>
 * The function will run provided the dictionary file has been loaded in the hashmap.<br>
 * The format of the file is in accordance with the openFST open source project.</br>
 * For more details of the openFST library, visit the site
 * <a href = "http://www.openfst.org/twiki/bin/view/FST/WebHome"> OpenFST </a> <br> <br>
 * The mapping for &lts>, &lt/s> and &ltsil> has to be 'SIL'. This is hardcoded in the function.<br>
 * Start silence has to be &lts>, end silence &lt/s>. <br>
 * Start markers (#1) and end word markers (#0) have been added to get the start and end triphone.
 *
 * The mapping from filler to filler has to be same.
 *  e.g. +bn+	+bn+
 *
 * @param transcription
 *
 */

int Fst::getFstTxtFiles(
		string transcription,
		Dictionary dict) {

	lexiconPos = 2;
	symsPos = 0;
//	FstCreateFiles fg = new FstCreateFiles();
	hashMapSyms = map<string, string>();
	wordGraph.open((m_fstFolderPath + "//wordGraph.stxt").c_str());
	lexicon.open((m_fstFolderPath + "//lexicon.stxt").c_str());
	syms.open((m_fstFolderPath + "//lexicon.syms").c_str());

	//write the silences, epsilon and fillers in the symbol file
	StringArray fillerWords = dict.getFillerWords();
	printSilencesAndFillers(fillerWords);

	//write start and end word markers in the symbol file
	printSyms("#1");
	printSyms("#0");

	string phoneme;
	StringArray words = removeAllFillers(Dictionary::string2words(transcription), fillerWords);
//	StringArray words=Dictionary::string2words(transcription);

	int totalWords = words.size();

	int pos = 0;
//	print start silence in wordGraph.
	printWordGraph(pos, "<s>");
	pos++;
	printFillerInWordGraph(pos, fillerWords);
	for (int i = 0; i < totalWords; i++) {
		printWordGraph(pos, words[i]);
		phoneme = dict.getPhonemes(words[i]);
		if (phoneme.compare("NULL") == 0) {
			cout << words[i] << " not present is the dictionary" << endl;
			exit(0);
		}
		createHash(words[i], phoneme);
		writeLexicon(phoneme, words[i]);

		//alternate pronunciations are marked as (2), (3), ...
		//e.g. hello(2), hello(3),..
		int no = 1;
		std::ostringstream tempStr;
		tempStr << words[i] << "(" << no << ")";
		string altPronunciation = tempStr.str();

		while ((phoneme = dict.getPhonemes(altPronunciation)).compare("NULL") != 0) {
			//printing the alternate pronunciation
			printWordGraph(pos, altPronunciation);
			no++;
			createHash(altPronunciation, phoneme);
			writeLexicon(phoneme, altPronunciation);

			tempStr.str("");
			tempStr.clear(); // Clear state flags.
			tempStr << words[i] << "(" << no << ")";
			altPronunciation = tempStr.str();
		}
		pos++;
		if (i < (totalWords - 1)) {
			//printing intermediate silence/ epsilon
			printWordGraph(pos, "<eps>");
			printWordGraph(pos, "<sil>");
			pos++;
			printFillerInWordGraph(pos, fillerWords);

		}

	}
//	print end fillers in wordGraph
	printFillerInWordGraph(pos, fillerWords);

//	print end silence in wordGraph.
	printWordGraph(pos, "<eps>");
	printWordGraph(pos, "</s>");
	pos++;

//	print silence and fillers in the lexicon file.
	printSilAndFillersInLexicon(fillerWords);

//	printing the end position in lexicon and wordGraph file in accordance with the OpenFst notations.
	lexicon << 1;
	wordGraph << pos;

	if (wordGraph != NULL)
		wordGraph.close();
	if (lexicon != NULL)
		lexicon.close();
	if (syms != NULL)
		syms.close();

	return 0;
}

int Fst::getFstTxtFilesMultipleTrans(
		string wfstFile,
		string transcription,
		Dictionary dict) {

	ifstream wfstFileID(wfstFile.c_str());
	if (!wfstFileID.is_open()) {
		cout << "ERROR: Unable to open : " << wfstFile << endl;
		exit(0);
	}

	StringArray words;
	IntArray sInd;
	IntArray eInd;
	int maxInd = -1;
	string tempLine;
	while (getline(wfstFileID, tempLine)) {
		StringArray wordList = Dictionary::string2words(tempLine);
		if (wordList.size() > 1) {
			sInd.push_back(atoi(wordList[0].c_str()));
			eInd.push_back(atoi(wordList[1].c_str()));
			if (wordList[2].compare("SIL") == 0) {
				words.push_back("<sil>");
			} else {
				words.push_back(wordList[2]);
			}
		} else {
			maxInd = atoi(wordList[0].c_str());
		}
	}

	lexiconPos = 2;
	symsPos = 0;
//	FstCreateFiles fg = new FstCreateFiles();
	hashMapSyms = map<string, string>();
	wordGraph.open((m_fstFolderPath + "//wordGraph.stxt").c_str());
	lexicon.open((m_fstFolderPath + "//lexicon.stxt").c_str());
	syms.open((m_fstFolderPath + "//lexicon.syms").c_str());

	//write the silences, epsilon and fillers in the symbol file
	StringArray fillerWords = dict.getFillerWords();
	printSilencesAndFillers(fillerWords);

	//write start and end word markers in the symbol file
	printSyms("#1");
	printSyms("#0");

	string phoneme;

	int totalWords = words.size();
	for (int i = 0; i < totalWords; i++) {
		if ((words[i].compare("<sil>") == 0) && sInd[i] != 0) // && eInd[i]!=maxInd - add this to force SIL at the end
				{
			printWordGraphFromTo(sInd[i], eInd[i], "<eps>");
		}
		printWordGraphFromTo(sInd[i], eInd[i], words[i]);

		if ((words[i].compare("<sil>") != 0) && (words[i].compare("<eps>") != 0)) {

			phoneme = dict.getPhonemes(words[i]);
			if (phoneme.compare("NULL") == 0) {
				cout << words[i] << " not present is the dictionary" << endl;
				exit(0);
			}
			createHash(words[i], phoneme);
			writeLexicon(phoneme, words[i]);

			//alternate pronunciations are marked as (1), (2), ...
			//e.g. hello(1), hello(2),..
			int no = 1;
			std::ostringstream tempStr;
			tempStr << words[i] << "(" << no << ")";
			string altPronunciation = tempStr.str();

			while ((phoneme = dict.getPhonemes(altPronunciation)).compare("NULL") != 0) {
				//printing the alternate pronunciation
				printWordGraphFromTo(sInd[i], eInd[i], altPronunciation);
				no++;
				createHash(altPronunciation, phoneme);
				writeLexicon(phoneme, altPronunciation);

				tempStr.str("");
				tempStr.clear(); // Clear state flags.
				tempStr << words[i] << "(" << no << ")";
				altPronunciation = tempStr.str();
			}
		}
	}
//	print silence and fillers in the lexicon file.
	printSilAndFillersInLexicon(fillerWords);

	lexicon << 1;
	if (maxInd == -1) {
		cout << "ERROR: Parsing JSGF failed " << endl;
		exit(0);
	}
	wordGraph << maxInd;

	if (wordGraph != NULL)
		wordGraph.close();
	if (lexicon != NULL)
		lexicon.close();
	if (syms != NULL)
		syms.close();

	return 0;

}

/**
 * This function will write the silences, epsilon and fillers in the symbol file
 *
 */
void Fst::printSilencesAndFillers(
		StringArray fillerWords) {
	printSyms("<eps>");
	printSyms("SIL");

	for (int i = 0; i < (int) fillerWords.size(); i++) {
		printSyms(fillerWords[i]);
	}
}

void Fst::printSyms(
		string c) {
	syms << c << " " << symsPos << endl;
//	cout<<c <<" " <<symsPos <<endl;
	symsPos++;
}

void Fst::printWordGraph(
		int pos,
		string c) {
	wordGraph << pos << " " << (pos + 1) << " " << c << " " << c << endl;
//	cout<<pos << " " << (pos+1) << " " << c << " " << c <<endl;
}

void Fst::printWordGraphFromTo(
		int pos,
		int nextPos,
		string c) {
//					wordGraph.write(pos + " " + (pos+1) + " " + c + " " + c + "\n");
	wordGraph << pos << " " << nextPos << " " << c << " " << c << endl;
}

/**
 * This function will print the fillers in the wordGraph file
 * in addition to the <eps> as fillers are not compulsory.
 */
void Fst::printFillerInWordGraph(
		int & pos,
		StringArray fillerWords) {

	//0,1,2 reserved for silences
	if (fillerWords.size() > 3) {
		for (int i = 3; i < (int) fillerWords.size(); i++)
			printWordGraph(pos, fillerWords[i]);

		printWordGraph(pos, "<eps>");
		pos++;
	}
}

void Fst::createHash(
		string word,
		string phoneme) {

	if (hashMapSyms.find(word) == hashMapSyms.end()) {
		hashMapSyms[word] = "1";
		printSyms(word);
	}
	StringArray monophone = Dictionary::string2words(phoneme);
	for (StringArray::iterator it = monophone.begin(); it != monophone.end(); ++it) {
		if (hashMapSyms.find(*it) == hashMapSyms.end()) {
			hashMapSyms[*it] = "1";
			printSyms(*it);
		}
	}
}

void Fst::writeLexicon(
		string phoneme,
		string word) {
	StringArray monophone = Dictionary::string2words(phoneme);
	int length = monophone.size();
//			#1 is the start of word marker
	printLexicon(0, lexiconPos, "#1", word);

	for (int i = 0; i < length; i++) {
		printLexicon(lexiconPos, lexiconPos + 1, monophone[i], "<eps>");
		lexiconPos++;
	}
//			#0 is end of word marker
	printLexicon(lexiconPos, 1, "#0", "<eps>");
	lexiconPos++;
}

void Fst::printLexicon(
		int pos,
		int nextPos,
		string c,
		string d) {
	lexicon << pos << " " << nextPos << " " << c << " " << d << endl;
//	cout<<pos << " " << nextPos<< " " << c << " " << d << endl;
}

void Fst::printSilAndFillersInLexicon(
		StringArray fillerWords) {
	//printing the silences in the lexicon
	printLexicon(0, 1, "SIL", "<s>");
	printLexicon(0, 1, "SIL", "</s>");
	printLexicon(0, 1, "SIL", "<sil>");
	printLexicon(0, 1, "#0", "#0");
	//0,1,2 reserved for silences
	for (int i = 3; i < (int) fillerWords.size(); i++)
		printLexicon(0, 1, fillerWords[i], fillerWords[i]);
}

bool Fst::isFiller(
		string word,
		StringArray fillerWords) {
	if (word.compare("SIL") == 0) {
		return true;
	}
	for (int i = 0; i < (int) fillerWords.size(); i++)
		if (word.compare(fillerWords[i]) == 0) {
			return true;
			break;
		}
	return false;
}

StringArray Fst::removeAllFillers(
		StringArray words,
		StringArray fillerWords) {
	StringArray noFillerWordList;
	for (unsigned int i = 0; i < words.size(); i++) {
		if (!isFiller(words[i], fillerWords)) {
			noFillerWordList.push_back(words[i]);
		}
	}
	return noFillerWordList;
}

} /* namespace std */
