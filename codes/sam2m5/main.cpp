
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;

typedef struct
{
	string qName;
	int qLen;
	int qStart;
	int qEnd;
	char qStrand;
	string rName;
	int64_t rLen;
	int64_t rStart;
	int64_t rEnd;
	char rStrand;
	float score;
	int nMatch;
	int nMisMatch;
	int nIns;
	int nDel;
	int mapQ;
	string qTrans;
	string mTrans;
	string rTrans;
} m5_t;

typedef struct
{
	string qName;
	int flag;
	string rName;
	int64_t rStart;
	int mapQ;
	string cigar;
	int tLen;
	string seq;
} sam_t;

char lowVal[128] = {
   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
  64, 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o',
 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',  91,  92,  93,  94,  95,
  96, 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o',
 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 123, 124, 125, 126, 127
};

char upVal[128] = {
   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
  64, 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',  91,  92,  93,  94,  95,
  96, 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 123, 124, 125, 126, 127
};

char revVal[128] = {
   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
  64, 'T', 'B', 'G', 'D', 'E', 'F', 'C', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
 'P', 'Q', 'R', 'S', 'A', 'U', 'V', 'W', 'X', 'Y', 'Z',  91,  92,  93,  94,  95,
  96, 't', 'b', 'g', 'd', 'e', 'f', 'c', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o',
 'p', 'q', 'r', 's', 'a', 'u', 'v', 'w', 'x', 'y', 'z', 123, 124, 125, 126, 127
};

template <typename T> T str2type(string str);
string str2Lower(string str);
string str2Upper(string str);
bool   getNextSamEntry(ifstream &fin, sam_t &m);
void getAlignmentStat(string cigar, string &qSeq, string &rSeq, int64_t rStartIn, int &qStart, int &qEnd, int64_t &rStart, int64_t &rEnd, 
int &nMatch, int &nMisMatch, int &nIns, int &nDel, string &transQ, string &transA, string &transR);
string revComplement(string &str);
void   getAllReads(string path, map<string, string> &readList);
bool convertSam2M5(sam_t &mapS, m5_t &mapM, map<string, string> &chrList);
void printM5(ostream &out, m5_t &m);

int main(int argc, char *argv[])
{
	if(argc != 4)
	{
		cerr << "USAGE: ./sam2mx <ref.fa> <input.sam> <mode>" << endl;
		cerr << endl;
		cerr << "mode   0: print primary" << endl;
		cerr << "       1: print primary and supplementary" << endl;
		cerr << "       2: print all" << endl;
		return EXIT_FAILURE;
	}
	map<string, string> refList;
	getAllReads(argv[1], refList);

	//
	ifstream finSam(argv[2]);
	if(finSam.is_open() == false)
	{
		cerr<< "Could not open file: " << argv[2] << endl;
		exit(EXIT_FAILURE);
	}

	int mode = str2type<int>(argv[3]);

	vector<string> field;
	sam_t mapSam;
	m5_t mapM5;

	string currentLrId = "";
	string currentLrSeq = "";
	string currentLrSeq_rev = "";

	while(getNextSamEntry(finSam, mapSam))
	{
		if(mapSam.qName != currentLrId)
		{
			currentLrId = mapSam.qName;
			if(mapSam.flag & 16)
			{
				currentLrSeq_rev = mapSam.seq;
				currentLrSeq = revComplement(mapSam.seq);
			}
			else
			{
				currentLrSeq = mapSam.seq;
				currentLrSeq_rev = revComplement(mapSam.seq);
			}
			//
			if(convertSam2M5(mapSam, mapM5, refList))
				printM5(cout, mapM5);
		}
		else // another mapping of the same query
		{
			if( (mode==1 && (mapSam.flag & 2048)) || mode==2 )
			{
				// fix SEQ if not stored
				// if(mapSam.seq == "*")
				// {
					if(mapSam.flag & 16)
					{
						mapSam.seq = currentLrSeq_rev;
					}
					else
					{
						mapSam.seq = currentLrSeq;
					}
				// }
				//
				if(convertSam2M5(mapSam, mapM5, refList))
					printM5(cout, mapM5);
			}
		}
		// break;
	}

	return EXIT_SUCCESS;
}

template <typename T>
T str2type(string str)
{
	T n;
	istringstream sin(str);
	sin >> n;
	return n;
}

// template <typename T>
// string type2str(T v)
// {
// 	ostringstream sout;
// 	sout << v;
// 	return sout.str();
// }

string str2Lower(string str)
{
	string lowerStr(str.size(), 'n');
	for(int i=0; i<str.size(); i++)
	{
		lowerStr[i] = lowVal[str[i]];
	}
	return lowerStr;
}

string str2Upper(string str)
{
	string upperStr(str.size(), 'N');
	for(int i=0; i<str.size(); i++)
	{
		upperStr[i] = upVal[str[i]];
	}
	return upperStr;
}

string revComplement(string &str)
{
	int n = str.size();
	string revStr(n, 'N');
	for(int i=0; i<n; i++)
		revStr[i] = revVal[str[n-i-1]];
	return revStr;
}

void getAllReads(string path, map<string, string> &readList)
{
	readList.clear();
	// cerr << "[getAllReads] NOTE: loading sequences from file " << path << endl;
	gzFile readFile;
	kseq_t *readObj;
	readFile = gzopen(path.c_str(), "r");
	if(readFile==NULL)
	{
		cerr<< "[getAllReads] ERROR: could not open file: " << path << endl;
		exit(0);
	}

	int64_t cnt = 0;
	readObj = kseq_init(readFile);
	while(kseq_read(readObj) >= 0)
	{
		cnt++;
		readList.insert(pair<string, string>(readObj->name.s, str2Upper(readObj->seq.s)));
	}
	// cerr << "[getAllReads] NOTE: " << cnt << " sequences are loaded" << endl;
	gzclose(readFile);
}

bool getNextSamEntry(ifstream &fin, sam_t &m)
{
	bool entryFound = false;
	string line;
	while(getline(fin, line))
	{
		if(line.empty() == false && line[0] != '@')
		{
			entryFound = true;
			break;
		}
	}
	if(entryFound==false)
		return false;
	//
	vector<string> field;
	istringstream sin(line);
	string s;
	while(sin >> s)
	{
		field.push_back(s);
	}
	//
	m.qName = field[0];
	m.flag = str2type<int>(field[1]);
	m.rName = field[2];
	m.rStart = str2type<int64_t>(field[3]);
	m.mapQ = str2type<int>(field[4]);
	m.cigar = field[5];
	m.tLen = str2type<int>(field[8]);
	m.seq = str2Upper(field[9]);
	//
	return true;
}

bool convertSam2M5(sam_t &mapS, m5_t &mapM, map<string, string> &chrList)
{
	if( (mapS.flag & 4) && mapS.rName == "*" )
		return false;

	mapM.qName = mapS.qName;
	mapM.qLen = mapS.seq.size();
	mapM.qStrand = '+';// ( mapS.flag & 16 ? '-' : '+' );
	mapM.rName = mapS.rName;
	mapM.rLen = chrList[mapS.rName].size();
	mapM.rStrand = ( mapS.flag & 16 ? '-' : '+' );
	mapM.mapQ = mapS.mapQ;
	getAlignmentStat(mapS.cigar, mapS.seq, chrList[mapS.rName], mapS.rStart - 1, mapM.qStart, mapM.qEnd, mapM.rStart, mapM.rEnd, 
		mapM.nMatch, mapM.nMisMatch, mapM.nIns, mapM.nDel, mapM.qTrans, mapM.mTrans, mapM.rTrans);
	if(mapM.mTrans.size() == 0)
		return false;
	mapM.score = (float)mapM.nMatch / mapM.mTrans.size() * 100;
	if(mapM.rStrand == '-')
	{
		int qStart_tmp = mapM.qStart;
		mapM.qStart = mapM.qLen - mapM.qEnd;
		mapM.qEnd = mapM.qLen - qStart_tmp;
		mapM.qTrans = revComplement(mapM.qTrans);
		mapM.mTrans = revComplement(mapM.mTrans);
		mapM.rTrans = revComplement(mapM.rTrans);
	}

	return true;
}

void printM5(ostream &out, m5_t &m)
{
	out << std::setprecision(2);
	out << std::fixed;
	out << m.qName << " " << m.qLen << " " << m.qStart << " " << m.qEnd << " " << m.qStrand << " "
	    << m.rName << " " << m.rLen << " " << m.rStart << " " << m.rEnd << " " << m.rStrand << " "
	    << m.score << " " << m.nMatch << " " << m.nMisMatch << " " << m.nIns << " " << m.nDel << " " << m.mapQ << " "
	    << m.qTrans << " " << m.mTrans << " " << m.rTrans << "\n";
	    // << "\n" << m.qTrans << "\n" << m.mTrans << "\n" << m.rTrans << "\n";
}

void getAlignmentStat(string cigar, string &qSeq, string &rSeq, int64_t rStartIn, int &qStart, int &qEnd, int64_t &rStart, int64_t &rEnd, 
	int &nMatch, int &nMisMatch, int &nIns, int &nDel, string &transQ, string &transA, string &transR)
{
	int cnt = 0;
	int qIndex = 0;
	int64_t rIndex = rStartIn;
	int i;
	istringstream sin(cigar);
	int n;
	char c;

	qStart = 0;
	qEnd = 0;
	rStart = rIndex;
	rEnd = 0;
	nMatch = 0;
	nMisMatch = 0;
	nIns = 0;
	nDel = 0;
	transQ = "";
	transA = "";
	transR = "";

	while(sin>> n >> c)
	{
		switch (c)
		{
			case 'H':
			//	break;
			case 'S':
				if(cnt==0)
				{
					qIndex += n;
					qStart = qIndex;
				}
				break;
			case 'M':
				for(i=0; i<n; i++)
				{
					transQ += qSeq[qIndex];
					transR += rSeq[rIndex];
					if(qSeq[qIndex]==rSeq[rIndex])
					{
						transA += "|";
						nMatch++;
					}
					else
					{
						transA += "*";
						nMisMatch++;
					}
					qIndex++;
					rIndex++;
				}
				break;
			case '=':
				for(i=0; i<n; i++)
				{
					transQ += qSeq[qIndex];
					transR += rSeq[rIndex];
					transA += "|";
					nMatch++;
					qIndex++;
					rIndex++;
				}
				break;
			case 'X':
				for(i=0; i<n; i++)
				{
					transQ += qSeq[qIndex];
					transR += rSeq[rIndex];
					transA += "*";
					nMisMatch++;
					qIndex++;
					rIndex++;
				}
				break;
			case 'I':
				for(i=0; i<n; i++)
				{
					transQ += qSeq[qIndex];
					transA += '*';
					transR += '-';

					nIns++;
					qIndex++;
				}
				break;
			case 'D':
				for(i=0; i<n; i++)
				{
					transQ += '-';
					transA += '*';
					transR += rSeq[rIndex];

					nDel++;
					rIndex++;
				}
				break;
			default:
				cerr<< "[getAlignmentStat] ERROR: cannot parse cigar element " << n << c << endl;
				exit(EXIT_FAILURE);
		}
		cnt++;
	}
	qEnd = qIndex;
	rEnd = rIndex;
}
