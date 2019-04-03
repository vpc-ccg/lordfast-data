#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <set>
#include <algorithm>
#include <zlib.h>
#include <dirent.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

using namespace std;

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

class info_t
{
public:
	string qName;
	int64_t qLen;
	string tName;
	int64_t  tStart;
	int64_t  tEnd;
	int64_t  tLen;
	char   is_rev;
	float score;
	int nMatch;
	int nMismatch;
	int nIns;
	int nDel;
	string transQ;
	string transA;
	string transR;
	string read_seq;
};

void get_fields(string &line, vector<string> &v)
{
	v.clear();
	string s;
	istringstream sin(line);
	while(sin >> s)
		v.push_back(s);
}

template <typename T>
T str2type(string str)
{
	T n;
	istringstream sin(str);
	sin >> n;
	return n;
}

template <typename T>
string type2str(T v)
{
	ostringstream sout;
	sout << v;
	return sout.str();
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

int countN(string seq)
{
	int cnt = 0;
	for(int i=0; i<seq.size(); i++)
		cnt += (seq[i]=='N' || seq[i]=='n');
	return cnt;
}

void calcInfo(info_t &m)
{
	m.nMatch = 0;
	m.nMismatch = 0;
	m.nIns = 0;
	m.nDel = 0;
	m.transA = "";
	for(int i = 0; i < m.transQ.size(); i++)
	{
		if(m.transQ[i] == m.transR[i])
		{
			m.transA += '|';
			m.nMatch++;
		}
		else
		{
			m.transA += '*';
			if(m.transR[i] == '-')
			{
				m.nIns++;
			}
			else if(m.transQ[i] == '-')
			{
				m.nDel++;
			}
			else // mis-match
			{
				m.nMismatch++;
			}
		}
	}
	//
	if(m.is_rev == '-')
	{
		m.transQ = revComplement(m.transQ);
		m.transA = revComplement(m.transA);
		m.transR = revComplement(m.transR);
	}
	m.score = (float)m.nMatch / m.transA.size() * 100;
}

int main(int argc, char* argv[])
{
	if(argc < 5)
	{
		cerr<< "USAGE: ./extract pbsim_dir numOfReads minLen maxLen 1> reads.fasta 2> mapping.txt" << endl;
		return EXIT_FAILURE;
	}

	int64_t maxNumOfReads = str2type<int64_t>(argv[2]);
	int minLen = str2type<int>(argv[3]);
	int maxLen = str2type<int>(argv[4]);
	
	DIR *dir = opendir(argv[1]);
	if(dir == NULL)
	{
		cerr<< "ERROR: could not open directory: " << argv[1] << endl;
		return EXIT_FAILURE;
	}
	vector<string> lst_prefix;
	string prefix;
	dirent *ent;
	while ((ent = readdir (dir)) != NULL)
	{
		prefix = argv[1];
		prefix += "/";
		prefix += ent->d_name;
		if(prefix.find(".maf") != string::npos)
		{
			lst_prefix.push_back( prefix.substr(0, prefix.find(".maf")) );
		}
	}
	closedir(dir);
	// to make the experiment reproducible
	sort(lst_prefix.begin(), lst_prefix.end());


	vector<string> all_ids;
	vector<int> all_lens;
	set<string> selected_ids;
	int64_t numOfReads;
	ostringstream sout, serr;

	for(int i=0; i<lst_prefix.size(); i++)
	{
		// go through each file and store the read ids
		{
			gzFile file_read = gzopen((lst_prefix[i]+".fastq").c_str(), "r");
			kseq_t *seq = kseq_init(file_read);
			string read_id;
			while (kseq_read(seq) >= 0)
			{
				if(countN(seq->seq.s)==0 && seq->seq.l >= minLen && seq->seq.l <= maxLen)
				{
					all_ids.push_back(seq->name.s);
					all_lens.push_back(seq->seq.l);
				}
			}
			kseq_destroy(seq);
			gzclose(file_read);
		}
	}

	numOfReads = all_ids.size();
	if(maxNumOfReads > numOfReads)
	{
		maxNumOfReads = numOfReads;
	}

	int64_t totalBases = 0;
	// srand (time(NULL));
	srand (0); // to make the experiment reproducible
	while(selected_ids.size() < maxNumOfReads)
	{
		int r = rand() % numOfReads;
		if(selected_ids.count(all_ids[r]) == 0)
		{
			selected_ids.insert(all_ids[r]);
			totalBases += all_lens[r];
		}
	}
	// cerr<< totalBases << endl;
	// return 0;

	int64_t zmwId = 1;
	for(int i=0; i<lst_prefix.size(); i++)
	{
		string chr_name;
		int64_t chr_len;
		// read reference chromosome name
		{
			gzFile file_read = gzopen( (lst_prefix[i]+".ref").c_str(), "r");
			kseq_t *seq = kseq_init(file_read);
			if(kseq_read(seq) >= 0)
			{
				chr_name = seq->name.s;
				chr_len = seq->seq.l;
			}
			else
			{
				cerr<< "ERROR: no sequence in the reference file: " << (lst_prefix[i]+".ref") << endl;
				return EXIT_FAILURE;
			}
			kseq_destroy(seq);
			gzclose(file_read);
		}
		// read all the read information
		{
			gzFile file_read = gzopen((lst_prefix[i]+".fastq").c_str(), "r");
			ifstream file_maf((lst_prefix[i]+".maf").c_str());
			kseq_t *seq = kseq_init(file_read);
			string read_id;
			vector<string> fields;
			info_t tmp_info;
			string line;
			while (kseq_read(seq) >= 0)
			{
				read_id = seq->name.s;
				if(selected_ids.count(read_id) > 0) // The read id is selected for output. Print this read
				{
					tmp_info.read_seq = seq->seq.s;
					tmp_info.qName = seq->name.s;
					tmp_info.qLen = seq->seq.l;
					tmp_info.qName += "/" + type2str<int64_t>(zmwId) + "/0_" + type2str<int64_t>(tmp_info.qLen); zmwId++; // pacbio fasta header format
					tmp_info.tName = chr_name;
					tmp_info.tLen = chr_len;
					getline(file_maf, line); // a
					getline(file_maf, line); // ref info
					get_fields(line, fields);
					tmp_info.tStart = str2type<int64_t>(fields[2]);
					tmp_info.tEnd = tmp_info.tStart + str2type<int64_t>(fields[3]);
					tmp_info.transR = str2Upper(fields[6]);
					getline(file_maf, line); // read info
					get_fields(line, fields);
					tmp_info.is_rev = fields[4][0];
					tmp_info.transQ = str2Upper(fields[6]);
					getline(file_maf, line); // empty line
					calcInfo(tmp_info);

					// print the sequence in the standard output
					sout<< ">" << tmp_info.qName << "\n"
						<< tmp_info.read_seq << "\n";
					// print the alignment in the standard error
					serr<< std::setprecision(2);
					serr<< std::fixed;
					serr<< tmp_info.qName << " " << tmp_info.qLen << " " 
						<< "0 " << tmp_info.qLen << " + " 
						<< tmp_info.tName << " " << tmp_info.tLen << " "
						<< tmp_info.tStart << " " << tmp_info.tEnd << " "
						<< tmp_info.is_rev << " " << tmp_info.score << " "
						<< tmp_info.nMatch << " " << tmp_info.nMismatch << " "
						<< tmp_info.nIns << " " << tmp_info.nDel << " 255 "
						<< tmp_info.transQ << " " << tmp_info.transA << " "
						<< tmp_info.transR << "\n";

					if(sout.str().size() > 1000000)
					{
						cout<< sout.str();
						cout.flush();
						sout.str("");
						sout.clear();

						cerr<< serr.str();
						cerr.flush();
						serr.str("");
						serr.clear();
					}
				}
				else // Ignore this read
				{
					getline(file_maf, line); // a
					getline(file_maf, line); // ref info
					getline(file_maf, line); // read info
					getline(file_maf, line); // empty line
				}
			}
			kseq_destroy(seq);
			gzclose(file_read);
		}
	}

	if(sout.str().size() > 0)
	{
		cout<< sout.str();
		cout.flush();
		sout.str("");
		sout.clear();

		cerr<< serr.str();
		cerr.flush();
		serr.str("");
		serr.clear();
	}

	return 0;
}
