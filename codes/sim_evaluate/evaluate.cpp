#include "kseq.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <zlib.h>
#include <stdint.h>

using namespace std;

// #define DODEBUG

#ifdef DODEBUG 
#define DEBUG(cmd) cmd
#else 
#define DEBUG(cmd)
#endif

int mappedRankUB = 1;
int correctRankUB = 1;
float minOverlap = 0.95;
int corDistLimit = 100;

class m5_t
{
public:
    string qName;
    int32_t qLen;
    int32_t qStart;
    int32_t qEnd;
    string tName;
    int64_t tLen;
    int64_t tStart;
    int64_t tEnd;
    bool isRev;
    double alnScore;
    int32_t nMatch;
    int32_t nMismatch;
    int32_t nInsertion;
    int32_t nDeletion;
    string qTranscript;
    string aTranscript;
    string tTranscript;
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

bool getNextM5(ifstream &fin, m5_t &map)
{
    string line;
    if(getline(fin, line))
    {
        vector<string> fields;
        get_fields(line, fields);
        map.qName = fields[0];
        map.qLen = str2type<int32_t>(fields[1]);
        map.qStart = str2type<int32_t>(fields[2]);
        map.qEnd = str2type<int32_t>(fields[3]);
        map.tName = fields[5];
        map.tLen = str2type<int64_t>(fields[6]);
        map.tStart = str2type<int64_t>(fields[7]);
        map.tEnd = str2type<int64_t>(fields[8]);
        map.isRev = (fields[9][0] == '-' ? true : false);
        map.alnScore = str2type<double>(fields[10]);
        map.nMatch = str2type<int32_t>(fields[11]);
        map.nMismatch = str2type<int32_t>(fields[12]);
        map.nInsertion = str2type<int32_t>(fields[13]);
        map.nDeletion = str2type<int32_t>(fields[14]);
        map.qTranscript = fields[16];
        map.aTranscript = fields[17];
        map.tTranscript = fields[18];
        return true;
    }
    else
    {
        return false;
    }
}

int calcRank(m5_t &mTrue, vector<m5_t> &mList)
{
    int i;
    int r = -1;
    for(i = 0; i < mList.size(); i++)
    {
        if(mList[i].isRev == mTrue.isRev && mList[i].tName == mTrue.tName)
        {
            if((mList[i].tStart < mTrue.tStart && mList[i].tEnd >= mTrue.tStart && (float)(mList[i].tEnd - mTrue.tStart)/(mTrue.tEnd - mTrue.tStart) > minOverlap)
            || (mList[i].tStart >= mTrue.tStart && mList[i].tStart <= mTrue.tEnd) && (float)(min(mList[i].tEnd, mTrue.tEnd) - mTrue.tStart)/(mTrue.tEnd - mTrue.tStart) > minOverlap)
            {
                r = i+1;
                break;
            }   
        }
    }
    return r;
}

bool isCorrectlyMapped(m5_t &mTrue, vector<m5_t> &mList)
{
    int i;
    int ovlp = 0;
    for(i = 0; i < mList.size(); i++)
    {
        if(mList[i].isRev == mTrue.isRev && mList[i].tName == mTrue.tName)
        {
            if(mList[i].tStart < mTrue.tStart)
            {
                if(mList[i].tEnd >= mTrue.tStart)
                {
                    ovlp += min(mList[i].tEnd, mTrue.tEnd) - mTrue.tStart;
                }
            }
            else // mList[i].tStart >= mTrue.tStart
            {
                if(mList[i].tStart <= mTrue.tEnd)
                {
                    ovlp += min(mList[i].tEnd, mTrue.tEnd) - mList[i].tStart;
                }
            }
        }
    }

    return ovlp/(float)(mTrue.tEnd - mTrue.tStart) > minOverlap;
}

void calcIntervals(m5_t &mTrue, vector<pair<int64_t, int64_t> > &intvList)
{
    int32_t i;
    int32_t qIdx = 0;
    int32_t tIdx = 0;
    int32_t lenTrans = mTrue.qTranscript.size();

    for(i = 0; i < lenTrans; i++)
    {
        if(mTrue.aTranscript[i] == '|') // match
        {
            intvList.push_back(make_pair(mTrue.tStart + tIdx - corDistLimit, mTrue.tStart + tIdx + corDistLimit));
            qIdx++;
            tIdx++;
            if(qIdx != intvList.size())
            {
                cerr<< "ERROR: incorrect index for the list of intervals!" << endl;
                exit(EXIT_FAILURE);
            }
        }
        else // mis-match, insertion, or deletion
        {
            if(mTrue.qTranscript[i] != '-' && mTrue.tTranscript[i] != '-') // mis-match
            {
                intvList.push_back(make_pair(mTrue.tStart + tIdx - corDistLimit, mTrue.tStart + tIdx + corDistLimit));
                qIdx++;
                tIdx++;
                if(qIdx != intvList.size())
                {
                    cerr<< "ERROR: incorrect index for the list of intervals!" << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else if(mTrue.tTranscript[i] == '-') // insertion
            {
                intvList.push_back(make_pair(mTrue.tStart + tIdx - 1 - corDistLimit, mTrue.tStart + tIdx + corDistLimit));
                qIdx++;
                if(qIdx != intvList.size())
                {
                    cerr<< "ERROR: incorrect index for the list of intervals!" << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else // deletion
            {
                tIdx++;
            }
        }
    }
}

void calcBaseStat(m5_t &mRep, vector<pair<int64_t, int64_t> > &intvList, int64_t &baseCorrect, int64_t &baseIncorrect, int64_t &baseUnmapped)
{
    // cout<< mRep.qTranscript << endl
    //  << mRep.aTranscript << endl
    //  << mRep.tTranscript << endl;

    int32_t i;
    int32_t qIdx = mRep.qStart;
    int32_t tIdx = 0;
    int32_t lenTrans = mRep.qTranscript.size();

    baseCorrect = 0;
    baseIncorrect = 0;
    baseUnmapped = 0;

    for(i = 0; i < lenTrans; i++)
    {
        if(mRep.aTranscript[i] == '|') // match
        {
            // cout<< qIdx << " " << tIdx << "\t\t" << mRep.tStart + tIdx << "\t\t" << intvList[qIdx].first << ", " << intvList[qIdx].second << endl;
            if(mRep.tStart + tIdx >= intvList[qIdx].first && mRep.tStart + tIdx <= intvList[qIdx].second)
                baseCorrect++;
            else
                baseIncorrect++;
            qIdx++;
            tIdx++;
            // cout<< baseCorrect << " " << baseIncorrect << endl;
            // exit(0);
        }
        else // mis-match, insertion, or deletion
        {
            if(mRep.qTranscript[i] != '-' && mRep.tTranscript[i] != '-') // mis-match
            {
                if(mRep.tStart + tIdx >= intvList[qIdx].first && mRep.tStart + tIdx <= intvList[qIdx].second)
                    baseCorrect++;
                else
                    baseIncorrect++;
                qIdx++;
                tIdx++;
            }
            else if(mRep.tTranscript[i] == '-') // insertion
            {
                if(mRep.tStart + tIdx >= intvList[qIdx].first && mRep.tStart + tIdx <= intvList[qIdx].second)
                    baseCorrect++;
                else
                    baseIncorrect++;
                qIdx++;
            }
            else // deletion
            {
                tIdx++;
            }
        }
    }

    // calculate the unmapped bases
    baseUnmapped += mRep.qStart;
    baseUnmapped += mRep.qLen - mRep.qEnd;
}

KSEQ_INIT(gzFile, gzread)

int main(int argc, char* argv[])
{
    if(argc < 6)
    {
        cerr<< "USAGE: ./evaluate reads.fasta true.map reported.map corDistLimit minOverlap" << endl;
        return EXIT_FAILURE;
    }

    gzFile fileRead = gzopen(argv[1], "r");
    if(fileRead == NULL)
    {
        cerr<< "ERROR: could not open file: " << argv[1] << endl;
        return EXIT_FAILURE;
    }

    ifstream fileTrue(argv[2]);
    if(fileTrue.is_open() == false)
    {
        cerr<< "ERROR: could not open file: " << argv[2] << endl;
        return EXIT_FAILURE;
    }

    ifstream fileRep(argv[3]);
    if(fileRep.is_open() == false)
    {
        cerr<< "ERROR: could not open file: " << argv[3] << endl;
        return EXIT_FAILURE;
    }

    corDistLimit = str2type<int>(argv[4]);
    if(corDistLimit <= 0)
    {
        cerr<< "ERROR: corDistLimit has to be a positive integer." << endl;
        return EXIT_FAILURE;
    }

    minOverlap = str2type<float>(argv[5]);
    if(minOverlap < 0 || minOverlap > 1)
    {
        cerr<< "ERROR: minOverlap has to be a float number [0..1]." << endl;
        return EXIT_FAILURE;
    }

    m5_t m5True;
    m5_t m5Rep;
    string currentId;
    vector<string> fields;

    int64_t resReadTotal = 0;
    int64_t resReadMapped = 0;
    int64_t resReadCorrect = 0;
    int64_t resBaseTotal = 0;
    int64_t resBaseCorrect = 0;
    int64_t resBaseIncorrect = 0;
    int64_t resBaseUnmapped = 0;

    // get the first m5 entry from the reported mapping file
    getNextM5(fileRep, m5Rep);
    // for each long read in the fileRead
    kseq_t *seq = kseq_init(fileRead);
    while (kseq_read(seq) >= 0)
    {
        DEBUG({
            fprintf(stderr, ">%s length:%u\n", seq->name.s, seq->seq.l);
        });
        resReadTotal++;
        resBaseTotal += seq->seq.l;
        currentId = seq->name.s;
        getNextM5(fileTrue, m5True);
        if(m5True.qName != currentId)
        {
            cerr<< "ERROR: could not find the correct mapping for this read! Should not happen!" << endl;
            exit(EXIT_FAILURE);
        }
        vector<m5_t> mapList;
        while(m5Rep.qName == currentId)
        {
            mapList.push_back(m5Rep);
            if(getNextM5(fileRep, m5Rep) == false)
                break;
        }

        if(mapList.size() > 0)
        {
            resReadMapped++;
            if(isCorrectlyMapped(m5True, mapList))
            {
                // fprintf(stderr, "%s\n", seq->name.s);
                resReadCorrect++;
                //
                int i;
                int64_t baseCorrect, baseIncorrect, baseUnmapped;
                int64_t baseCorrect_sum = 0, baseIncorrect_sum = 0, baseUnmapped_sum = 0;
                vector<pair<int64_t, int64_t> > intervals;
                //
                calcIntervals(m5True, intervals);
                //
                for(i = 0; i < mapList.size(); i++)
                {
                    // FIXME: if supplementary alignments are overlapping, counts the overlapping bases twice
                    calcBaseStat(mapList[i], intervals, baseCorrect, baseIncorrect, baseUnmapped);
                    baseCorrect_sum += baseCorrect;
                    baseIncorrect_sum += baseIncorrect;
                }
                // FIXME: what if supplementary alignments overlap?
                baseUnmapped_sum = ((int64_t)seq->seq.l - baseCorrect_sum - baseIncorrect_sum > 0 ? (int64_t)seq->seq.l - baseCorrect_sum - baseIncorrect_sum : 0);
                //
                resBaseCorrect += baseCorrect_sum;
                resBaseIncorrect += baseIncorrect_sum;
                resBaseUnmapped += baseUnmapped_sum;
            }
            else
            {
                resBaseIncorrect += seq->seq.l;
            }
        }
        else
        {
            resBaseUnmapped += seq->seq.l;
        }

        // int rank;
        // int64_t baseCorrect;
        // int64_t baseIncorrect;
        // int64_t baseUnmapped;
        // vector<pair<int64_t, int64_t> > intervals;

        // rank = calcRank(m5True, mapList);
        // if(rank <= 0 || rank > mappedRankUB)
        // {
        //  resBaseUnmapped += seq->seq.l;
        // }
        // else // [1 .. mappedRandUB]
        // {
        //  resReadMapped++;
        //  if(rank <= correctRankUB)
        //  {
        //      resReadCorrect++;
        //      calcIntervals(m5True, intervals);
        //      calcBaseStat(mapList[rank - 1], intervals, baseCorrect, baseIncorrect, baseUnmapped);
        //      resBaseCorrect += baseCorrect;
        //      resBaseIncorrect += baseIncorrect;
        //      resBaseUnmapped += baseUnmapped;
        //  }
        //  else // (correctRankUB, mappedRankUB]
        //  {
        //      resBaseIncorrect += seq->seq.l;
        //  }
        // }
    }

    double resSensitivity = (double)resBaseCorrect / (resBaseCorrect + resBaseIncorrect + resBaseUnmapped) * 100;
    double resPrecision = (double)resBaseCorrect / (resBaseCorrect + resBaseIncorrect) * 100;

    // print results
    // fprintf(stdout, "%10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "rTot", "rMap", "rCor", "bTot", "bCor", "bInc", "bUnmap", "sens", "prec");
    // fprintf(stdout, "%10lld %10lld %10lld %10lld %10lld %10lld %10lld %10.2lf %10.2lf\n", resReadTotal, resReadMapped, resReadCorrect, resBaseTotal, resBaseCorrect, resBaseIncorrect, resBaseUnmapped, resSensitivity, resPrecision);
    // cout<< resReadTotal << " & "
    //  << resReadMapped << " & "
    //  << resReadCorrect << " & "
    //  << resBaseTotal << " & "
    //  << resBaseCorrect << " & "
    //  << resBaseIncorrect << " & "
    //  << resBaseUnmapped << " & "
    //  << resSensitivity << " & "
    //  << resPrecision << endl;

    fprintf(stdout, "%7lld & %7lld & %7lld & %10lld & %10lld & %10lld & %10lld & %6.2lf & %6.2lf\n", 
        resReadTotal, resReadMapped, resReadCorrect, resBaseTotal, resBaseCorrect, resBaseIncorrect, resBaseUnmapped, resSensitivity, resPrecision);

    kseq_destroy(seq);
    gzclose(fileRead);

    return EXIT_SUCCESS;
}
