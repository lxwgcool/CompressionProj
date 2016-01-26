#ifndef CLSMULTITHREAD_H
#define CLSMULTITHREAD_H
#include <map>
#include <string>
#include "KmerUtils.h"
using namespace std;


struct St_KmerInfo
{
    string strSeq;
    unsigned int iID;
    int iNum;
    St_KmerInfo():strSeq(""), iID(0), iNum(0)
    {}
};

struct St_SuperKmer
{
    St_KmerInfo stInfo;
    vector<St_KmerInfo> vSimilar;
    map<KmerTypeShort, int> mpSimilar;
    map<string, int> mpOrgSimilar;

    int GetFreqOrg()
    {
        return GetSizeOrg() + stInfo.iNum;
    }

    void BuildMapOrg()
    {
        mpOrgSimilar.clear();
        if(stInfo.strSeq == "")
            return;
        for(int iPos = 0; iPos < stInfo.strSeq.length(); iPos++)
        {
            string strPos = stInfo.strSeq.substr(iPos, 1);

            //For "A"
            mpOrgSimilar[stInfo.strSeq.replace(iPos, 1, "A")] = 0;
            //For "T"
            mpOrgSimilar[stInfo.strSeq.replace(iPos, 1, "T")] = 0;
            //For "G"
            mpOrgSimilar[stInfo.strSeq.replace(iPos, 1, "G")] = 0;
            //For "C"
            mpOrgSimilar[stInfo.strSeq.replace(iPos, 1, "C")] = 0;

            //Recover the original value
            stInfo.strSeq.replace(iPos, 1, strPos.c_str());
        }
    }

    void BuildMap()
    {
        mpSimilar.clear();
        if(stInfo.strSeq == "")
            return;

        for(int iPos = 0; iPos < stInfo.strSeq.length(); iPos++)
        {
            string strPos = stInfo.strSeq.substr(iPos, 1);
            unsigned int iLen = stInfo.strSeq.length();
            //for "A"
            string strCurSimilarSeq = stInfo.strSeq.replace(iPos, 1, "A");
            int iKmerPosStart = 0;
            KmerTypeShort vValue = 0;
            FormKmerTypeShortSeg(strCurSimilarSeq.c_str(), iKmerPosStart, iLen, vValue);
            mpSimilar[vValue] = 0;
            //iKmerPosStart++;

            //For "T"
            strCurSimilarSeq = stInfo.strSeq.replace(iPos, 1, "T");
            FormKmerTypeShortSeg(strCurSimilarSeq.c_str(), iKmerPosStart, iLen, vValue);
            mpSimilar[vValue] = 0;
            //iKmerPosStart++;

            //For "G"
            strCurSimilarSeq = stInfo.strSeq.replace(iPos, 1, "G");
            FormKmerTypeShortSeg(strCurSimilarSeq.c_str(), iKmerPosStart, iLen, vValue);
            mpSimilar[vValue] = 0;
            //iKmerPosStart++;

            //For "C"
            strCurSimilarSeq = stInfo.strSeq.replace(iPos, 1, "C");
            FormKmerTypeShortSeg(strCurSimilarSeq.c_str(), iKmerPosStart, iLen, vValue);
            mpSimilar[vValue] = 0;

            //Back to the original value
            stInfo.strSeq.replace(iPos, 1, strPos.c_str());
        }
    }

    int GetSize()
    {
        int iSize = 0;
        for(map<KmerTypeShort, int>::iterator itr = mpSimilar.begin(); itr != mpSimilar.end(); itr++)
        {
            if(itr->second != 0)
                iSize++;
        }
        return iSize;
    }

    int GetSizeOrg()
    {
        int iSize = 0;
        for(map<string, int>::iterator itr = mpOrgSimilar.begin(); itr != mpOrgSimilar.end(); itr++)
        {
            if(itr->second != 0)
                iSize++;
        }
        return iSize;
    }
};

struct St_SimilarKmer
{
    vector<St_KmerInfo> vRegKmer;
    map<string, int> vSuperKmer;
};

class  ClsThreadFunc
{
public:
    static void* SimilarDetect(void* args);
};

class ClsMultiThread
{
public:
    ClsMultiThread();

public:
    void FindSimilarKmer(map<string, int>& mpCombine, vector<St_KmerInfo>& vKmerRegular);
};

#endif // CLSMULTITHREAD_H
