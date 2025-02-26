#ifndef CLSBWA_H
#define CLSBWA_H
#include "string"
using namespace std;

class ClsBWA
{
public:
    static ClsBWA& GetInstance()
    {
        static ClsBWA instance;
        return instance;
    }
private:
    ClsBWA(){}
    ClsBWA(const ClsBWA&);
    ClsBWA& operator = (const ClsBWA&);

public:
    //strBamFolder
    string CreateBamBySingleReads(string& strRefPath, string& strReadsPath,
                                  string strBamFileName = "", string strBamFolder="",
                                  bool bMapAll = false, bool bLooseMatch = false,
                                  bool bSortByName = false, int iThreadNum = 2);

    string CreateBamByPEReads(string& strRefPath, string& strReads1Path, string& strReads2Path,
                              bool bSortByName = false, string strBamFolder = "",
                              bool bMapAll = false, bool bLooseMatch = false,
                              string strBamFileNameCustmize = "", int iThreadNum=2);
    //the long reads could be took into consideration
    //bwa bwasw ref.fa long_read.fq > aln.sam
};

#endif // CLSBWA_H
