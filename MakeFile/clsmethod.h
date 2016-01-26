#ifndef CLSMETHOLD_H
#define CLSMETHOLD_H

#include <string>
#include "clsfastareader.h"
#include "clsbasealgorithm.h"
#include "api/BamReader.h"
#include "clsmultithread.h"
using namespace BamTools;
using namespace std;

//Method -->Encode All of Reads
//Method -->Decode All of Reads

struct St_Fastq
{
    string strName;
    string strComments;
    string strSeq;
    string strQuality;

    St_Fastq()
    {
        Init();
    }

    void Init()
    {
        strName = "";
        strComments = "";
        strSeq = "";
        strQuality = "";
    }
};

class ClsReadsCoding
{
public:
    ClsReadsCoding();
    ~ClsReadsCoding();

private:


public:
    void EnCodeReads();
    void DeCodeReads();
};

struct St_KmerPos
{
    string::size_type iPos;
    bool bRC;
    St_KmerPos():iPos(0), bRC(false)
    {}

    St_KmerPos(string::size_type iV1, bool bV2)
    {
        iPos = iV1;
        bRC = bV2;
    }
};

class ClsMethod
{
private:
    ClsReadsCoding* m_pReadsCoding;    
    vector<St_Fastq> m_vHighQualityFastq; //A Group Fastq Sequence Code  --> No reads Contain N

public:
    vector<St_Fastq> m_vFastq; //A Group Fastq Sequence Code

public:
    ClsMethod();
    ~ClsMethod();
    void EncodeFastqFile(string strReadsFaPath);

public:
    void Init();
    void ReadFastqFile(string strPath, vector<St_Fastq>& vFastq);
    void FilterNItems(); // filter the reads which contain N --> Go!!!
    void SaveHighQualityFastqToFile(string strFqPath, string strFaPath); // Do not Contain Comments

    string GetSolidKmer(string strReadsFaPath, int iKmerLen,
                        string strAnchorKmerFileName = "AnchorKmer.fa",
                        bool bUseFixThreshould = false, int iCoverageThreshold=10); // by jellyfish

    void HuffmanEncoding(string strAnchorKmerPath, string strSeqFaPath, int iKmerNum);
    void HuffmanDecoding();

    //Transfer the case from non-reference based to reference based
    void CodingByKmerAlignment(string strAnchorKmerPath, string strSeqFaPath, int iKmerLen);

private:
    void PrepareData(string strAnchorKmerPath, string strReadsFaPath,
                     string& strContigFa, string& strMapReadsBam, string& strUnmapReadsFa);
    void ReferenceBasedEncode(string strContigFa, string strMapReadsBam);
    void RegularEncode();

    void GetvGens(vector<St_Fasta>& vKmer, vector<string>& vGens,
                  vector<string::size_type>& vDelta, vector<bool>& vRCStatus,
                  string strSeqFaPath, int iKmerNum);

private:
    string GetUnMappedReads(string strReadsFa, string strRefFa, BamReader* pBamReader, string strRoot);

    void FirstRoundContigCreate(string strRoot, ClsFastaReader* pFastaReader, string strAnchorKmerFa,
                                BamReader* pBamReader, int iVelvetKmerLen,
                                vector<St_Fasta>& vContigFasta,
                                int& iTotalNumKmer, int& iNumContigFirstTime,
                                string& strUnassembledKmerFa, string& strUnassembledKmerFq);

    void SecondRoundContigCreate(string strRoot, ClsFastaReader* pFastaReader, string strAnchorKmerFa,
                                 BamReader* pBamReader, int iVelvetKmerLen,
                                 vector<St_Fasta>& vContigFasta, int& iNumContigSecondTime,
                                 string& strUnassembledKmerFa, string& strUnassembledKmerFq);

    void ThirdRoundContigCreate(string strRoot, ClsFastaReader* pFastaReader, string strAnchorKmerFa,
                                BamReader* pBamReader, int iVelvetKmerLen,
                                vector<St_Fasta>& vContigFasta, int& iNumContigThirdTime,
                                string& strUnassembledKmerFa, string& strUnassembledKmerFq);

    void MapByOrgReads(string strRoot, string strReadsFa, string strRef,
                       ClsFastaReader* pFastaReader, BamReader* pBamReader,
                       vector<St_Fasta>& vUnmappedReads, map<int, string>& mpType,
                       int& iFirstTypeMap_OrgReads, int& iSecondTypeMap_OrgReads,
                       int& iThirdTypeMap_OrgReads, int& iMockTypeMap_OrgReads);

    void CreateMiniReads(string strRoot, vector<St_Fasta>& vOrgReads,
                         string& strMiniReadsPathFa, string& strMiniReadsPathFq);

    void MapByMiniReads(string strRoot, string strRef, string strMiniReadsFq,
                        BamReader* pBamReader, map<int, string> mpType, int iTotalMinReadsNum,
                        int& iFirstTypeMap, int& iSecondTypeMap, int& iThirdTypeMap,
                        int& iMockTypeMap, vector<St_Fasta>& vUnmappedMiniReads,
                        vector<St_Fasta>& vMappedMiniReads);

    void MappRemainingMiniReads(string strRoot, vector<St_Fasta>& vUnmappedMiniReads,
                                BamReader* pBamReader);

    void MultiThredMatch(); //--> 明天搞起 //使用多线程match

public:
    //This is for the remaining part of mini reads compression, --> which fail to be aligned by reference file
    void UnAlignedMiniReadsCompression(string strFasta, string strFastq, int iKmerLen);
    void UnAlignedMiniReadsComprByHuffmanCoding(string strFasta, int iKmerLen=12);
};

#endif // CLSMETHOLD_H
