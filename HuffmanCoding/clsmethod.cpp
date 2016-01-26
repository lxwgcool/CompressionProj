#include "clsmethod.h"
#include <iostream>
#include <fstream>

#include "clsvelvet.h"
#include "clsbwa.h"

#include <zlib.h>
#include <stdio.h>
#include "./FastqFileParse/kseq.h"
KSEQ_INIT(gzFile, gzread)

#include <algorithm> //For function: "sort"

//For Bit File Write --->
#include "./DNAZip/output.h"
#include "./DNAZip/huffman.h"
//<---

//For Decompress the Bit File
#include "./DNAZip/input.h"

#include "KmerUtils.h"
#include "math.h"

#define TWOBITPERNEU

//Reads Coding
ClsReadsCoding::ClsReadsCoding()
{}

ClsReadsCoding::~ClsReadsCoding()
{}

//*************************************
//************ClsMethod****************
//*************************************
ClsMethod::ClsMethod()
{
    m_pReadsCoding = NULL;
}

ClsMethod::~ClsMethod()
{
    if(m_pReadsCoding)
    {
        delete m_pReadsCoding;
        m_pReadsCoding = NULL;
    }

    m_vFastq.clear();
}

void ClsMethod::Init()
{
    if(m_pReadsCoding)
    {
        m_pReadsCoding = new ClsReadsCoding();
    }
}

void ClsMethod::ReadFastqFile(string strPath, vector<St_Fastq>& vFastq)
{
    vFastq.clear();
    St_Fastq stFastq;

    gzFile fp;
    kseq_t *seq;
    int l;

    if(access(strPath.c_str(), 0) != 0) //File do not existed
    {
        cout << "Error: File does not existed!" << endl;
        return;
    }

    fp = gzopen(strPath.c_str(), "r");
    seq = kseq_init(fp);
    while ((l = kseq_read(seq)) >= 0)
    {
        stFastq.Init();
        //Record Name
        stFastq.strName = seq->name.s;
        //Record Comments
        if (seq->comment.l)
            stFastq.strComments = seq->comment.s;
        //Record Sequence
        stFastq.strSeq = seq->seq.s;
        //Record Quality
        if (seq->qual.l)
            stFastq.strQuality = seq->qual.s;
        vFastq.push_back(stFastq);
    }
    cout << IntToStr(vFastq.size()) << endl;
    kseq_destroy(seq);
    gzclose(fp);
    return;
}

void ClsMethod::FilterNItems() // filter the reads which contain N
{
    if(m_vFastq.empty())
        return;
    m_vHighQualityFastq.clear();
    for(vector<St_Fastq>::iterator itr = m_vFastq.begin(); itr != m_vFastq.end(); itr++)
    {
        if(itr->strSeq.find('N') != string::npos ||
           itr->strSeq.find('n') != string::npos)
            continue;
        m_vHighQualityFastq.push_back(*itr);
    }
    cout << IntToStr(m_vHighQualityFastq.size()) << endl;
}

void ClsMethod::SaveHighQualityFastqToFile(string strFqPath, string strFaPath) // Do not Contain Comments
{
    bool bFqExist = false;
    if(access(strFqPath.c_str(), 0) == 0)
    {
        cout << "High Quality Fq File existed!" << endl;
        bFqExist = true;
    }

    bool bFaExist = false;
    if(access(strFaPath.c_str(), 0) == 0)
    {
        cout << "High Quality Fa File existed!" << endl;
        bFaExist = true;
    }

    if(bFqExist && bFaExist)
    {
        cout << "Both Exsited" << endl;
        return;
    }

    //Save Fq
    ofstream ofsFq;
    if(!bFqExist)
        ofsFq.open(strFqPath.c_str());

    ofstream ofsFa;
    if(!bFaExist)
        ofsFa.open(strFaPath.c_str());

    for(vector<St_Fastq>::iterator itr = m_vHighQualityFastq.begin();
        itr != m_vHighQualityFastq.end(); itr++)
    {
        if(!bFqExist)
        {
            //write name
            ofsFq << "@" << itr->strName << endl;
            //write sequence
            ofsFq << itr->strSeq << endl;
            //write quality
            ofsFq <<"+" << endl;
            ofsFq << itr->strQuality << endl;
        }

        if(!bFaExist)
        {
            //write name
            ofsFa << ">" << itr->strName << endl;
            //write sequence
            ofsFa << itr->strSeq << endl;
        }
    }

    if(!bFqExist)
        ofsFq.close();

    if(!bFaExist)
        ofsFa.close();
}

void ClsMethod::EncodeFastqFile(string strReadsFaPath)
{
    //1: Get K-mers and Calculate its frequency // by jellyfish  (kmc2 is not work)
    int iKmerNum = 12;
    string strAnchorKmerFa = GetSolidKmer(strReadsFaPath, iKmerNum);

    //2: Divide it to two seperate parts: one is successful mapping, another is non-mapping
    string strContigFa, strMapReadsBam, strUnmapReadsFa;
    PrepareData(strAnchorKmerFa, strReadsFaPath,
                strContigFa, strMapReadsBam, strUnmapReadsFa);

    //3: Compress the first part
    ReferenceBasedEncode(strContigFa, strMapReadsBam);

    //4: Compress the second part
    RegularEncode();
}

string ClsMethod::GetSolidKmer(string strReadsFaPath, int iKmerLen,
                               string strAnchorKmerFileName,
                               bool bUseFixThreshould, int iCoverageThreshold)
{
    string strRoot = GetHigherFolderPath(GetCurExeFolderPath());

    //Target is get the AnchorKmer.fa
    string strAnchorKmerFasta = strRoot + "ThirdPartyTools/Jellyfish/data/" + strAnchorKmerFileName;
    if(access(strAnchorKmerFasta.c_str(), 0) == 0)
    {
        cout << "Anchor Kmer Exsited!" << endl;
        return strAnchorKmerFasta;
    }

    //1: Run JellyFish, create mer_count.jf file -->
    string strJellyFishPath = strRoot + "ThirdPartyTools/Jellyfish/bin/jellyfish";
    string strKmerLen = IntToStr(iKmerLen);
    string strThreadsNum = IntToStr(4);
    string strHushElementNum = "150M";
    string strMerCountFile = strRoot + "ThirdPartyTools/Jellyfish/data/mer_counts.jf";


    string strCmd = strJellyFishPath + " count -m " + strKmerLen + " -s " + strHushElementNum +
                    " -o " + strMerCountFile +
                    " -t " + strThreadsNum + " -C " + strReadsFaPath;
    if(access(strMerCountFile.c_str(), 0) != 0) // do not exist
        system(strCmd.c_str());

    //2: dump the kmer information out of mer_count.jf
    string strDumpFile = strRoot + "ThirdPartyTools/Jellyfish/data/mer_counts_dumps.fa";
    strCmd = strJellyFishPath + " dump " + strMerCountFile +
             " > " + strDumpFile;
    if(access(strDumpFile.c_str(), 0) != 0) // do not exist
        system(strCmd.c_str());

    //3: Create Histogram
    string strHistoPath = strRoot + "ThirdPartyTools/Jellyfish/data/histo.ini";
    strCmd = strJellyFishPath + " histo " + strMerCountFile + " > " + strHistoPath;
    system(strCmd.c_str());

    //Read Histogram file
    ifstream ifs;
    ifs.open(strHistoPath.c_str());
    string strLine;
    int iSumFreq = 0;
    int iSumKmerNum = 0;
    getline(ifs, strLine);
    while(!ifs.eof())
    {
        getline(ifs, strLine);
        string::size_type stpzPos = strLine.find(" ");
        string strFreq = strLine.substr(0, stpzPos);
        string strKmerLen = strLine.substr(stpzPos + 1, strLine.length() - stpzPos);
        iSumKmerNum += atoi(strKmerLen.c_str());
        iSumFreq += atoi(strFreq.c_str()) * atoi(strKmerLen.c_str());
    }
    float fAverage = (float)iSumFreq / iSumKmerNum;
    cout << "Average Coverage: " << FloatToStr(fAverage);
    ifs.close();

    //4: Load Fa File and Calculate the sum of kmer and the sum of appeared times.
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    vector<St_Fasta> vKmerFasta;
    pFastaReader->ReadFastaRegular(strDumpFile, vKmerFasta); // collect the kmer which frequency larger than fAverage
    cout << IntToStr(vKmerFasta.size()) << endl;
    vector<St_Fasta> vHighCovKmerFasta;
    for(vector<St_Fasta>::iterator itr = vKmerFasta.begin(); itr != vKmerFasta.end(); itr++)
    {
        if(itr->strSeq == "")
            continue;
        if(atoi(itr->strName.c_str()) > fAverage)
            vHighCovKmerFasta.push_back(*itr);
    }
    vKmerFasta.clear();
    cout << IntToStr(vHighCovKmerFasta.size()) << endl;
    delete pFastaReader;
    pFastaReader = NULL;

    //Save those new Anchor Kmer as a new fasta file
    ofstream ofsAnchorKmer;
    ofsAnchorKmer.open(strAnchorKmerFasta.c_str());
    for(vector<St_Fasta>::iterator itr = vHighCovKmerFasta.begin();
        itr != vHighCovKmerFasta.end(); itr++)
    {
        ofsAnchorKmer << ">" << itr->strName << endl;
        ofsAnchorKmer << itr->strSeq << endl;
    }
    ofsAnchorKmer.close();
    vHighCovKmerFasta.clear();
    return strAnchorKmerFasta;
}

void ClsMethod::PrepareData(string strAnchorKmerPath, string strReadsFaPath,
                            string& strContigFa, string& strMapReadsBam, string& strUnmapReadsFa)
{
    string strRoot = GetHigherFolderPath(GetCurExeFolderPath());

    //0: transfer the fasta reads to fastq --> Then reason: we need to keep the quality vlaue is
    //as the same as the mapped reads and unmapped reads
    string strReadsFastqPath = strRoot + "Output/FullReadsFastq/fullreads.fq";
    string strCmd = "perl " + strRoot + "/ThirdPartyTools/fasta_to_fastq.pl " +
                    strReadsFaPath + " > " + strReadsFastqPath;
    system(strCmd.c_str());

    //1: velvet (create the contigs)
    strContigFa = ClsVelvet::GetInstance().LocalAssembly(strAnchorKmerPath, 40, 21);

    //2: align those reads to those contigs, and create the bam file
    //(1)Transfer from fasta to fastq
    string strAnchorKmerFqPath = strRoot + "ThirdPartyTools/Jellyfish/data/AnchorKmer.fq";
    strCmd = "perl " + strRoot + "/ThirdPartyTools/fasta_to_fastq.pl " +
             strAnchorKmerPath + " > " + strAnchorKmerFqPath;
    system(strCmd.c_str());
    //(2): align those file back to the original reference geno to create the bam file
    string strBamFileFolder = "Output/BamFile/";
    string strBamFilePath = ClsBWA::GetInstance().CreateBamBySingleReads(strContigFa,
                                                                         strReadsFaPath,
                                                                         "",
                                                                         strBamFileFolder,
                                                                         false, true);
    //3: just pick up the successfully mapping one
    BamReader* pBamReader = new BamReader();
    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;
    //解析sort by name 之后的文件
    vector<St_Fasta> vMappedReads;
    vector<St_Fasta> vUnMappedReads;
    St_Fasta stFasta;
    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        if(al.QueryBases == "")
            continue;

        stFasta.strName = al.Name;
        stFasta.strSeq = al.QueryBases;
        if(al.IsMapped()) // becase it has been sorted        
            vMappedReads.push_back(stFasta);
        else
            vUnMappedReads.push_back(stFasta);
    }
    delete pBamReader;
    pBamReader = NULL;

    //4: regenerate the bam file again by the good reads --> Create Bam File
    //(1) Save Mapped Reads to Fastq File
    string strMappedReadsFa = strRoot + "Output/MappedReadsFastq/MappedReads.fa";
    string strMappedReadsFq = strRoot + "Output/MappedReadsFastq/MappedReads.fq";
    ofstream ofs;
    ofs.open(strMappedReadsFa.c_str());
    for(vector<St_Fasta>::iterator itr = vMappedReads.begin(); itr != vMappedReads.end(); itr++)
    {
        ofs << ">" << itr->strName << endl;
        ofs << itr->strSeq << endl;
    }
    ofs.close();
    strCmd = "perl " + strRoot + "/ThirdPartyTools/fasta_to_fastq.pl " +
             strMappedReadsFa + " > " + strMappedReadsFq;
    system(strCmd.c_str());

    //(2) Align this fastq file to reference to create bam file
    strBamFileFolder = "Output/MappedBamFile/";
    strMapReadsBam = ClsBWA::GetInstance().CreateBamBySingleReads(strContigFa,
                                                                  strMappedReadsFq,
                                                                  "",
                                                                  strBamFileFolder,
                                                                  false, true);

    //5: collect the unmapping reads and create a special single file -->Create Fasta File
    //(1) Save as fasta file
    strUnmapReadsFa = strRoot + "Output/UnmappedReadsFasta/UnmappedReads.fa";
    string strUnmapReadsFq = strRoot + "Output/UnmappedReadsFasta/UnmappedReads.fq";
    ofs.open(strUnmapReadsFa.c_str());
    for(vector<St_Fasta>::iterator itr = vUnMappedReads.begin();
        itr != vUnMappedReads.end(); itr++)
    {
        //Name
        ofs << ">" << itr->strName << endl;
        //Sequence
        ofs << itr->strSeq << endl;
    }
    //(2) Save as fastq file
    strCmd = "perl " + strRoot + "/ThirdPartyTools/fasta_to_fastq.pl " +
             strUnmapReadsFa + " > " + strUnmapReadsFq;
    system(strCmd.c_str());

    ofs.close();
}

void ClsMethod::ReferenceBasedEncode(string strContigFa, string strMapReadsBam)
{
    //Use the third party tool to compress the reference based
    string strRoot = GetHigherFolderPath(GetCurExeFolderPath());
    string strCramToolPath = strRoot + "ThirdPartyTools/cramtools-3.0/cramtools-3.0.jar";
    string strOutputCramPath = strRoot + "Output/CompressResult/MapReadsCompress.Cram";
    string strCmd = "java -jar " + strCramToolPath +
                    " cram --input-bam-file " + strMapReadsBam +
                    " --reference-fasta-file " + strContigFa +
                    " --output-cram-file " + strOutputCramPath;
    system(strCmd.c_str());
}

void ClsMethod::RegularEncode()
{
    //compress the remaning data use gzip or bzip2
}

//Some additional idea: create the reference and use the reference based encode to compress all of those data directly.

////////////////////////////////////////////////////////////////////////////////////////////////////////

bool sortfunctionfasta(St_Fasta stFa1, St_Fasta stFa2)
{
    return atoi(stFa1.strName.c_str()) > atoi(stFa2.strName.c_str());
}
//Make Huffman Code
void ClsMethod::HuffmanEncoding(string strAnchorKmerPath, string strSeqFaPath, int iKmerNum)
{
    //1: Read Anchor Kmer (Fasta File)
    vector<St_Fasta> vKmerFasta;
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    pFastaReader->ReadFastaRegular(strAnchorKmerPath, vKmerFasta);
    sort(vKmerFasta.begin(), vKmerFasta.end(), sortfunctionfasta); //我们在这里是从大到小进行排序的

    //2: 再对这些Kmer筛选，不要这么多 我们只要前1024个
    int iMaxNum = 512;
    if(vKmerFasta.size() > iMaxNum)
    {
        vKmerFasta.erase(vKmerFasta.begin() + iMaxNum, vKmerFasta.end());
    }

    //3:计算这些Kmer在原始序列中的位置，并且得到剩下的ACTG出现的次数
    //先走一步看一步的存存看看，看看村后的结果跟存储之前的差异多大
    vector<string> vGens;
    vector<string::size_type> vDelta; //Record All of Delta Pos of Kmer
    vector<bool> vRCStatus;
    GetvGens(vKmerFasta, vGens, vDelta, vRCStatus,
             strSeqFaPath, iKmerNum); // 这里的Kmer Fasta已经是排序好了的了 --> 前256 个
                                               // 走过这个函数，vKmerFasta 和 vGens 都发生了变化

    //4.1:Save Delta File
    string strRootPath =  GetHigherFolderPath(GetCurExeFolderPath());
    bit_file_c bf;
    string strDeltaPath = strRootPath + "Test/CompressedFile/Delta.bit";
    bf.Open(strDeltaPath.c_str(), BF_WRITE );
    for(vector<string::size_type>::iterator itr = vDelta.begin(); itr != vDelta.end(); itr++)
    {
        writeBitVINT(bf, *itr);
    }
    bf.Close();

    //4.2 Save RC status
    string strDestFile = strRootPath + "Test/CompressedFile/RCStatus.bit";
    bf.Open( strDestFile.c_str(), BF_WRITE );
    writeBitVINT( bf, vRCStatus.size() );
    writeBits( bf, &vRCStatus );
    bf.ByteAlign();
    bf.Close();

    //我们最后将Anchor Kmer save到bit file中
    //5: Save Sequence Info and Frequency Info into the File
    // Question: Why the size been decreased ??!!!         
    strDestFile = strRootPath + "Test/CompressedFile/AnchorKmer.bit";
    bf.Open(strDestFile.c_str(), BF_WRITE );

    writeBitVINT(bf, vKmerFasta.size());
    map<string, int> mpKmerFreq;
    cout << "Anchor Kmer Numer: " << IntToStr(vKmerFasta.size()) << endl;
    //int iCount = 0;
    for(vector<St_Fasta>::iterator itr = vKmerFasta.begin();
        itr != vKmerFasta.end(); itr++)
    {
        //if(iCount > 1023)
        //    break;
        //set the value for mpKmerFreq;
        mpKmerFreq[itr->strSeq] = atoi(itr->strName.c_str());

        //record result
        writeString(bf,  itr->strSeq);
        writeBitVINT(bf, atoi(itr->strName.c_str()));

        //iCount++;
    }
    bf.Close();
    delete pFastaReader;
    pFastaReader = NULL;

    //6:Huffman Tree
    //(1) Init Huffman Tree -->
    Hufftree<string, int>* hufftree = NULL;
    hufftree = new Hufftree<string, int>(mpKmerFreq.begin(), mpKmerFreq.end());

#ifdef TWOBITPERNEU
    vector<bool> vCode;
    vCode.resize(2);
    //For A //00
    vCode[0] = false;
    vCode[1] = false;
    hufftree->UpdateEncoding("A", vCode);

    //For C //01
    vCode[0] = false;
    vCode[1] = true;
    hufftree->UpdateEncoding("C", vCode);

    //For G //10
    vCode[0] = true;
    vCode[1] = false;
    hufftree->UpdateEncoding("G", vCode);

    //For T //11
    vCode[0] = true;
    vCode[1] = true;
    hufftree->UpdateEncoding("T", vCode);
#endif

    //(2) Code each part -->
    //输出几个相应的编码试试看如何-->
    vector<bool> vValue1 = hufftree->encode("A");
    for(std::vector<bool>::iterator itr = vValue1.begin(); itr != vValue1.end(); itr++)
    {
        if(*itr)
            cout << "1";
        else
            cout << "0";
    }
    cout << endl;
    vector<bool> vValue2 = hufftree->encode("T");
    for(std::vector<bool>::iterator itr = vValue2.begin(); itr != vValue2.end(); itr++)
    {
        if(*itr)
            cout << "1";
        else
            cout << "0";
    }
    cout << endl;
    vector<bool> vValue3 = hufftree->encode("G");
    for(std::vector<bool>::iterator itr = vValue3.begin(); itr != vValue3.end(); itr++)
    {
        if(*itr)
            cout << "1";
        else
            cout << "0";
    }
    cout << endl;
    vector<bool> vValue4 = hufftree->encode("TTAATAAA");
    for(std::vector<bool>::iterator itr = vValue4.begin(); itr != vValue4.end(); itr++)
    {
        if(*itr)
            cout << "1";
        else
            cout << "0";
    }
    cout << endl;
    vector<bool> vValue5 = hufftree->encode("AAAATAAA");
    for(std::vector<bool>::iterator itr = vValue5.begin(); itr != vValue5.end(); itr++)
    {
        if(*itr)
            cout << "1";
        else
            cout << "0";
    }
    cout << endl;

    //(3) Encoding
    vector<bool> encodedLongGens = hufftree->encode(vGens.begin(), vGens.end());
    cout << "Bool Size: " << IntToStr(encodedLongGens.size()) << endl;

    //(4) Save File
    strDestFile = strRootPath + "Test/CompressedFile/OrgFilCompress.bit";
    bit_file_c destBf;
    destBf.Open( strDestFile.c_str(), BF_WRITE );
    writeBitVINT( destBf, encodedLongGens.size() );
    writeBits( destBf, &encodedLongGens );
    destBf.ByteAlign();

    destBf.Close();

    //(5) Release memory
    vGens.clear();

    delete hufftree;
    hufftree = NULL;
}

void ClsMethod::HuffmanDecoding()  // 我现在相当于之前的简化版本 --> 可以类比于只要考虑一个insertion
{
    string strRootPath = ::GetHigherFolderPath(get_current_dir_name());
    string strAncherMerBitPath = strRootPath + "Test/CompressedFile/AnchorKmer.bit";
    string strCompressSeqBitPath = strRootPath + "Test/CompressedFile/OrgFilCompress.bit";

    //1: Load AncherMerBit to Build Huffman Tree -->
    map<string, int> frequencies;

    bit_file_c bf;
    bf.Open(strAncherMerBitPath.c_str(), BF_READ);

    unsigned num = readBitVINT(bf);
    for( unsigned i = 0; i < num; i++ )
    {
        string gens = readString(bf);
        int freq = readBitVINT(bf);
        frequencies[ gens ] = freq;
    }
    bf.Close();

    //2: Create Huffman Tree
    Hufftree<string, int>* hufftree = NULL;
    hufftree = new Hufftree<string, int>(frequencies.begin(), frequencies.end());

    //3: Load the compression File
    bf.Open(strCompressSeqBitPath.c_str(), BF_READ);

    //4: Read Bit Length and Bit Content
    vector<bool> genBits;
    unsigned bitCnt = readBitVINT(bf);
    readBits(bf, bitCnt, &genBits);

    //5: Huffman Decode to get the original sequence value
    vector<string> vSeq;
    hufftree->decode(genBits, std::back_inserter(vSeq));
    bf.Close();

    //6: Recover this vector<string> to one long string, and save it to a single file
    string strSeq = "";
    for(vector<string>::iterator itr = vSeq.begin(); itr != vSeq.end(); itr++)
    {
        strSeq += *itr;
    }
    string strRecoverFilePath = strRootPath + "Test/CompressedFile/Recover.txt";
    ofstream ofs;
    ofs.open(strRecoverFilePath.c_str());
    ::DisplayString(ofs, strSeq, 70);
    ofs.close();
}

bool sortfunctionkmerpos(St_KmerPos st1, St_KmerPos st2) //From small to large
{
    return st1.iPos < st2.iPos;
}

void ClsMethod::GetvGens(vector<St_Fasta>& vKmer, vector<string>& vGens,
                         vector<string::size_type>& vDelta, vector<bool>& vRCStatus,
                         string strSeqFaPath, int iKmerNum) // 直接存储用于之后编码的容器
{
    //1: read the fasta file
    vector<St_Fasta> vFasta;
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    pFastaReader->ReadFastaRegular(strSeqFaPath, vFasta); //我们在这里只有一个fasta文件的一个单元项要读，因此直接后面用引用对它进行表示
    string strRefSeq = vFasta[0].strSeq;

    // 先修改里面的内容，然后再进行逐步的解析
    // 在这里怎么去取呢？ 我们记住position就好了 ->
    string strReplace = "";
    for(int i=0; i<iKmerNum; i++)
        strReplace += "$";
    const char* czReplace = strReplace.c_str();
    const int ciNum = iKmerNum;

    vector<St_KmerPos> vPos;
    for(vector<St_Fasta>::iterator itr = vKmer.begin(); itr != vKmer.end(); itr++)
    {
        string& strCurKmer = itr->strSeq;
        string strReverseCmp = ::GetReverseCompelement(strCurKmer);
        string::size_type sztpCurPos = 0;
        //正向序列
        sztpCurPos = strRefSeq.find(strCurKmer.c_str(), sztpCurPos);
        while(sztpCurPos != string::npos)
        {
            vPos.push_back(St_KmerPos(sztpCurPos, false));
            strRefSeq.replace(sztpCurPos, ciNum, czReplace);
            sztpCurPos += ciNum;
            sztpCurPos = strRefSeq.find(strCurKmer.c_str(), sztpCurPos);
        }

        //反向互补序列
        sztpCurPos = 0;
        sztpCurPos = strRefSeq.find(strReverseCmp.c_str(), sztpCurPos);
        while(sztpCurPos != string::npos)
        {
            vPos.push_back(St_KmerPos(sztpCurPos, true));
            strRefSeq.replace(sztpCurPos, ciNum, czReplace);
            sztpCurPos += ciNum;
            sztpCurPos = strRefSeq.find(strReverseCmp.c_str(), sztpCurPos);
        }
    }
    //我们现在已经拿到相应的pos了，队得到的pos进行排序
    sort(vPos.begin(), vPos.end(), sortfunctionkmerpos); // 从小到大排序了
    //(1)Record the Delta Pos -->
    vDelta.clear();
    string::size_type sztpPos = 0;
    for(vector<St_KmerPos>::iterator itr = vPos.begin();
        itr != vPos.end(); itr++)
    {
        vDelta.push_back(itr->iPos - sztpPos);
        vRCStatus.push_back(itr->bRC);
        sztpPos = itr->iPos + iKmerNum; // jump the length of kmer
    }

    strRefSeq = vFasta[0].strSeq;
    //(2)我们还是按照最简单的方法来做
    for(vector<St_KmerPos>::iterator itr = vPos.end()-1;
        itr >= vPos.begin(); itr--) // 我们从后往前去遍历 --->添加特殊符号进去 "$"
    {
        strRefSeq.insert(itr->iPos, "$");
    }

    //遍历所有的每个字符，如果碰到"$"那么就读取从其后面开始的8个字符，反之就直接将这些字符存入容器中
    //统计ATGC的相应的出现的频率 --> 纯次数，不做任何归一化
    vGens.clear();
    map<string, int> mpNucleotide;
    for(string::size_type sztpCurPos = 0; sztpCurPos < strRefSeq.length(); sztpCurPos++)
    {
        if(strRefSeq[sztpCurPos] != '$')
        {
            vGens.push_back(strRefSeq.substr(sztpCurPos, 1));
            mpNucleotide[strRefSeq.substr(sztpCurPos, 1)]++;
        }
        else
        {
            vGens.push_back(strRefSeq.substr(sztpCurPos + 1, ciNum));
            sztpCurPos += ciNum;
        }
    }

    cout << "Size: " << vGens.size() << endl;

#ifndef  TWOBITPERNEU
    //Add Those 4 Regular Nucletide into Kmer -->
    St_Fasta stFasta;
    for(map<string, int>::iterator itr = mpNucleotide.begin(); itr != mpNucleotide.end(); itr++)
    {
        stFasta.strSeq = itr->first;
        stFasta.strName = IntToStr(itr->second);
        vKmer.push_back(stFasta);
    }
#endif
    delete pFastaReader;
    pFastaReader = NULL;
}

/**************************************************************************************
 Compressed the file based on the idea: transfer the case of non-reference to reference
//*************************************************************************************/
//Some fundermentery Function
string ClsMethod::GetUnMappedReads(string strReadsFa, string strRefFa, BamReader* pBamReader,
                                   string strRoot)
{
    string strReadsFq = strRoot + "ThirdPartyTools/Jellyfish/data/AnchorKmer.fq";
    string strCmd = "perl " + strRoot + "/ThirdPartyTools/fasta_to_fastq.pl " +
             strReadsFa + " > " + strReadsFq;
    system(strCmd.c_str());

    string strBamFileFolder = "Output/BamFile/";
    string strBamFilePath = ClsBWA::GetInstance().CreateBamBySingleReads(strRefFa,
                                                                         strReadsFq,
                                                                         "",
                                                                         strBamFileFolder);

    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;
    //解析sort by name 之后的文件
    vector<St_Fasta> vUnMappedReads; // only record the unmapped kmer
    St_Fasta stFasta;
    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        if(al.QueryBases == "")
            continue;

        if(!al.IsMapped())
        {
            stFasta.strName = al.Name;
            stFasta.strSeq = al.QueryBases;
            vUnMappedReads.push_back(stFasta);
        }
    }
    pBamReader->Rewind();

    //Create Fasta File
    string strUnmapReadsFa = strRoot + "Output/UnmappedReadsFasta/UnmappedKmer.fa"; // for velvet
    ofstream ofs;
    ofs.open(strUnmapReadsFa.c_str());
    for(vector<St_Fasta>::iterator itr = vUnMappedReads.begin();
        itr != vUnMappedReads.end(); itr++)
    {
        //Name
        ofs << ">" << itr->strName << endl;
        //Sequence
        ofs << itr->strSeq << endl;
    }
    ofs.close();
    vUnMappedReads.clear();
    return strUnmapReadsFa;
}

//First Round fasta file
void ClsMethod::FirstRoundContigCreate(string strRoot, ClsFastaReader* pFastaReader,
                                       string strAnchorKmerFa,
                                       BamReader* pBamReader, int iVelvetKmerLen,
                                       vector<St_Fasta>& vContigFasta,
                                       int& iTotalNumKmer, int& iNumContigFirstTime,
                                       string& strUnassembledKmerFa, string& strUnassembledKmerFq)
{
    //1: Read Solid(Anchor) Kmer
    vector<St_Fasta> vKmerFasta;
    pFastaReader->ReadFastaRegular(strAnchorKmerFa, vKmerFasta);
    iTotalNumKmer = vKmerFasta.size();

    //2:Make Local Assembly    
    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << "Run Valvet" << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    //erase all of old files -->
    string strCmd = "rm " + strRoot + "/Output/LocalAssembly/*";
    system(strCmd.c_str());

    //<--
    string strContigFa = ClsVelvet::GetInstance().LocalAssembly(strAnchorKmerFa, 60, iVelvetKmerLen);
    pFastaReader->ReadFastaRegular(strContigFa, vContigFasta);
    iNumContigFirstTime = vContigFasta.size();

    //3: change the name of contigs created by this phase
    for(vector<St_Fasta>::iterator itr = vContigFasta.begin(); itr != vContigFasta.end(); itr++)
    {
        itr->strName += "_first";
    }

    //4: Get Un-Assembled Kmer
    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << "Run BWA" << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
    strUnassembledKmerFa = GetUnMappedReads(strAnchorKmerFa, strContigFa, pBamReader, strRoot);// For Velvet

    //5: Transfer Fa to Fq
    strUnassembledKmerFq = strRoot + "Output/UnmappedReadsFasta/UnmappedKmer.fq"; // for bwa
    strCmd = "perl " + strRoot + "/ThirdPartyTools/fasta_to_fastq.pl " +
             strUnassembledKmerFa + " > " + strUnassembledKmerFq;
    system(strCmd.c_str());
}

//Second Round
void ClsMethod::SecondRoundContigCreate( string strRoot, ClsFastaReader* pFastaReader,
                                         string strAnchorKmerFa,
                                         BamReader* pBamReader, int iVelvetKmerLen,
                                         vector<St_Fasta>& vContigFasta, int& iNumContigSecondTime,
                                         string& strUnassembledKmerFa, string& strUnassembledKmerFq)
{
    //1: Local Assembly by Velvet and record the relevant contigs
    //erase all of old files -->
    string strCmd = "rm " + strRoot + "/Output/LocalAssembly/*";
    system(strCmd.c_str());
    //<--

    string strContigFa = ClsVelvet::GetInstance().LocalAssembly(strAnchorKmerFa, 60, iVelvetKmerLen);
    vector<St_Fasta> vContigTemp;
    pFastaReader->ReadFastaRegular(strContigFa, vContigTemp);
    iNumContigSecondTime = vContigTemp.size();
    for(vector<St_Fasta>::iterator itr = vContigTemp.begin(); itr != vContigTemp.end(); itr++)
    {
        itr->strName += "_second";
    }
    vContigFasta.insert(vContigFasta.end(), vContigTemp.begin(), vContigTemp.end());
    vContigTemp.clear();

    //2: Get Un-Assembled Kmer Fasta File
    // For Velvet
    strUnassembledKmerFa = GetUnMappedReads(strAnchorKmerFa, strContigFa, pBamReader, strRoot);

    //3: Transfer Fasta to Fastq
    strUnassembledKmerFq = strRoot + "Output/UnmappedReadsFasta/UnmappedKmer.fq"; // for bwa
    strCmd = "perl " + strRoot + "/ThirdPartyTools/fasta_to_fastq.pl " +
             strUnassembledKmerFa + " > " + strUnassembledKmerFq;
    system(strCmd.c_str());
}

//Third Round
void ClsMethod::ThirdRoundContigCreate( string strRoot, ClsFastaReader* pFastaReader,
                                        string strAnchorKmerFa,
                                        BamReader* pBamReader, int iVelvetKmerLen,
                                        vector<St_Fasta>& vContigFasta, int& iNumContigThirdTime,
                                        string& strUnassembledKmerFa, string& strUnassembledKmerFq)
{
    //1: Local Assembly by Velvet
    //erase all of old files -->
    string strCmd = "rm " + strRoot + "/Output/LocalAssembly/*";
    system(strCmd.c_str());
    //<--
    string strContigFa = ClsVelvet::GetInstance().LocalAssembly(strAnchorKmerFa, 60, iVelvetKmerLen);
    vector<St_Fasta> vContigTemp;
    pFastaReader->ReadFastaRegular(strContigFa, vContigTemp);
    iNumContigThirdTime = vContigTemp.size();
    for(vector<St_Fasta>::iterator itr = vContigTemp.begin(); itr != vContigTemp.end(); itr++)
    {
        itr->strName += "_third";
    }
    vContigFasta.insert(vContigFasta.end(), vContigTemp.begin(), vContigTemp.end());
    vContigTemp.clear();

    //2: Get Un-Assembled Kmer Fasta File
    // For Velvet
    strUnassembledKmerFa = GetUnMappedReads(strAnchorKmerFa, strContigFa, pBamReader, strRoot);

    //3: Transfer Fasta to Fastq
    strUnassembledKmerFq = strRoot + "Output/UnmappedReadsFasta/UnmappedKmer.fq"; // for bwa
    strCmd = "perl " + strRoot + "/ThirdPartyTools/fasta_to_fastq.pl " +
             strUnassembledKmerFa + " > " + strUnassembledKmerFq;
    system(strCmd.c_str());
}

//Forth Type of Contig
void ClsMethod::MapByOrgReads( string strRoot, string strReadsFa, string strRef,
                               ClsFastaReader* pFastaReader, BamReader* pBamReader,
                               vector<St_Fasta>& vUnmappedReads, map<int, string>& mpType,
                               int& iFirstTypeMap_OrgReads, int& iSecondTypeMap_OrgReads,
                               int& iThirdTypeMap_OrgReads, int& iMockTypeMap_OrgReads)
{    
    //Read Reads File -->
    vector<St_Fasta> vReads;
    pFastaReader->ReadFastaRegular(strReadsFa, vReads, false);

    //1: Transfer Fasta to Fq
    string strReadsFq = strRoot + "Output/UnmappedReadsFasta/Reads.fq"; // for bwa
    string strCmd = "perl " + strRoot + "/ThirdPartyTools/fasta_to_fastq.pl " +
                    strReadsFa + " > " + strReadsFq;
    system(strCmd.c_str());

    int i10UnmapPart = 0; int i1020UnmapPart = 0; int i2050UnmapPart = 0; int i50UnmapPart = 0;
    //2:Use BWA for Mapping
    string strBamFileFolder = "Output/BamFile/";
    string strBamFilePath = ClsBWA::GetInstance().CreateBamBySingleReads( strRef,
                                                                          strReadsFq,
                                                                          "",
                                                                          strBamFileFolder);
    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    BamAlignment al;
    St_Fasta stFasta;
    int iZeroBases = 0;
    int iBadRefID = 0;    
    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        if(al.Length == 0)
        {
            iZeroBases++;
            continue;
        }

        if(!al.IsMapped())
        {
            stFasta.strName = al.Name;
            stFasta.strSeq = al.QueryBases;
            vUnmappedReads.push_back(stFasta);
        }
        else // Mapping successfully // 原因可能是在于在mapping的时候，map到了多个位置
        {
            //Mapping Type
            if(mpType[al.RefID] == "first")
                iFirstTypeMap_OrgReads++;
            else if(mpType[al.RefID] == "second")
                iSecondTypeMap_OrgReads++;
            else if(mpType[al.RefID] == "third")
                iThirdTypeMap_OrgReads++;
            else if(mpType[al.RefID] == "mock")
                iMockTypeMap_OrgReads++;
            else
                iBadRefID++;

            //Mapping Ratio
            float fUnmapRatio = (float)(al.Length - ::GetMatchPartLen(al)) / al.Length;
            if(fUnmapRatio < .1)
            {
                i10UnmapPart++;
            }
            else if(fUnmapRatio >= .1 && fUnmapRatio < .2)
            {
                i1020UnmapPart++;
            }
            else if(fUnmapRatio >= .2 && fUnmapRatio < .5)
            {
                i2050UnmapPart++;
            }
            else // this is >= 0.5
            {
                i50UnmapPart++;
            }            
        }
    }

    //percentage of each type of mapping reads
    int iTotalMappedNumber = vReads.size() - vUnmappedReads.size();
    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    cout << "iZeroBases: " << IntToStr(iZeroBases) << endl;
    cout << "iBadRefID: " << IntToStr(iBadRefID) << endl;
    cout << "Total Mapped Number: " << IntToStr(iTotalMappedNumber) << endl;
    cout << "First Type: " << IntToStr(iFirstTypeMap_OrgReads) << " --> " << GetRatio((float)iFirstTypeMap_OrgReads / iTotalMappedNumber) << endl;
    cout << "First Type: " << IntToStr(iSecondTypeMap_OrgReads) << " --> " << GetRatio((float)iSecondTypeMap_OrgReads / iTotalMappedNumber) << endl;
    cout << "First Type: " << IntToStr(iThirdTypeMap_OrgReads) << " --> " << GetRatio((float)iThirdTypeMap_OrgReads / iTotalMappedNumber) << endl;
    cout << "First Type: " << IntToStr(iMockTypeMap_OrgReads) << " --> " << GetRatio((float)iMockTypeMap_OrgReads / iTotalMappedNumber) << endl;
    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;

    //Print Unmap Part
    cout << endl;
    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    cout << "Totoal Number of Mapping : " << IntToStr(iTotalMappedNumber) << endl;
    cout << " < 10% Unmap Part        : " << IntToStr(i10UnmapPart) << " --> " << GetRatio((float)i10UnmapPart / iTotalMappedNumber) << endl;
    cout << " >=10% && <20% Unmap Part: " << IntToStr(i1020UnmapPart) << " --> " << GetRatio((float)i1020UnmapPart / iTotalMappedNumber) << endl;
    cout << " >=20% && <50% Unmap Part: " << IntToStr(i2050UnmapPart) << " --> " << GetRatio((float)i2050UnmapPart / iTotalMappedNumber) << endl;
    cout << " >=50% Unmap Part        : " << IntToStr(i50UnmapPart) << " --> " << GetRatio((float)i50UnmapPart / iTotalMappedNumber) << endl;
    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;

    pBamReader->Rewind();
}

void ClsMethod::CreateMiniReads(string strRoot, vector<St_Fasta>& vOrgReads,
                                string& strMiniReadsPathFa, string& strMiniReadsPathFq)
{
    if(vOrgReads.empty())
        return;

    ///1: Divide Unmapped Reads to small pices (3 pices per reads) and Create the relevant Fasta File
    ofstream ofs;
    strMiniReadsPathFa = strRoot + "Test/MiniReads.fa";
    strMiniReadsPathFq = strRoot + "Test/MiniReads.fq";
    ofs.open(strMiniReadsPathFa.c_str());

    int iMinLen = vOrgReads[0].strSeq.length() / 3;
    for(vector<St_Fasta>::iterator itr = vOrgReads.begin();
        itr != vOrgReads.end(); itr++)
    {
        ofs << ">" << itr->strName + ":1" << endl;
        ofs << itr->strSeq.substr(0, iMinLen) << endl;

        ofs << ">" << itr->strName + ":2" << endl;
        ofs << itr->strSeq.substr(iMinLen, iMinLen) << endl;

        ofs << ">" << itr->strName + ":3" << endl;
        ofs << itr->strSeq.substr(iMinLen*2, itr->strSeq.length() - iMinLen*2) << endl;
    }
    ofs.close();

    ///2: Transfer fa to fq -->
    string strCmd = "perl " + strRoot + "/ThirdPartyTools/fasta_to_fastq.pl " +
                    strMiniReadsPathFa + " > " + strMiniReadsPathFq;
    system(strCmd.c_str());
}

void ClsMethod::MapByMiniReads( string strRoot, string strRef, string strMiniReadsFq,
                                BamReader* pBamReader, map<int, string> mpType, int iTotalMinReadsNum,
                                int& iFirstTypeMap, int& iSecondTypeMap, int& iThirdTypeMap,
                                int& iMockTypeMap, vector<St_Fasta>& vUnmappedMiniReads,
                                vector<St_Fasta>& vMappedMiniReads)
{   
    vUnmappedMiniReads.clear();
    vMappedMiniReads.clear();

    ///3: Mapping those reads file back to reference
    string strBamFileFolder = "Output/BamFile/";
    string strBamFilePath = ClsBWA::GetInstance().CreateBamBySingleReads(strRef,
                                                                  strMiniReadsFq,
                                                                  "",
                                                                  strBamFileFolder,
                                                                  false, true);
    int i10UnmapPart = 0; int i1020UnmapPart = 0; int i2050UnmapPart = 0; int i50UnmapPart = 0;
    int iTotalMappedNumber = 0;
    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    St_Fasta stFasta;
    BamAlignment al;
    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        if(al.QueryBases == "")
            continue;

        if(!al.IsMapped())
        {
            stFasta.strName = al.Name;
            stFasta.strSeq = al.QueryBases;
            vUnmappedMiniReads.push_back(stFasta);
        }
        else //for the mapped case
        {
            if(al.RefID == -1)
                continue;

            stFasta.strName = al.Name;
            stFasta.strSeq = al.QueryBases;
            vMappedMiniReads.push_back(stFasta);

            //cout << al.Name << endl;
            if(mpType[al.RefID] == "first")
                iFirstTypeMap++;
            else if(mpType[al.RefID] == "second")
                iSecondTypeMap++;
            else if(mpType[al.RefID] == "third")
                iThirdTypeMap++;
            else if(mpType[al.RefID] == "mock")
                iMockTypeMap++;

            //Mapping Ratio
            float fUnmapRatio = (float)(al.Length - ::GetMatchPartLen(al)) / al.Length;
            if(fUnmapRatio < .1)
            {
                i10UnmapPart++;
            }
            else if(fUnmapRatio >= .1 && fUnmapRatio < .2)
            {
                i1020UnmapPart++;
            }
            else if(fUnmapRatio >= .2 && fUnmapRatio < .5)
            {
                i2050UnmapPart++;
            }
            else // this is >= 0.5
            {
                i50UnmapPart++;
            }

            iTotalMappedNumber++;
        }
    }

    //Print Unmap Part
    cout << endl;
    cout << "*******************Mini Reads******************" << endl;
    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
    cout << "Totoal Number of Mapping : " << IntToStr(iTotalMappedNumber) << endl;
    cout << "Mapping Percent ------>  : " << GetRatio((float)iTotalMappedNumber / iTotalMinReadsNum) << endl;
    cout << " < 10% Unmap Part        : " << IntToStr(i10UnmapPart) << " --> " << GetRatio((float)i10UnmapPart / iTotalMappedNumber) << endl;
    cout << " >=10% && <20% Unmap Part: " << IntToStr(i1020UnmapPart) << " --> " << GetRatio((float)i1020UnmapPart / iTotalMappedNumber) << endl;
    cout << " >=20% && <50% Unmap Part: " << IntToStr(i2050UnmapPart) << " --> " << GetRatio((float)i2050UnmapPart / iTotalMappedNumber) << endl;
    cout << " >=50% Unmap Part        : " << IntToStr(i50UnmapPart) << " --> " << GetRatio((float)i50UnmapPart / iTotalMappedNumber) << endl;
    cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;

    pBamReader->Rewind();
}

void ClsMethod::MappRemainingMiniReads(string strRoot, vector<St_Fasta>& vReads, BamReader* pBamReader)
{

    //1: Create fasta file for unmapped reads
    int iTotalReadsNum = vReads.size();
    string strReadsFa = strRoot + "Output/UnmappedReadsFasta/UnMappedMiniReads.fa";
    string strReadsFq = strRoot + "Output/UnmappedReadsFasta/UnMappedMiniReads.fq";
    ofstream ofs;
    ofs.open(strReadsFa.c_str());
    for(vector<St_Fasta>::iterator itr = vReads.begin();
        itr != vReads.end(); itr++)
    {
        ofs << ">" << itr->strName << endl;
        ofs << itr->strSeq << endl;
    }
    ofs.close();

    //2: translate fasta file to fastq file
    string strCmd = "perl " + strRoot + "/ThirdPartyTools/fasta_to_fastq.pl " +
                    strReadsFa + " > " + strReadsFq;
    system(strCmd.c_str());

    //3: Velvet --> Create Contigs
    strCmd = "rm " + strRoot + "/Output/LocalAssembly/*";
    system(strCmd.c_str());

    int iVelvetKmerLen = 21;
    string strContigFa = ClsVelvet::GetInstance().LocalAssembly(strReadsFa, 50, iVelvetKmerLen);

    //4: See how many remaining reads could be mapped back to those contigs
    string strBamFileFolder = "Output/BamFile/";
    string strBamFilePath = ClsBWA::GetInstance().CreateBamBySingleReads( strContigFa,
                                                                          strReadsFq,
                                                                          "",
                                                                          strBamFileFolder);
    pBamReader->Open(strBamFilePath);
    pBamReader->OpenIndex(strBamFilePath + ".bai");
    vector<St_Fasta>& vUnmappedReads = vReads;
    vUnmappedReads.clear();
    St_Fasta stFasta;
    BamAlignment al;
    while(pBamReader->GetNextAlignment(al)) //Get each alignment result
    {
        if(al.QueryBases == "")
            continue;

        if(!al.IsMapped())
        {
            stFasta.strName = al.Name;
            stFasta.strSeq = al.QueryBases;
            vUnmappedReads.push_back(stFasta);
        }
    }
    cout << "**********************" << endl;
    cout << "Mapping Percentage: " << GetRatio( 1- (float)vUnmappedReads.size() / iTotalReadsNum) << endl;
    cout << "**********************" << endl;
}

/////--->Go!!!
void ClsMethod::CodingByKmerAlignment(string strAnchorKmerPath, string strSeqFaPath, int iKmerLen)
{
    int iTotalNumKmer  = 0;
    int iNumContigFirstTime = 0;
    int iNumContigSecondTime = 0;
    int iNumContigThirdTime = 0;
    int iNumKmerForMock = 0;

    int iVelvetKmerLen = 21;
    string strRoot = GetHigherFolderPath(GetCurExeFolderPath());
    // Create Fasta Reader
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    // Create Bam Reader
    BamReader* pBamReader = new BamReader();
    // Use it to record all of valid contigs
    vector<St_Fasta> vContigFasta; // this is the container

    //******1: Use velvet create the relevant contigs (use velvet 3 times repeatly)******    
    cout << "*****************" << "I am the first time" << "*****************" << endl;
    string strUnAssembledKmerFa, strUnAssembledKmerFq;
    FirstRoundContigCreate(strRoot, pFastaReader, strAnchorKmerPath, pBamReader,
                           iVelvetKmerLen, vContigFasta, iTotalNumKmer, iNumContigFirstTime,
                           strUnAssembledKmerFa, strUnAssembledKmerFq);

    /// (2)--> The second time
    cout << "*****************" << "I am the second time" << "*****************" << endl;
    SecondRoundContigCreate(strRoot, pFastaReader, strUnAssembledKmerFa, pBamReader,
                            iVelvetKmerLen, vContigFasta, iNumContigSecondTime,
                            strUnAssembledKmerFa, strUnAssembledKmerFq);

    /// (3)--> The third time
    cout << "*****************" << "I am the third time" << "*****************" << endl;
    ThirdRoundContigCreate( strRoot, pFastaReader, strUnAssembledKmerFa, pBamReader,
                            iVelvetKmerLen, vContigFasta, iNumContigThirdTime,
                            strUnAssembledKmerFa, strUnAssembledKmerFq);

    /// (4): the remaing unmapping part should be concatenated to create the contig -->
    /*cout << "*****************" << "I am the remaining part" << "*****************" << endl;
    vector<St_Fasta> vUnAssembledKmer;
    pFastaReader->ReadFastaRegular(strUnAssembledKmerFa, vUnAssembledKmer);
    string strMockContig = "";
    for(vector<St_Fasta>::iterator itr = vUnAssembledKmer.begin();
        itr != vUnAssembledKmer.end(); itr++)
    {
        strMockContig += itr->strSeq;
    }
    iNumKmerForMock = vUnAssembledKmer.size();
    vContigFasta.push_back(St_Fasta("MockContig", strMockContig));
    vUnAssembledKmer.clear();*/

    //************ 2: Save those contig to file --> Create ***********************        
    //-->Erase the old file
    string strCmd = "rm " + strRoot + "Output/FinalContig/*";
    system(strCmd.c_str());
    //<--
    string strFinalContigFa = strRoot + "Output/FinalContig/FinalContig.fa"; // for velvet

    ofstream ofs;
    ofs.open(strFinalContigFa.c_str());
    map<int, string> mpType;
    int iIndex = 0;
    for(vector<St_Fasta>::iterator itr = vContigFasta.begin();
        itr != vContigFasta.end(); itr++)
    {
        //Name
        ofs << ">" << itr->strName << endl;
        //Sequence
        ofs << itr->strSeq << endl;

        if(itr->strName.find("first") != string::npos)
            mpType[iIndex] = "first";
        else if(itr->strName.find("second") != string::npos)
            mpType[iIndex] = "second";
        else if(itr->strName.find("third") != string::npos)
            mpType[iIndex] = "third";
        else if(itr->strName.find("Mock") != string::npos)
            mpType[iIndex] = "mock";

        iIndex++;
    }
    ofs.close();
    int iTotalNumContig = vContigFasta.size();

    //******3: Mapping Reads Back to contigs, to see how many reads could be mapped back successfully
    vector<St_Fasta> vHighQualityReads;
    pFastaReader->ReadFastaRegular(strSeqFaPath, vHighQualityReads);
    int iTotalMinReadsNum = vHighQualityReads.size() * 3; // Now we divid the reads into 3 parts

    ///(0) Get Mini Reads
    string strMiniReadsPathFa = ""; string strMiniReadsPathFq = "";
    CreateMiniReads(strRoot, vHighQualityReads, strMiniReadsPathFa, strMiniReadsPathFq);

    ///(1) 我们不这么去考虑，因为先将整体考虑去mapping (对现有的reads不进行进一步细分)，
    /// 这样即使map上了，我们同样存在大量的softclip使得我们无法进行有效的compression
    //1: SoftClipNum in Org Reads -->
    //int i10SoftClip = 0; int i1020SoftClip = 0; int i2050SoftClip = 0; int i50SoftClip = 0;
    //vector<St_Fasta> vUnmappedReads;
    //int iFirstTypeMap_OrgReads = 0; int iSecondTypeMap_OrgReads = 0;
    //int iThirdTypeMap_OrgReads = 0; int iMockTypeMap_OrgReads = 0;
    //MapByOrgReads(strRoot, strSeqFaPath, strFinalContigFa,
    //              pFastaReader, pBamReader,
    //              vUnmappedReads, mpType,
    //              iFirstTypeMap_OrgReads, iSecondTypeMap_OrgReads,
    //              iThirdTypeMap_OrgReads, iMockTypeMap_OrgReads);
    //int iRemainingTotal = vUnmappedReads.size() * 3;
    //int iHitOrgSuccucess = iTotalMinReadsNum - iRemainingTotal;

    ///(2)Use Mini Reads For Mapping
    int iFirstTypeMap = 0;
    int iSecondTypeMap = 0;
    int iThirdTypeMap = 0;
    int iMockTypeMap = 0;
    vector<St_Fasta> vUnmappedMiniReads;
    vector<St_Fasta> vMappedMiniReads;
    MapByMiniReads(strRoot, strFinalContigFa, strMiniReadsPathFq,
                   pBamReader, mpType, iTotalMinReadsNum,
                   iFirstTypeMap, iSecondTypeMap, iThirdTypeMap, iMockTypeMap,
                   vUnmappedMiniReads, vMappedMiniReads);

    //Save mapped mini Reads --> to see how large they are
    string strMappedMiniReads = strRoot + "Output/UnmappedReadsFasta/MappedMiniReads.fa"; // for velvet
    ofs.open(strMappedMiniReads.c_str());
    for(vector<St_Fasta>::iterator itr = vMappedMiniReads.begin();
        itr != vMappedMiniReads.end(); itr++)
    {
        ofs << ">" << itr->strSeq <<endl;
    }
    ofs.close();

    //Out put unmapped mini reads --> this is the whole remaining thing -->  The first time
    MappRemainingMiniReads(strRoot, vUnmappedMiniReads, pBamReader);

    //This is the second time mapping the remaining part
    MappRemainingMiniReads(strRoot, vUnmappedMiniReads, pBamReader);

    //Till this step: the un-mapped-mini-reads we get is the input reads that need to be compressed by
    //Huffman-coding
    //01/17/2016 -> 1: Save those unmapped reads to single file
    //(1) 正常的Fsata文件 --> 用于后面的压缩
    string strReadsFa = strRoot + "Output/UnmappedReadsFasta/UnMappedMiniReads.fa";
    ofs.open(strReadsFa.c_str());
    for(vector<St_Fasta>::iterator itr = vUnmappedMiniReads.begin();
        itr != vUnmappedMiniReads.end(); itr++)
    {
        ofs << ">" << itr->strName << endl;
        if(itr->strSeq.length() > 33)
            ofs << itr->strSeq.substr(0,33) << endl; // 在这里为了保证能够简化后面的运算，我们只取前面33个字符（对于mini reads而言）
        else
            ofs << itr->strSeq << endl;
    }
    ofs.close();
    //(2) 纯sequence文件，用于压缩size的比对
    string strPureReadsSeq = strRoot + "Output/UnmappedReadsFasta/PureUnMappedMiniReadsSeq.fa";
    ofs.open(strPureReadsSeq.c_str());
    for(vector<St_Fasta>::iterator itr = vUnmappedMiniReads.begin();
        itr != vUnmappedMiniReads.end(); itr++)
    {
        ofs << ">" << itr->strSeq << endl;
    }
    ofs.close();

    //***********************Print Result***********************************************
    cout << "Total Number of Kmer                          : " << IntToStr(iTotalNumKmer) << endl;
    cout << "The Contig We get by the frist time of velvet : " << IntToStr(iNumContigFirstTime) << endl;
    cout << "The Contig We get by the second time of velvet: " << IntToStr(iNumContigSecondTime) << endl;
    cout << "The Contig We get by the third time of velvet : " << IntToStr(iNumContigThirdTime) << endl;
    cout << "The remaining kmer to mock contig             : " << IntToStr(iNumKmerForMock) << endl;
    cout << endl;
    cout << "Total Number of Contig is                     : " << IntToStr(iTotalNumContig) << endl;
    cout << endl;

    int iUnmappedNum = vUnmappedMiniReads.size();
    int iMappedNum = iTotalMinReadsNum - iUnmappedNum; //iHitOrgSuccucess + (iRemainingTotal - iUnmappedNum);
    float fMapRatio = (float)iMappedNum / iTotalMinReadsNum;
    //float fRegularReadsMap = (float)iHitOrgSuccucess / iMappedNum;
    float fFirstTypeMap = (float)iFirstTypeMap/iMappedNum;
    float fSecondTypeMap = (float)iSecondTypeMap/iMappedNum;
    float fThirdTypeMap = (float)iThirdTypeMap/iMappedNum;
    float fMockTypeMap = (float)iMockTypeMap/iMappedNum;

    cout << "Total Mini Reads Num         : " << IntToStr(iTotalMinReadsNum) << endl;
    cout << "Total Mapped Mini Reads Num  : " << IntToStr(iMappedNum) << endl;
    cout << "Total Unmapped Mini Reads Num: " << IntToStr(iUnmappedNum) << endl;
    cout << "Mapped Percentage            : " << GetRatio(fMapRatio) << endl;
    cout << endl;

    /*
    cout << "Regular Reads Map            : " << IntToStr(iHitOrgSuccucess) << " --> " << GetRatio(fRegularReadsMap) << endl;
    //amplify 3 times for each mapping type by original reads
    cout << ">>>>>>>>>>>>>>>>>>>>iFirstTypeMap_OrgReads : " << IntToStr(iFirstTypeMap_OrgReads) <<
            " --> " << GetRatio((float)(iFirstTypeMap_OrgReads*3) / iHitOrgSuccucess) << endl;
    cout << ">>>>>>>>>>>>>>>>>>>>iSecondTypeMap_OrgReads: " << IntToStr(iSecondTypeMap_OrgReads) <<
            " --> " << GetRatio((float)(iSecondTypeMap_OrgReads*3) / iHitOrgSuccucess) << endl;
    cout << ">>>>>>>>>>>>>>>>>>>>iThirdTypeMap_OrgReads : " << IntToStr(iThirdTypeMap_OrgReads) <<
            " --> " << GetRatio((float)(iThirdTypeMap_OrgReads*3) / iHitOrgSuccucess) << endl;
    cout << ">>>>>>>>>>>>>>>>>>>>iMockTypeMap_OrgReads  : " << IntToStr(iMockTypeMap_OrgReads) <<
            " --> " << GetRatio((float)(iMockTypeMap_OrgReads*3) / iHitOrgSuccucess) << endl;*/

    cout << endl;
    cout << "First Type Map               : " << IntToStr(iFirstTypeMap) << " --> " << GetRatio(fFirstTypeMap) << endl;
    cout << "Second Type Map              : " << IntToStr(iSecondTypeMap) << " --> " << GetRatio(fSecondTypeMap) << endl;
    cout << "Third Type Map               : " << IntToStr(iThirdTypeMap) << " --> " << GetRatio(fThirdTypeMap) << endl;
    cout << "Mock Type Map                : " << IntToStr(iMockTypeMap) << " --> " << GetRatio(fMockTypeMap) << endl;

    //


    //Release Memory
    delete pBamReader;
    pBamReader = NULL;

    delete pFastaReader;
    pFastaReader = NULL;
}

/* We try to align those un-aligned mini read
 *
 * 1: Feature of those reaming unassembled & unaligned mini reads:
 * (1) Those reads do not share long overlap range.
 * (2) The length of those mini rads is 33 or 34
 *
 * 2: Strategy of Compression
 *
 * (1)8 bps fragment for Huffman Coding
 *    a) Use Kmer Counting to find out 12 Kmer (12 would bring us the best result)
 *    b) Choose the high frequency one
 *    c) Merge the kmers which just one character mismatch -> Merge the frequency number
 *    d) Use Huffman Tree to encode those kmers
 *    e) Code the corresponding 8-mer part in original mini reads by the old strategy
 *       * Find the most frequency one first and than try to find out the remaning part and so on
 *       * Need to compare with run length coding ->  Use RLC encoding the whole length directly if
 *         the compression size is less than huffman coding
 *       * with this comparison, all the reads which should use RLC for compression obviously will be
 *         picked out.
 *
 * (2) Use the remaning part run 8-mer huffman coding again
 *    * This is because the concatenated part may contain those 8 -mer as well.
 *
 * (3) 4 bps kmer fragment for huffman coding
 *    a) Divide reads into 4 parts and use the huffman coding
 *    b) Huffman Coding or RLC depends on which one could give us the small size --> Now we have done
 */
void ClsMethod::UnAlignedMiniReadsCompression(string strFasta, string strFastq, int iKmerLen)
{
    if(access(strFasta.c_str(), 0) != 0)
    {
        cout << "Fail to find Un Aligned Mini Reads" << endl;
        return;
    }

    string strRoot = GetHigherFolderPath(GetCurExeFolderPath());
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    BamReader* pBamReader = new BamReader();

    ///1: 12-mer fragment for Huffman Coding --->
    string strAnchorKMerPath = this->GetSolidKmer(strFasta, iKmerLen, //12 is the best one
                                            "MiniReadsAnchorKmer.fa");//Here, we still use average value for selection
    //(1)Select the high frequency 12 mer --> use velvelt to assemble those mini reads by k=13 不靠谱！！！
    //这说明我的想法是对的，对于短的kmer或者reads，即使很小的kmer能够存在overlap，但是基本上也是不能很好组装成理想的contig的
    /*string strCmd = "rm " + strRoot + "/Output/LocalAssembly/*";
    system(strCmd.c_str());
    int iVelvetKmerLen = 11;
    string strContigFa = ClsVelvet::GetInstance().LocalAssembly(strAnchorKMerPath, 50, iVelvetKmerLen);
    string strUnassembledKmerFa = GetUnMappedReads(strFasta, strContigFa, pBamReader, strRoot);*/
    vector<St_Fasta> vHighFreqKmer;
    pFastaReader->ReadFastaRegular(strAnchorKMerPath, vHighFreqKmer);
    vector< vector<St_Fasta> > mpSimilar;
    //1: Try to record & erase the similar one (just one different than others)
    int iLen = vHighFreqKmer[0].strSeq.length();
    int iOffSet = 0;
    vector<St_Fasta> vOneRoundSimilar;
    for(vector<St_Fasta>::iterator itr = vHighFreqKmer.begin(); itr != vHighFreqKmer.end(); itr++)
    {
        vOneRoundSimilar.clear();
        //for(int iIndex = vHighFreqKmer.size() -1; iIndex > iOffSet; iIndex--)
        for(vector<St_Fasta>::iterator subItr = vHighFreqKmer.end() - 1;
            subItr > itr; subItr--)

        {
            //Check if the difference is <= 1
            int iDiff = 0;
            bool bSimilar = true;
            for(int i=0; i<iLen; i++)
            {
                if(itr->strSeq[i] != subItr->strSeq[i])
                    iDiff++;
                if(iDiff > 1)
                {
                    bSimilar = false;
                    break;
                }
            }
            if(bSimilar)
            {
                vOneRoundSimilar.push_back(*subItr);
                vHighFreqKmer.erase(subItr);
            }
        }
        if(!vOneRoundSimilar.empty())
            mpSimilar.push_back(vOneRoundSimilar);
        iOffSet++;

        if(mpSimilar.size() > 6000)
            break;
    }

    cout << IntToStr(mpSimilar.size()) << endl;
    string strSimilarFile = strRoot + "/Output/LocalAssembly/Similar.txt";
    ofstream ofs;
    ofs.open(strSimilarFile.c_str());
    for(vector< vector<St_Fasta> >::iterator itr = mpSimilar.begin(); itr != mpSimilar.end(); itr++)
    {
        ofs << IntToStr(itr->size()) << endl;
    }
    ofs.close();


    delete pFastaReader;
    pFastaReader = NULL;
    delete pBamReader;
    pBamReader = NULL;
    cout << strAnchorKMerPath << endl;
}

/*
 * Scheme: --> Go!!!
 * 1: Load those unmapped mini-reads
 * 2: Collect the 12-mer of those reads
 *    (1) Those 12-mer is NOT the full set
 *    (2) For example: 33 mini reads: 1-12, 12-24, 24-33(this could be expressed by
 *            any 12-mer share the same first 9 characters with this remaining part)
 *    (3) Record the frequency while collecting those kmers
 * 3: Get the super k-mer
 *    (1) 1 12-mer = 4*12 = 48 similar mers -> the max number of items in one group is 48
 *    (2) Sort the kmer by frequence --> number of super k-mer = the number of unique kmer / 48
 * 4: Divde the kmer in two big group:
 *    groupA: Super Kmer Candidate   GroupB: Low frequency kmer
 * 5: Build the super-kmer clusters
 *    (1) Travrsal "groupA" and collect its similar kmer in both groupA and groupB
 *        *Do the traversal until the groupA is empty
 *    (2) See the remaning part in groupB as super-kmer
 * 6: Make huffman tree
 *    (1) Build Huffman tree by super-kmer
 *    (2) Build Huffman tree by group of k-mers in superkmer
 * 7: Encode those unmapped mini-reads by the huffman coding obtained by step 6
 *
 * Question:
 * 1: Too slow:
 *   (1) Find Super Kmer and its similar group
 *
 * 2: Coding
 *   (1) Search which part from Super Kmer and which part from which similar group
 *
 * 3: The feature similar do not be used
 *
 * 1: In fact: here we could make sure the huffman coding is better than normal ways, however,
 *    we do not know the space consuming of huffman tree
 */
bool sortfunction_kmer(St_KmerInfo stKmer1, St_KmerInfo stKmer2) // from large to small (frequency)
{
    return stKmer1.iNum > stKmer2.iNum;
}

void ClsMethod::UnAlignedMiniReadsComprByHuffmanCoding(string strFasta, int iKmerLen)
{
    //Step 1: Load those unmapped mini-reads (those reads are unmapped )
    ClsFastaReader* pFastaReader = new ClsFastaReader();
    vector<St_Fasta> vReads;
    pFastaReader->ReadFastaRegular(strFasta, vReads);
    //int iKeptSize = vReads.size() * .5;
    //vReads.erase(vReads.begin() + iKeptSize, vReads.end());

    //Save those pure sequence to a single file
    ofstream ofs;
    string strRootPath =  GetHigherFolderPath(GetCurExeFolderPath());
    string strFilePath = strRootPath + "Test/CompressedFile/PureSeqFile.txt";
    ofs.open(strFilePath.c_str());
    for(vector<St_Fasta>::iterator itr = vReads.begin(); itr != vReads.end(); itr++)
    {
        ofs << itr->strSeq << endl;
    }
    ofs.close();

    delete pFastaReader;
    pFastaReader = NULL;

    //Step 2: Collect 12-mer of those reads
    vector<string> vGens; // 用于顺序记录reads的组成，因为reads是等长的，
                          // 所以在解析的时候可以将他们先组合到一起，然后再进行等长的分割就可以了
    map<string, int> mpKmer; //"string" is the kmer seq, "int" is frequency
    for(vector<St_Fasta>::iterator itr = vReads.begin(); itr != vReads.end(); itr++)
    {
        string& strSeq = itr->strSeq;
        int iPos = 0;
        while(iPos < strSeq.length())
        {
            string strCurKmer = "";
            if(iPos + iKmerLen < strSeq.length())
                strCurKmer = strSeq.substr(iPos, iKmerLen);
            else
                strCurKmer = strSeq.substr(iPos, strSeq.length()-iPos);
            vGens.push_back(strCurKmer);
            mpKmer[strCurKmer]++;
            iPos += iKmerLen;
        }
    }

    //--> 看看在这个阶段是不是都能够匹配上
    int iii = 0;
    int zeroItem = 0;
    for(vector<string>::iterator itr = vGens.begin(); itr != vGens.end(); itr++)
    {
        if(mpKmer.find(*itr) == mpKmer.end())
            iii++;
        else
        {
            if(mpKmer.find(*itr)->second == 0)
                zeroItem++;
        }
    }


    //Chose Half Frequency -->small value of sequence could bring us the small size of file
    for(map<string, int>::iterator itr = mpKmer.begin(); itr != mpKmer.end(); itr++)
    {
        itr->second = ceil((float)itr->second / 2);
    }

    //Step 3: Get super-k-mer
    //a)Move Map to Vector
    vector<St_KmerInfo> vKmer;
    vKmer.resize(mpKmer.size());
    map<string, int>::iterator itrMpKmer = mpKmer.begin();
    for(vector<St_KmerInfo>::iterator itr = vKmer.begin(); itr != vKmer.end(); itr++)
    {
        if(itrMpKmer->first == "ACCTTCATTTA")
        {
            int i = 0;
            i++;
        }
        itr->strSeq = itrMpKmer->first;
        itr->iNum = itrMpKmer->second;
        itrMpKmer++;
    }

    //b)Sort Vector by Frequency
    sort(vKmer.begin(), vKmer.end(), sortfunction_kmer);

    vector<St_SuperKmer> vSuperKmer;
    UpdateSuperKmer(vKmer, vSuperKmer, vGens);

    /*
    //c)The number of super kmer is: the top totalnum/(4*iKmerLen) high frequency kmer
    int iNumSuperKmer = vKmer.size() / (4*iKmerLen);
    //d)Divide k-mer into two groups(candidate, and regular)
    vector<St_KmerInfo> vKmerCandidate(vKmer.begin(), vKmer.begin() + iNumSuperKmer);
    vector<St_KmerInfo> vKmerRegular(vKmer.begin() + iNumSuperKmer, vKmer.end());
    vKmer.clear();
    //e)Collect the similar group for those kmers    

    St_SuperKmer stSuperKmer;
    //  i)优先遍历candidate中的
    for(vector<St_KmerInfo>::iterator itrCandi = vKmerCandidate.begin();
        itrCandi != vKmerCandidate.end(); itrCandi++)
    {        
        //init new superKmer -->
        bool bFind = false;
        for(vector<St_SuperKmer>::iterator itr = vSuperKmer.begin(); itr != vSuperKmer.end(); itr++)
        {
            //KmerTypeShort vValue = 0;
            //FormKmerTypeShortSeg(itrCandi->strSeq.c_str(), 0, itrCandi->strSeq.length(), vValue);
            if(itr->stInfo.strSeq == "ACCTTCATTTA")
            {
                int i = 0;
                i++;
            }

            if(itr->mpOrgSimilar.find(itrCandi->strSeq) !=  itr->mpOrgSimilar.end()) //find it
            {
                itr->mpOrgSimilar[itrCandi->strSeq].iCount += itrCandi->iNum;
                bFind = true;
                break;
            }
        }
        if(!bFind)
        {
            if(itrCandi->strSeq == "ACCTTCATTTA")
            {
                int i = 0;
                i++;
            }

            stSuperKmer.stInfo = *itrCandi;
            stSuperKmer.vSimilar.clear();
            stSuperKmer.BuildMapOrg(vSuperKmer.size());
            vSuperKmer.push_back(stSuperKmer);
        }
    }

    //这里读的速度太慢，我们从程序上采取这个方法来做: 所有备选的similar item 串到一个map中，这样我们就可以一次遍历所有的item了
    map<string, St_KmerIndex> mpCombine;
    for(vector<St_SuperKmer>::iterator itr = vSuperKmer.begin(); itr != vSuperKmer.end(); itr++)
    {
        if(itr->mpOrgSimilar.find("ACCTTCATTTA") != itr->mpOrgSimilar.end())
        {
            int i=0;
            i++;
        }
        mpCombine.insert(itr->mpOrgSimilar.begin(), itr->mpOrgSimilar.end());
    }

    ClsMultiThread *pMT = new ClsMultiThread();
    pMT->FindSimilarKmer(mpCombine, vKmerRegular);
    delete pMT;
    pMT = NULL;

    //将mpCombine 中的更新的值跟主的vSuperKmer进行相应的同步
    for(vector<St_SuperKmer>::iterator itr = vSuperKmer.begin(); itr != vSuperKmer.end(); itr++)
    {
        //对于里面的每一个元素
        for(map<string, St_KmerIndex>::iterator subItr = itr->mpOrgSimilar.begin();
            subItr != itr->mpOrgSimilar.end(); subItr++)
        {
            subItr->second.iCount = mpCombine[subItr->first].iCount;
            if(subItr->first == "ACCTTCATTTA")
            {
                int rere = 0;
                rere++;
            }

            if(subItr->second.iCount == 0)
            {
                map<string, St_KmerIndex>::iterator itrDelete = subItr;
                subItr++;
                itr->mpOrgSimilar.erase(itrDelete);
                subItr--;
            }
        }
    }

    //-------->check if there are some items with abnormal 0 value
    for(vector<St_SuperKmer>::iterator itr = vSuperKmer.begin(); itr != vSuperKmer.end(); itr++)
    {
        if(itr->GetSizeOrg() == 0)
        {
            int i = 0;
            i++;
            break;
        }
    }

    //  j) 将Regular中剩下的也看作是SuperKmer加入到容器中
    map<string, St_KmerIndex> mpSuperKmerSilimar;
    vector<St_SuperKmer> vSuperRegKmer;
    int i = 0;
    for(vector<St_KmerInfo>::iterator itrReg = vKmerRegular.begin();
        itrReg != vKmerRegular.end(); itrReg++)
    {
        if(itrReg->strSeq == "ACCTTCATTTA")
        {
            int qwe = 0;
            qwe++;
        }

        if(mpSuperKmerSilimar.find(itrReg->strSeq) != mpSuperKmerSilimar.end()) //如果遍历到
        {            
            mpSuperKmerSilimar[itrReg->strSeq].iCount += itrReg->iNum;
        }
        else // 如果遍历不到
        {
            if(itrReg->strSeq == "ACCTTCATTTA")
            {
                int qwe = 0;
                qwe++;
            }

            St_SuperKmer stSuperKmer;
            stSuperKmer.stInfo = *itrReg;
            vSuperRegKmer.push_back(stSuperKmer);
            //---> InSert 相应的 变量到map中去
            St_KmerIndex stKmerIndex;
            stKmerIndex.iGroupIndex = vSuperRegKmer.size()-1;
            stKmerIndex.iCount = 0;
            for(int iPos = 0; iPos < itrReg->strSeq.length(); iPos++)
            {
                string strPos = itrReg->strSeq.substr(iPos, 1);

                //For "A"
                mpSuperKmerSilimar[itrReg->strSeq.replace(iPos, 1, "A")] = stKmerIndex;
                //For "T"
                mpSuperKmerSilimar[itrReg->strSeq.replace(iPos, 1, "T")] = stKmerIndex;
                //For "G"
                mpSuperKmerSilimar[itrReg->strSeq.replace(iPos, 1, "G")] = stKmerIndex;
                //For "C"
                mpSuperKmerSilimar[itrReg->strSeq.replace(iPos, 1, "C")] = stKmerIndex;

                //Recover the original value
                itrReg->strSeq.replace(iPos, 1, strPos.c_str());
            }
            mpSuperKmerSilimar[itrReg->strSeq].iCount = itrReg->iNum;
        }

        if(mpSuperKmerSilimar.find("ACCTTCATTTA") != mpSuperKmerSilimar.end())
        {
            if(mpSuperKmerSilimar["ACCTTCATTTA"].iCount == 0)
            {
                int dsds = 0;
                dsds++;
            }
        }
        i++;
    }
    //更新新的容器中的记录的map的值
    for(map<string, St_KmerIndex>::iterator itr = mpSuperKmerSilimar.begin();
        itr != mpSuperKmerSilimar.end(); itr++)
    {        
        if(itr->first == "ACCTTCATTTA")
        {
            int i = 0;
            i++;
        }

        if(itr->second.iCount == 0)
            continue;
        else
            vSuperRegKmer[itr->second.iGroupIndex].mpOrgSimilar.insert(pair<string,St_KmerIndex>(itr->first, itr->second));
    }

    for(vector<St_SuperKmer>::iterator itr = vSuperRegKmer.begin(); itr != vSuperRegKmer.end(); itr++)
    {
        if(itr->GetSizeOrg() == 0)
        {
            int i = 0;
            i++;
            break;
        }
    }

    vSuperKmer.insert(vSuperKmer.end(), vSuperRegKmer.begin(), vSuperRegKmer.end());*/
    //-->更新super kmer节点的count值
    //for(vector<St_SuperKmer>::iterator itr = vSuperKmer.begin(); itr < vSuperKmer.end(); itr++)
    //{
    //    itr->stInfo.iNum = itr->GetSizeOrg();
    //}


    cout << "*********************Huffman Code Info***********************" << endl;
    cout << "Number of SuperKmer is: " << IntToStr(vSuperKmer.size()) << endl;
    cout << "The first Similar Group size is: " << IntToStr(vSuperKmer[0].GetSizeOrg()) << endl;
    cout << "The second Similar Group size is: " << IntToStr(vSuperKmer[1].GetSizeOrg()) << endl;
    cout << "The third Similar Group size is: " << IntToStr(vSuperKmer[2].GetSizeOrg()) << endl;
    cout << "The third Similar Group size is: " << IntToStr(vSuperKmer[3].GetSizeOrg()) << endl;
    cout << "The third Similar Group size is: " << IntToStr(vSuperKmer[4].GetSizeOrg()) << endl;
    cout << "The third Similar Group size is: " << IntToStr(vSuperKmer[5].GetSizeOrg()) << endl;
    cout << "*************************************************************" << endl;

    //Step4: Build Huffman Tree for SuperKmer --> Double Level --> Go!! -->终于到了建huffman tree的时候了
    /******scheme
    /* 1： 我们首先建最外层的
     * 2： 然后针对每一个单独的superkmer，我们都进行相应的huffman树的建立
     * 3： 合并主huffman树和子huffman树，将结果作为最后节点的huffman code
     * 4： 通过map进行相应的查找，然后路分段编码，直至最后获得最终的压缩结果。
     */

    //1：  我们首先建最外层的
    map<string, int> mpKmerFreq;
    for(vector<St_SuperKmer>::iterator itr = vSuperKmer.begin();
        itr != vSuperKmer.end(); itr++)
    {
        mpKmerFreq[itr->stInfo.strSeq] = itr->GetSizeOrg();
    }

    Hufftree<string, int>* hufftree = NULL;
    hufftree = new Hufftree<string, int>(mpKmerFreq.begin(), mpKmerFreq.end());  //这里出错，醒来看看是怎么回事儿

    //Compress the sequence by huffman coding Huffman Coding
    map<string, vector<bool> > mpKmerCode;
    for(vector<St_SuperKmer>::iterator itr = vSuperKmer.begin();
        itr != vSuperKmer.end(); itr++)
    {
        //1: Get super kmer code
        vector<bool> vSuperKmerCode = hufftree->encode(itr->stInfo.strSeq);

        //2: Get current kmer code
        map<string, int> mpSubKmerFreq;
        for(map<string, St_KmerIndex>::iterator subItr = itr->mpOrgSimilar.begin();
            subItr != itr->mpOrgSimilar.end(); subItr++)
        {
            mpSubKmerFreq[subItr->first] = subItr->second.iCount;
        }
        Hufftree<string, int>* subHufftree = NULL;
        subHufftree = new Hufftree<string, int>(mpSubKmerFreq.begin(), mpSubKmerFreq.end());
        for(map<string, St_KmerIndex>::iterator subItr = itr->mpOrgSimilar.begin();
            subItr != itr->mpOrgSimilar.end(); subItr++)
        {
            vector<bool> vSubKmerCode = subHufftree->encode(subItr->first);
            vector<bool> vSum;
            vSum.insert(vSum.end(), vSuperKmerCode.begin(), vSuperKmerCode.end());
            vSum.insert(vSum.end(), vSubKmerCode.begin(), vSubKmerCode.end());
            mpKmerCode[subItr->first] = vSum;
        }
        delete subHufftree;
        subHufftree = NULL;
    }

    //取得相应的Kmer code map后，我们根据Gen中的东西，对相应的值统一存储
    vector<bool> vSum;
    int i = 0;
    for(vector<string>::iterator itr = vGens.begin(); itr != vGens.end(); itr++)
    {
        if(mpKmerCode.find(*itr) == mpKmerCode.end()) //这个证明找不到-->the case is unexpected!!!
        {
            i++;
            continue;
        }
        vector<bool>& vCode = mpKmerCode[*itr];
        vSum.insert(vSum.end(), vCode.begin(), vCode.end());
    }

    //(4) Save File
    string strDestFile = strRootPath + "Test/CompressedFile/OrgFilCompress.bit";
    bit_file_c destBf;
    destBf.Open( strDestFile.c_str(), BF_WRITE );
    //writeBitVINT( destBf, vEncodedLongGens.size() ); // Here we do not need to write bits to it
    writeBits(destBf, &vSum );
    destBf.ByteAlign();

    destBf.Close();

    //(5) Release memory
    vGens.clear();

    delete hufftree;
    hufftree = NULL;
}

struct St_KmerLabel
{
    int iGroupIndex;
    int iFreq;
    St_KmerLabel(): iGroupIndex(-1), iFreq(0)
    {}
};

void ClsMethod::UpdateSuperKmer(vector<St_KmerInfo>& vKmer,
                                vector<St_SuperKmer>& vSuperKmer,
                                vector<string>& vGens)
{
    map<string, int> mpTest;
    for(vector<St_KmerInfo>::iterator itr = vKmer.begin(); itr != vKmer.end(); itr++)
        mpTest[itr->strSeq] = itr->iNum;

    int i = 0;
    int iZero = 0;
    for(vector<string>::iterator itr = vGens.begin(); itr != vGens.end(); itr++)
    {
        if(mpTest.find(*itr) == mpTest.end())
            i++;
        else
        {
            if(mpTest.find(*itr)->second == 0)
                iZero++;
        }
    } // 都存在！！！

    vSuperKmer.clear();
    map<string, St_KmerLabel> mpLabel;
    St_SuperKmer stSuperKmer;
    //使用之前的方法，--> add label in to the map and use the binary search by API "find"
    for(vector<St_KmerInfo>::iterator itr = vKmer.begin(); itr != vKmer.end(); itr++)
    {
        if(itr->strSeq == "ACCTTCATTTA")
        {
            int fdfdf = 0;
            fdfdf++;
        }

        map<string, St_KmerLabel>::iterator mpItr = mpLabel.find(itr->strSeq);
        if(mpItr != mpLabel.end())
        {
            mpItr->second.iFreq += itr->iNum;
            int fdfdf = 0;
            fdfdf++;
        }
        else
        {
            //Add the relevant label into Map
            St_KmerLabel stKmerLabel;
            stKmerLabel.iGroupIndex = vSuperKmer.size(); // 因为是要相加的，所以应该是直接用size进行定位
            stSuperKmer.stInfo = *itr;
            vSuperKmer.push_back(stSuperKmer);
            //current seq
            stKmerLabel.iFreq = itr->iNum;
            mpLabel[itr->strSeq] = stKmerLabel;
            //similar seq
            stKmerLabel.iFreq = 0;
            for(unsigned int iPos = 0; iPos < itr->strSeq.length(); iPos++)
            {
                string strPos = itr->strSeq.substr(iPos, 1);
                string strOrg = itr->strSeq;
                //For "A"
                string strNewLabel = "";
                if(strPos != "A")
                {
                    strNewLabel = itr->strSeq.replace(iPos, 1, "A");
                    if(mpLabel.find(strNewLabel) == mpLabel.end())
                        mpLabel[strNewLabel] = stKmerLabel;

                    if(strNewLabel == "ACCTTCATTTA")
                    {
                        int fdfdf = 0;
                        fdfdf++;
                    }

                }
                //For "T"
                if(strPos != "T")
                {
                    strNewLabel = itr->strSeq.replace(iPos, 1, "T");
                    if(mpLabel.find(strNewLabel) == mpLabel.end())
                        mpLabel[strNewLabel] = stKmerLabel;

                    if(strNewLabel == "ACCTTCATTTA")
                    {
                        int fdfdf = 0;
                        fdfdf++;
                    }
                }
                //For "G"
                if(strPos != "G")
                {
                    strNewLabel = itr->strSeq.replace(iPos, 1, "G");
                    if(mpLabel.find(strNewLabel) == mpLabel.end())
                        mpLabel[strNewLabel] = stKmerLabel;

                    if(strNewLabel == "ACCTTCATTTA")
                    {
                        int fdfdf = 0;
                        fdfdf++;
                    }
                }
                //For "C"
                if(strPos != "C")
                {
                    strNewLabel = itr->strSeq.replace(iPos, 1, "C");
                    if(mpLabel.find(strNewLabel) == mpLabel.end())
                        mpLabel[strNewLabel] = stKmerLabel;

                    if(strNewLabel == "ACCTTCATTTA")
                    {
                        int fdfdf = 0;
                        fdfdf++;
                    }
                }

                //Recover the original value
                itr->strSeq.replace(iPos, 1, strPos.c_str());
            }
        }
    }

    i = 0;
    iZero = 0;
    for(vector<string>::iterator itr = vGens.begin(); itr != vGens.end(); itr++)
    {
        if(*itr == "ACCTTCATTTA")
        {
            int fdfdf = 0;
            fdfdf++;
        }

        map<string, St_KmerLabel> ::iterator itrLabel = mpLabel.find(*itr);
        if(itrLabel == mpLabel.end())
            i++;
        else
        {
            if(itrLabel->second.iFreq == 0)
                iZero++;
        }
    }

    //Set the additional value of vSuperKmer from mpLabel;
    mpTest.clear();
    for(map<string, St_KmerLabel>::iterator itr = mpLabel.begin(); itr != mpLabel.end(); itr++)
    {
        if(itr->second.iFreq == 0)
            continue;
        St_SuperKmer& stCurSuperKmer = vSuperKmer[itr->second.iGroupIndex];
        St_KmerIndex stKmerIndex;
        stKmerIndex.iCount = itr->second.iFreq;
        stKmerIndex.iGroupIndex = itr->second.iGroupIndex;
        stCurSuperKmer.mpOrgSimilar[itr->first] = stKmerIndex;
        mpTest[itr->first] = itr->second.iFreq;
    }

    //看看 vGens -->是不是全部能够存在于 map org similar 中 ！！！  --->Go!!!
    int iUnHit = 0;
    for(vector<string>::iterator itr = vGens.begin(); itr != vGens.end(); itr++)
    {
        if(mpTest.find(*itr) == mpTest.end())
            iUnHit++;
    }
    cout << IntToStr(iUnHit) << endl; // Do not Correct !!!! Check Why!!!
}
