#include "clsmethod.h"
#include <iostream>

int main(int argc, char * argv[]) //1: Kmer length, 2: Do HuffmanCoding
                                  //-> We find Kmer = 38 could give us the best result
{    
    string strTemp = "AAAAAAAA";
    string dsds = strTemp.replace(0,1,"B");
    int i = 0;
    /*
    string strFqPath = "/home/lq/lxwg/Software/GAGE_DATA_Scaf_Gap_Filling_Testing/Data/Stapylococcus_aureus/Data/original/frag_1.fastq";
    ClsMethod* pMethod = new ClsMethod();

    string strHighQualityFastqPath = "/home/lq/lxwg/WorkStudio/Prototype/Scaf_Gap_Filling_Tools/Huffman_Coding/Test/HighQuality.fq";
    string strHighQualityFastaPath = "/home/lq/lxwg/WorkStudio/Prototype/Scaf_Gap_Filling_Tools/Huffman_Coding/Test/HighQuality.fa";

    //--> Encode the fastq file
    pMethod->ReadFastqFile(strFqPath, pMethod->m_vFastq);
    pMethod->FilterNItems();
    pMethod->SaveHighQualityFastqToFile(strHighQualityFastqPath,
                                        strHighQualityFastaPath);
    pMethod->EncodeFastqFile(strHighQualityFastaPath);
    //<--*/

    ClsMethod* pMethod = new ClsMethod();

    /*/The method based on Kmer Coding + regular Coding ******
    string strSeqFaPath = "/home/lq/lxwg/WorkStudio/Prototype/Scaf_Gap_Filling_Tools/MinEntropy/Test/OrgSeq.fa";//genome.fasta";
    //1: Get Solid Kmer -->
    int iKmerLen = 18;
    string strAnchorKMerPath = pMethod->GetSolidKmer(strSeqFaPath, iKmerLen);

    //2: 将这些Kmer建成Huffman树 --> 去建吧
    pMethod->HuffmanEncoding(strAnchorKMerPath, strSeqFaPath, iKmerLen);

    //3: Decode the Compressed File
    //pMethod->HuffmanDecoding();
    */

    //The new method coding those two parts
    //1: Transfer fq file to fa --> Prepare the original "fasta" and "fastq" file

    string strRoot = GetHigherFolderPath(GetCurExeFolderPath());

    /*
    string strFqPath = strRoot + "Test/frag_1.fastq";
    pMethod->ReadFastqFile(strFqPath, pMethod->m_vFastq);
    pMethod->FilterNItems();
    string strHighQualityFastqPath = strRoot + "Test/HighQuality.fq";
    string strHighQualityFastaPath = strRoot + "Test/HighQuality.fa";
    pMethod->SaveHighQualityFastqToFile(strHighQualityFastqPath, strHighQualityFastaPath);

    //2: Get the soild kmer
    int iKmerLen = atoi(argv[1]); //40;
    int iCoverageThreshold = 10; //在这里我们还是使用之前的 average value去进行solid kmer的筛选
    string strAnchorKMerPath = pMethod->GetSolidKmer(strHighQualityFastaPath, iKmerLen);
    cout << strAnchorKMerPath << endl;

    //3: Assemble those kmers. and coding by the method of reference based -->
    if(atoi(argv[2]) == 0) // 0: do not do coding, 1: do coding
        return 0;
    pMethod->CodingByKmerAlignment(strAnchorKMerPath, strHighQualityFastaPath, iKmerLen); */

    //--->Coding the remaining part by Huffman Code
    string strUnAlignedMiniReadsFa = strRoot + "Output/UnmappedReadsFasta/UnMappedMiniReads.fa";
    //string strUnAlignedMiniReadsFq = strRoot + "Output/UnmappedReadsFasta/UnMappedMiniReads.fq";
    //int iKmerLen = atoi(argv[1]); //40;
    //pMethod->UnAlignedMiniReadsCompression(strUnAlignedMiniReadsFa, strUnAlignedMiniReadsFq, iKmerLen);
    pMethod->UnAlignedMiniReadsComprByHuffmanCoding(strUnAlignedMiniReadsFa, 11);
    //<---

    delete pMethod;
    pMethod = NULL;
    return 0;
}

