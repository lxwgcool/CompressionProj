//
//  KmerUtils.h
//  
//
//  Created by Yufeng Wu on 8/31/14.
//
//  Process kmer related stuff
//

#ifndef ____KmerUtils__
#define ____KmerUtils__

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include "stdint.h"
using namespace std;


//////////////////////////////////////////////////////////////////////
// define kmer that is short than 32 nts

typedef char ReadNtType;
typedef uint64_t KmerTypeShort;

// define the map for k-mer frequency
typedef map<KmerTypeShort, double> MapShortKmerFreq;

//////////////////////////////////////////////////////////////////////
// init kmer by a pointer to an array, and a position
void FormKmerTypeShortSeg( const ReadNtType *pArrayIn, int posKmer, int kmerLen,KmerTypeShort &kmer );
// or create by another kmer and then left-shift by one nt
void FormKmerTypeShortShift( const KmerTypeShort &kmerTypePrev, int kmerLen, ReadNtType ntRightMost, KmerTypeShort &kmerNew);
// get all kmers from a short sequence
void GetAllKmersFromSeq( const ReadNtType *pSeq, int seqLen, int kmerLen, vector<KmerTypeShort> &listKmers);
// dump kmer out
void DumpKmer( const KmerTypeShort &kmer, int kmerLen );
// convert integer-based kmer to string representation
void ConvKmerToString( const KmerTypeShort &kmer, int szKmer, char *strKmer );
// read a set of kmer from file and store in a hash table; return kmer length
int ReadInKmerToHashMap( const char *fileName, MapShortKmerFreq &mapKmerFreq );
// test whether a new seq (i.e. read) contain at least certain number
// of frequent k-mer as in the map
bool IsReadContainingFreqKmers(  const ReadNtType *pSeq, int seqLen, int kmerLen, int minOccurKmerThres, const MapShortKmerFreq &mapKmerFreqIn);
// write out a read
void OutputReadFastq(ostream &outputStream, const ReadNtType *pSeq, int seqLen, const char *idReadDesc, const char *qualityRead);

#endif /* defined(____KmerUtils__) */
