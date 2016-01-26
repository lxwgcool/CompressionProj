#include "clsmultithread.h"
#include "clsbasealgorithm.h"
#include "pthread.h"
#include "ctime"
#include "time.h"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <map>
#include <stdlib.h>
#include <tr1/memory>
#include <cstdio>

//we do multi thread here !! -->
std::string exec(const char* cmd)
{
    std::tr1::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while (!feof(pipe.get())) {
        memset(buffer, 0, 128);
        if (fgets(buffer, 128, pipe.get()) != NULL)
        {
            result += buffer;
        }
    }
    return result;
}

void* ClsThreadFunc::SimilarDetect(void* args)
{
    St_SimilarKmer* pStSimilarKmer = (St_SimilarKmer*) args;
    // Search in Regular Container
    int iTemp = 0;
    for(vector<St_KmerInfo>::iterator itrReg = pStSimilarKmer->vRegKmer.begin();
        itrReg != pStSimilarKmer->vRegKmer.end(); itrReg++)
    {
        bool bFind = false;
        if(pStSimilarKmer->vSuperKmer.find(itrReg->strSeq) != pStSimilarKmer->vSuperKmer.end()) //find it
        {
            pStSimilarKmer->vSuperKmer[itrReg->strSeq]++;
            bFind = true;
        }
        if(bFind)
        {
            pStSimilarKmer->vRegKmer.erase(itrReg);
            itrReg--;
        }
        iTemp++;
        /*
        for(vector<St_SuperKmer>::iterator itr = vSuperKmer.begin(); itr != vSuperKmer.end(); itr++)
        {
            //KmerTypeShort vValue = 0;
            //FormKmerTypeShortSeg(itrReg->strSeq.c_str(), 0, itrReg->strSeq.length(), vValue);
            if(itr->mpOrgSimilar.find(itrReg->strSeq) !=  itr->mpOrgSimilar.end()) //find it
            {
                itr->mpOrgSimilar[itrReg->strSeq]++;
                bFind = true;
                break;
            }
        }
        if(bFind)
        {
            vKmerRegular.erase(itrReg);
            itrReg--;
        }*/
        //iTemp++;
    }
}

ClsMultiThread::ClsMultiThread()
{
}

void ClsMultiThread::FindSimilarKmer(map<string, int>& mpCombine, vector<St_KmerInfo>& vKmerRegular)
{
    //1: Get the number of CPU
    string strNum = exec("nproc");
    if(strNum.substr(strNum.length()-1,1) == "\n")
        strNum = strNum.substr(0, strNum.length()-1);
    int iCpuNum = atoi(strNum.c_str());

    //2: Divid KmerRegular Table evenly by the number of cpu
    int iLen = vKmerRegular.size() / iCpuNum;
    vector<St_SimilarKmer> vSubGroupKmer;
    vSubGroupKmer.resize(iCpuNum);
    for(int i=0; i<iCpuNum; i++)
    {
        if((i+1)*iLen < vKmerRegular.size())
            vSubGroupKmer[i].vRegKmer.insert(vSubGroupKmer[i].vRegKmer.begin(),
                                             vKmerRegular.begin() + i*iLen, vKmerRegular.begin() + (i+1)*iLen);
        else
            vSubGroupKmer[i].vRegKmer.insert(vSubGroupKmer[i].vRegKmer.begin(),
                                             vKmerRegular.begin() + i*iLen,
                                             vKmerRegular.begin() + (vKmerRegular.size() - i*iLen));
        vSubGroupKmer[i].vSuperKmer.insert(mpCombine.begin(), mpCombine.end());
    }

    //Create Threads
    int start_s = time(NULL);
    pthread_t tids[iCpuNum];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for(int i=0; i<iCpuNum; i++)
    {
        int ret = pthread_create(&tids[i], &attr, ClsThreadFunc::SimilarDetect, (void*)&(vSubGroupKmer[i]));
        if(ret != 0)
            cout << "Thread error: " << ret << endl;
    }
    pthread_attr_destroy(&attr);

    void* status;
    for(int i=0; i<iCpuNum; i++)
    {
        int ret = pthread_join(tids[i], &status);
        if(ret != 0)
            cout << "Thread Error: " << ret << endl;
        else
            cout << "Thread Status: " << (long)status << endl;
    }
    pthread_exit(NULL);
    int stop_s = time(NULL);
    double dRuningTime = difftime(stop_s, start_s);
    string strTimeFormat = ::GetHMSTimeFormat(dRuningTime);
    cout << "Running Time: " << strTimeFormat << endl;

}
