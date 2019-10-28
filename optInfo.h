#ifndef __OPT_INFO_H_
#define __OPT_INFO_H_

#include<iomanip>
#include<vector>

class OptInfo
{
    public:
        OptInfo(int it, double ps, double ds, double pd, double dd): iter(it), primalS(ps), dualS(ds), primalD(pd), dualD(dd), jDistN(0), jDistD(0) {}
        OptInfo(int it, double ps, double ds, double pd, double dd, int jdn, int jdd): iter(it), primalS(ps), dualS(ds), primalD(pd), dualD(dd), jDistN(jdn), jDistD(jdd){}
        int iter;
        double primalS;
        double dualS;
        double primalD;
        double dualD;
        size_t nSFM;
        int    jDistN;
        int    jDistD;
};

inline void print(std::vector<OptInfo> & stats)
{
    int count = 0;
    for(std::vector<OptInfo>::iterator it = stats.begin(); it!=stats.end(); it++)
    {
        std::cout << std::fixed << std::setprecision(3) << "Iter :" << it->iter <<  " " << it->primalS  << " " << it->dualS << " " << it->primalD << " " << it->dualD << std::endl;
    }

}

inline void printAll(std::vector<OptInfo> & stats)
{
    int count = 0;
    for(std::vector<OptInfo>::iterator it = stats.begin(); it!=stats.end(); it++)
    {
        std::cout << std::fixed << std::setprecision(3) << "Iter :" << it->iter <<  " " << it->primalS  << " " << it->dualS << " " << it->primalD << " " << it->dualD << " " << it->jDistN << " " << it->jDistD << std::endl;
    }

}
#endif
