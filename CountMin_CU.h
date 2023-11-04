#ifndef _COUNTMINCU_H_
#define _COUNTMINCU_H_

#include <iostream>
#include <vector>
#include <functional>
#include <limits.h>

#include <sdsl/wm_int.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/wt_huff.hpp>

using namespace std;
using namespace sdsl;

class CountMinCU {
private:
    int depth;
    int width;
    vector<int> sketch;

public:
    CountMinCU(int width, int depth);
    void insert(uint32_t element, int delta = 1);
    int estimate(uint32_t element);
    int estimate_wt_huff(uint32_t element,wt_huff<rrr_vector<15>> CMCU_wt_huff);
    int estimate_wm_int(uint32_t element,wm_int<rrr_vector<15>> CMCU_wm_int);
    wm_int<rrr_vector<15>> compress_wm_int();
    wt_huff<rrr_vector<15>> compress_wt_huff();
    size_t sizeInBytes();
};

#endif

