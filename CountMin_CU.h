#ifndef _COUNTMINCU_H_
#define _COUNTMINCU_H_

#include <iostream>
#include <vector>
#include <functional>
#include <limits.h>

#include <sdsl/wm_int.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/wt_huff.hpp>

class CountMinCU {
private:
    int num_hashes;
    int table_size;
    std::vector<std::vector<int>> sketch;

public:
    CountMinCU(int width, int depth);
    void insert(uint32_t element, int delta = 1);
    int estimate(uint32_t element);
    int estimate_compress(uint32_t element,sdsl::wm_int<sdsl::rrr_vector<15>> CMCU_wm_int);
    size_t sizeInBytes();
    sdsl::wm_int<sdsl::rrr_vector<15>> compress_wm_int();
    //sdsl::wt_huff<sdsl::rrr_vector<15>> compress_wt_huff();
};

#endif

