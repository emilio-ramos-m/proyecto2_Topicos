#ifndef _COUNTMINCU_H_
#define _COUNTMINCU_H_

#include <iostream>
#include <vector>
#include <functional>
#include <limits.h>


class CountMinCU {
private:
    int num_hashes;
    int table_size;
    std::vector<std::vector<int>> sketch;

public:
    CountMinCU(int width, int depth);
    void insert(uint32_t element, int delta = 1);
    int estimate(uint32_t element);
    //sdsl::wt_huff<sdsl::sd_vector<>> compress_wt_huff();
    //sdsl::wm_int<sdsl::sd_vector<>> compress_wm_int();
};

#endif

