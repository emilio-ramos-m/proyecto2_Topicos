#ifndef _HYPERLOGLOG_H_
#define _HYPERLOGLOG_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#include <sdsl/wm_int.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/wt_huff.hpp>

using namespace sdsl;
using namespace std;


class HyperLogLog {
private:
    int precision;
    int M;
    double alpha;
    std::vector<uint8_t> registers;

public:
    HyperLogLog(int precision = 14);
    void insert(const std::string& kmer);
    double estimateCardinality();
    void Union(const HyperLogLog& other);
    void Union_wm_int(wm_int<rrr_vector<15>> hll1, wm_int<rrr_vector<15>> hll2);
    void Union_wt_huff(wt_huff<rrr_vector<15>> hll1, wt_huff<rrr_vector<15>> hll2);
    wm_int<rrr_vector<15>> compress_wm_int();
    wt_huff<rrr_vector<15>> compress_wt_huff();
    size_t sizeInBytes();
};

#endif