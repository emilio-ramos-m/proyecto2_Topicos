#ifndef _HYPERLOGLOG_H_
#define _HYPERLOGLOG_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#include <sdsl/wm_int.hpp>
#include <sdsl/csa_wt.hpp>
#include <sdsl/wt_huff.hpp>


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
    size_t sizeInBytes();
    uint32_t compress_wm_int();
    uint32_t compress_wt_huff();
};

#endif