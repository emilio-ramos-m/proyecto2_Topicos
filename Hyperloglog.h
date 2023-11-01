#ifndef _HYPERLOGLOG_H_
#define _HYPERLOGLOG_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>



class HyperLogLog {
private:
    int precision;
    int M;
    double alpha;
    std::vector<int> registers;

public:
    HyperLogLog(int precision = 14);
    void insert(const std::string& kmer);
    double estimateCardinality();
    void merge(const HyperLogLog& other);
    //sdsl::wt_huff<sdsl::csa_wt<sdsl::wt_huff<>>> compress_wt_huff();
    //sdsl::wm_int<sdsl::sd_vector<>> compress_wm_int();
    //sdsl::wm_int<sdsl::rrr_vector<15>> HyperLogLog::compress_wm_int();
};

#endif