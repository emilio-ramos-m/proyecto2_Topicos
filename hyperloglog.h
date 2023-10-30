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
    HyperLogLog(int precision);
    void insert(const std::string& kmer);
    double estimateCardinality();
    void merge(const HyperLogLog& other);
};

#endif