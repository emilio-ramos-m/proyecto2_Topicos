#include "hyperloglog.h"
#include "MurmurHash3.h"


HyperLogLog::HyperLogLog(int precision) {
    precision = precision;
    M = pow(2, precision);
    if (M == 16) alpha = 0.673;
    else if (M == 32) alpha = 0.697;
    else if (M == 64) alpha = 0.709;
    else alpha = 0.7213 / (1 + 1.079 / M);
    registers = std::vector<int>(M, 0);
}

void HyperLogLog::insert(const std::string& kmer) {
    uint32_t hash;
    MurmurHash3_x86_32(kmer.c_str(), kmer.size(), 0, &hash);
    unsigned int p = hash >> (32 - precision);
    unsigned int b = hash << precision;
   
    int w = 1;
    int aux = 1 << precision;
    while ((b & aux) == 0 && aux != 0) {
        w++;
        aux >>= 1;
    }
    registers[p] = std::max(registers[p], w);
}

double HyperLogLog::estimateCardinality() {
    // Computar el estimador
    double sum = 0;
    for (int i = 0; i < M; i++) {
        sum += pow(2, -registers[i]);
    }
    double estimate = alpha * M * M / sum;

    if (estimate <= 5.0 / 2.0 * M) {
        int zeroRegisters = std::count(registers.begin(), registers.end(), 0);
        if (zeroRegisters != 0) {
            // Corrección de pequeños valores
            estimate = M * log(static_cast<double>(M) / zeroRegisters);
        }
    } else if (estimate > 1.0 / 30.0 * pow(2, 32)) {
        // Corrección de grandes valores
        estimate = -pow(2, 32) * log(1.0 - estimate / pow(2, 32));
    }

    return estimate;
}

void HyperLogLog::merge(const HyperLogLog& other) {
    for (int i = 0; i < M; i++) {
        registers[i] = std::max(registers[i], other.registers[i]);
    }
}