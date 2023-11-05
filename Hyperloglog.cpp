#include "Hyperloglog.h"
#include "hash_functions/MurmurHash3.h"



HyperLogLog::HyperLogLog(int precision) {
    this->precision = precision;
    M = pow(2, precision);
    switch (M) {
        case 16: alpha = 0.673; break;
        case 32: alpha = 0.697; break;
        case 64: alpha = 0.709; break;
        default: alpha = 0.7213 / (1 + 1.079 / M); break;
    }
    registers = vector<uint8_t>(M, 0);
}

void HyperLogLog::insert(string kmer) {
    uint32_t hash;
    MurmurHash3_x86_32(kmer.c_str(), kmer.size(), 1, &hash);
    unsigned int p = hash >> 32 - precision;
    unsigned int b = hash << precision;
   
    uint8_t w = 1;
    unsigned int aux = 1 << 31;
    while ((b & aux) == 0 && aux != 0) {
        w++;
        aux = aux >> 1;
    }
    registers[p] = max(registers[p], w);
}

double HyperLogLog::estimateCardinality() {
    // Computar el estimador
    double sum = 0;
    for (int i = 0; i < M; i++) {
        sum += 1/pow(2, registers[i]);
    }
    double estimate = alpha * M * M / sum;

    if (estimate <= 5.0 / 2.0 * M) {
        int zeroRegisters = count(registers.begin(), registers.end(), 0);
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

void HyperLogLog::Union(HyperLogLog other) {
    for (int i = 0; i < M; i++) {
        registers[i] = max(registers[i], other.registers[i]);
    }
}

void HyperLogLog::Union_wm_int(wm_int<rrr_vector<15>> hll1, wm_int<rrr_vector<15>> hll2) {
    if(hll1.size() != hll2.size()){
        cout<<"Los tamaños de los hll no son iguales"<<endl;
        exit(0);
    }
    int correcion = 8;
    for(int i = 0; i < registers.size(); i++){
        registers[i] = max(hll1[i+correcion], hll2[i+correcion]);
    }
}

void HyperLogLog::Union_wt_huff(wt_huff<rrr_vector<15>> hll1, wt_huff<rrr_vector<15>> hll2) {
    if(hll1.size() != hll2.size()){
        cout<<"Los tamaños de los hll no son iguales"<<endl;
        exit(0);
    }
    int correcion = 8;
    for(int i = 0; i < registers.size(); i++){
        registers[i] = max(hll1[i+correcion], hll2[i+correcion]);
    }
}

wm_int<rrr_vector<15>> HyperLogLog::compress_wm_int() {
    wm_int<rrr_vector<15>> wm_int;
    construct_im(wm_int, registers, 1);
    return wm_int;
}

wt_huff<rrr_vector<15>> HyperLogLog::compress_wt_huff() {
    wt_huff<rrr_vector<15>> wt_huff;
    construct_im(wt_huff, registers, 1);
    return wt_huff;
} 

size_t HyperLogLog::sizeInBytes() {
    return registers.size() * sizeof(uint8_t);
}

