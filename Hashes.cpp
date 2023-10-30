#include "Hashes.h"


Hashes::Hashes() {
}

int Hashes::hash(uint32_t element, int i, int table_size) {
    uint32_t hash_value;
    switch (i%5) {
    case 0:
        uint32_t hash;
        MurmurHash3_x86_32(&element, sizeof(int), i, &hash);
        hash_value = hash % table_size;
        break;
    case 1:
        hash_value = MurmurHash2(&element, sizeof(int), i) % table_size;
        break;
    case 2:
        hash_value = MurmurHash1(&element, sizeof(int), i) % table_size;
        break;
    case 3:
        uint32_t hash1;
        sha1_32a(&element, sizeof(int), i, &hash1);
        hash_value = hash1 % table_size;
        break;
    case 4:
        hash_value = lookup3(&element, sizeof(int), i) % table_size;
        break;           
    default:
        break;
    }
    return hash_value;
}