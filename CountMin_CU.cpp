#include "CountMin_CU.h"
#include "Hashes.h"

using namespace std;

CountMinCU::CountMinCU(int width, int depth) {
    num_hashes = depth;
    table_size = width;
    sketch = vector<vector<int>>(depth, vector<int>(width, 0));
}

void CountMinCU::insert(uint32_t element, int delta) {
    Hashes hashes;
    int j, min_frec = estimate(element);
    for (int i = 0; i < num_hashes; i++) {
        j = hashes.hash(element, i,table_size);
        if(sketch[i][j] == min_frec){
            sketch[i][j] += delta;
        }
    }
}

int CountMinCU::estimate(uint32_t element) {
    int min_count = INT_MAX;
    Hashes hashes;  
    for (int i = 0; i < num_hashes; i++) {
        min_count = min(min_count, sketch[i][hashes.hash(element, i,table_size)]);
    }
    return min_count;
}

/*

sdsl::wt_huff<sdsl::sd_vector<>> compress_wt_huff() {
        // Representar 'sketch' como una secuencia de bits
        sdsl::bit_vector bv(num_hashes * table_size, 0);
        for (int i = 0; i < num_hashes; ++i) {
            for (int j = 0; j < table_size; ++j) {
                bv[i * table_size + j] = (sketch[i][j] != 0);
            }
        }

        // Crear estructuras de datos de sdsl-lite
        sdsl::wt_huff<sdsl::sd_vector<>> compressed(bv);

        return compressed;
}

sdsl::wm_int<sdsl::sd_vector<>> compress_wm_int() {
        // Representar 'sketch' como una secuencia de bits
        sdsl::bit_vector bv(num_hashes * table_size, 0);
        for (int i = 0; i < num_hashes; ++i) {
            for (int j = 0; j < table_size; ++j) {
                bv[i * table_size + j] = (sketch[i][j] != 0);
            }
        }

        // Crear estructuras de datos de sdsl-lite
        sdsl::wm_int<sdsl::sd_vector<>> compressed(bv);

        return compressed;
}

*/