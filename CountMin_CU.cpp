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

int CountMinCU::estimate_compress(uint32_t element, sdsl::wm_int<sdsl::rrr_vector<15>> CMCU_wm_int){
    int min_count = INT_MAX;
    Hashes hashes;
    for (int i=0 ; i<num_hashes; i++) {
        min_count = min(min_count, (int)CMCU_wm_int[hashes.hash(element, i,table_size) + num_hashes*i]);
    }
    return min_count;
}

size_t CountMinCU::sizeInBytes() {
    return sizeof(int) * num_hashes * table_size;
}

sdsl::wm_int<sdsl::rrr_vector<15>> CountMinCU::compress_wm_int(){
    sdsl::int_vector<> sketch_int(num_hashes*table_size);
    for(int i=0;i<num_hashes;i++){
       for(int j=0;j<table_size;j++){
           sketch_int[i*table_size+j]=sketch[i][j];
       }
    }
    sdsl::wm_int<sdsl::rrr_vector<15>> wm_int;
    sdsl::construct_im(wm_int, sketch_int, 4);
    return wm_int;
}
/*
sdsl::wt_huff<sdsl::rrr_vector<15>> CountMinCU::compress_wt_huff(){
    sdsl::int_vector<> sketch_int(num_hashes*table_size);
    for(int i=0;i<num_hashes;i++){
       for(int j=0;j<table_size;j++){
           sketch_int[i*table_size+j]=sketch[i][j];
       }
    }
    sdsl::wt_huff<sdsl::rrr_vector<15>> wt_huff;
    sdsl::construct_im(wt_huff, sketch_int, 4);
    return wt_huff;
}
*/