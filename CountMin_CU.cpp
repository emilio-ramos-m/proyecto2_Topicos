#include "CountMin_CU.h"
#include "Hashes.h"


CountMinCU::CountMinCU(int width, int depth) {
    this->depth = depth;
    this->width = width;
    sketch = vector<int>(depth*width, 0);
}

void CountMinCU::insert(uint32_t element, int delta) {
    Hashes hashes;
    int j, min_frec = estimate(element);
    for (int i = 0; i < depth; i++) {
        j = hashes.hash(element, i,width);
        if(sketch[i*width+j] == min_frec){
            sketch[i*width+j] += delta;
        }
    }
}

int CountMinCU::estimate(uint32_t element) {
    int min_count = INT_MAX;
    Hashes hashes;  
    for (int i = 0; i < depth; i++) {
        int j = hashes.hash(element, i,width);
        min_count = min(min_count, sketch[i*width+j]);
    }
    return min_count;
}

int CountMinCU::estimate_wt_huff(uint32_t element,wt_huff<rrr_vector<15>> CMCU_wt_huff){
    int min_count = INT_MAX, correcion = 2;
    Hashes hashes;
    for (int i=0 ; i<depth; i++) {
        int j = hashes.hash(element, i,width);
        min_count = min(min_count, (int)CMCU_wt_huff[i*width+(j+correcion)]);
    }
    return min_count;
}

int CountMinCU::estimate_wm_int(uint32_t element, wm_int<rrr_vector<15>> CMCU_wm_int){
    int min_count = INT_MAX, correcion = 2;
    Hashes hashes;
    for (int i=0 ; i<depth; i++) {
        int j = hashes.hash(element, i,width);
        min_count = min(min_count, (int)CMCU_wm_int[i*width+(j+correcion)]);
    }
    return min_count;
}


wm_int<rrr_vector<15>> CountMinCU::compress_wm_int(){
    wm_int<rrr_vector<15>> wm_int;
    construct_im(wm_int, sketch, 4);
    return wm_int;
}

wt_huff<rrr_vector<15>> CountMinCU::compress_wt_huff(){
    wt_huff<rrr_vector<15>> wt_huff;
    construct_im(wt_huff, sketch, 4);
    return wt_huff;
}

size_t CountMinCU::sizeInBytes() {
    return sizeof(int) * depth * width;
}