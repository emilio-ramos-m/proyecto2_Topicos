#ifndef _HASHES_H_
#define _HASHES_H_

#include "hash_functions/MurmurHash3.h"
#include "hash_functions/MurmurHash2.h"
#include "hash_functions/MurmurHash1.h"
#include "hash_functions/Platform.h"
#include "hash_functions/sha1.h"

class Hashes {
    public:
        Hashes();
        int hash(uint32_t element, int i, int table_size);
};

#endif