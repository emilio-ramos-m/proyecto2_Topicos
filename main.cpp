#include "Hyperloglog.h"
#include "CountMin_CU.h"
#include <fstream> 

std::vector<std::string> generarKmers(std::string line, int k){
    std::vector<std::string> kmers;
    std::string kmer_aux;
    int i = 0;
    while(i < line.size()){
        if(i+k <= line.size()){
            kmer_aux = line.substr(i, k);
            kmers.push_back(kmer_aux);
            i++;
        }else{
            break;
        }
    }
    return kmers;
}

int main(){
    HyperLogLog hll;

    // Abre el archivo de k-mers
    std::string filename = "GCF_000717965.1_ASM71796v1_genomic.fna";
    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        std::cerr << "Error al abrir el archivo: " << filename << std::endl;
        return 0;
    }

    // Inserta k-mers en el HyperLogLog
    std::string line;
    while (std::getline(inputFile, line)) {
        if(line[0]=='>') {
            continue;
        }
        // Extrae e inserta los k-mers de la línea 
        std::vector<std::string> kmers = generarKmers(line, 4);
        for(int i = 0; i < kmers.size(); i++){
            hll.insert(kmers[i]);
        }
    }
    
    // Cierra el archivo
    inputFile.close();

    // Calcula la cardinalidad estimada
    double estimate = hll.estimateCardinality();
    std::cout << "Cardinalidad estimada: " << estimate << std::endl;
    std::cout<<"Bytes HLL: "<<hll.sizeInBytes()<<std::endl;
    std::cout<<"Bytes HLL comprimido en wm_int: "<<hll.compress_wm_int()<<std::endl;
    std::cout<<"Bytes HLL comprimido en wt_huff: "<<hll.compress_wt_huff()<<std::endl;

    return 0;
}