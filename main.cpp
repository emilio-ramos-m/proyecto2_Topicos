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
        std::vector<std::string> kmers = generarKmers(line, 5);
        for(int i = 0; i < kmers.size(); i++){
            hll.insert(kmers[i]);
        }
    }
    
    // Cierra el archivo
    inputFile.close();

    // Calcula la cardinalidad estimada
    double estimate = hll.estimateCardinality();
    std::cout << "Cardinalidad estimada: " << estimate << std::endl;

    return 0;
}