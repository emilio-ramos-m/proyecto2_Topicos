#include "Hyperloglog.h"
#include "CountMin_CU.h"
#include <fstream> 


int main(){
    HyperLogLog hll;

    // Abre el archivo de k-mers
    std::string filename = "GCF_000717965.1_ASM71796v1_genomic.fna";
    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        std::cerr << "Error al abrir el archivo: " << filename << std::endl;
        return 0;
    }

    // Elimina la primera línea
    std::string line;
    std::getline(inputFile, line);
    // Inserta el k-mer en el HyperLogLog
    while (std::getline(inputFile, line)) { 
        hll.insert(line);
    }
    
    // Cierra el archivo
    inputFile.close();

    // Calcula la cardinalidad estimada
    double estimate = hll.estimateCardinality();
    std::cout << "Cardinalidad estimada: " << estimate << std::endl;

    return 0;
}