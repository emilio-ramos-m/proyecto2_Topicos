#include "Hyperloglog.h"
#include "CountMin_CU.h"
#include <fstream> 

using namespace std;
using namespace sdsl;

vector<string> generarKmers(string line, int k);

int main(int argc, char* argv[]){ 
    string input = argv[1];
    
    if (input == "hll") {
        HyperLogLog hll;

        // Abre el archivo de k-mers
        string filename = "GCF_000717965.1_ASM71796v1_genomic.fna";
        ifstream inputFile(filename);

        if (!inputFile.is_open()) {
            cerr << "Error al abrir el archivo: " << filename << endl;
            return 0;
        }

        // Inserta k-mers en el HyperLogLog
        string line;
        while (getline(inputFile, line)) {
            if(line[0]=='>') {
                continue;
            }
            // Extrae e inserta los k-mers de la línea 
            vector<string> kmers = generarKmers(line, 4);
            for(int i = 0; i < kmers.size(); i++){
                hll.insert(kmers[i]);
            }
        }
        
        // Cierra el archivo
        inputFile.close();

        // Calcula la cardinalidad estimada
        double estimate = hll.estimateCardinality();
        cout << "Cardinalidad estimada: " << estimate << endl;
        cout<<"Tamaño de HLL "<<hll.sizeInBytes()<<endl;

        //Razones de compresion     
        wm_int<rrr_vector<15>> HLL_wm_int = hll.compress_wm_int();
        cout<<"Tamaño de compresion wm_int "<<size_in_bytes(HLL_wm_int)<<" bytes"<<endl;
        cout<<"Razon de compresion wm_int "<<size_in_bytes(HLL_wm_int)/(double)hll.sizeInBytes()<<" bytes"<<endl;
        cout<<endl;
        wt_huff<rrr_vector<15>> HLL_wt_huff = hll.compress_wt_huff();
        cout<<"Tamaño de compresion wt_huff "<<size_in_bytes(HLL_wt_huff)<<" bytes"<<endl;
        cout<<"Razon de compresion wt_huff "<<size_in_bytes(HLL_wt_huff)/(double)hll.sizeInBytes()<<" bytes"<<endl;
        cout<<endl;   

    } else if (input == "cmcu") {    
        CountMinCU countMinCu(1000, 4);
        unordered_map<uint32_t, int> frec;
        
        //Abrir el archivo
        string line, filename = "Chicago-20080319.txt";
        ifstream file(filename);
        if(!file){
            cerr << "No se pudo abrir el archivo " << filename << endl;
            return 0;
        }

        //Leer el archivo
        while (getline(file, line)) {
            istringstream iss(line);
            uint32_t number;
            if (iss >> number) {
                countMinCu.insert(number);
                frec[number] = frec.find(number) == frec.end()? 1 : frec[number]+1;
            }
        }
        file.close();

        //Razones de compresion
        wm_int<rrr_vector<15>> CMCU_wm_int = countMinCu.compress_wm_int();
        cout<<"Tamaño de CMCU "<<countMinCu.sizeInBytes()<<" bytes"<<endl;
        cout<<"Tamaño de compresion wm_int "<<size_in_bytes(CMCU_wm_int)<<" bytes"<<endl;
        cout<<"Razon de compresion wm_int "<<size_in_bytes(CMCU_wm_int)/(double)countMinCu.sizeInBytes()<<" bytes"<<endl;
        cout<<endl;
        /*
        wt_huff<rrr_vector<15>> CMCU_wt_huff = countMinCu.compress_wt_huff();
        cout<<"Tamaño de compresion wt_huff "<<size_in_bytes(CMCU_wt_huff)<<" bytes"<<endl;
        cout<<"Razon de compresion wt_huff "<<size_in_bytes(CMCU_wt_huff)/(double)countMinCu.sizeInBytes()<<" bytes"<<endl;
        cout<<endl;
        */
    }else{
        cout<<"troleo mano"<<endl;
    }

    return 0;
}


    vector<string> generarKmers(string line, int k){
        vector<string> kmers;
        string kmer_aux;
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