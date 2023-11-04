#include "Hyperloglog.h"
#include "CountMin_CU.h"
#include <fstream> 

using namespace std;
using namespace sdsl;

vector<string> generarKmers(string line, int k);
vector<string> loadFile(string filename);

int main(int argc, char* argv[]){ 
    string input = argv[1];
    
    if (input == "hll") {
        HyperLogLog hll1,hll2;

        // Cargar primer archivo
        vector<string> file1 = loadFile("GCF_000717965.1_ASM71796v1_genomic.fna");
        // Insertar los kmers
        for(int i = 0; i < file1.size(); i++){
            vector<string> kmers = generarKmers(file1[i], 4);
            for(int j = 0; j < kmers.size(); j++){
                hll1.insert(kmers[j]);
            }
        }

        // Cargar segundo archivo
        vector<string> file2 = loadFile("GCF_000006945.2_ASM694v2_genomic.fna");
        // Insertar los kmers
        for(int i = 0; i < file2.size(); i++){
            vector<string> kmers = generarKmers(file2[i], 4);
            for(int j = 0; j < kmers.size(); j++){
                hll2.insert(kmers[j]);
            }
        }
        
        // Calcula la cardinalidad estimada
        double estimate = hll1.estimateCardinality();
        cout << "Cardinalidad estimada: " << estimate << endl;
        cout<<"Tamaño de HLL "<<hll1.sizeInBytes()<<endl;

        // Razones de compresion     
        wm_int<rrr_vector<15>> HLL_wm_int1 = hll1.compress_wm_int();
        cout<<"Tamaño de compresion wm_int "<<size_in_bytes(HLL_wm_int1)<<" bytes"<<endl;
        cout<<"Razon de compresion wm_int "<<size_in_bytes(HLL_wm_int1)/(double)hll1.sizeInBytes()<<" bytes"<<endl;
        cout<<endl;
        wt_huff<rrr_vector<15>> HLL_wt_huff1 = hll1.compress_wt_huff();
        cout<<"Tamaño de compresion wt_huff "<<size_in_bytes(HLL_wt_huff1)<<" bytes"<<endl;
        cout<<"Razon de compresion wt_huff "<<size_in_bytes(HLL_wt_huff1)/(double)hll1.sizeInBytes()<<" bytes"<<endl;
        cout<<endl;   

        // Union de los dos HLL
        hll1.Union(hll2);

        // Compresion de los HLL
        wm_int<rrr_vector<15>> HLL_wm_int2 = hll2.compress_wm_int();
        hll1.Union_wm_int(HLL_wm_int1,HLL_wm_int2);

        wt_huff<rrr_vector<15>> HLL_wt_huff2 = hll2.compress_wt_huff();
        hll1.Union_wt_huff(HLL_wt_huff1,HLL_wt_huff2);



    } else if (input == "cmcu") {    
        CountMinCU countMinCu(100000, 4);
        unordered_map<uint32_t, int> frec;

        // Cargar el archivo
        vector<string> file = loadFile("Chicago-20080319.txt");
        // Insertar los datos
        uint32_t num;
        for(int i = 0; i < file.size(); i++){
            string aux = file[i];
            stringstream ss(aux);
            ss>>num;
            countMinCu.insert(num);
            frec[num] = frec.find(num) == frec.end() ? 1 : frec[num] + 1;
        }
        
        //Razones de compresion
        wm_int<rrr_vector<15>> CMCU_wm_int = countMinCu.compress_wm_int();
        cout<<"Tamaño de CMCU "<<countMinCu.sizeInBytes()<<" bytes"<<endl;
        cout<<"Tamaño de compresion wm_int "<<size_in_bytes(CMCU_wm_int)<<" bytes"<<endl;
        cout<<"Razon de compresion wm_int "<<size_in_bytes(CMCU_wm_int)/(double)countMinCu.sizeInBytes()<<" bytes"<<endl;
        cout<<endl;

        wt_huff<rrr_vector<15>> CMCU_wt_huff = countMinCu.compress_wt_huff();
        cout<<"Tamaño de compresion wt_huff "<<size_in_bytes(CMCU_wt_huff)<<" bytes"<<endl;
        cout<<"Razon de compresion wt_huff "<<size_in_bytes(CMCU_wt_huff)/(double)countMinCu.sizeInBytes()<<" bytes"<<endl;
        cout<<endl;
        
        // Estimacion de frecuencias
        for(auto it : frec){
            int estimacion_CMCU = countMinCu.estimate(it.first);
            int estimacion_wt_huff = countMinCu.estimate_wt_huff(it.first,CMCU_wt_huff);
            int estimacion_wm_int = countMinCu.estimate_wm_int(it.first,CMCU_wm_int);
        }

        
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

vector<string> loadFile(string filename){
    vector<string> file;
    string line;
    ifstream inputFile(filename);
    if(!inputFile){
        cerr << "No se pudo abrir el archivo " << filename << endl;
        return file;
    }
    while (getline(inputFile, line)) {
        if(line[0] == '>') continue;
        file.push_back(line);
    }
    inputFile.close();
    return file;
}