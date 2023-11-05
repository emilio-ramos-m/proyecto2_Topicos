#include "Hyperloglog.h"
#include "CountMin_CU.h"
#include <fstream> 

using namespace std;
using namespace sdsl;
using namespace chrono;

void insertKmer(HyperLogLog *h, vector<string> file,int k);
vector<string> loadFile(string filename);


int main(){ 
    cout<<"\n### Evaluacion Experimental ###\n"<<endl;

    cout<<"## HyperLogLog ##"<<endl;
    HyperLogLog hll1,hll2;
    
    // Cargar los archivos
    cout<<"Cargando genoma 1..."<<endl;
    vector<string> file1 = loadFile("GCF_000717965.1_ASM71796v1_genomic.fna");
    cout<<"Cargando genoma 2..."<<endl;
    vector<string> file2 = loadFile("GCF_001182945.1_P44_Wales_1_VIM_2_11_12_genomic.fna");
    
    cout<<"Insertando k-mers..." <<endl;
    insertKmer(&hll1,file1,10);
    insertKmer(&hll2,file2,10);
    cout<<endl;

    // Calcula la cardinalidad estimada
    cout<<"# Cardinalidad estimada #"<<endl;
    cout<<"Genoma 1: " << hll1.estimateCardinality() << endl;
    cout<<"Genoma 2: " << hll2.estimateCardinality() << endl;
    cout<<endl;

    // Razones de compresion
    cout<<"# Razones de compresion #"<<endl;
    wm_int<rrr_vector<15>> HLL_wm_int1 = hll1.compress_wm_int();
    wt_huff<rrr_vector<15>> HLL_wt_huff1 = hll1.compress_wt_huff();
    wm_int<rrr_vector<15>> HLL_wm_int2 = hll2.compress_wm_int();
    wt_huff<rrr_vector<15>> HLL_wt_huff2 = hll2.compress_wt_huff();
    int size_hll1 = hll1.sizeInBytes();
    int size_hll2 = hll2.sizeInBytes();
    int size_HLL_wm_int1 = size_in_bytes(HLL_wm_int1);
    int size_HLL_wt_huff1 = size_in_bytes(HLL_wt_huff1);
    int size_HLL_wm_int2 = size_in_bytes(HLL_wm_int2);
    int size_HLL_wt_huff2 = size_in_bytes(HLL_wt_huff2);
    

    cout<<"Tamaño de HLL1 "<<hll1.sizeInBytes()<<" bytes"<<endl;
    cout<<"Tamaño de compresion wm_int "<<size_HLL_wm_int1<<" bytes"<<endl;
    cout<<"Razon de compresion wm_int "<<size_HLL_wm_int1/(double)size_hll1<<" bytes"<<endl;
    cout<<"Tamaño de compresion wt_huff "<<size_HLL_wt_huff1<<" bytes"<<endl;
    cout<<"Razon de compresion wt_huff "<<size_HLL_wt_huff1/(double)size_hll1<<" bytes"<<endl;
    cout<<endl;

    cout<<"Tamaño de HLL2 "<<hll2.sizeInBytes()<<" bytes"<<endl;
    cout<<"Tamaño de compresion wm_int "<<size_HLL_wm_int2<<" bytes"<<endl;
    cout<<"Razon de compresion wm_int "<<size_HLL_wm_int2/(double)size_hll2<<" bytes"<<endl;
    cout<<"Tamaño de compresion wt_huff "<<size_HLL_wt_huff2<<" bytes"<<endl;
    cout<<"Razon de compresion wt_huff "<<size_HLL_wt_huff2/(double)size_hll2<<" bytes"<<endl;
    cout<<endl;

    // Union de los dos HLL
    cout<<"# Tiempos Union de HLL1 y HLL2 #"<<endl;
    HyperLogLog hll11 = hll1; //original
    HyperLogLog hll12 = hll1; //comp. wm_int
    HyperLogLog hll13 = hll1; //comp. wt_huff
    hll12.Union_wm_int(HLL_wm_int1, HLL_wm_int2);
    hll13.Union_wt_huff(HLL_wt_huff1, HLL_wt_huff2);
    hll1.Union(hll2);

    //Comparacion de tiempos de computos con/sin comprimir Union
    int time_Union = 0, time_Union_wm_int = 0, time_Union_wt_huff = 0;
    int count = 20, n = 20;
    while(count > 0){
        //Union sin compresion
        HyperLogLog aux = hll11;
        auto start = chrono::high_resolution_clock::now();
        aux.Union(hll2);
        auto end = chrono::high_resolution_clock::now();
        time_Union += duration_cast<nanoseconds> (end - start).count();

        //Union con wm_int
        aux = hll11;
        start = chrono::high_resolution_clock::now();
        aux.Union_wm_int(HLL_wm_int1, HLL_wm_int2);
        end = chrono::high_resolution_clock::now();
        time_Union_wm_int += duration_cast<nanoseconds> (end - start).count();

        //Union con wt_huff
        aux = hll11;
        start = chrono::high_resolution_clock::now();
        aux.Union_wt_huff(HLL_wt_huff1, HLL_wt_huff2);
        end = chrono::high_resolution_clock::now();
        time_Union_wt_huff += duration_cast<nanoseconds> (end - start).count();

        count--;
    }
    
    cout<<"Union HLL1-HLL2: "<<time_Union/n<<" ns"<<endl;
    cout<<"Union HLL1-HLL2 con wm_int: "<<time_Union_wm_int/n<<" ns"<<endl;
    cout<<"Union HLL1-HLL2 con wt_huff: "<<time_Union_wt_huff/n<<" ns"<<endl;
    cout<<endl;
    cout<<endl;

    // CountMinCU
    cout<<"## CountMinCU ##"<<endl;
    CountMinCU countMinCu(100000, 4);
    unordered_map<uint32_t, int> frec;
    
    // Cargar el archivo
    cout<<"Cargando dataset..."<<endl;
    vector<string> file = loadFile("Chicago-20080319.txt");
    // Insertar los datos
    cout<<"Insertando datos..."<<endl;
    uint32_t num;
    for(int i = 0; i < file.size(); i++){
        string aux = file[i];
        stringstream ss(aux);
        ss>>num;
        countMinCu.insert(num);
        frec[num] = frec.find(num) == frec.end() ? 1 : frec[num] + 1;
    }
    cout<<endl;
    
    //Razones de compresion
    wm_int<rrr_vector<15>> CMCU_wm_int = countMinCu.compress_wm_int();
    wt_huff<rrr_vector<15>> CMCU_wt_huff = countMinCu.compress_wt_huff();
    int size_cmcu = countMinCu.sizeInBytes();
    int size_CMCU_wm_int = size_in_bytes(CMCU_wm_int);
    int size_CMCU_wt_huff = size_in_bytes(CMCU_wt_huff);
    
    cout<<"# Razones de compresion #"<<endl;
    cout<<"Tamaño de CMCU "<<size_cmcu<<" bytes"<<endl;
    cout<<"Tamaño de compresion wm_int "<<size_CMCU_wm_int<<" bytes"<<endl;
    cout<<"Razon de compresion wm_int "<<size_CMCU_wm_int/(double)size_cmcu<<" bytes"<<endl;
    cout<<"Tamaño de compresion wt_huff "<<size_CMCU_wt_huff<<" bytes"<<endl;
    cout<<"Razon de compresion wt_huff "<<size_CMCU_wt_huff/(double)size_cmcu<<" bytes"<<endl;
    cout<<endl;
    
    // Tiempos de estimacion de frecuencia
    cout<<"# Tiempos de estimacion de frecuencia #"<<endl;
    int time_cmcu = 0, time_cmcu_wm_int = 0, time_cmcu_wt_huff = 0;
    
    // Tiempo countMinCu
    auto start = high_resolution_clock::now();
    for(auto it : frec){
        countMinCu.estimate(it.first);
    }
    auto end = high_resolution_clock::now();
    time_cmcu += duration_cast<nanoseconds> (end - start).count();

    // Tiempo countMinCu con wm_int
    start = high_resolution_clock::now();
    for(auto it : frec){
        countMinCu.estimate_wm_int(it.first,CMCU_wm_int);
    }
    end = high_resolution_clock::now();
    time_cmcu_wm_int += duration_cast<nanoseconds> (end - start).count();

    // Tiempo countMinCu con wt_huff
    start = high_resolution_clock::now();
    for(auto it : frec){
        countMinCu.estimate_wt_huff(it.first,CMCU_wt_huff);
    }
    end = high_resolution_clock::now();
    time_cmcu_wt_huff += duration_cast<nanoseconds> (end - start).count();
    
    cout<<"Tiempo de estimacion de frecuencia de CMCU: "<<time_cmcu/frec.size()<<" ns"<<endl;
    cout<<"Tiempo de estimacion de frecuencia de CMCU con wm_int: "<<time_cmcu_wm_int/frec.size()<<" ns"<<endl;
    cout<<"Tiempo de estimacion de frecuencia de CMCU con wt_huff: "<<time_cmcu_wt_huff/frec.size()<<" ns"<<endl;

    return 0;
}

void insertKmer(HyperLogLog *h, vector<string> file,int k){
    string kmer;
    for(int i = 0; i < file.size(); i++){
        if(file[i].size() < k) continue; 
        for(int j = 0; j < file[i].size()-k+1; j++){
            kmer = file[i].substr(j,k);
            h->insert(kmer);
        }
    }
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