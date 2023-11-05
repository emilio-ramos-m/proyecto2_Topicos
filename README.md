# Proyecto2_Topicos
## Integrantes:
- **Emilio Ramos Montesino**
- **Esteban Chandia Cifuentes**

## Descripción:

### Bibliotecas utilizadas:
- **sdsl-lite:** Para la comprimir los sketches con las estructuras wm_int y wt_huff.
#### Instalación:
```
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
sudo ./install.sh /usr/local/
```

### Compilación y Ejecución:
```
g++ -std=c++11 *.cpp hash_functions/*.cpp  -O3 -DNDEBUG -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64 && ./a.out
```