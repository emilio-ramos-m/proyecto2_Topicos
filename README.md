# Proyecto2_Topicos
## Integrantes:
- **Emilio Ramos Montesino**
- **Esteban Chandia Cifuentes**

## Descripción:
En este proyecto se implementan los algoritmos de sketching HyperLogLog y CountMin-CU, los cuales soportarán las estructuras de compresión wm_int y wt_huff de la biblioteca [SDSL-Lite](https://github.com/simongog/sdsl-lite) para reducir el espacio de almacenamiento de los sketches. Adicionalmente se realizarán experimentos para comparar el desempeño de los algoritmos  con y sin compresión. Evaluando el tiempo de ejecución y el espacio de almacenamiento de los sketches. Las funciones hash utilizadas fueron extraídas del repositorio [Smhasher](https://github.com/rurban/smhasher).

### Bibliotecas utilizadas:
- **sdsl-lite:** Para la comprimir los sketches con las estructuras wm_int y wt_huff.
#### Instalación:
```
git clone https://github.com/simongog/sdsl-lite.git
cd sdsl-lite
sudo ./install.sh /usr/local/
```

### Compilación y Ejecución:
Para compilar y ejecutar el proyecto se debe ejecutar el siguiente comando:
```
g++ -std=c++11 *.cpp hash_functions/*.cpp  -O3 -DNDEBUG -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64 && ./a.out
```