## FastCM+
This is the source code for the paper "Xin Sun, Xin Huang, and Di Jin. Fast Algorithms for Core Maximization on Large Graphs" published at Proceedings of the VLDB Endowment (2022). 
## Execute
```
cd src/
mkdir build
cd build
cmake ..
make
./fastcm ../../dataset/xxx.txt "k" "b"
```
Example:
```
./fastcm ../../dataset/Brightkite_edges.txt 20 200
```
The results will be written into results.txt.

## Experiment Datasets
As the graph datasets are large, we do not upload them. You can download the graphs used in our paper by following the instructions in Section 7.

## Note
If you have any further questions, please feel free to contact me.
Please cite our paper, if you use our source code.
