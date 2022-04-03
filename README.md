## FastCM+
This is the source code for the paper "Xin Sun, Xin Huang, Di Jin. Fast Algorithms for Core Maximization on Large Graphs" published at PVLDB 2022. 
## Execute
```
cd FastCMs/
mkdir build
cd build
cmake ..
make
./fastcm ../../dataset/xxx.txt 
```
Example:
```
./fastcm ../../dataset/Brightkite_edges.txt 20 200
```
The results will be written into results.txt.

## Experiment Datasets
You can download the graphs used in our paper by following the instructions in Section 7. We list some of the datasets.

## Acknowledgement
If you have any further questions, please feel free to contact me.
If you used this code, please kindly cite the paper.
