graph_name=".txt"
k=20
budget="50 100 150 200 250"
for b in budget
    do
        mkdir build
        cd FastCMs/build
        ./fastcm ../../datasets/${graph_name} ${k} ${b}