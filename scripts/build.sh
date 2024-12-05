cd ..
g++ -fopenmp -O3 -o ./hcnng hcnng.cpp -I ./common.h

datasets=('sift' 'gist' )

min_cluster=1000
num_cluster=20

for data in "${datasets[@]}"
do  
    echo "Indexing - ${data}"

    data_path=/data/vector_datasets/${data}
    index_path=/data/tmp/hcnng/${data}

    if [ ! -d "$index_path" ]; then 
        mkdir -p "$index_path"
    fi

    data_file="${data_path}/${data}_base.fvecs"
    
    index_file="${index_path}/hcnng_${data}.ivecs"

    ./hcnng $data_file $min_cluster $num_cluster $index_file
done