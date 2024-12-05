
cd ..
export USE_OPENMP=0
g++ -fopenmp -O3 -o ./search ./search.cpp -I ./common.h -DOPEN_MP=$USE_OPENMP

result_path=./results/
datasets=('sift')
K=100
num_calc=-1

for data in "${datasets[@]}"
do

    data_path=/data/vector_datasets/${data}
    index_path=/data/tmp/hcnng/${data}
    data_file="${data_path}/${data}_base.fvecs"
    index_file="${index_path}/hcnng_${data}.ivecs"
    query="${data_path}/${data}_query.fvecs"
    gnd="${data_path}/${data}_groundtruth.ivecs"
    result_file="${result_path}/${data}_K${K}_num_calc${num_calc}.log"

    # K: number of neighbors
    # num_calc: maximum number of calculations
    ./search $data_file $query $gnd $index_file $K $num_calc $result_file
done


