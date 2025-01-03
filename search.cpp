// compile with: g++ search.cpp -o search -std=c++11 -fopenmp -O3

#if defined (_MSC_VER)  // Visual studio
    #define thread_local __declspec( thread )
#elif defined (__GCC__) // GCC
    #define thread_local __thread
#endif

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <map>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <tuple>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "common.h"
#include "utils.h"
#include <omp.h>
#include <random>
#include <time.h>
#include <thread>

#include <chrono>

using namespace std;

// Because of the random construction of the MST and combine them together
// If we use furthest_dist to prune during the `beam search phase`, 
// we cannot reach high recall (i.e., higher than 92%)

// However, if we only need not near-perfect recall, we can add pruning, it can improve the QPS

tuple<vector<int>, vector<float>> search_KNN(float *query, int K, AdjList &graph, Matrix<float> &points, int start, int max_calc){
	int N = points.rows;
	int calc_left = max_calc-1;
	unordered_set<int> visited; visited.insert(start);
	priority_queue<tuple<float, int> > q, knn;
	float furthest_dist = dist_L2(points[start], query, points.cols);
	q.push(make_tuple(-furthest_dist, start));
	knn.push(make_tuple(furthest_dist, start));
	while(!q.empty() && calc_left > 0){
		float d; int v;
		tie(d, v) = q.top();
		q.pop();
		// if(-d > furthest_dist) continue;
		for(int u : graph[v]){
			if(calc_left <= 0) break;
			if(in_set(u, visited))
				continue;
			visited.insert(u);
			d = dist_L2(points[u], query, points.cols);
			calc_left--;
			q.push(make_tuple(-d, u));
			knn.push(make_tuple(d, u));
			if(knn.size() > K)
				knn.pop();
			// if(knn.size() == K) {
			// 	tie(d, v) = knn.top();
			// 	furthest_dist = std::min(d, furthest_dist);
			// }
		}
	}

	vector<int> nearests;
	vector<float> dists;
	while(!knn.empty()){
		float x; int y;
		tie(x, y) = knn.top();
		nearests.push_back(y);
		dists.push_back(x);
		knn.pop();
	}
	reverse(nearests.begin(), nearests.end());
	reverse(dists.begin(), dists.end());
	return make_tuple(nearests, dists);
}


void run_on_testset(Matrix<float> &queries, int K, Matrix<float> &points, vector<vector<int> > &GT, AdjList &graph, int max_calc){
	float recall = 0;
	int N = points.rows;
	int num_queries = queries.rows;
	
	unsigned long long time = 0;

	StopW stopw = StopW();
#ifdef USE_OPENMP
	#pragma omp parallel for
#endif
	for(int i = 0; i < num_queries; i++){
		// choose start point randomly
		int start = rand_int(0, N-1);
		auto knn = search_KNN(queries[i], K, graph, points, start, max_calc);
#ifdef USE_OPENMP
		#pragma omp critical
#endif
		{
			recall += get_recall(GT[i], get<0>(knn), K);
		}
	}
	time += stopw.getElapsedTimeMicro();

	float time_us_per_query = time / num_queries;

	cout << "------------------------------------------------" << endl;
	cout << "K = " << K << " max_clac= " << max_calc <<  endl;
	cout << "Recall = " << recall * 100.000 / num_queries << endl;
	cout << "Time = " << time_us_per_query << " us \t QPS = " << 1e6 / (time_us_per_query) << " query/s" << endl;
}


int main(int argc, char** argv){

	// argv[1]: dataset file
	// argv[2]: queries file
	// argv[3]: ground-truth file
	// argv[4]: graph file
	// argv[5]: K, number of neighbors
	// argv[6]: maximum number of calculations

	int N, Dim, num_queries, nn_gt;

	string file_dataset(argv[1]);
	string file_queries(argv[2]);
	string file_gt(argv[3]);
	string file_graph(argv[4]);
	int K = atoi(argv[5]);
	int max_calc = atoi(argv[6]);
	const char* result_file(argv[7]);

	Matrix<float> points = read_fvecs(file_dataset, N, Dim);
	printf("base read (%d,%d) ...\n", N, Dim);
	Matrix<float>  queries = read_fvecs(file_queries, num_queries, Dim);
	printf("queries read (%d,%d)...\n", num_queries, Dim);
	vector<vector<int> > gt = read_ivecs(file_gt, num_queries, nn_gt);
	printf("groundtruth read...\n");

	freopen(result_file, "a", stdout);

	AdjList graph = read_adjlist(file_graph, points, true);
	
	if(max_calc>0)
		run_on_testset(queries, K, points, gt, graph, max_calc);
	else{
		float p=0;
		// max_calc means the maximum calculation time (i.e., the number of explored vertex)
		for(p=4.0;p>=1.99; p-=0.10){
			max_calc = (int) ((float)N / pow(10.0, p));
			run_on_testset(queries, K, points, gt, graph, max_calc);
		}
	}

	fclose(stdout);
	return 0;
}

