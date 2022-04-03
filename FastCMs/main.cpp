#include "graph.h"
#include "core.h"
#include <iostream>
#include <chrono>
using namespace std;

int K, b, record_b;
int lambda = 1;
int main(int argc, char* argv[]){
	string graph_name = argv[1];
	K = atoi(argv[2]);
	b = atoi(argv[3]);
	record_b = b;
	string out_file = "results.txt";
	int num_followers;
	Readin_Graph(graph_name.c_str());
	std::unique_ptr<FastCMAlgorithm> algo(new FastCMAlgorithm());
	core_decompostion();

	// *****  FastCM procedure:
	auto start = std::chrono::high_resolution_clock::now();
	num_followers = algo->FastCM(b);
	auto end = std::chrono::high_resolution_clock::now();
	double total_t = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
	algo->output_data(out_file, total_t, string("FastCM"), num_followers, graph_name);
	algo->clear_everthing();	
	new_edges.clear();
	// *****  FastCM+ procedure: 
	start = std::chrono::high_resolution_clock::now();
	num_followers = 0;
	cout << "bbb:" << b << endl;
	
	while (b > 0) {
		algo->clear_everthing();
		if (algo->cannot_insert) break;
		num_followers += algo->FastCM_plus(b);
		new_edges.insert(new_edges.end(), shell_new_edges.begin(), shell_new_edges.end());
		b -= shell_new_edges.size();
		lambda++;
	}
	end = std::chrono::high_resolution_clock::now();
	total_t = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

	algo->output_data(out_file, total_t, string("FastCM+"), num_followers, graph_name);

	return 0;
}

