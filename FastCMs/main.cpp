#include "graph_IO.h"
#include <iostream>
#include <chrono>
using namespace std;

int K, b;

int main(int argc, char* argv[]){
	//Tsize=5;
	double timeuse = 0;
	//ofstream sout(argv[3],ios::out|ios::app|ios::binary);
	//int numProcs = omp_get_num_procs();
	//cout << "omp_get_num_procs() = " << numProcs << endl;
	int k;
	cout << "Please input the coreness k and budget b" << endl;
	cin >> k >> b;
	K = k;
	//Readin_Graph(argv[1]);
	//Readin_Graph("facebook_combined.txt");
	//Readin_Graph("Gowalla_edges.txt");
	//Readin_Graph("twitter_combined.txt");
	//Readin_Graph("quchong_route_network.txt");
	//Readin_Graph("live.txt");
	//Readin_Graph("Email-Enron.txt");
	//Readin_Graph("as-skitter.txt");
	//Readin_Graph("socfb-konect.txt");
	//Readin_Graph("Brightkite_edges.txt");
	//Readin_Graph("Russia_route.txt");
	//Readin_Graph("Australia_route.txt");
	//Readin_Graph("facebook_combined.txt");
	//Readin_Graph("sample_web-Stanford.txt");
	//truss_decompostion_light(13507);
	//Readin_Graph("baidu.txt");
	string graph_name = argv[1];
	string out_file = "results.txt";
	int num_followers;
	Readin_Graph(graph_name.c_str());
	core_decompostion();

	// *****  FastCM+ procedure: 
	auto start = std::chrono::high_resolution_clock::now();
	num_followers = FastCM_plus(b);
	auto end = std::chrono::high_resolution_clock::now();
	double total_t = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
	output_data(out_file, total_t, string("FastCM+"), num_followers, graph_name);
	clear_everthing();

	// *****  FastCM procedure:
	start = std::chrono::high_resolution_clock::now();
	num_followers = FastCM(b);
	end = std::chrono::high_resolution_clock::now();
	total_t = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
	output_data(out_file, total_t, string("FastCM"), num_followers, graph_name);

	delete [] eg; 
	return 0;
}

