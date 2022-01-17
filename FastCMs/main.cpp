#include "graph_IO.h"
#include <iostream>
#include <chrono>
using namespace std;

int K, b, record_b;
int lambda = 1;
int main(int argc, char* argv[]){
	double timeuse = 0;
	int k;
	//cout << "Please input the coreness k and budget b" << endl;
	K = atoi(argv[2]);
	b = atoi(argv[3]);
	record_b = b;
	//cin >> k >> b;
	//Readin_Graph("facebook_combined.txt");
	//Readin_Graph("socfb-konect.txt");
	//Readin_Graph("sample_web-Stanford.txt");
	//Readin_Graph("quchong_route_network.txt");
	//Readin_Graph("sample_soc-dogster.txt");
	//Readin_Graph("Email-Enron.txt");
	//Readin_Graph("as-skitter.txt");
	//Readin_Graph("socfb-konect.txt");
	//Readin_Graph("Gowalla_edges.txt");
	//Readin_Graph("Russia_route.txt");
	//Readin_Graph("Brightkite_edges.txt");
	//Readin_Graph("facebook_combined.txt");
	//Readin_Graph("socfb-konect.txt");
	//truss_decompostion_light(13507);
	//Readin_Graph("baidu.txt");
	//Readin_Graph("sample_twitter_combined.txt");
	string out_file = "results.txt";
	string graph_name = argv[1];
	int num_followers;
	Readin_Graph(graph_name.c_str());
	core_decompostion();

	// *****  FastCM procedure:
	auto start = std::chrono::high_resolution_clock::now();
	num_followers = FastCM(b);
	auto end = std::chrono::high_resolution_clock::now();
	double total_t = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;
	output_data(out_file, total_t, string("FastCM"), num_followers, graph_name);
	clear_everthing();	
	new_edges.clear();
	// *****  FastCM+ procedure: 
	start = std::chrono::high_resolution_clock::now();
	num_followers = 0;
	cout << "bbb:" << b << endl;
	
	while (b > 0) {
		clear_everthing();
		if (cannot_insert) break;
		num_followers += FastCM_plus(b);
		new_edges.insert(new_edges.end(), shell_new_edges.begin(), shell_new_edges.end());
		b -= shell_new_edges.size();
		//cout << "****" << new_edges.size() << endl;
		lambda++;
		//insert_new_edges();
	}
	end = std::chrono::high_resolution_clock::now();
	total_t = double(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()) / 1000;

	output_data(out_file, total_t, string("FastCM+"), num_followers, graph_name);


	delete [] eg; 
	return 0;
}

