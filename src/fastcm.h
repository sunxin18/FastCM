#include <fstream>
#include <unordered_set>
#include <queue> 
#include <set>
#include <unordered_map>
#include <vector>
#include <map>
using namespace std;

inline bool cmp(const pair<int, int>& p1, const pair<int, int>& p2) {
	return p1.second < p2.second;
}

extern vector<pair<int, int>> shell_new_edges;
extern vector<pair<int, int>> new_edges;

class FastCMAlgorithm {
public:
    FastCMAlgorithm() {
        cannot_insert = false;
    }
    void output_data(const string &out_file, const double &total_t, const string algorithm, int &num_followers, const string &graph_name);
    void clear_everthing();
    void print_shell();
    int solution_selection(const int &b);
    int greedy_solution_selection(int b);
    void parition_shell();
    void partial_conversion(vector<int> &ver);
    void complete_conversion(vector<int> ver, int flag);
    int FastCM_plus(int b);
    int FastCM(int b);
    bool cannot_insert;
private:
    map<int, int> is_collapse;
    vector<vector<pair<int, int>>> solution_condidates;
    vector<int> followers_number;
    vector<vector<int>> followers_items;
    map<int, int> k_degree;
    vector<int> single_component;
    unordered_map<int, vector<int>> cc;
};