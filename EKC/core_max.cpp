#include <bits/stdc++.h>
#define rep(i, j, k) for (int i = j; i < k; ++i)
#define repr(i, j, k) for (int i = j; i <= k; ++i)
#define dep(i, j, k) for (int i = j; i > k; --i)
#define depr(i, j, k) for (int i = j; i >= k; --i)
#define getnum(x) scanf("%d", &x);
#define _print(x,y) cout << #x << " : " << y << endl;
#define mp make_pair
#define pb push_back
#define fi first
#define se second
using namespace std;
typedef long long LL;
typedef pair<int, int> PII;
typedef map<int, int> MII;
typedef vector<int> VI;
typedef vector<VI> VVI;
//timer
map<string, double> map_timer_;
inline void resetTimer(string str)
{
	map_timer_[str] = (double)clock() / CLOCKS_PER_SEC;
}
inline double getTimer(string str)
{
	return (double)clock() / CLOCKS_PER_SEC - map_timer_[str];
}
// --------------------
//#define DEBUG
// variable
map<int, int> vert_set_index_;
VI deg_;
VVI vert_set_;

map<int, int> shell_vert_set_index_;
VI shell_vert_set_;
set<int> core_set_;
VVI shell_down_nbr_, shell_all_nbr_, core_nbr_;
VI shell_ub_deg_, shell_status_, shell_layer_num_;

set<PII> exist_edge_;
set<PII> pruned_edge_;
int pruned_count = 0, pruned_to_core = 0;

inline bool isExist(int u, int v) {
	if (u > v) swap(u, v);
	if (exist_edge_.find({u, v}) != exist_edge_.end()) {
		return true;
	}
	return false;
}

inline void addEdge(int u, int v) {
	if (u > v) swap(u, v);
	exist_edge_.insert({u, v});
}

bool compareShellDeg(PII e1, PII e2) {
	const int& idx_u1 = shell_vert_set_index_[e1.first];
	const int& idx_v1 = shell_vert_set_index_[e1.second];
	const int& idx_u2 = shell_vert_set_index_[e2.first];
	const int& idx_v2 = shell_vert_set_index_[e2.second];
	return shell_all_nbr_[idx_u1].size()+shell_all_nbr_[idx_v1].size()
			> shell_all_nbr_[idx_u2].size()+shell_all_nbr_[idx_v2].size();
}

//bool compareLayerSum()
// status
int k_, b_, edge_num_, node_num_;
int max_deg_;

// constant
string FILE_IN("sample_soc-dogster.txt");
//Brightkite_edges
//twitter_combined
// for core maintenance
VI follower_num_;
vector<PII> anchor_edges_;

vector<PII> candidate_edge_;
int kcore_num_, km1_shell_num_;

VI candidate_node_;

VI getKmxCoreTag(const VVI& vert_set, const MII& vert_set_index, VI& vert_deg, int x) {
	int N = vert_set.size();
	VI vert_tag(N, 0);
	
	queue<int> to_delete;
	rep(i, 0, N) {
		if (vert_deg[i] < k_-x) {
			vert_tag[i] = -1;
			to_delete.push(i);
		}
	}
	while (!to_delete.empty()) {
		const int& cur_idx = to_delete.front();
		to_delete.pop();
		rep(j, 1, vert_set[cur_idx].size()) {
			const int& nbr_id = vert_set[cur_idx][j];
			const int& nbr_idx = vert_set_index.at(nbr_id);
			if (vert_deg[nbr_idx] >= k_-x) {
				--vert_deg[nbr_idx];
				if (vert_deg[nbr_idx] < k_-x) {
					to_delete.push(nbr_idx);
					vert_tag[nbr_idx] = -1;
				}
			}
		}
	}
	return vert_tag;
}

void reformVertSet(VVI& vert_set, MII& vert_set_index, VI& vert_deg, const VI& vert_tag) {
	VVI vert_set_temp;
	rep(idx, 0, vert_set.size()) {
		if (!vert_tag[idx]) {
			const int &cur_id = vert_set[idx][0];
			VI vert_set_insert;
			vert_set_insert.push_back(cur_id);
			rep(j, 1, vert_set[idx].size()) {
				int nbr_id = vert_set[idx][j];
				int nbr_idx = vert_set_index[nbr_id];
				if (!vert_tag[nbr_idx]) {
					vert_set_insert.push_back(nbr_id);
					if (!isExist(cur_id, nbr_id)) {
						addEdge(cur_id, nbr_id);
					}
				}
			}
			vert_set_temp.push_back(vert_set_insert);
		}
	}
	vert_set.swap(vert_set_temp);
	
	// reform vert_set_index, vert_deg
	const int &N = vert_set.size();
	vert_set_index.clear();
	rep(i, 0, N) {
		vert_set_index[vert_set[i][0]] = i;
		vert_deg[i] = vert_set[i].size() - 1;
	}
}

/*	shell_all_nbr_, shell_vert_set_index_
	shell_ub_deg_
	shell_status_ : -1 deleted, 0 unexplored, 1 survived
*/
void shrink(int u_id) {
	const int &u_idx = shell_vert_set_index_[u_id];
	shell_status_[u_idx] = -1;
	for (auto nbr_vert : shell_all_nbr_[u_idx]) {
		const int &nbr_idx = shell_vert_set_index_[nbr_vert];
		if (shell_status_[nbr_idx] == 1) {
			--shell_ub_deg_[nbr_idx];
			if (shell_ub_deg_[nbr_idx] < k_) {
				shrink(nbr_vert);
			}
		}
	}
}

/*	[use]
	core_nbr_, shell_all_nbr_, shell_down_nbr_, shell_vert_set_, shell_vert_set_index_, shell_layer_num_
	[init]
	shell_ub_deg_
	shell_status_: -1 deleted, 0 unexplored, 1 survived
	[parameter]
	PII edge: first.layer_number < second.layer_number
*/
int findFollower(PII edge) {
	// local variable
	priority_queue<PII, vector<PII>, greater<PII> > que;	// layer_num, id
	set<int> shell_in_que;
	// variable
	const int &N = shell_vert_set_.size();
	int u = edge.first, v = edge.second;
	int u_idx = shell_vert_set_index_[u], v_idx = shell_vert_set_index_[v];
	int l_u = shell_layer_num_[u_idx], l_v = shell_layer_num_[v_idx];
	// init shell_ub_deg_
	shell_ub_deg_.assign(N, 0);
	// init shell_status_
	shell_status_.assign(N, 0);
	rep(idx, 0, N) {
		int id = shell_vert_set_[idx];
		int l_id = shell_layer_num_[idx];
		if (l_id < l_u) {
			shell_status_[idx] = -1;
		} else if (l_id == l_u && id != u && id != v) {
			shell_status_[idx] = -1;
		}
	}
	// priority_queue
	que.push({l_u, u});
	shell_in_que.insert(u);
	if (l_v == l_u) {
		que.push({l_v, v});
		shell_in_que.insert(v);
	}
	// algorithm
	while (!que.empty()) {
		PII cur_node = que.top();
		que.pop();
		const int &cur_l = cur_node.first;
		const int &cur_id = cur_node.second;
		const int &cur_idx = shell_vert_set_index_[cur_id];
		// computer shell_ub_deg_
		int ub_deg = core_nbr_[cur_idx].size();
		for (auto nbr_vert : shell_all_nbr_[cur_idx]) {
			const int &nbr_idx = shell_vert_set_index_[nbr_vert];
			const int &nbr_l = shell_layer_num_[nbr_idx];
			if (shell_status_[nbr_idx] == 1 || (!shell_status_[nbr_idx] && (nbr_l > cur_l || shell_in_que.find(nbr_vert)!= shell_in_que.end()))) {
				++ub_deg;
			}
		}
		if (ub_deg >= k_) {
			shell_ub_deg_[cur_idx] = ub_deg;
			shell_status_[cur_idx] = 1;
			for (auto nbr_vert : shell_down_nbr_[cur_idx]) {
				if (shell_in_que.find(nbr_vert) == shell_in_que.end()) {
					que.push({ shell_layer_num_[shell_vert_set_index_[nbr_vert]], nbr_vert });
					shell_in_que.insert(nbr_vert);
				}
			}
		} else {
			// [to-do] check whether anchor_node nbr
			#ifdef EARLY_TERM_IN_EDGE_NBR
			/*for (auto nbr_vert : shell_all_nbr_[u_idx]) {
				
			}*/
			#endif
			shrink(cur_id);
		}
	}
	
	pruned_to_core = 0;
	pruned_edge_.clear();
	VI follower;
	int ret = 0;
	rep(idx, 0, shell_status_.size()) {
		if (shell_status_[idx] == 1) {
			++ret;
			follower.push_back(shell_vert_set_[idx]);

			pruned_to_core += core_set_.size() - core_nbr_[idx].size();
		}
	}
	rep(i, 0, follower.size()) {
		rep(j, i+1, follower.size()) {
			int u = follower[i], v = follower[j];
			if (u > v) {
				swap(u, v);
			}
			pruned_edge_.insert({u, v});
		}
	}
	return ret;
}

void getCandidateEdge()
{
	candidate_edge_.clear();
	rep(u_idx, 0, shell_vert_set_.size()) {
		const int& u = shell_vert_set_[u_idx];
		int d_u = shell_all_nbr_[u_idx].size() + core_nbr_[u_idx].size();
		int l_u = shell_layer_num_[u_idx];
		
		// shell to core
		if (d_u == k_-1) {
			if (core_set_.size() > core_nbr_[u_idx].size()) {
				candidate_node_.pb(u);
			}
		}
		
		// shell to shell
		rep(v_idx, u_idx+1, shell_vert_set_.size()) {
			const int& v = shell_vert_set_[v_idx];
			int d_v = shell_all_nbr_[v_idx].size() + core_nbr_[v_idx].size();
			int l_v = shell_layer_num_[v_idx];
			
			if (isExist(u, v)) continue;
			if (l_u == l_v) {
				if (d_u < k_-1 || d_v < k_-1) {
					continue;
				}
				candidate_edge_.pb({u, v});
			} else {
				int outer = l_u < l_v ? u : v;
				int inner = u + v - outer;
				int outer_deg = (outer == u ? d_u : d_v);
				if (outer_deg < k_-1) {
					continue;
				}
				candidate_edge_.pb({outer, inner});
			}
		}
	}
	sort(candidate_edge_.begin(), candidate_edge_.end(), compareShellDeg);
	_print(candidate_edge, candidate_edge_.size());
	_print(candidate_node, candidate_node_.size());
}

/*	core_nbr_, shell_all_nbr_, shell_down_nbr_, shell_vert_set_, shell_vert_set_index_, shell_layer_num_
	shell_ub_deg_
*/
PII anchorOneEdge()
{
	getCandidateEdge();
	cout << "the size of candidate edges" << candidate_edge_.size(); 
	if (candidate_edge_.size() == 0) {
		puts("No candidate!");
		return mp(0, 0);
	}
	int max_follower_num = 0, max_pruned_to_core = 0;
	PII anchor = candidate_edge_[0];

	int _cnt = 0;
	resetTimer("batch round");
	for (auto edge : candidate_edge_) {
		const int& u = edge.first;
		const int& v = edge.second;
		const int& u_idx = shell_vert_set_index_[u];
		const int& v_idx = shell_vert_set_index_[v];
		
		// check follower pruning
		int u1 = u, v1 = v;
		if (u1 > v1) {
			swap(u1, v1);
		}
		if (pruned_edge_.find({u1, v1}) != pruned_edge_.end()){
			++pruned_count;
			continue;
		}

		shell_all_nbr_[u_idx].push_back(v);
		shell_all_nbr_[v_idx].push_back(u);
		if (shell_layer_num_[u_idx] != shell_layer_num_[v_idx]) {
			shell_down_nbr_[u_idx].push_back(v);
		}
		
		int cur_follower_num = findFollower(edge);
		if (cur_follower_num > max_follower_num) {
			max_follower_num = cur_follower_num;
			anchor = edge;
			max_pruned_to_core = pruned_to_core;
		}
		
		shell_all_nbr_[u_idx].pop_back();
		shell_all_nbr_[v_idx].pop_back();
		if (shell_layer_num_[u_idx] != shell_layer_num_[v_idx]) {
			shell_down_nbr_[u_idx].pop_back();
		}
		++_cnt;
		if (_cnt % 100000 == 0) {
			printf("candidate time: %lf\n", getTimer("batch round"));
			resetTimer("batch round");
		}
	}
	
	int u = anchor.first;
	int v = anchor.second;
	const int &u_idx = vert_set_index_[u];
	const int &v_idx = vert_set_index_[v];
	vert_set_[u_idx].push_back(v);
	vert_set_[v_idx].push_back(u);
	++deg_[u_idx];
	++deg_[v_idx];
	addEdge(u, v);
printf("shell_vert: %d, max_follower_num: %d, pruned_count: %d, pruned_to_core: %d\n", shell_vert_set_.size(), max_follower_num, pruned_count, max_pruned_to_core);
	follower_num_.pb(max_follower_num);
	return anchor;
}

void constrcutKm1Shell()
{
	shell_vert_set_.clear();
	shell_vert_set_index_.clear();
	core_nbr_.clear();
	shell_down_nbr_.clear();
	shell_all_nbr_.clear();
	shell_layer_num_.clear();
	shell_ub_deg_.clear();
	
	// compute km1core	
	VI shell_deg(deg_.begin(), deg_.end());
	VI shell_tag = getKmxCoreTag(vert_set_, vert_set_index_, shell_deg, 1);
	
	// compute k-1 shell
	queue<int> to_delete;
	int cnt = 1;
	rep(i, 0, shell_tag.size()) {
		if (!shell_tag[i] && shell_deg[i] < k_) {
			to_delete.push(i);
			shell_tag[i] = cnt;
		}
	}
	do {
		queue<int> tmp_que;
		++cnt;
		while (!to_delete.empty()) {
			const int& idx = to_delete.front();
			to_delete.pop();
			rep(i, 1, vert_set_[idx].size()) {
				const int &v = vert_set_[idx][i];
				const int &v_idx = vert_set_index_[v];
				if (!shell_tag[v_idx]) {
					--shell_deg[v_idx];
					if (shell_deg[v_idx] < k_) {
						tmp_que.push(v_idx);
						shell_tag[v_idx] = cnt;
					}
				}
			}
		}
		to_delete.swap(tmp_que);
	} while (!to_delete.empty());
	
	// construct shell_sturcture(shell_layer_num_, shell_vert_set_, core_set_, shell_vert_set_index_)
	rep(i, 0, shell_tag.size()) {
		if (shell_tag[i] > 0) {
			shell_vert_set_.push_back(vert_set_[i][0]);
		} else if (shell_tag[i] == 0) {
			core_set_.insert(vert_set_[i][0]);
		}
	}
	const int &N = shell_vert_set_.size();
	rep(i, 0, N) {
		shell_vert_set_index_[shell_vert_set_[i]] = i;
	}
	shell_layer_num_.resize(N);
	rep(i, 0, N) {
		const int& u_id = shell_vert_set_[i];
		const int& u_idx = vert_set_index_[u_id];
		shell_layer_num_[i] = shell_tag[u_idx];		
		const int &l_u = shell_layer_num_[i];
		const VI& nbr = vert_set_[u_idx];
		
		VI all_nbr_insert, one_nbr_insert, up_nbr_insert, down_nbr_insert, core_nbr_insert;
		
		rep(j, 1, nbr.size()) {
			const int& v_id = nbr[j];
			const int& v_idx = vert_set_index_[v_id];
			int l_v = shell_tag[v_idx];
			if (l_v < 0) {
				continue;
			}
			if (l_v == 0) {
				core_nbr_insert.push_back(v_id);
			} else {
				all_nbr_insert.push_back(v_id);
				if (l_v == l_u) {
//					one_layer_nbr_insert.push_back(v_id);
				} else if (l_v < l_u) {
//					up_layer_nbr_insert.push_back(v_id);
				} else {
					down_nbr_insert.push_back(v_id);
				}
			}
		}
		shell_all_nbr_.push_back(all_nbr_insert);
		core_nbr_.push_back(core_nbr_insert);
//		one_layer_nbr_.push_back(one_layer_nbr_insert);
//		up_layer_nbr_.push_back(up_layer_nbr_insert);
		shell_down_nbr_.push_back(down_nbr_insert);
	}
	#ifdef DEBUG
	int _l_add = 1;
	bool flag;
	do {
		cout << "layer :" << _l_add<<endl;
		flag = false;
		rep(idx, 0, shell_vert_set_.size()) {
			if (shell_layer_num_[idx] == _l_add) {
				flag = true;
				cout << "id: "<<shell_vert_set_[idx] << endl;
				for (auto nbr_vert : shell_down_nbr_[idx]) {
					cout << nbr_vert << ' ';
				}
				cout << endl;
			}
		}
		++_l_add;
	} while (flag);
	#endif
}

void printShell() {
	int l_num = 0, cnt = 0;
	do {
		cnt = 0;
		++l_num;
		rep(i, 0, shell_vert_set_.size()) {
			if (shell_layer_num_[i] == l_num) {
				cout << shell_vert_set_[i] << ' ';
				++cnt;
			}
		}
		cout << endl;
	} while (cnt != 0);
}

void algorithm()
{
	// anchoring
	rep(i, 0, b_) {
		// data struction, get km1_tag, km1_deg
		resetTimer("one round");
		constrcutKm1Shell();
		cout << shell_vert_set_.size() << ' ' << core_set_.size()<< endl;

		// anchoring algorithm
		PII anchor_edge = anchorOneEdge();
		anchor_edges_.push_back(anchor_edge);
		printf("anchor one round time: %lf\n", getTimer("one round"));
	}
}

/*	[init]
	VI deg_
	VVI vert_set_
	map<int, int> vert_set_index_
*/
void readData()
{
	VI vert_set_insert;

	scanf("%d%d", &k_, &b_);
	FILE *fin = fopen(FILE_IN.c_str(), "r");
	if (fin == NULL) {
		puts("Error: Can not open dataset");
		return;
	}
	int u, v;
	while (~fscanf(fin, "%d%d", &u, &v)) {
		if (vert_set_insert.empty()) {
			vert_set_insert.push_back(u);
			vert_set_insert.push_back(v);
		} else if (vert_set_insert[0] == u) {
			vert_set_insert.push_back(v);
		} else {
			vert_set_.push_back(vert_set_insert);
			vert_set_insert.clear();
			vert_set_insert.push_back(u);
			vert_set_insert.push_back(v);
		}
	}
	if (!vert_set_insert.empty()) {
		vert_set_.push_back(vert_set_insert);
	}
	// init vert_set_index_, deg_
	const int &N = vert_set_.size();
	deg_.resize(N);
	rep(i, 0, N) {
		vert_set_index_[vert_set_[i][0]] = i;
		deg_[i] = vert_set_[i].size() - 1;
	}

	fclose(fin);
}

/*	[change]
	VI deg_
	VVI vert_set_
	map<int, int> vert_set_index_
*/
void init()
{
	readData();
	// compute km1core	
	VI kmb_tag = getKmxCoreTag(vert_set_, vert_set_index_, deg_, 1);
	puts("kmb");
	reformVertSet(vert_set_, vert_set_index_, deg_, kmb_tag);
_print(init vert_set size, vert_set_.size())
}

void outputResult()
{
	int total_fol = 0;
	int cnt = 0;
	for (auto edge : anchor_edges_) {
		if (edge.first == 0 && edge.second == 0)
			break; 
		if (follower_num_[cnt] == 0) break;
		printf("(%d %d)\n", edge.first, edge.second);
		cout << follower_num_[cnt++] << endl;
		total_fol+=follower_num_[cnt-1];
	}
	cout << "the used budget" << cnt << endl;
	cout << total_fol << endl;
	std::string str1 = "vary_k_youtube20.txt";
	std::fstream ff(str1.c_str(), std::ios::out|ios::app);
	ff << "k: " << k_ << "followers: " << total_fol << "time " <<  getTimer("algo");
	ff.close();
}

int main(int argc, char *argv[])
{
	//k_ = atoi(argv[1]);
	//b_ = atoi(argv[2]);
	printf("xxxx"); 
	resetTimer("init");
	init();
	printf("init time: %lf\n", getTimer("init"));
	
	resetTimer("algo");
	algorithm();
	printf("algo time: %lf\n", getTimer("algo"));

	outputResult();
	return 0;
}
