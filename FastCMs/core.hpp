#include <fstream>
#include <unordered_set>
#include "graph.h"
typedef struct endpoint
{
	int u;
	int eid;
} endpoint;
map<int, vector<endpoint>> g;
int n, m, dmax;
int M;

typedef struct Edgetype
{
	int x, y,truss;
}Edgetype;
Edgetype *eg;

map<pair<int, int>, int> eset;
map<int, int> coreness;

unordered_map<int, vector<int>> cc;

extern int K, b;
map<int, int> is_collapse;
vector<vector<pair<int, int>>> solution_condidates;
vector<int> followers_number;
vector<vector<int>> followers_items;
vector<int> k_core_vertices;
map<int, int> k_degree;
vector<int> single_component;
vector<pair<int, int>> new_edges;

inline bool cmp(const pair<int, int>& p1, const pair<int, int>& p2) {
	return p1.second < p2.second;
}

void output_data(const string &out_file, const double &total_t, const string algorithm, int &num_followers, const string &graph_name) {
	std::fstream ff(out_file.c_str(), std::ios::out|std::ios::app);
	ff << "Graph: " <<  graph_name << endl;
	ff << algorithm << "+: k = " << K << " b:" << b << " time:" << total_t << " the number of followers:" << num_followers << endl;
	ff << "The new edges identified by:" << algorithm << endl;
	for (const auto &edge: new_edges) {
	ff << edge.first << " " << edge.second << endl;
	}
	ff.close();
}

void clear_everthing() {
	cc.clear();
	solution_condidates.clear();
	followers_items.clear();
	followers_number.clear();
	new_edges.clear();
	k_degree.clear();
	single_component.clear();
}

void core_decompostion()
{
	std::unique_ptr<vector<int>[]> L(new vector <int>[dmax + 1]);
	std::unique_ptr<int[]> degg(new int[n + 1]);
	std::unique_ptr<bool[]> del(new bool[n + 1]);
	int maxdeg = 0;
	int maxcore = 0;
	map<int, vector <endpoint>>::iterator p;
	for (p = g.begin(); p != g.end(); p++)
	{
		int x = p ->first;
		degg[x] = g[x].size();
		if (maxdeg < degg[x])
		{
			maxdeg = degg[x];
		}
		L[degg[x]].push_back(x);
	}
	for (int i = 0; i < n + 1; i++)
	{
		del[i] = false;
	}
	int kmax = 0;
	int md = maxdeg;
	for (int k = 0; k <= maxdeg; k++)
	{
		for (int i = 0; i < L[k].size(); i++)
		{
			int x = L[k][i];
			if (del[x] == true)
			{
				continue;
			}
			coreness[x] = degg[x];
			if (maxcore < coreness[x])
			{
				maxcore = coreness[x];
			}
			del[x] = true;
			for (int j = 0; j < g[x].size(); ++j) {
				int y = g[x][j].u;
				if (del[y] == true)
				{
					continue;
				}
				int temp = y;
				if (degg[temp] > k)
				{
					degg[temp]--;
					L[degg[temp]].push_back(temp);
				}
			}
			//printf("%d %d id=%d, truss=%d\n", x, y, eid, truss[eid]);
		}
	}
;
	printf("maxdegree=%d, maxcore=%d \n", md, maxcore);
	int count = 0;
	int count1 = 0;
	float sumcore = 0;
	int shell = 0;
	for (p = g.begin(); p != g.end(); p++)
	{
		int i = p->first;
		sumcore += coreness[i];
		if (coreness[i] >= K)
			count++;
		if (coreness[i] >= K) {
			k_core_vertices.push_back(i);
			count1++;
		}
		if (coreness[i] == K - 1)
			shell++;
	}
	int sum = 0;
	map<int, int> mp;
	for (int i = 0; i < n; i++)
	{
		if (coreness.count(i)) {
			mp[coreness[i]]++;
			sum += coreness[i];
		}
	}
	vector<pair<int, int>> vpr;
	for (map<int, int>::iterator it = mp.begin(); it != mp.end(); it++) {
		vpr.emplace_back(make_pair(it->first, it->second));
	}
	// Print the size of k-shells
	// sort(vpr.begin(), vpr.end(), cmp);
	// for (vector<pair<int, int> >::iterator it = vpr.begin(); it != vpr.end(); it++) {
	// 	cout << it->first << ":" << it->second << endl;
	// }
	cout << "coreness_mean" << sumcore / n << endl;
	cout << "the number of" << K << "core is " << count << endl;
	cout << "the number of k-1-shell:" << shell << endl;
}

int solution_selection(const int &b) {
	vector<int> w;
	vector<vector<pair<int, int>>>::iterator x;
	vector<int> followers;
	for (const auto &x: solution_condidates) {
		w.push_back(x.size());
	}
	vector<vector<int>>::iterator p;
	int m = w.size();
	int n = b;
	vector<vector<int>> dp(m + 1, vector<int>(n + 1));
	for (int i = 1; i <= m; i++) {
		for (int j = 1; j <= n; j++) {
			if (j < w[i - 1]) {
				dp[i][j] = dp[i - 1][j];
			}
			else
			{
				dp[i][j] = max(dp[i - 1][j], dp[i - 1][j - w[i - 1]] + followers_number[i - 1]);
			}
		}
	}
	int n1 = n;
	for (int i = m; i >= 1; i--) {
		if (n1 >= w[i - 1] && dp[i][n1] == dp[i - 1][n1 - w[i - 1]] + followers_number[i - 1]) {
			if (new_edges.empty()) {
				new_edges = solution_condidates[i - 1];
				followers = followers_items[i - 1];
			}
			else
			{
				new_edges.insert(new_edges.end(), solution_condidates[i - 1].begin(), solution_condidates[i - 1].end());
				followers.insert(followers.end(), followers_items[i - 1].begin(), followers_items[i - 1].end());
			}
			n1 = n1 - w[i - 1];
		}
	}
	if (new_edges.size() == b) {
		return dp[m][n];
	} else {
		int rem = b - new_edges.size();
		int addi_2 = 0;
		int addi_1 = 0;
		map<int, int> del;
		vector<int> new_single(single_component);
		for (int i = 0; i < new_single.size(); i++) {
			if (rem == 0) break;
			int u = new_single[i];
			if (del[u] == 1) continue;
			for (int j = i + 1; j < new_single.size(); j++) {
				if (i == j) continue;
				int v = new_single[j];
				if (del[v] == 1) continue;
				if (eset.find(pair<int, int>(u, v)) == eset.end()) {
					followers.push_back(u);
					followers.push_back(v);
					pair<int, int>p;
					p = make_pair(u, v);
					new_edges.emplace_back(p);
					rem--;
					addi_2++;
					del[u] = 1;
					del[v] = 1;
					break;
				}
			}
		}
		for (int i = 0; i < new_single.size(); i++) {
			if (rem == 0) break;
			if (del[new_single[i]] == 1) continue;
			int u = new_single[i];
			while(1) {
				int a = int(rand() % (k_core_vertices.size()));
				int v = k_core_vertices[a];
				if (eset.find(pair<int, int>(u, v)) == eset.end())
				{
					followers.push_back(u);
					pair<int, int>p;
					p = make_pair(u, v);
					new_edges.push_back(p);
					addi_1++;
					k_core_vertices.erase(remove(k_core_vertices.begin(), k_core_vertices.end(), v), k_core_vertices.end());
					rem--;
					break;
				}
			}
		}
		for (vector<pair<int, int>>::iterator p = new_edges.begin(); p != new_edges.end(); p++) {
			int u = p->first;
			int v = p->second;
			//cout << u << " " << v << endl;
		}
		//cout << "the size of used budget" << new_edges.size();
		return dp[m][n] + 2 * addi_2 + addi_1;	
	}
}

int greedy_solution_selection(int b) {
	int n = solution_condidates.size();
	map<double, vector<int>, greater<double>> index;
	vector<double> ratio(n);
	for (int i = 0; i < n; i++) {
		ratio[i] = followers_number[i] / solution_condidates[i].size();
		index[ratio[i]].push_back(i);
	}
	int C = 0;
	int followers = 0;
	for (auto& item : index) {
		auto vec = item.second;
		for (auto pos : vec) {
			C += solution_condidates[pos].size();
			if (C >= b) {
				return followers;
			}
			followers += followers_number[pos];
			new_edges.insert(new_edges.end(), solution_condidates[pos].begin(), solution_condidates[pos].end());
		}
	}
	int rem = b - C;
	int addi = 0;
	map<int, int> del;
	vector<int> new_single(single_component);
	for (int i = 0; i < new_single.size(); i++) {
		if (rem == 0) break;
		int u = new_single[i];
		if (del[u] == 1) continue;
		for (int j = i + 1; j < new_single.size(); j++) {
			if (i == j) continue;
			int v = new_single[j];
			if (del[v] == 1) continue;
			if (eset.find(pair<int, int>(u, v)) == eset.end()) {
				rem--;
				addi++;
				del[u] = 1;
				del[v] = 1;
				break;
			}
		}
	}
	return followers + addi * 2;
}


void parition_shell() {                    //partiton the k-l-shell
	vector<int> visit(n);
	queue<int> q;
	int index = 0;
	map<int, vector <endpoint>>::iterator p;
	for (p = g.begin(); p != g.end(); p++) {
		int i = p->first;
		if (g[i].size() == 0) continue;
		if (coreness[i] == K - 1 && visit[i] == 0) {
			q.push(i);
			visit[i] = 1;
			cc[index].push_back(i);
			while (!q.empty()) {
				int u = q.front();
				q.pop();
				int count = 0;
				for (int j = 0; j < g[u].size(); j++) {
					int v = g[u][j].u;
					if (coreness[v] >= K - 1) {
						count++;
						if (visit[v] == 1)continue;
						if (coreness[v] == K - 1) {
							q.push(v);
							visit[v] = 1;
							cc[index].push_back(v);
						}
					}
				}
				k_degree[u] = count;
				if (count < K) {
					is_collapse[u] = 1;
				}
			}
			sort(cc[index].begin(), cc[index].end());        //??
			index++;
		}
	}
}
void partial_conversion(vector<int> &ver) {
	map<int, vector<int>> layer;
	map<int, int> layer_degree;
	map<int, int> get_layer;
	map<int, int> del;
	int l = 0;

	// generate layer strcuture
	while (!ver.empty()) {
		for (int i = 0; i < ver.size(); i++) {
			if (k_degree[ver[i]] < K) {
				del[ver[i]] = 1;
			}
		}
		for (vector<int>::iterator it = ver.begin(); it != ver.end();) {  
			int v = *it;
			if (del[v] == 1) {
				get_layer[v] = l;
				layer[l].push_back(v);
				layer_degree[v] = k_degree[v];
				for (int j = 0; j < ver.size(); j++) {
					int u = ver[j];
					if (u == v) continue;
					if (eset.find(pair<int, int>(u, v)) != eset.end()) {
						if (del[u] != 1) {
							k_degree[u]--;
						}
					}
				}
				it = ver.erase(it);
			}
			else {
				it++;
			}
		}
		l++;
	}
	vector<int> anchored_vertices;
	for (int la = 1; la < l; la++) {
		map<int, int> layer_degree_c(layer_degree);
		map<int, vector<int>> act;
		map<int, int> valid_degree;
		anchored_vertices.clear();
		
		// Compute Follower Gain 
		for (int j = 0; j < la; j++) {           
			for (int t = 0; t < layer[j].size(); t++) {
				int IncDeg = 0, valid = 0;
				int u = layer[j][t];
				for (auto t : g[u]) {
					int v = t.u;
					if (get_layer[v] >= la || coreness[v] >= K) {
						valid++;
					}
					if (get_layer[v] == la) {
						IncDeg++;
					}
				}
				int Cost = K - valid;
				if (IncDeg - Cost >= 0) {          
					for (auto t : g[u]) {
						int v = t.u;
						if (coreness[v] != K - 1) continue;
						if (get_layer[v] == la) {
							act[u].push_back(v);
						}
					}
					valid_degree[u] = valid;
				}
			}
		}

		// Anchor low-layered vertices
		map<int, vector<int>>::iterator it, p;
		while (!act.empty()) {			
			int size = 0;
			bool flag = false;
			for (it = act.begin(); it != act.end(); it++) {
				if (it->second.size() > size) {
					size = it->second.size();
					p = it;
					flag = true;
				}
			}
			if (flag == false) break;
			anchored_vertices.push_back(p->first);
			vector<int>::iterator p1;
			for (const auto &u : p->second) {
				layer_degree_c[u]++;
				if (layer_degree_c[u] == K) {
					for (it = act.begin(); it != act.end(); it++) {
						if (it->first == p->first) continue;
						for (p1 = it->second.begin(); p1 != it->second.end(); ) {
							if (*p1 == u) {
								p1 = it->second.erase(p1);
								break;
							}
							else
								p1++;
						}
					}
				}
			}
			act.erase(p->first);
		}
		sort(anchored_vertices.begin(), anchored_vertices.end());
		if (anchored_vertices.size() > 2 * b) continue;
		int rem = 0;
		for (int t = 0; t < layer[la].size(); t++) {
			int u = layer[la][t];
			if (layer_degree_c[u] < K) {
				rem += K - layer_degree_c[u];
			}
		}
		if (anchored_vertices.size() != 0) {
			for (int i = 0; i < anchored_vertices.size() - 1; i++) {
				for (int j = i + 1; j < anchored_vertices.size(); j++) {
					if (eset.find(pair<int, int>(anchored_vertices[i], anchored_vertices[j])) != eset.end()) {
						valid_degree[anchored_vertices[i]]++;
						valid_degree[anchored_vertices[j]]++;
					}
				}
			}
		}
		if (rem < 2 * b) {
			int followers = anchored_vertices.size();
			map<int, int> rem_degree;
			int anchored_size = anchored_vertices.size();
			for (int i = 0; i < anchored_size; i++) {
				int u = anchored_vertices[i];
				if (valid_degree[u] >= K) continue;
				rem_degree[u] = K - valid_degree[u];
			}
			for (int t = 0; t < layer[la].size(); t++) {
				int u = layer[la][t];
				if (layer_degree_c[u] < K) {
					int require = K - layer_degree_c[u];
					rem_degree[u] = require;
					while (require > 0) {
						anchored_vertices.push_back(u);
						require--;
					}
				}
			}
			vector<pair<int, int>> res;
			map<int, int>::iterator it, it2;
			for (it = rem_degree.begin(); it != rem_degree.end(); it++) {
				if (it->second == 0) continue;
				for (it2 = it; it2 != rem_degree.end(); it2++) {
					if (*it == *it2) continue;
					if (it2->second == 0) continue;
					int u = it->first;
					int v = it2->first;
					if (eset.find(pair<int, int>(u, v)) == eset.end()) {
						pair<int, int>pp;
						if (u < v)
							pp = make_pair(u, v);
						else
							pp = make_pair(v, u);
						if (find(res.begin(), res.end(), pp) != res.end()) {
							continue;
						}
						res.emplace_back(pp);
						rem_degree[u]--;
						rem_degree[v]--;
						if(rem_degree[u] == 0) {
							break;
						}
					}
				}
			}
			for (it = rem_degree.begin(); it != rem_degree.end(); it++) {
				int u = it->first;
				while (rem_degree[u] > 0) {
					int a = rand() % (k_core_vertices.size());
					int v = k_core_vertices[a];
					if (eset.find(pair<int, int>(u, v)) == eset.end())
					{
						pair<int, int>p;
						p = make_pair(u, v);
						res.emplace_back(p);
						k_core_vertices.erase(remove(k_core_vertices.begin(), k_core_vertices.end(), v), k_core_vertices.end());
					}
					rem_degree[u]--;
				}
			}
			vector<int> record_follower;
			for (int r = la; r < l; r++) {
				 followers += layer[r].size();
				 record_follower.insert(record_follower.end(), layer[r].begin(), layer[r].end());
			}
			record_follower.insert(record_follower.end(), anchored_vertices.begin(), anchored_vertices.end());
			if (res.size() > b) {
				anchored_vertices.clear();
				continue;
			}
			solution_condidates.push_back(res);
			followers_number.push_back(followers);
			followers_items.push_back(record_follower);
			for (vector < pair<int, int>>::iterator p = res.begin(); p != res.end(); p++) {
				int u = p->first;
				int v = p->second;
				//cout << u << " " << v << endl;
			}
			//cout << "inserted size:" << res.size();
			//cout << "followers size:" << followers << endl;
			return;
		}
		else
			anchored_vertices.clear();
	}
}



void complete_conversion(vector<int> ver, int flag) {
	vector<pair<int, int>> res;
	vector<int> weak;
	vector<int> collapse;            
	for (int i = 0; i < ver.size(); i++ ) {
		if (is_collapse[ver[i]] == 1) {
			collapse.push_back(ver[i]);           // k-1-collpase
		}
	}
	if (collapse.size() > 2 * b) {
		//cout << "can't convert all the k-1-collapse vertices" << endl;
		if (flag == 1) {
			partial_conversion(ver);
		}
		return;
	}
	int rem_degree;
	map<int, int> del;
	for (int i = 0; i < collapse.size(); i++) {
		int u = collapse[i];
		if (del[u] == 1) continue;
		for (int j = i + 1; j < collapse.size(); j++) {
			int v = collapse[j];
			if (del[v] == 1) continue;
			if (eset.find(pair<int, int>(u, v)) == eset.end()) {
				pair<int, int>p;
				p = make_pair(u, v);
				res.emplace_back(p);
				del[u] = 1;
				del[v] = 1;
				break;
			}
		}
		if (del[u] == 0) {
			weak.push_back(u);
		}
	}
	for (int i = 0; i < weak.size(); i++) {
		int u = weak[i];
		int flag = 0;
		for (int j = 0; j < ver.size(); j++) {
			int v = ver[j];
			if (u == v) continue;
			if (eset.find(pair<int, int>(u, v)) == eset.end()) {
				pair<int, int>p;
				p = make_pair(u, v);
				res.emplace_back(p);
				flag = 1;
				break;
			}
		}
		if (flag == 0) {
			//cout << "error";
			int a = rand() % (k_core_vertices.size());
			int v = k_core_vertices[a];
			pair<int, int>p;
			p = make_pair(u, v);
			res.emplace_back(p);
		}
	}
	if (res.size() == 0) return;
	solution_condidates.push_back(res);
	followers_number.push_back(ver.size());
	followers_items.push_back(ver);
	for (vector < pair<int, int>>::iterator p = res.begin(); p != res.end(); p++) {
		int u = p->first;
		int v = p->second;
		//cout << u << " " << v << endl;
	}
	//cout << "inserted size:" << res.size();
	//cout << "followers size:" << ver.size() << endl;
}
int FastCM_plus(int b) {
	parition_shell();
	int co = 0;
	for (int nn = 0; nn < cc.size(); nn++) {
		co += cc[nn].size();
		if (cc[nn].size() == 1) {
			single_component.push_back(cc[nn][0]);
			continue;
		}
		complete_conversion(cc[nn], 1);
	}
	//cout << "The number of k-1 shell:" << co << endl;
	//cout << "The size of components:" << cc.size() << endl;
	return solution_selection(b);
}

int FastCM(int b) {
	parition_shell();
	int co = 0;
	for (int nn = 0; nn < cc.size(); nn++) {
		co += cc[nn].size();
		if (cc[nn].size() == 1) {
			single_component.push_back(cc[nn][0]);
			continue;
		}
		complete_conversion(cc[nn], 0);
	}
	//cout << "The number of k-1 shell:" << co << endl;
	//cout << "The size of components:" << cc.size() << endl;
	return greedy_solution_selection(b);
}

