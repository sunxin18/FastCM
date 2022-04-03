#include "core.h"
#include "graph.h"

int k_core_size;

extern int K, b, lambda, record_b;
vector<pair<int, int>> shell_new_edges;
vector<pair<int, int>> new_edges;


void FastCMAlgorithm::output_data(const string &out_file, const double &total_t, const string algorithm, int &num_followers, const string &graph_name) {
	std::fstream ff(out_file.c_str(), std::ios::out|std::ios::app);
	ff << "Graph: " <<  graph_name << endl;
	ff << algorithm << "+: k = " << K << " b:" << record_b << " time:" << total_t << " the number of followers:" << num_followers << endl;
	ff << "The new edges identified by:" << algorithm << endl;
	for (const auto &edge: new_edges) {
		ff << edge.first << " " << edge.second << endl;
	}
	ff.close();
}

void FastCMAlgorithm::clear_everthing() {
	cc.clear();
	shell_new_edges.clear();
	solution_condidates.clear();
	followers_items.clear();
	followers_number.clear();
	k_degree.clear();
	single_component.clear();
}

void FastCMAlgorithm::print_shell() {
	map<int, vector <endpoint>>::iterator p;
	int count = 0;
	float sumcore = 0;
	int shell = 0;
	int sum = 0;
	map<int, int> mp;
	for (p = g.begin(); p != g.end(); p++)
	{
		int i = p->first;
		sumcore += coreness[i];
		if (coreness[i] >= K)
			count++;
		if (coreness[i] == K - 1) {
			shell++;
		}
	}
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
	sort(vpr.begin(), vpr.end(), cmp);
	for (vector<pair<int, int> >::iterator it = vpr.begin(); it != vpr.end(); it++) {
		cout << it->first << ":" << it->second << endl;
	}
	cout << "coreness_mean" << sumcore / n << endl;
	k_core_size = count;
	cout << "the number of" << K << "core is " << count << endl;
	cout << "the number of k-1-shell:" << shell << endl;
}


int FastCMAlgorithm::solution_selection(const int &b) {
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
			if (shell_new_edges.empty()) {
				shell_new_edges = solution_condidates[i - 1];
				followers = followers_items[i - 1];
			}
			else
			{
				shell_new_edges.insert(shell_new_edges.end(), solution_condidates[i - 1].begin(), solution_condidates[i - 1].end());
				followers.insert(followers.end(), followers_items[i - 1].begin(), followers_items[i - 1].end());
			}
			n1 = n1 - w[i - 1];
		}
	}
	if (shell_new_edges.size() == b) {
		return dp[m][n];
	} else {
		if (lambda == 1) {
			int rem = b - shell_new_edges.size();
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
						shell_new_edges.emplace_back(p);
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
						shell_new_edges.push_back(p);
						addi_1++;
						k_core_vertices.erase(remove(k_core_vertices.begin(), k_core_vertices.end(), v), k_core_vertices.end());
						rem--;
						break;
					}
				}
			}
			return dp[m][n] + 2 * addi_2 + addi_1;
		} 
		else {
			int rem = b - shell_new_edges.size();
			map<int, int> single_count;
			for (const auto v: single_component) {
				single_count[v] = lambda;
			}
			if (rem > lambda * single_component.size()) {
				for (int i = 0; i < single_component.size(); i++) {
					int u = single_component[i];
					if (single_count[u] == 0) continue;
					for (int j = i + 1; j < single_component.size(); j++) {
						if (single_count[u] == 0) break;
						int v = single_component[j];
						if (single_count[v] == 0) continue;
						if (eset.find(pair<int, int>(u, v)) == eset.end()) {
							pair<int, int>p;
							if (u < v) p = make_pair(u, v);
							else p = make_pair(v, u);
							shell_new_edges.emplace_back(p);
							single_count[u]--;
							single_count[v]--;
						}
					}
					while (single_count[u] > 0) {
						int a = rand() % (k_core_vertices.size());
						int v = k_core_vertices[a];
						pair<int, int>p;
						if (eset.find(pair<int, int>(u, v)) != eset.end()) continue;
						if (u < v) p = make_pair(u, v);
						else p = make_pair(v, u);
						if (find(shell_new_edges.begin(), shell_new_edges.end(), p) != shell_new_edges.end()) {
							continue;
						}
						single_count[u]--;
						shell_new_edges.emplace_back(p);
					}
				}
				return dp[m][n] + single_component.size();
			}
			int addi = 0;
			for (const auto &item : single_count) {
				int u = item.first, d = item.second;
				if (rem < d) {
					cannot_insert = true;
					while (rem > 0) {
						int a = rand() % (k_core_vertices.size());
						int v = k_core_vertices[a];
						pair<int, int>p;
						if (eset.find(pair<int, int>(u, v)) != eset.end()) continue;
						if (u < v) p = make_pair(u, v);
						else p = make_pair(v, u);
						if (find(shell_new_edges.begin(), shell_new_edges.end(), p) != shell_new_edges.end()) {
							continue;
						}
						rem--;
					}
					break;
				}
				if (d == 0) continue;
				while (d > 0) {
					int a = rand() % (k_core_vertices.size());
					int v = k_core_vertices[a];
					pair<int, int>p;
					if (eset.find(pair<int, int>(u, v)) != eset.end()) continue;
					if (u < v) p = make_pair(u, v);
					else p = make_pair(v, u);
					if (find(shell_new_edges.begin(), shell_new_edges.end(), p) != shell_new_edges.end()) {
						continue;
					}
					d--;
					rem--;
					shell_new_edges.emplace_back(p);
				}
				addi++;
			}
			return dp[m][n] + addi;
			for (vector<pair<int, int>>::iterator p = shell_new_edges.begin(); p != shell_new_edges.end(); p++) {
				int u = p->first;
				int v = p->second;
				//cout << u << " " << v << endl;
			}
			//cout << "the size of used budget" << shell_new_edges.size() << endl;
		}
	}
}

int FastCMAlgorithm::greedy_solution_selection(int b) {
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
			followers += followers_number[pos];
			new_edges.insert(new_edges.end(), solution_condidates[pos].begin(), solution_condidates[pos].end());
			if (C >= b) {
				return followers;
			}
		}
	}
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
				followers += 2;
				pair<int, int>p;
				p = make_pair(u, v);
				new_edges.emplace_back(p);
				rem--;
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
				pair<int, int>p;
				p = make_pair(u, v);
				new_edges.push_back(p);
				followers++;
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
	return followers;
}


void FastCMAlgorithm::parition_shell() {                    //partiton the k-l-shell
	vector<int> visit(n);
	queue<int> q;
	int index = 0;
	map<int, vector <endpoint>>::iterator p;
	for (p = g.begin(); p != g.end(); p++) {
		int i = p->first;
		if (g[i].size() == 0) continue;
		if (coreness[i] == K - lambda && visit[i] == 0) {
			q.push(i);
			visit[i] = 1;
			cc[index].push_back(i);
			while (!q.empty()) {
				int u = q.front();
				q.pop();
				int count = 0;
				for (int j = 0; j < g[u].size(); j++) {
					int v = g[u][j].u;
					if (coreness[v] >= K - lambda) {
						count++;
						if (visit[v] == 1) continue;
						if (coreness[v] == K - lambda) {
							q.push(v);
							visit[v] = 1;
							cc[index].push_back(v);
						}
					}
				}
				k_degree[u] = count;
				if (count == K - lambda) {
					is_collapse[u] = 1;
				}
			}
			sort(cc[index].begin(), cc[index].end());        //??
			index++;
		}
	}
}
void FastCMAlgorithm::partial_conversion(vector<int> &ver) {
	map<int, vector<int>> layer;
	map<int, int> layer_degree;
	map<int, int> get_layer;
	map<int, int> del;
	priority_queue<pair<int, int>, vector<pair<int, int>>, less<pair<int, int>>> heap;
	int l = 0;
	// generate layer strcuture
	while (!ver.empty()) {
		for (int i = 0; i < ver.size(); i++) {
			if (k_degree[ver[i]] <= K - lambda) {
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

	if (lambda > 1) {
		for (int la = 1; la < l; la++) {      
			unordered_map<int, int> scores;
			map<int, int> valid_degree;
			vector<int> anchored_vertices;
			vector<pair<int, int>> res;
			vector<int> weak;
			vector<int> partial_vertices;
			anchored_vertices.clear();
			for (int y = la; y < l; y++) {
				for (int t = 0; t < layer[y].size(); t++) {
					partial_vertices.push_back(layer[y][t]);
				}
			}
			map<int, int> count;
			int total = 0;
			for (auto &u: partial_vertices) {
				int c = 0;
				for (auto t : g[u]) {
					int v = t.u;
					if (coreness[v] > K - lambda || (coreness[v] == K - lambda && get_layer[v] >= l)) {
						c++;
					} 
				}
				total += K - c;
				count[u] = c;
			}

			// Compute Follower Gain 
			for (int j = 0; j < la; j++) {           
				for (int t = 0; t < layer[j].size(); t++) {
					int IncDeg = 0, valid = 0;
					int u = layer[j][t];
					for (auto t : g[u]) {
						int v = t.u;
						if (get_layer[v] >= la || coreness[v] > K - lambda) {
							valid++;
						}
						if (get_layer[v] >= la) {
							IncDeg++;
						}
					}
					int Cost = K - valid;
					int score = IncDeg - Cost;
					if (score >= 0) {          
						scores[u] = score;  
						heap.emplace(score, u);
						valid_degree[u] = valid;
					}
				}
			}

			// Anchor low-layered vertices
			unordered_map<int, int> if_anchored;
			while (!heap.empty()) {			
				while(!heap.empty()) {
					auto p = heap.top();
					if (scores[p.second] != p.first) {
						heap.pop();
					} else {
						break;
					}
				}
				if (heap.empty()) break;
				auto p = heap.top();
				heap.pop();
 				int anchored_vertex = p.second;
				anchored_vertices.push_back(anchored_vertex);
				if_anchored[anchored_vertex] = 1;
				for (auto t : g[anchored_vertex]) {
					int v = t.u;
					if (count[v] == 0) continue;
						count[v]--;
					if (count[v] == 0) {
						for (auto t : g[v]) {
							int w = t.u;
							if (coreness[w] != K - 1) {
								continue;
							}
							if (if_anchored[w] == 0 && scores.find(w) != scores.end()) {
								if (scores[w] >= 0) {
									scores[w]--;
									if (scores[w] >= 0) {
										heap.emplace(scores[w], w);
									}
								}
							}
						}
					}
				}
			}
			//sort(anchored_vertices.begin(), anchored_vertices.end());
			if (anchored_vertices.size() > 2 * b) continue;
			int rem = 0;
			// for (int t = 0; t < layer[la].size(); t++) {
			// 	int u = layer[la][t];
			// 	if (layer_degree_c[u] < K) {
			// 		rem += K - layer_degree_c[u];
			// 	}
			// }
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
			for (int i = 0; i < anchored_vertices.size(); i++) {
				int u = anchored_vertices[i];
				partial_vertices.push_back(u);
				if (valid_degree[u] >= K) continue;
				count[u] = K - valid_degree[u];
			}


			//if (total > 2 * b) continue;
			for (int i = 0; i < partial_vertices.size(); i++) {
				int u = partial_vertices[i];
				if (count[u] == 0) continue;
				for (int j = i + 1; j < partial_vertices.size(); j++) {
					if (count[u] == 0) break;
					int v = partial_vertices[j];
					if (count[v] == 0) continue;
					if (eset.find(pair<int, int>(u, v)) == eset.end()) {
						pair<int, int>p;
						if (u < v) p = make_pair(u, v);
						else p = make_pair(v, u);
						res.emplace_back(p);
						count[u]--;
						count[v]--;
					}
				}
				if (count[u] > 0) {
					weak.push_back(u);
				}
			}
			for (int i = 0; i < weak.size(); i++) {
				int u = weak[i];
				for (int j = 0; j < ver.size(); j++) {
					if (count[u] == 0) break;
					int v = ver[j];
					if (u == v) continue;
					if (eset.find(pair<int, int>(u, v)) == eset.end()) {
						pair<int, int>p;
						if (u < v) p = make_pair(u, v);
						else p = make_pair(v, u);
						if (find(res.begin(), res.end(), p) != res.end()) {
							continue;
						}
						res.emplace_back(p);
						count[u]--;
					}
				}
				while (count[u] > 0) {
					int a = rand() % (k_core_vertices.size());
					int v = k_core_vertices[a];
					if (eset.find(pair<int, int>(u, v)) == eset.end()) {
						pair<int, int>p;
						if (u < v) p = make_pair(u, v);
						else p = make_pair(v, u);
						if (find(res.begin(), res.end(), p) != res.end()) {
							continue;
						}
						count[u]--;
						res.emplace_back(p);	
					}
				}
			}
			if (res.size() > b) {
				anchored_vertices.clear();
				valid_degree.clear();
				scores.clear();
				continue;
			}
			solution_condidates.push_back(res);
			followers_number.push_back(partial_vertices.size());
			followers_items.push_back(partial_vertices);	
			return;
		}	
	}


	for (int la = 1; la < l; la++) {
		map<int, int> layer_degree_c(layer_degree);
		unordered_map<int, int> scores;
		map<int, int> valid_degree;
		vector<int> anchored_vertices;
		
		// Compute Follower Gain 
		for (int j = 0; j < la; j++) {           
			for (int t = 0; t < layer[j].size(); t++) {
				int IncDeg = 0, valid = 0;
				int u = layer[j][t];
				for (auto t : g[u]) {
					int v = t.u;
					if (get_layer[v] >= la || coreness[v] > K - 1) {
						valid++;
					}
					if (get_layer[v] == la) {
						IncDeg++;
					}
				}
				int Cost = K - valid;
				int score = IncDeg - Cost;
				if (score >= 0) {       
					scores[u] = score;  
					heap.emplace(score, u);
					valid_degree[u] = valid;
				}
			}
		}

		// Anchor low-layered vertices
		unordered_map<int, int> if_anchored;
		while (!heap.empty()) {	
			while(!heap.empty()) {
				auto p = heap.top();
				if (scores[p.second] != p.first) {
					heap.pop();
				} else {
					break;
				}
			}
			if (heap.empty()) break;
			auto p = heap.top();
			cout << p.second << " " << p.first << endl;
			heap.pop();
			int anchored_vertex = p.second;
			anchored_vertices.push_back(anchored_vertex);
			if_anchored[anchored_vertex] = 1;
			for (auto t : g[anchored_vertex]) {
				int v = t.u;
				if (get_layer[v] != la) continue;
				if (layer_degree_c[v] == K) continue;
				layer_degree_c[v]++;
				if (layer_degree_c[v] == K) {
					for (auto t : g[v]) {
						int w = t.u;
						if (coreness[w] != K - 1) {
							continue;
						}
						if (scores.find(w) != scores.end() && if_anchored[w] == 0) {
							if (scores[w] >= 0) {
								scores[w]--;
								if (scores[w] >= 0) {
									heap.emplace(scores[w], w);
								}
							}
						}
					}
				}
			}
		}
		if (anchored_vertices.size() > 2 * b) continue;
		int rem = 0;
		map<int, int> rem_degree;
		for (int t = 0; t < layer[la].size(); t++) {
			int u = layer[la][t];
			if (layer_degree_c[u] < K) {
				int require = K - layer_degree_c[u];
				rem += require;
				rem_degree[u] = require;
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
		for (int i = 0; i < anchored_vertices.size(); i++) {
			int u = anchored_vertices[i];
			if (valid_degree[u] >= K) continue;
			rem_degree[u] = K - valid_degree[u];
			rem += K - valid_degree[u];
		}
		if (rem < 2 * b) {
			int followers = anchored_vertices.size();
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
						if (u < v) p = make_pair(u, v);
						else p = make_pair(v, u);
						res.emplace_back(p);
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
				continue;
			}
			solution_condidates.push_back(res);
			followers_number.push_back(followers);
			followers_items.push_back(record_follower);
			for (vector < pair<int, int>>::iterator p = res.begin(); p != res.end(); p++) {
				int u = p->first;
				int v = p->second;
				cout << u << " " << v << endl;
			}
			 cout << "inserted size:" << res.size();
			 cout << "followers size:" << followers << endl;
			// for (auto f : record_follower) {
			// 	cout << f << " ";
			// }
			return;
		}
		else
			anchored_vertices.clear();
	}
}

void FastCMAlgorithm::complete_conversion(vector<int> ver, int flag) {
	vector<pair<int, int>> res;
	vector<int> weak;
	vector<int> collapse;            
	for (int i = 0; i < ver.size(); i++) {
		if (is_collapse[ver[i]] == 1) {
			collapse.push_back(ver[i]);           // k-1-collpase
		}
	}
	if (collapse.size() * lambda > 2 * b) {
		cout << "can't convert all the k-1-collapse vertices" << endl;
		cout << "current size" << ver.size() << endl;
		if (flag == 1) {
			partial_conversion(ver);
		}
		return;
	}
	int rem_degree;
	map<int, int> del;
	map<int, int> count;
	if (lambda > 1) collapse = ver;
	for (const auto v: collapse) {
		count[v] = max(0, K - k_degree[v]);
	}
	for (int i = 0; i < collapse.size(); i++) {
		int u = collapse[i];
		if (count[u] == 0) continue;
		for (int j = i + 1; j < collapse.size(); j++) {
			if (count[u] == 0) break;
			int v = collapse[j];
			if (count[v] == 0) continue;
			if (eset.find(pair<int, int>(u, v)) == eset.end()) {
				pair<int, int>p;
				if (u < v) p = make_pair(u, v);
				else p = make_pair(v, u);
				res.emplace_back(p);
				count[u]--;
				count[v]--;
			}
		}
		if (count[u] > 0) {
			weak.push_back(u);
		}
	}
	for (int i = 0; i < weak.size(); i++) {
		int u = weak[i];
		for (int j = 0; j < ver.size(); j++) {
			if (count[u] == 0) break;
			int v = ver[j];
			if (u == v) continue;
			if (eset.find(pair<int, int>(u, v)) == eset.end()) {
				pair<int, int>p;
				if (u < v) p = make_pair(u, v);
				else p = make_pair(v, u);
				if (find(res.begin(), res.end(), p) != res.end()) {
					continue;
				}
				res.emplace_back(p);
				count[u]--;
			}
		}
		while (count[u] > 0) {
			//cout << "error";
			int a = rand() % (k_core_vertices.size());
			int v = k_core_vertices[a];
			if (eset.find(pair<int, int>(u, v)) == eset.end()) {
				pair<int, int>p;
				if (u < v) p = make_pair(u, v);
				else p = make_pair(v, u);
				if (find(res.begin(), res.end(), p) != res.end()) {
					continue;
				}
				count[u]--;
				res.emplace_back(p);	
			}
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
int FastCMAlgorithm::FastCM_plus(int b) {
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

int FastCMAlgorithm::FastCM(int b) {
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

