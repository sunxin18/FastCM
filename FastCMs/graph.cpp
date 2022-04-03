#include "graph.h"
// Read Input Graph Dataset  and store in g.
extern int K;
map<pair<int, int>, int> eset;
map<int, int> coreness;
map<int, vector<endpoint>> g;
vector<int> k_core_vertices;
int n, m, dmax, M;

void Readin_Graph(const char *str)
{
	FILE *fin = fopen(str, "r");
	if (fin == nullptr) {
		printf("Fail to open file\n");
		exit(0);
	}
	int x, y, z;
	
	n = m = 0;
	dmax = 0;
	while (fscanf(fin, "%d %d", &x, &y) != EOF)
	{
		n = max(n, x);
		n = max(n, y);
	}
	n++;
	fclose(fin);
	fin = fopen(str, "r");
	endpoint a, b;
	int eid = 0;

	while (fscanf(fin, "%d %d", &x, &y)!=EOF)
	{
		if (x == y) continue;
		if (y < x) swap(x, y);
		if (x < y)
		{
			if (eset.find(pair<int, int>(x, y)) == eset.end())
			{
				eset[pair<int, int>(x, y)] = eid;
				eset[pair<int, int>(y, x)] = eid;
				a.u = y;
				a.eid = eid;
				g[x].push_back(a);
				b.u = x;
				b.eid = eid;
				g[y].push_back(b);
				eid++;
				m += 2;
			}
		}
		n = max(y, n);
	}
	map<int, vector <endpoint>>:: iterator p;
	for (p = g.begin(); p != g.end(); p++)
	{
		int i = p->first;
		dmax = max(dmax, (int)g[i].size());
	}
	/*for (int i = 0; i < n; i++)
	{
		dmax = max(dmax, (int)g[i].size());
	}*/
	cout << "read over" << endl;
	printf("n = %d, m = %d, dmax = %d, eid=%d\n", g.size(), eid, dmax, eid);
	fclose(fin);
}


void core_decompostion() {
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
	for (p = g.begin(); p != g.end(); p++)
	{
		int i = p->first;
		if (coreness[i] >= K) {
			k_core_vertices.push_back(i);
		}
	}
}


