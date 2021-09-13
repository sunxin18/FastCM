#include "core.hpp"
// Read Input Graph Dataset  and store in g.
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
	Edgetype e;
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
	eg = new Edgetype [eid];
	map<int, vector <endpoint>>:: iterator p;
	for (p = g.begin(); p != g.end(); p++) {
		for (int j = 0; j < p->second.size(); j++)
		{
			int i = p->first;
			int y = g[i][j].u;
			int z = g[i][j].eid;
			if (x < y)
			{
				e.x = x;
				e.y = y;
				eg[z] = e;
				//printf("%d %d id=%d\n", x, y, z);
			}

		}
	}
	
	M = eid;
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


