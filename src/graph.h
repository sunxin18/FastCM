#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <time.h>
#include <iostream>
#include <string>
#include <memory>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <map>
#include <vector>
using namespace std;


// graph 
struct endpoint
{
	int u;
	int eid;
};

extern map<pair<int, int>, int> eset;
extern map<int, int> coreness;
extern map<int, vector<endpoint>> g;
extern vector<int> k_core_vertices;
extern int n, m, dmax;

extern void Readin_Graph(const char *str);
extern void core_decompostion();

#endif  /* _GRAPH_H_ */