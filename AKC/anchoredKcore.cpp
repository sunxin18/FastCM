/*
Copyright (c) of 2016 by Fan Zhang <fanzhang.cs@gmail.com>
*/
#pragma warning(disable:4996)
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<vector>
#include<time.h>
#include<iostream>
#include<algorithm>
#include<cstdlib>
#include<ctime>

using namespace std;

long inputK, inputB;
double runTime;
double algStartTime, anchorAndfollowerNumber;
string infile, outfile;
long anchorNumber, kcoreNumber, maxKm1coreNumber = 0, maxKm1NeighborNumber = 0, kmbNeighborNumber = 0, kmbcoreNumber;
long kcoreNum, km1ShellSize, km1ShellCanSize, potentialFollowersIn1 = 0, potentialFollowersIn4 = 0, potentialFollowersIn6 = 0;
unsigned long maxLayerNum;
long anchorTries = 0, anchorTriesFE = 0, anchorTriesUB = 0, anchorTriesRA = 0, noFollowerTries = 0, totalTriedFollowers = 0, lastLayersNumber;
long long anchoringTries = 0, activatedCandidates = 0;
long numberAfterLastLayers;
vector<long> subgraphTag, subgraphDegree, verDegree, anchorTag, followTag;
vector<long> kmbTag, kmbDegree, kmbIndex, km1Tag, km1Degree, kTag, kDegree, km1DeletTag, kDeletTag, kcoreDeletTag, km1ShellTag, km1ShellCandidates, km1ShellCandidateIds, lastCandidates, km1ShellLayerRecord, triedFollowTag;
vector<long> followers, anchoredTag, followersRecord;
vector<long> verSetIndex, setTagIDs;
vector<long> km1ShellUpperBounds;
vector<long> km1Neighbors;
vector<long> anchorIDs;
vector<vector<long> > verSet, km1IdSubTag, km1Nei, kmbNei, km1ShellDownLayerNei, km1ShellOneLayerNei, km1ShellUpLayerNei, km1ShellLayerCand;

bool compareUB(const long& a, const long& b)//sort comparison 
{
	return km1ShellUpperBounds[a] > km1ShellUpperBounds[b];
}

void dataInput()
{
	double time0 = (double)clock() / CLOCKS_PER_SEC;

	//read edges, build verSet //need first row ordered data
	long vertexID, neighborID;
	vector<long> verSetInsertion;
	FILE* fe = NULL;
	fe = fopen(infile.c_str(), "r");
	if (fe == NULL)
	{
		printf("ERROR!");
		exit(1);
	}
	else
	{
		char fech = 's';
		while (fech != '\377')
		{
			vertexID = 0;
			fech = getc(fe);
			if (fech == '\377') break;
			while (fech != '\t')
			{
				vertexID = vertexID * 10 + (int)fech - 48;
				fech = getc(fe);
			}
			neighborID = 0;
			fech = getc(fe);
			while (fech != '\n' && fech != '\377')
			{
				neighborID = neighborID * 10 + (int)fech - 48;
				fech = getc(fe);
			}

			if (verSet.size() == 0)
			{
				verSetInsertion.clear();
				verSetInsertion.push_back(vertexID);
				verSetInsertion.push_back(neighborID);
				verSet.push_back(verSetInsertion);
			}
			else
			{
				if (vertexID == verSet[verSet.size() - 1][0])
				{
					verSet[verSet.size() - 1].push_back(neighborID);
				}
				else
				{
					verSetInsertion.clear();
					verSetInsertion.push_back(vertexID);
					verSetInsertion.push_back(neighborID);
					verSet.push_back(verSetInsertion);
				}
			}
		}
	}
	fclose(fe);

	long degree;
	vector<long> neiTemp;
	neiTemp.clear();
	for (unsigned long i = 0; i < verSet.size(); i++)
	{
		kTag.push_back(1);
		anchoredTag.push_back(0);
		anchorTag.push_back(0);
		followTag.push_back(0);
		subgraphTag.push_back(0);
		kcoreDeletTag.push_back(0);
		kmbNei.push_back(neiTemp);
		km1ShellDownLayerNei.push_back(neiTemp);
		km1ShellOneLayerNei.push_back(neiTemp);
		km1ShellUpLayerNei.push_back(neiTemp);
		km1ShellTag.push_back(0);
		km1ShellCandidateIds.push_back(0);
		km1ShellUpperBounds.push_back(0);
	}
	km1Tag = kTag;
	kmbTag = kTag;
	long maxDegree = 0;
	double averageDegree = 0;
	setTagIDs.clear();
	for (unsigned long i = 0; i < verSet.size(); i++)
	{
		degree = verSet[i].size() - 1;
		if (degree < inputK)
		{
			//setTagIDs.push_back(i);
			kTag[i] = 2;
		}
		if (degree < inputK - inputB)
		{
			setTagIDs.push_back(i);
			kmbTag[i] = 2;
		}
		verDegree.push_back(degree);
		if (degree > maxDegree)
		{
			maxDegree = degree;
		}
		averageDegree += degree;
	}

	kDegree = verDegree;
	km1Degree = verDegree;
	kmbDegree = verDegree;

	subgraphDegree = subgraphTag;

	averageDegree /= (double)(verSet.size());
	printf("max Degree: %ld\naverage Degree: %lf\n", maxDegree, averageDegree);

	//Index verSet
	for (unsigned long i = 0; i < verSet.size(); i++)
	{
		verSetIndex[verSet[i][0]] = i;
	}

	double time1 = (double)clock() / CLOCKS_PER_SEC;
	printf("read time: %lf\n", time1 - time0);

}

void computeKmbKcore()
{
	for (unsigned long i = 0; i < setTagIDs.size(); i++)
	{
		long setTagID = setTagIDs[i];
		kmbTag[setTagID] = 0; //cur vertex deleted
		kmbDegree[setTagID] = 0;
		kmbcoreNumber--;
		for (unsigned long j = 1; j < verSet[setTagID].size(); j++)
		{
			long neiVertexID = verSet[setTagID][j];
			long neiverSetID = verSetIndex[neiVertexID];
			if (verSet[neiverSetID][0] == neiVertexID && kmbTag[neiverSetID] != 0)
			{
				kmbDegree[neiverSetID]--; //neighbor degree - 1
				if (kmbTag[neiverSetID] == 1 && kmbDegree[neiverSetID] < inputK - inputB) //new candidate for computing
				{
					setTagIDs.push_back(neiverSetID);
					kmbTag[neiverSetID] = 2;
				}
			}
		}
	}

	//kTag: k-core tag towards verSet
	//kDegree : k-core degree
	kTag = kmbTag;
	kDegree = kmbDegree;
	setTagIDs.clear();
	kcoreNumber = kmbcoreNumber;
	for (unsigned long i = 0; i < kDegree.size(); i++)
	{
		if (kDegree[i] < inputK && kTag[i])
		{
			setTagIDs.push_back(i);
			kTag[i] = 2;
		}
	}
	for (unsigned long i = 0; i < setTagIDs.size(); i++)
	{
		long setTagID = setTagIDs[i];
		kTag[setTagID] = 0; //cur vertex deleted
		kDegree[setTagID] = 0;
		kcoreNumber--;
		for (unsigned long j = 1; j < verSet[setTagID].size(); j++)
		{
			long neiVertexID = verSet[setTagID][j];
			long neiverSetID = verSetIndex[neiVertexID];
			if (verSet[neiverSetID][0] == neiVertexID && kTag[neiverSetID] != 0)
			{
				kDegree[neiverSetID]--; //neighbor degree - 1
				if (kTag[neiverSetID] == 1 && kDegree[neiverSetID] < inputK) //new candidate for computing
				{
					setTagIDs.push_back(neiverSetID);
					kTag[neiverSetID] = 2;
				}
			}
		}
	}
}

void constrcutKm1SubsV6() //km1Tag, km1Degree, km1ShellTag, km1Nei, km1ShellDownLayerNei, km1ShellOneLayerNei, km1ShellUpLayerNei, km1ShellCandidates, lastCandidates
{
	km1Tag = kmbTag;
	km1Degree = kmbDegree;
	long km1coreNumber = kmbcoreNumber + kmbNeighborNumber;
	vector<long> temp, temp1, temp2;

	//get candidates for (k-1)-core deleting
	setTagIDs.clear();
	for (unsigned long i = 0; i < kmbDegree.size(); i++)
	{
		if (kmbDegree[i] < inputK - 1 && kmbTag[i] && !anchorTag[i] && !followTag[i])
		{
			setTagIDs.push_back(i);
			km1Tag[i] = 2;//tag for first layer of (k-1)-shell 
		}
	}

	//find (k-1)-core
	for (unsigned long i = 0; i < setTagIDs.size(); i++)
	{
		long setTagID = setTagIDs[i];
		km1Tag[setTagID] = 0; //cur vertex deleted
		km1Degree[setTagID] = 0;

		km1coreNumber--;
		for (unsigned long j = 1; j < verSet[setTagID].size(); j++)
		{
			long neiVertexID = verSet[setTagID][j];
			long neiverSetID = verSetIndex[neiVertexID];
			if (verSet[neiverSetID][0] == neiVertexID && km1Tag[neiverSetID] != 0)
			{
				km1Degree[neiverSetID]--; //neighbor degree - 1
				if ((km1Tag[neiverSetID] == 1 || km1Tag[neiverSetID] == -1) && !anchorTag[neiverSetID] && !followTag[neiverSetID] && km1Degree[neiverSetID] < inputK - 1) //new candidate for computing
				{
					setTagIDs.push_back(neiverSetID);
					km1Tag[neiverSetID] = 2;
				}
			}
		}
	}

	if (km1coreNumber > maxKm1coreNumber)
	{
		maxKm1coreNumber = km1coreNumber;
	}

	km1ShellTag = km1Tag;
	vector<long> kDegree1 = km1Degree;

	//find neighbors of (k-1)-core
	km1Neighbors.clear();
	for (unsigned long i = 0; i < km1Tag.size(); i++)
	{
		if (km1Tag[i] == 1)
		{
			for (unsigned long j = 1; j < verSet[i].size(); j++)
			{
				long neiVertexID = verSet[i][j];
				long neiverSetID = verSetIndex[neiVertexID];
				if (verSet[neiverSetID][0] == neiVertexID && !km1Tag[neiverSetID] && kmbTag[neiverSetID])
				{
					km1Neighbors.push_back(neiverSetID);
					km1Tag[neiverSetID] = -1;
					////km1Degree[i]++;
				}
			}
		}
	}
	for (unsigned long i = 0; i < km1Neighbors.size(); i++)
	{
		long id = km1Neighbors[i];
		long degree = 0;
		for (unsigned long j = 1; j < verSet[id].size(); j++)
		{
			long neiVertexID = verSet[id][j];
			long neiverSetID = verSetIndex[neiVertexID];
			if (verSet[neiverSetID][0] == neiVertexID && km1Tag[neiverSetID])
			{
				if (km1Tag[neiverSetID] == 1)
				{
					km1Degree[neiverSetID]++;
				}
				degree++;
			}
		}
		km1Degree[id] = degree;
	}
	//sort(km1Neighbors.begin(), km1Neighbors.end(), compare2);
	long km1NeighborNumber = km1Neighbors.size();
	if (km1NeighborNumber > maxKm1NeighborNumber)
	{
		maxKm1NeighborNumber = km1NeighborNumber;
	}

	//build neighbor set //km1Tag = -1: neghbor of k-1 core
	km1Nei = kmbNei;
	for (unsigned long i = 0; i < km1Tag.size(); i++)
	{
		if (km1Tag[i])
		{
			for (unsigned long j = 1; j < verSet[i].size(); j++)
			{
				long neiVertexID = verSet[i][j];
				long neiverSetID = verSetIndex[neiVertexID];
				if (verSet[neiverSetID][0] == neiVertexID && km1Tag[neiverSetID])
				{
					km1Nei[i].push_back(neiverSetID);
				}
			}
		}
	}

	//build (k-1)-shell: km1ShellCandidateTag, km1ShellDownLayerNei

	//get candidates for k-core deleting
	setTagIDs.clear();
	for (unsigned long i = 0; i < km1Degree.size(); i++)
	{
		if (kDegree1[i] < inputK && km1ShellTag[i] && !anchorTag[i] && !followTag[i])
		{
			setTagIDs.push_back(i);
			km1ShellTag[i] = 2;//tag for first layer of (k-1)-shell 
		}
	}

	//find k-core
	unsigned long layerNum = 3;
	unsigned long layerTag = setTagIDs.size() - 1;

	for (unsigned long i = 0; i < setTagIDs.size(); i++)
	{
		if (i > layerTag)
		{
			layerTag = setTagIDs.size() - 1;
			layerNum++;
		}
		long setTagID = setTagIDs[i];

		for (unsigned long j = 0; j < km1Nei[setTagID].size(); j++)
		{
			long nei = km1Nei[setTagID][j];
			if (km1ShellTag[nei] == 1)
			{
				kDegree1[nei]--; //neighbor degree - 1
				if (kDegree1[nei] < inputK && !anchorTag[nei] && !followTag[nei]) //new candidate for computing
				{
					setTagIDs.push_back(nei);
					km1ShellTag[nei] = layerNum;
				}
			}
		}
	}
	maxLayerNum = layerNum;

	//tag and record: km1ShellDownLayerNei, km1ShellOneLayerNei, km1ShellUpLayerNei, km1ShellCandidates, lastCandidates
	km1ShellSize = setTagIDs.size();
	km1ShellLayerRecord.clear();
	km1ShellLayerCand.clear();
	//renew km1ShellDownLayerNei, km1ShellOneLayerNei, km1ShellUpLayerNei
	temp.clear();
	for (unsigned long i = 0; i < lastCandidates.size(); i++)
	{
		long id = lastCandidates[i];
		km1ShellDownLayerNei[id] = temp;
		km1ShellOneLayerNei[id] = temp;
		km1ShellUpLayerNei[id] = temp;
	}

	//build km1ShellLayerRecord, km1ShellLayerCand; 
	for (unsigned long i = 0; i < layerNum; i++)
	{
		km1ShellLayerRecord.push_back(0);
		km1ShellLayerCand.push_back(temp);
	}

	for (unsigned long i = 0; i < setTagIDs.size(); i++)
	{
		long setTagID = setTagIDs[i];
		long layerNum = km1ShellTag[setTagID];
		temp.clear();
		temp1.clear();
		temp2.clear();
		for (unsigned long j = 0; j < km1Nei[setTagID].size(); j++)
		{
			long nei = km1Nei[setTagID][j];
			long num = km1ShellTag[nei];
			if (num > layerNum || num == 1)
			{
				temp.push_back(nei);
			}
			else if (num == layerNum)
			{
				temp1.push_back(nei);
			}
			else if (num)
			{
				temp2.push_back(nei);
			}
		}
		km1ShellDownLayerNei[setTagID] = temp;
		km1ShellOneLayerNei[setTagID] = temp1;
		km1ShellUpLayerNei[setTagID] = temp2;
	}

	//km1ShellCandidates = setTagIDs;

	for (unsigned long i = 0; i < km1Neighbors.size(); i++)
	{
		long id = km1Neighbors[i];
		temp.clear();
		//km1ShellCandidates.push_back(id);
		for (unsigned long j = 0; j < km1Nei[id].size(); j++)
		{
			long nei = km1Nei[id][j];
			if (km1ShellTag[nei])
			{
				temp.push_back(nei);
				km1ShellUpLayerNei[nei].push_back(id);
			}
		}
		km1ShellDownLayerNei[id] = temp;
	}

	//order by km1shell neighbors and layers
	km1ShellCandidates = km1Neighbors;
	for (unsigned long i = 0; i < setTagIDs.size(); i++)
	{
		km1ShellCandidates.push_back(setTagIDs[i]);
	}

	lastCandidates = km1ShellCandidates;

	//build upperbounds
	if (1)
	{
		vector<vector<long> > km1ShellLayers;
		temp.clear();
		for (unsigned long i = 0; i < maxLayerNum; i++)
		{
			km1ShellLayers.push_back(temp);
		}
		for (unsigned long i = 0; i < km1ShellCandidates.size(); i++)
		{
			long id = km1ShellCandidates[i];
			if (km1Tag[id] == -1)
			{
				km1ShellLayers[1].push_back(id);
			}
			else
			{
				long layer = km1ShellTag[id];
				km1ShellLayers[layer].push_back(id);
			}
		}
		for (unsigned long i = maxLayerNum - 1; i > 0; i--)
		{
			for (unsigned long j = 0; j < km1ShellLayers[i].size(); j++)
			{
				long id = km1ShellLayers[i][j];
				for (unsigned long k = 0; k < km1ShellUpLayerNei[id].size(); k++)
				{
					long nid = km1ShellUpLayerNei[id][k];
					km1ShellUpperBounds[nid] = km1ShellUpperBounds[nid] + km1ShellUpperBounds[id] + 1;
				}
			}
		}
		sort(km1ShellCandidates.begin(), km1ShellCandidates.end(), compareUB);
	}
}

void backtrackingEarlyStopToComputeKcore(long ver, vector<long>& kcoreFailers)
{
	kcoreFailers.push_back(ver);
	subgraphTag[ver] = -1;
	subgraphDegree[ver] = 0;
	km1ShellCanSize--;
	for (unsigned long i = 0; i < km1ShellUpLayerNei[ver].size(); i++)
	{
		long id = km1ShellUpLayerNei[ver][i];
		if (subgraphTag[id] > 0)
		{
			subgraphDegree[id]--;
			if (subgraphDegree[id] < inputK)
			{
				backtrackingEarlyStopToComputeKcore(id, kcoreFailers);
			}
		}
	}
	for (unsigned long i = 0; i < km1ShellOneLayerNei[ver].size(); i++)
	{
		long id = km1ShellOneLayerNei[ver][i];

		if (subgraphTag[id] == 1)
		{
			subgraphDegree[id]--;
			if (subgraphDegree[id] < inputK)
			{
				backtrackingEarlyStopToComputeKcore(id, kcoreFailers);
			}
		}
	}
	for (unsigned long i = 0; i < km1ShellDownLayerNei[ver].size(); i++)
	{
		long id = km1ShellDownLayerNei[ver][i];
		if (subgraphTag[id] > 0)
		{
			subgraphDegree[id]--;
			if (subgraphDegree[id] < inputK)
			{
				backtrackingEarlyStopToComputeKcore(id, kcoreFailers);
			}
		}
	}
}

long activatedNumInV62S;

void buildKcoreCandidatesV62S(long lay, unsigned long layeri, vector<long>& kcoreCandidates, vector<long>& kcoreFailers) //bfs, sugraphTag, subgraphDegree, kcoreCandidates, 2:bfs, 3:bfs+early stop
{
	long anc = km1ShellLayerCand[lay][layeri];
	unsigned long downNum = km1ShellDownLayerNei[anc].size();

	subgraphDegree[anc] += downNum;

	for (unsigned long i = 0; i < km1ShellOneLayerNei[anc].size(); i++)
	{
		long id = km1ShellOneLayerNei[anc][i];
		if (subgraphTag[id] > 0)
		{
			subgraphDegree[anc]++;
		}
	}
	activatedNumInV62S++;
	if (subgraphDegree[anc] >= inputK)
	{
		subgraphTag[anc] = 1;
		kcoreCandidates.push_back(anc);
		if (!anchorTag[anc] && !followTag[anc])
		{
			km1ShellCanSize++;
		}
		for (unsigned long i = 0; i < downNum; i++)
		{
			long layer;
			long id = km1ShellDownLayerNei[anc][i];
			subgraphDegree[id]++;
			if (!subgraphTag[id] && !anchorTag[id] && !followTag[id])
			{
				layer = km1ShellTag[id];

				if (layer != 1)
				{
					if (!km1ShellLayerRecord[layer])
					{
						km1ShellLayerRecord[layer] = 1;
					}
					km1ShellLayerCand[layer].push_back(id);
					subgraphTag[id] = 2;
				}
			}
		}
	}
	else
	{
		kcoreFailers.push_back(anc);
		subgraphTag[anc] = -1;
		subgraphDegree[anc] = 0;
		for (unsigned long i = 0; i < km1ShellUpLayerNei[anc].size(); i++)
		{
			long id = km1ShellUpLayerNei[anc][i];
			if (subgraphTag[id] > 0)
			{
				subgraphDegree[id]--;
				if (subgraphDegree[id] < inputK)
				{
					backtrackingEarlyStopToComputeKcore(id, kcoreFailers);
				}
			}
		}
		for (unsigned long i = 0; i < km1ShellOneLayerNei[anc].size(); i++)
		{
			long id = km1ShellOneLayerNei[anc][i];
			if (subgraphTag[id] == 1)
			{
				subgraphDegree[id]--;
				if (subgraphDegree[id] < inputK)
				{
					backtrackingEarlyStopToComputeKcore(id, kcoreFailers);
				}
			}
		}
	}
	layeri++;
	if (layeri < km1ShellLayerCand[lay].size())
	{
		buildKcoreCandidatesV62S(lay, layeri, kcoreCandidates, kcoreFailers);
	}
	else
	{
		lay++;
		for (unsigned long i = lay; i < km1ShellLayerCand.size(); i++)
		{
			if (km1ShellLayerCand[i].size())
			{
				buildKcoreCandidatesV62S(i, 0, kcoreCandidates, kcoreFailers);
				break;
			}
		}
	}
}

void anchorOneV6SRecord(long anc, long maxSize)
{
	long bestAnchor = -1;
	vector<long> kcoreCandidates, kcoreFailers, deletSet;

	anchorTag[anc] = 1;
	//tag candidates
	kcoreCandidates.clear();
	kcoreFailers.clear();

	for (unsigned long j = 0; j < km1ShellLayerRecord.size(); j++)
	{
		if (km1ShellLayerRecord[j])
		{
			km1ShellLayerCand[j].clear();
			km1ShellLayerRecord[j] = 0;
		}
	}
	long layer;
	if (km1Tag[anc] == -1)
	{
		layer = 1;
	}
	else
	{
		layer = km1ShellTag[anc];
	}
	subgraphDegree[anc] += inputK;
	km1ShellLayerRecord[layer] = 1;
	km1ShellLayerCand[layer].push_back(anc);
	buildKcoreCandidatesV62S(layer, 0, kcoreCandidates, kcoreFailers);
	for (unsigned long j = 0; j < kcoreFailers.size(); j++)
	{
		long id = kcoreFailers[j];
		subgraphTag[id] = 0;
		//subgraphDegree[id] = 0;
	}

	//count number, //recover subgraphTag, subgraphDegree
	long incNum = 1;
	for (unsigned long j = 0; j < kcoreCandidates.size(); j++)
	{
		long id = kcoreCandidates[j];
		if (subgraphTag[id])
		{
			if (subgraphTag[id] == 1)
			{
				if (!kTag[id] && !anchorTag[id] && !followTag[id])
				{
					incNum++;
					followTag[id] = 1;
					followersRecord.push_back(id);
				}
			}
			subgraphTag[id] = 0;
		}
		if (subgraphDegree[id])
		{
			subgraphDegree[id] = 0;
		}
	}

	anchorNumber++;
	anchorAndfollowerNumber += incNum;

}

long anchorOneV6S() //shrink early stop in bfs
{
	long bestAnchor = -1, maxSize = 0;
	vector<long> kcoreCandidates, kcoreFailers, deletSet;
	//vector<long> kcoreDeletInitalCans, kcoreDeletCandidates, kcoreTagRecord;
	vector<long> triedFollowTag = followTag;
	vector<long> updatedUpperBounds = km1ShellUpperBounds;
	unsigned long km1shellcansize = km1ShellCandidates.size();

	for (unsigned long i = 0; i < km1shellcansize; i++)
	{
		long anc = km1ShellCandidates[i];
		anchorTries++;
		anchorTriesRA++;

		if (!triedFollowTag[anc])
		{
			anchorTriesFE++;
			if (updatedUpperBounds[anc] >= maxSize)//|| (km1shellCandidatesOrder == 6 && i < lastLayersNumber))
			{
				anchorTriesUB++;
				anchorTag[anc] = 1;
				//tag candidates
				km1ShellCanSize = 1;
				kcoreCandidates.clear();
				kcoreFailers.clear();
				activatedNumInV62S = 0;
				for (unsigned long j = 0; j < km1ShellLayerRecord.size(); j++)
				{
					if (km1ShellLayerRecord[j])
					{
						km1ShellLayerCand[j].clear();
						km1ShellLayerRecord[j] = 0;
					}
				}
				long layer;
				if (km1Tag[anc] == -1)
				{
					layer = 1;
				}
				else
				{
					layer = km1ShellTag[anc];
				}
				subgraphDegree[anc] += inputK;
				km1ShellLayerRecord[layer] = 1;
				km1ShellLayerCand[layer].push_back(anc);
				buildKcoreCandidatesV62S(layer, 0, kcoreCandidates, kcoreFailers);
				for (unsigned long j = 0; j < kcoreFailers.size(); j++)
				{
					long id = kcoreFailers[j];
					subgraphTag[id] = 0;
					//subgraphDegree[id] = 0;
				}

				anchoringTries++;
				activatedCandidates += activatedNumInV62S;
				//k-core cumputation
				long incNum;
				if (1)
				{
					//kcore check

					//count number, //recover subgraphTag, subgraphDegree
					incNum = 1;
					for (unsigned long j = 0; j < kcoreCandidates.size(); j++)
					{
						long id = kcoreCandidates[j];
						if (subgraphTag[id])
						{
							if (subgraphTag[id] == 1 && !kTag[id] && !anchorTag[id] && !followTag[id])
							{
								incNum++;
								totalTriedFollowers++;
								if (!triedFollowTag[id])
								{
									triedFollowTag[id] = 1;
								}
							}
							subgraphTag[id] = 0;
						}
						if (subgraphDegree[id])
						{
							subgraphDegree[id] = 0;
						}
					}
					anchorTag[anc] = 0;

					if (!incNum)
					{
						incNum = 0;
					}


					if (incNum > maxSize)
					{
						maxSize = incNum;
						bestAnchor = anc;
					}
				}
				else //recover
				{
					anchorTag[anc] = 0;
					for (unsigned long j = 0; j < kcoreCandidates.size(); j++)
					{
						long id = kcoreCandidates[j];
						if (subgraphTag[id])
						{
							subgraphTag[id] = 0;
						}
						if (subgraphDegree[id])
						{
							subgraphDegree[id] = 0;
						}
					}
				}
				updatedUpperBounds[anc] = incNum - 1;
			}
			else // all upperbounds < maxSize
			{
				break;
			}
		}
		if (i == lastLayersNumber)
		{
			numberAfterLastLayers += maxSize;
		}
	}

	//record anchor and followers
	if (bestAnchor != -1)
	{
		long layer;
		if (km1Tag[bestAnchor] == -1)
		{
			layer = 1;
		}
		else
		{
			layer = km1ShellTag[bestAnchor];
		}
		anchorOneV6SRecord(bestAnchor, maxSize);
	}

	return bestAnchor;
}

void algorithm()
{
	algStartTime = (double)clock() / CLOCKS_PER_SEC;
	long anchorID;

	anchorNumber = 0;
	anchorAndfollowerNumber = 0;

	//(k-b)-core pre-process
	kmbcoreNumber = kmbTag.size();
	//get (k-b)-core and k-core
	computeKmbKcore();
	//get neighbors for anchoring
	vector<long> kmbNeighbors;
	kmbNeighbors.clear();
	for (unsigned long i = 0; i < kmbTag.size(); i++)
	{
		if (kmbTag[i] == 1)
		{
			for (unsigned long j = 1; j < verSet[i].size(); j++)
			{
				long neiVertexID = verSet[i][j];
				long neiverSetID = verSetIndex[neiVertexID];
				if (verSet[neiverSetID][0] == neiVertexID && !kmbTag[neiverSetID])
				{
					kmbNeighbors.push_back(neiverSetID);
					kmbTag[neiverSetID] = -1;
					//kmbDegree[i]++;
				}
			}
		}
	}
	for (unsigned long i = 0; i < kmbNeighbors.size(); i++)
	{
		long id = kmbNeighbors[i];
		long degree = 0;
		for (unsigned long j = 1; j < verSet[id].size(); j++)
		{
			long neiVertexID = verSet[id][j];
			long neiverSetID = verSetIndex[neiVertexID];
			if (verSet[neiverSetID][0] == neiVertexID && kmbTag[neiverSetID])
			{
				if (kmbTag[neiverSetID] == 1)
				{
					kmbDegree[neiverSetID]++;
				}
				degree++;
			}
		}
		kmbDegree[id] = degree;
	}
	kmbNeighborNumber = kmbNeighbors.size();

	//anchoring
	for (long i = 0; i < inputB; i++)
	{
		//data struction, get km1Tag, km1Degree, km1Nei, km1ShellDownLayerNei
		constrcutKm1SubsV6();
		lastLayersNumber = 0;

		//anchoring algorithm
		anchorID = anchorOneV6S();
		anchorIDs.push_back(anchorID);
	}

	runTime = (double)clock() / CLOCKS_PER_SEC - algStartTime;
	printf("time: %lf\n", runTime);
}

void dataOutput()
{
	//write
	char record[100];
	FILE* fs;
	fs = fopen(outfile.c_str(), "a");
	char fsch;

	if (fs == NULL)
	{
		printf("ERROR!");
		exit(1);
	}
	else
	{
		sprintf(record, "k%ld", inputK);
		fwrite(record, sizeof(*record), strlen(record), fs);
		fsch = putc('\t', fs);
		sprintf(record, "b%ld", inputB);
		fwrite(record, sizeof(*record), strlen(record), fs);
		fsch = putc('\t', fs);
		sprintf(record, "%.2lf", runTime);
		fwrite(record, sizeof(*record), strlen(record), fs);
		fsch = putc('\t', fs);
		sprintf(record, "%.3lf", anchorAndfollowerNumber - inputB);
		fwrite(record, sizeof(*record), strlen(record), fs);
		fsch = putc('\t', fs);
		sprintf(record, "kc%ld", kcoreNumber);
		fwrite(record, sizeof(*record), strlen(record), fs);
		fsch = putc('\t', fs);
		for (unsigned long i = 0; i < anchorIDs.size(); i++)
		{
			long id = anchorIDs[i];
			if (id != -1)
			{
				fsch = putc('\t', fs);
				sprintf(record, "a%ld", verSet[id][0]);
				fwrite(record, sizeof(*record), strlen(record), fs);
			}
			else
			{
				fsch = putc('\t', fs);
				sprintf(record, "aX");
				fwrite(record, sizeof(*record), strlen(record), fs);
			}

		}
		fsch = putc('\n', fs);
	}

}

int main(int argc, char* argv[])
{
	//configuration
	long maxVerID = 1000000000; //max vertex id
	verSetIndex.resize(maxVerID);
	infile = "sample_youtube.txt";

	outfile = "as-res.txt";
	scanf("%ld %ld", &inputK, &inputB);

	//input data 
	dataInput();

	//algorithm start
	algorithm();

	//output data
	dataOutput();

	return 0;
}