#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

ifstream inputP0("random_scenario.txt",ios::in);
ifstream inputPs("merge_scenarios.txt",ios::in);

const int S  = 5;
const int NT = 24;

double distanceF(double p1[],double p2[],int n = NT);

int main()
{
	//输入
	double randomSceanrio[NT];
	if(!inputP0)
	{
		cerr<<"file random_scenario is not exist!"<<endl;
		return -1;
	}
	for(int t = 0;t < NT;++t)
	{
		inputP0>>randomSceanrio[t];
	}
	
	double scenarios[S][NT];
	if(!inputPs)
	{
		cerr<<"file merge_scenarios is not exist!"<<endl;
		return -1;
	}
	for(int s = 0;s < S;++s)
	{
		for(int t = 0;t < NT;++t)
		{
			inputPs>>scenarios[s][t];
		}
	}
	
	double tmp1 = 0;
	double tmp2 = 0;
	
	for(int s = 0;s < S;++s)
	{
		tmp1 += distanceF(randomSceanrio,scenarios[s]);
	}
	tmp1 = tmp1 / S;
	
	for(int s1 = 0;s1 < S;++s1)
	{
		for(int s2 = 0;s2 < S;++s2)
		{
			if(s1 != s2)
			{
				tmp2 += distanceF(scenarios[s1],scenarios[s2]);
			}
		}
	}
	tmp2 = tmp2 / (2*S*S);
	
	double es = tmp1 - tmp2;
	cout<<es<<endl;
	
	return 0;
}

double distanceF(double p1[],double p2[],int n)
{
	if(n < 0 || NULL == p1 || NULL == p2)
	{
		return 0;
	}
	
	double distance = 0;
	for(int t = 0;t < n;++t)
	{
		distance += fabs(p1[t] - p2[t]);
	}
	return distance;
}