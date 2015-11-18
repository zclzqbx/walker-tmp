/*function:场景树生成程序（C++）
 *author:朱传林
 *date:2015.11.18
 *source:Scenario reduction and scenario tree construction for power management problems
 *detail:
 *程序主要功能：给定一批场景及第个时段保留的节点数生成场景树
 *程序分为两步：第一步是在给定场景的基础之上进行缩减，缩减到最终保留的场景个数；
 *第二步，对上一步生成的场景进一步缩减形成场景树。
 */

#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<time.h>

using namespace std;

const int Pre_N=200;//the number of scenario before reduction
const int Aft_N=20;//the number of scenario after reduction
const int NT=24;//time

const int J[NT]={0,1,0,0,1,0,2,0,0,1,0,0,2,0,2,0,0,3,0,3,0,4,0,180};//各个时段删除的节点数

double c_distance(double s1[NT],double s2[NT],int T)
{
	if(s1==NULL || s2==NULL || T<0)
		return -1;
	double result=0;
	for(int t=0;t<T;++t)
	{
		result+=sqrt((s1[t]-s2[t])*(s1[t]-s2[t]));
	}
	return result;
}

//two step of scenario tree generation
bool scenarioTreeStepOne(double scenario_origin[Pre_N][NT],double scenario_new[Aft_N][NT],
		bool deleteFlag[Pre_N],double scenarioP[Pre_N],double scenarioP_new[Aft_N])
{
	if(scenario_origin==NULL || scenario_new==NULL || 
			 deleteFlag==NULL || scenarioP==NULL || scenarioP_new==NULL)
		return false;

	int deleteNum=0;	
	
	while(deleteNum < Pre_N - Aft_N)
	{
		vector<double> z_l(Pre_N,-1);
		vector< vector<double> > c_l(Pre_N,vector<double>(Pre_N,-1));
		
		for(int i=0;i<Pre_N;++i)//Ckj
		{
			if(!deleteFlag[i])
			{
				double tmp=10000000.0;
				for(int k=0;k<Pre_N;++k)
				{
					if(k==i || deleteFlag[k])
					{
						for(int h=0;h<Pre_N;++h)
						{
							if(h!=i && !deleteFlag[h])
							{
								double tmp1=c_distance(scenario_origin[k],scenario_origin[h],NT);
								if(tmp>tmp1)
									tmp=tmp1;
							}
						}
						c_l[k][i]=tmp;
					}
				}
			}
		}

		for(int i=0;i<Pre_N;++i)//Zl
		{
			double tmp=0;
			if(!deleteFlag[i])
			{
				for(int k=0;k<Pre_N;++k)
				{
					if(k==i || deleteFlag[k])
					{
						tmp+=scenarioP[k]*c_l[k][i];
					}
				}
				z_l[i]=tmp;
			}
		}

		int deleteScenarioIndex=-1;//the next scenario being deleted
		double deleteScenarioValue=10000000.0;
		for(int i=0;i<Pre_N;++i)
		{
			if(!deleteFlag[i])
			{
				if(z_l[i]<deleteScenarioValue)
				{
					deleteScenarioValue=z_l[i];
					deleteScenarioIndex=i;
				}
			}
		}
		
		int mergeScenariosIndex=-1;
		double minDistance=10000000.0;

		for(int i=0;i<Pre_N;++i)//the nearest scenario
		{
			if(!deleteFlag[i] && i!=deleteScenarioIndex)
			{
				double tmp=c_distance(scenario_origin[deleteScenarioIndex],scenario_origin[i],NT);
				if(minDistance>tmp)
				{
					minDistance=tmp;
					mergeScenariosIndex=i;
				}
			}
		}

		scenarioP[mergeScenariosIndex]+=scenarioP[deleteScenarioIndex];//add the probability to the nearest scenaro
		deleteFlag[deleteScenarioIndex]=true;

		deleteNum++;
	}

	int theKeepScenarioNum=0;
	for(int i=0;i<Pre_N;++i)
	{
		if(!deleteFlag[i])
		{
			if(theKeepScenarioNum < Aft_N)
			{
				for(int t=0;t<NT;++t)
				{
					scenario_new[theKeepScenarioNum][t]=scenario_origin[i][t];
				}
				scenarioP_new[theKeepScenarioNum]=scenarioP[i];
			}
			theKeepScenarioNum++;
		}
		if(theKeepScenarioNum>Aft_N)
			return false;
	}
	return true;
}


bool scenarioTreeStepTwo(double scenario_new[Aft_N][NT],double scenarioP_new[Pre_N])
{
	return true;
}

//interface for user
bool scenario_reduce(double scenario_origin[Pre_N][NT],double scenario_new[Aft_N][NT],
		bool deleteFlag[Pre_N],double scenarioP[Pre_N],double scenarioP_new[Pre_N])
{
	bool stepOne=scenarioTreeStepOne(scenario_origin,scenario_new,
			deleteFlag,scenarioP,scenarioP_new);//the same to the user interface
	bool stepTwo=scenarioTreeStepTwo(scenario_new,scenarioP_new);
	
	if(stepOne && stepTwo)
		return true;	
	return false;
}

int main()
{
	double scenario_origin[Pre_N][NT];
	ifstream input("scenario_generate/mc_scenarios.txt",ios::in);
	if(!input)
	{
		cout<<"mc_scenarios is not exist!"<<endl;
		return 1;
	}
	bool deleteFlag[Pre_N];
	double scenarioP[Pre_N];
	double scenarioP_new[Pre_N]={0.0};

	for(int i=0;i<Pre_N;++i)
		scenarioP[i]=1.0/Pre_N;

	for(int i=0;i<Pre_N;++i)
		deleteFlag[i]=false;//all not be deleted

	for(int i=0;i<Pre_N;++i)//input data
	{
		for(int t=0;t<24;++t)
		{
			input>>scenario_origin[i][t];
		}
	}

	double scenario_new[Aft_N][NT];
	ofstream output("new_scenario.txt",ios::ate);
	
	time_t nowTime=time(0);//print the time
	struct tm* nowTimeStruct=localtime(&nowTime);
	output<<"系统当前时间："<<1900+nowTimeStruct->tm_year<<"."//year
		<<nowTimeStruct->tm_mon<<"."//month
			<<nowTimeStruct->tm_mday<<"  "//day
				<<nowTimeStruct->tm_hour<<":"//hour
					<<nowTimeStruct->tm_min<<":"//min
						<<nowTimeStruct->tm_sec<<endl;//scecond

	//scenario reduction function
	bool flag=scenario_reduce(scenario_origin,scenario_new,
				deleteFlag,scenarioP,scenarioP_new);
	if(!flag)
	{
		cout<<"the process failed..."<<endl;
		return 1;
	}

	for(int i=0;i<Aft_N;++i)
	{
		for(int t=0;t<NT;++t)
		{
			output<<scenario_new[i][t]<<"   ";
		}
		output<<endl;
	}

	input.close();
	output.close();
	return 0;
}
