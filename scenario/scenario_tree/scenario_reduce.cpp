
/*function:场景缩减程序（C++）
 *author:朱传林
 *date:2015.11.18
 *source:Scenario reduction and scenario tree construction for power management problems
 *detail:
 *function：reduce the number of scenario to the fixed scenario
 */


#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<limits.h>
#include<time.h>

using namespace std;

const int Pre_N=200;//the number of scenario before reduction
const int Aft_N=20;//the number of scenario after reduction
const int NT=24;//time


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
	double scenarioP[Pre_N];//probability of every scenario
	double scenarioP_new[Aft_N];

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
	
	/*time_t nowTime=time(0);
	struct tm* nowTimeStruct=localtime(&nowTime);
	output<<"系统当前时间："<<1900+nowTimeStruct->tm_year<<"."<<nowTimeStruct->tm_mon<<"."<<
		nowTimeStruct->tm_mday<<"  "<<nowTimeStruct->tm_hour<<":"<<nowTimeStruct->tm_min<<":"<<nowTimeStruct->tm_sec<<endl;

*/
	bool flag=scenarioTreeStepOne(scenario_origin,scenario_new,deleteFlag,scenarioP,scenarioP_new);//scenario reduction function
	if(!flag)
	{
		cout<<"the process failed..."<<endl;
		return 1;
	}

	for(int i=0;i<Aft_N;++i)
	{
		output<<scenarioP_new[i]<<endl;
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
