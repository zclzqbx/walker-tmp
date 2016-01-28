/*function:多风电场场景树生成程序（C++）
 *author:朱传林
 *date:2015.11.30
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

const int Pre_N   = 500;//the number of scenario before reduction
const int Aft_N   = 20;//the number of scenario after reduction
const int NT      = 24;//time
const int WINDNUM = 3;

double wind1[Pre_N][NT],wind2[Pre_N][NT],wind3[Pre_N][NT];//数据量过大，容易引起溢出
const int J[NT] = {0,0,0,0,1,0,2,0,0,1,0,0,2,0,2,0,0,3,0,3,0,4,0,480};//各个时段删除的节点数
//const int J[NT]={0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,498};//2
//const int J[NT]={0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,2,0,0,0,0,495};//5
int fAftToPre[Aft_N] = {0};

double c_distance(double s1[NT],double s2[NT],int T)
{//求两个场景的距离
	if(s1 == NULL || s2 == NULL || T < 0)
	{
		return -1;
	}
	double result = 0;
	for(int t = 0;t < T;++t)
	{
		result += sqrt((s1[t] - s2[t]) * (s1[t] - s2[t]));
	}
	return result;
}

//two step of scenario tree generation
bool scenarioTreeStepOne(double scenario_origin[Pre_N][NT],double scenario_new[Aft_N][NT],
		bool deleteFlag[Pre_N],double scenarioP[Pre_N],double scenarioP_new[Aft_N])
{
	if(scenario_origin == NULL || scenario_new == NULL || deleteFlag == NULL || 
			scenarioP == NULL || scenarioP_new == NULL)
	{
		return false;
	}

	int deleteNum = 0;	
	
	while(deleteNum < Pre_N - Aft_N)
	{
		vector<double> z_l(Pre_N,-1);
		vector< vector<double> > c_l(Pre_N,vector<double>(Pre_N,-1));
		
		for(int i=0;i<Pre_N;++i)//Ckj
		{
			if(!deleteFlag[i])
			{
				double tmp = 10000000.0;
				for(int k = 0;k < Pre_N;++k)
				{
					if(k == i || deleteFlag[k])
					{
						for(int h = 0;h < Pre_N;++h)
						{
							if(h != i && !deleteFlag[h])
							{
								double tmp1 = c_distance(scenario_origin[k],scenario_origin[h],NT);
								if(tmp > tmp1)
								{
									tmp = tmp1;
								}
							}
						}
						c_l[k][i] = tmp;
					}
				}
			}
		}

		for(int i = 0;i < Pre_N;++i)//Zl
		{
			double tmp = 0;
			if(!deleteFlag[i])
			{
				for(int k = 0;k < Pre_N;++k)
				{
					if(k == i || deleteFlag[k])
					{
						tmp += scenarioP[k] * c_l[k][i];
					}
				}
				z_l[i] = tmp;
			}
		}

		int deleteScenarioIndex    = -1;//the next scenario being deleted
		double deleteScenarioValue = 10000000.0;
		for(int i=0;i<Pre_N;++i)
		{
			if(!deleteFlag[i])
			{
				if(z_l[i] < deleteScenarioValue)
				{
					deleteScenarioValue = z_l[i];
					deleteScenarioIndex = i;
				}
			}
		}
		
		int mergeScenariosIndex = -1;
		double minDistance      = 10000000.0;

		for(int i=0;i<Pre_N;++i)//the nearest scenario
		{
			if(!deleteFlag[i] && i!=deleteScenarioIndex)
			{
				double tmp=c_distance(scenario_origin[deleteScenarioIndex],scenario_origin[i],NT);
				if(minDistance > tmp)
				{
					minDistance = tmp;
					mergeScenariosIndex = i;
				}
			}
		}

		//add the probability to the nearest scenaro
		scenarioP[mergeScenariosIndex]  += scenarioP[deleteScenarioIndex];
		deleteFlag[deleteScenarioIndex]  = true;

		deleteNum++;
	}

	int theKeepScenarioNum = 0;
	for(int i = 0;i < Pre_N;++i)
	{//提取出保存的场景
		if(!deleteFlag[i])
		{
			if(theKeepScenarioNum < Aft_N)
			{
				for(int t=0;t<NT;++t)
				{
					scenario_new[theKeepScenarioNum][t] = scenario_origin[i][t];					
				}
				scenarioP_new[theKeepScenarioNum] = scenarioP[i];
				fAftToPre[theKeepScenarioNum] = i;
			}
			theKeepScenarioNum++;
		}
		
		if(theKeepScenarioNum > Aft_N)
		{
			return false;
		}
	}
	return true;
}

bool scenarioTreeStepTwo(double scenario_new[Aft_N][NT],double scenarioP_new[Aft_N])
{
	if(scenario_new == NULL || scenarioP_new == NULL)
	{
		return false;
	}
	
	double pt[Aft_N][NT];//每个时段的概率
	bool deFlag[Aft_N];//删除标志
	
	//每次合并场景时，背后需要更改的场景，如果要删除场景1，场景1背后可能还有其他场景
	vector< vector<int> > hide(Aft_N,vector<int>(Aft_N,-1));
	
	for(int h=0;h < Aft_N;++h)
	{//初始化
		pt[h][NT-1] = scenarioP_new[h];
		deFlag[h]   = false;
	}
	
	for(int t = NT-2;t >= 0;--t)
	{
		if(J[t] > 0)//如果本阶段有缩减需求
		{
			int deNum = 0;
			while(deNum < J[t])
			{
				vector<double> z_l(Aft_N,-1);
				vector< vector<double> > c_l(Aft_N,vector<double>(Aft_N,-1));
				
				for(int i = 0;i < Aft_N;++i)//Ckj
				{
					if(!deFlag[i])
					{
						double tmp = 10000000.0;
						for(int k = 0;k < Aft_N;++k)
						{
							if(k == i || deFlag[k])
							{
								for(int h = 0;h < Aft_N;++h)
								{
									if(h != i && !deFlag[h])
									{
										double tmp1=c_distance(scenario_new[k],scenario_new[h],t+1);
										if(tmp > tmp1)
										{
											tmp = tmp1;
										}
									}
								}
								c_l[k][i] = tmp;
							}
						}
					}
				}

				for(int i = 0;i < Aft_N;++i)//Zl
				{
					double tmp = 0;
					if(!deFlag[i])
					{
						for(int k = 0;k < Aft_N;++k)
						{
							if(k == i || deFlag[k])
							{
								tmp += pt[k][t+1]*c_l[k][i];
							}
						}
						z_l[i] = tmp;
					}
				}
				
				int deleteScenarioIndex    = -1;//the next scenario being deleted
				double deleteScenarioValue = 10000000.0;
				for(int i = 0;i < Aft_N;++i)
				{
					if(!deFlag[i])
					{
						if(z_l[i] < deleteScenarioValue)
						{
							deleteScenarioValue = z_l[i];
							deleteScenarioIndex = i;
						}
					}
				}
				
				int mergeScenariosIndex = -1;
				double minDistance      = 10000000.0;

				for(int i = 0;i < Aft_N;++i)//the nearest scenario
				{
					if(!deFlag[i] && i != deleteScenarioIndex)
					{
						double tmp = c_distance(scenario_new[deleteScenarioIndex],scenario_new[i],t+1);
						if(minDistance > tmp)
						{
							minDistance         = tmp;
							mergeScenariosIndex = i;
						}
					}
				}
				
				pt[mergeScenariosIndex][t+1] += pt[deleteScenarioIndex][t+1];//合并概率
				
				int hide_index     = 0;				
				int hide_index_new = 0;
				int theScenarioNeedToBeModify = deleteScenarioIndex;
				
				while(hide[mergeScenariosIndex][hide_index_new] != -1)
				{
					hide_index_new++;
				}
				
				hide[mergeScenariosIndex][hide_index_new] = deleteScenarioIndex;//添加背后的场景
				hide_index_new++;
				
				do{//生成场景树
					for(int t1 = 0;t1 <= t;++t1)
					{
						scenario_new[theScenarioNeedToBeModify][t1] = scenario_new[mergeScenariosIndex][t1];
						wind1[fAftToPre[theScenarioNeedToBeModify]][t1] = wind1[fAftToPre[mergeScenariosIndex]][t1];
						wind2[fAftToPre[theScenarioNeedToBeModify]][t1] = wind2[fAftToPre[mergeScenariosIndex]][t1];
						wind3[fAftToPre[theScenarioNeedToBeModify]][t1] = wind3[fAftToPre[mergeScenariosIndex]][t1];
					}
					
					if(hide_index >= Aft_N || hide[deleteScenarioIndex][hide_index] == -1)
					{
						break;
					}
					else
					{
						theScenarioNeedToBeModify = hide[deleteScenarioIndex][hide_index];
						hide[mergeScenariosIndex][hide_index_new] = theScenarioNeedToBeModify;//背后所隐藏的场景也要转移
						hide_index++;
						hide_index_new++;
					}	
					
				}while(1);
				
				deFlag[deleteScenarioIndex] = true;
				for(int i=0;i<Aft_N;++i)
				{
					pt[i][t] = pt[i][t+1];
				}				
				deNum++;
			}
		}
	}	
	return true;
}

//interface for user
bool scenario_reduce(double scenario_origin[Pre_N][NT],double scenario_new[Aft_N][NT],
		bool deleteFlag[Pre_N],double scenarioP[Pre_N],double scenarioP_new[Pre_N])
{
	bool stepOne = scenarioTreeStepOne(scenario_origin,scenario_new,deleteFlag,scenarioP,scenarioP_new);//the same to the user interface
	bool stepTwo = scenarioTreeStepTwo(scenario_new,scenarioP_new);
	
	if(stepOne && stepTwo)
	{
		return true;	
	}
	return false;
}

int main()
{
	double scenario_origin[Pre_N][NT];
	ifstream input1("scenarios1.txt",ios::in);
	ifstream input2("scenarios2.txt",ios::in);
	ifstream input3("scenarios3.txt",ios::in);
	
	if(!input1)
	{
		cerr<<"scenarios1.txt is not exist!"<<endl;
		return 1;
	}
	for(int i = 0;i < Pre_N;++i)
	{
		for(int t = 0;t < NT;++t)
		{
			input1>>wind1[i][t];
		}
	}
	
	if(!input2)
	{
		cerr<<"scenarios2.txt is not exist!"<<endl;
		return 1;
	}
	for(int i = 0;i < Pre_N;++i)
	{
		for(int t = 0;t < NT;++t)
		{
			input2>>wind2[i][t];
		}
	}
	
	if(!input3)
	{
		cerr<<"scenarios3.txt is not exist!"<<endl;
		return 1;
	}
	for(int i = 0;i < Pre_N;++i)
	{
		for(int t = 0;t < NT;++t)
		{
			input3>>wind3[i][t];
		}
	}

	bool deleteFlag[Pre_N];
	double scenarioP[Pre_N];//缩减前后的概率
	double scenarioP_new[Pre_N]={0.0};

	for(int i = 0;i < Pre_N;++i)
	{
		scenarioP[i]=1.0/Pre_N;
	}

	for(int i = 0;i < Pre_N;++i)
	{
		deleteFlag[i] = false;//all not be deleted
	}

	for(int i = 0;i < Pre_N;++i)//input data
	{
		for(int t = 0;t < NT;++t)
		{
			scenario_origin[i][t] = (wind1[i][t] + wind2[i][t] +wind3[i][t]);
		}
	}

	double scenario_new[Aft_N][NT];
	ofstream output("new_scenario.txt",ios::ate);
	ofstream output1("new_wind1.txt",ios::ate);
	ofstream output2("new_wind2.txt",ios::ate);
	ofstream output3("new_wind3.txt",ios::ate);
	
	time_t nowTime=time(0);//print the time
	struct tm* nowTimeStruct=localtime(&nowTime);
	output<<"系统当前时间："<<1900+nowTimeStruct->tm_year<<"."//year
		<<nowTimeStruct->tm_mon + 1<<"."//month
			<<nowTimeStruct->tm_mday<<"  "//day
				<<nowTimeStruct->tm_hour<<":"//hour
					<<nowTimeStruct->tm_min<<":"//min
						<<nowTimeStruct->tm_sec<<endl;//scecond

	//scenario reduction function
	bool flag = scenario_reduce(scenario_origin,scenario_new,deleteFlag,scenarioP,scenarioP_new);
	if(!flag)
	{
		cerr<<"the process failed..."<<endl;
		return 1;
	}
	
	for(int i = 0;i < Aft_N;++i)
	{
		for(int t = 0;t < NT;++t)
		{
			output<<scenario_new[i][t]<<"   ";
			output1<<wind1[fAftToPre[i]][t]<<"   ";
			output2<<wind2[fAftToPre[i]][t]<<"   ";
			output3<<wind3[fAftToPre[i]][t]<<"   ";
		}
		output<<endl;
		output1<<endl;
		output2<<endl;
		output3<<endl;
	}

	input1.close();
	input2.close();
	input3.close();
	output.close();
	output1.close();
	output2.close();
	output3.close();
	return 0;
}
