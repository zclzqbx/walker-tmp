#include<iostream>
#include<time.h>
#include<fstream>
#include<math.h>
#include<vector>

using namespace std;

const int Pre_N=200;//the number of scenario before reduction
const int Aft_N=20;//the number of scenario after reduction
const int NT=24;//time

const int J[NT]={0,1,0,0,1,0,2,0,0,1,0,0,2,0,2,0,0,3,0,3,0,4,0,180};//各个时段删除的节点数

double c_distance(double s1[NT],double s2[NT],int T)
{
	if(s1==NULL || s2==NULL || T < 0)
		return -1;
	double result=0;
	for(int t=0;t<T;++t)
	{
		result+=sqrt((s1[t]-s2[t])*(s1[t]-s2[t]));
	}
	return result;
}

bool scenarioTreeStepTwo(double scenario_new[Aft_N][NT],double scenarioP_new[Aft_N])
{
	if(scenario_new==NULL || scenarioP_new==NULL)
		return false;
	
	double pt[Aft_N][NT];//每个时段的概率
	bool deFlag[Aft_N];//删除标志
	vector< vector<int> > hide(Aft_N,vector<int>(Aft_N,-1));//背后需要更改的场景
	
	for(int h=0;h<Aft_N;++h)
	{//初始化
		pt[h][NT-1]=scenarioP_new[h];
		deFlag[h]=false;
	}
	
	for(int t=NT-2;t>=0;--t)
	{
		if(J[t]>0)//如果本阶段有缩减需求
		{
			int deNum=0;
			while(deNum < J[t])
			{
				vector<double> z_l(Aft_N,-1);
				vector< vector<double> > c_l(Aft_N,vector<double>(Aft_N,-1));
				
				for(int i=0;i<Aft_N;++i)//Ckj
				{
					if(!deFlag[i])
					{
						double tmp=10000000.0;
						for(int k=0;k<Aft_N;++k)
						{
							if(k==i || deFlag[k])
							{
								for(int h=0;h<Aft_N;++h)
								{
									if(h!=i && !deFlag[h])
									{
										double tmp1=c_distance(scenario_new[k],scenario_new[h],t+1);
										if(tmp > tmp1)
											tmp=tmp1;
									}
								}
								c_l[k][i]=tmp;
							}
						}
					}
				}

				for(int i=0;i<Aft_N;++i)//Zl
				{
					double tmp=0;
					if(!deFlag[i])
					{
						for(int k=0;k<Aft_N;++k)
						{
							if(k==i || deFlag[k])
							{
								tmp+=pt[k][t+1]*c_l[k][i];
							}
						}
						z_l[i]=tmp;
					}
				}
				
				int deleteScenarioIndex=-1;//the next scenario being deleted
				double deleteScenarioValue=10000000.0;
				for(int i=0;i<Aft_N;++i)
				{
					if(!deFlag[i])
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

				for(int i=0;i<Aft_N;++i)//the nearest scenario
				{
					if(!deFlag[i] && i!=deleteScenarioIndex)
					{
						double tmp=c_distance(scenario_new[deleteScenarioIndex],scenario_new[i],t+1);
						if(minDistance>tmp)
						{
							minDistance=tmp;
							mergeScenariosIndex=i;
						}
					}
				}
				
				pt[mergeScenariosIndex][t+1]+=pt[deleteScenarioIndex][t+1];//合并概率
				
				int hide_index=0;				
				int hide_index_new=0;
				int theScenarioNeedToBeModify=deleteScenarioIndex;
				
				while(hide[mergeScenariosIndex][hide_index_new]!=-1)
					hide_index_new++;
				
				hide[mergeScenariosIndex][hide_index_new]=deleteScenarioIndex;
				hide_index_new++;
				
				do{//生成场景树
					for(int t1=0;t1<=t;++t1)
					{
						scenario_new[theScenarioNeedToBeModify][t1]=scenario_new[mergeScenariosIndex][t1];
					}
					
					if(hide_index>=Aft_N || hide[deleteScenarioIndex][hide_index]==-1)
					{
						break;
					}
					else
					{
						theScenarioNeedToBeModify=hide[deleteScenarioIndex][hide_index];
						hide[mergeScenariosIndex][hide_index_new]=theScenarioNeedToBeModify;
						hide_index++;
						hide_index_new++;
					}	
					
				}while(1);
				
				deFlag[deleteScenarioIndex]=true;
				for(int i=0;i<Aft_N;++i)
				{
					pt[i][t]=pt[i][t+1];
				}				
				deNum++;
			}
		}
	}	
	return true;
}

int main()
{
	ifstream input("new_scenario.txt",ios::in);
	if(!input)
	{
		cout<<"new_scenario is not exist!"<<endl;
		return 1;
	}
	double scenario_new[Aft_N][NT];//input data
	double scenarioP_new[Aft_N];
	for(int i=0;i<Aft_N;++i)
	{
		input>>scenarioP_new[i];//先输入场景
		for(int t=0;t<NT;++t)
		{
			input>>scenario_new[i][t];
		}
	}
	
	ofstream output("tree.txt",ios::ate);	
	time_t nowTime=time(0);
	struct tm* nowTimeStruct=localtime(&nowTime);
	
	output<<"The System Time:"<<1900+nowTimeStruct->tm_year<<"."
		<<nowTimeStruct->tm_mon<<"."
			<<nowTimeStruct->tm_mday<<"  "//day
				<<nowTimeStruct->tm_hour<<":"//hour
					<<nowTimeStruct->tm_min<<":"//min
						<<nowTimeStruct->tm_sec<<endl;//scecond
						
	bool stepTwo=scenarioTreeStepTwo(scenario_new,scenarioP_new);//计算
	if(stepTwo==false)
		return -1;
	
	for(int i=0;i<Aft_N;++i)
	{
		for(int t=0;t<NT;++t)
		{
			output<<scenario_new[i][t]<<"   ";
		}
		output<<endl;
	}	
	return 0;
}