#include<iostream>
#include<fstream>

using namespace std;

const int Pre_N=200;//the number of scenario before reduction
const int Aft_N=20;//the number of scenario after reduction
const int NT=24;//time

const int J[NT]={0,1,0,0,1,0,2,0,0,1,0,0,2,0,2,0,0,3,0,3,0,4,0,180};//各个时段删除的节点数

bool scenario_reduce(double scenario_origin[Pre_N][NT],double scenario_new[Aft_N][NT])
{


	return true;	
}


int main()
{
	double scenario_origin[Pre_N][NT];
	ifstream input("oringin_scenario.txt",ios::in);

	for(int i=0;i<Pre_N;++i)//input data
	{
		for(int t=0;t<24;++t)
		{
			input>>scenario_origin[i][t];
		}
	}

	double scenario_new[Aft_N][NT];
	ofstream output("new_scenario.txt",ios::ate);

	bool flag=scenario_reduce(scenario_origin,scenario_new);//scenario reduction function
	if(!bool)
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
