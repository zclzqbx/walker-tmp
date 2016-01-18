/*
 *author:朱传林
 *date:2015.11.30
 *function：将多个风电场的场景组合成一个场景
 *
 */
#include <iostream>
#include <fstream>

using namespace std;
 
const int WINDNUM = 3;
const int NT      = 24;
const int SCENNUM = 20;

//之所以分开存放是考虑到系统内存的限制
ifstream input1("scenarios1.txt",ios::in);//风电场1
ifstream input2("scenarios2.txt",ios::in);//风电场2
ifstream input3("scenarios3.txt",ios::in);//风电场3

ofstream output("merge_scenarios.txt",ios::ate);//综合场景按序存放

int main()
{
	double wind1[NT],wind2[NT],wind3[NT];
	if(!input1)
	{
		cerr<<"scenarios1.txt is not exist!"<<endl;
	}
	if(!input2)
	{
		cerr<<"scenarios2.txt is not exist!"<<endl;
	}
	if(!input3)
	{
		cerr<<"scenarios3.txt is not exist!"<<endl;
	}
	
	for(int sce = 0;sce < SCENNUM;++sce)
	{
		for(int t = 0;t < NT;++t)
		{//读取风电场数据，从不同的文件中读取
			input1>>wind1[t];
			input2>>wind2[t];
			input3>>wind3[t];
		}
		
		for(int t = 0;t < NT;++t)
		{//输出风电场1数据
			output<<wind1[t]<<"  ";
		}
		output<<endl;
		
		for(int t = 0;t < NT;++t)
		{//输出风电场2数据
			output<<wind2[t]<<"  ";
		}
		output<<endl;
		
		for(int t = 0;t < NT;++t)
		{//输出风电场3数据
			output<<wind3[t]<<"  ";
		}
		output<<endl<<endl;
	}
	return 0;
}

