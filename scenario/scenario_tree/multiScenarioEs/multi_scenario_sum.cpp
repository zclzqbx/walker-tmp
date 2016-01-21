/*
 *author:朱传林
 *date:2015.12.3
 *function：将多个风电场的场景组合成一个场景(求和)
 *
 */
#include <iostream>
#include <fstream>

using namespace std;
 
const int WINDNUM = 3;
const int NT      = 24;
const int SCENNUM = 8;

//之所以分开存放是考虑到系统内存的限制
ifstream input("scenarios.txt",ios::in);//风电场1

ofstream output("merge_scenarios.txt",ios::ate);//综合场景按序存放

int main()
{
	double wind1[NT],wind2[NT],wind3[NT];
	if(!input)
	{
		cerr<<"scenarios.txt is not exist!"<<endl;
	}

	
	for(int sce = 0;sce < SCENNUM;++sce)
	{
		for(int t = 0;t < NT;++t)
		{
			input>>wind1[t];
		}
		
		for(int t = 0;t < NT;++t)
		{
			input>>wind2[t];
		}
		
		for(int t = 0;t < NT;++t)
		{
			input>>wind3[t];
		}
		
		for(int t = 0;t < NT;++t)
		{//输出风电场1数据
			output<<wind1[t]+wind2[t]+wind3[t]<<"  ";
		}
		output<<endl;
	}
	return 0;
}

