/*              对抽样结果进行统计分析
 *author:朱传林
 *date:2016.01.15以前
 *function:对正态分布随机数进行分段统计，统计其各个区间出现的次数，用柱状图模型拟其概率密度函数
 *
 */
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <cstdlib>
using namespace std;

const int nsample = 20;//正态分布随机数精度
const int NUM     = 3000;//场景数
const int SegNum  = 120;//区间数

double uniform_distribution()
{
	//产生服从（-1，1）(由（0，1）-->(-1,1))均匀分布的随机数
	return ((double)rand()/RAND_MAX) * 2 - 1;		
}

long int ran()//随机数种子
{
	static long int a = 203;
	a = (a * 10001 + 19) % 18307;
	return a;
}

double normal_distribution(const int N=20,const double b=0.0475,const double u=0)
{//产生1个正态分布随机数
	double normal(0);
	srand((unsigned)(time(NULL))+ran());
	for(int i=0;i<N;++i)
	{
		normal+=uniform_distribution();
	}
	return normal*b*sqrt(3.0/N)+u;
}

int main()
{
	double u=300,b=30;//均值、方差
	int seg[SegNum]={0};
	
	ofstream output("mc_scenarios.txt",ios::ate);
	ofstream output1("normal_statistics.txt",ios::ate);
	
	double normal_value(0);
	for(int k=0;k<NUM;++k)
	{
		normal_value=normal_distribution(nsample,b,u);//产生一个服从正态分布的随机数
		output<<normal_value<<"   ";
		int tmp = normal_value / 5;
		seg[tmp]++;//对应区间的个数增一
	}
	output<<endl;	
	output.close();

	for(int i=0;i<SegNum;++i)
	{
		output1<<seg[i]<<endl;
	}
	output1.close();
	return 0;
}
