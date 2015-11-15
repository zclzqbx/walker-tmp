
/*              模拟蒙特卡洛方法
 *1)产生1个正态分布随机数需要nsample个服从（-1，1）均匀
 *分布的随机数，nsample越大，精度越高，但同时也会带来附加的
 *计算量。参考文献“正态分布随机数的生成方法_孙为民”。
 *
 */
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <cstdlib>
using namespace std;

const int nsample=20;
const int NUM=3000;
const int SegNum=60;

double uniform_distribution()
{
	return ((double)rand()/RAND_MAX)*2-1;
	//产生服从（-1，1）(由（0，1）-->(-1,1))均匀分布的随机数	
}

long int ran()//随机数种子
{
	static long int a=203;
	a=(a*10001+19)%18307;
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
	double u=300,b=30;//预测场景、均值、方差
	int seg[SegNum]={0};
	
	ofstream output("mc_scenarios.txt",ios::ate);
	
	double normal_value(0);
	for(int k=0;k<NUM;++k)
	{
		normal_value=normal_distribution(nsample,b,u);
		output<<normal_value<<"   ";
		int tmp=normal_value/10;
		seg[tmp]++;
	}
	output<<endl;	
	output.close();

	for(int i=0;i<SegNum;++i)
	{
		cout<<seg[i]<<endl;
	}
	return 0;
}
