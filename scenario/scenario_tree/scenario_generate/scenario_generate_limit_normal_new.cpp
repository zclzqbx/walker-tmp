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

const int NT         = 24;//误差场景时段
const int Num        = 1000;//场景个数
const int nsample    = 20;//与正态分布随机数精度有关 
const double pref    = 300;//切割点
const double Wz      = 350;//总装机容量

double uniform_distribution()
{	//产生服从（-1，1）(由（0，1）-->(-1,1))均匀分布的随机数
	return ((double)rand() / RAND_MAX) * 2 - 1;		
}

long int ran()
{//用于产生一个随机数
	static long int a = 203;
	a = (a * 10001 + 19) % 18307;
	return a;
}

double normal_distribution(const int N = 20,const double b=0.0475,const double u=0)
{//产生1个正态分布随机数
	double normal(0);
	srand((unsigned)(time(NULL)) + ran());
	for(int i = 0;i < N;++i)
	{
		normal += uniform_distribution();
	}
	return normal * b * sqrt(3.0/N) + u;
}


int main()
{
	double Pre[NT],u = 0,b = 0.0475;//预测场景、均值、方差
		
	ifstream input("wind_power.txt",ios::in);
	if(!input)
	{
		cerr<<"lack of file wind_power.txt!!!"<<endl;
		return -1;
	}
	for(int i=0;i<NT;++i)
	{//读取风电出力预测值
		input>>Pre[i];
	}
	
	ofstream output("mc_scenarios.txt",ios::ate);	
	
	for(int j=0;j<Num;++j)
	{//需要生成Num个场景序列
		for(int k=0;k<NT;++k)
		{
			b = Pre[k] / 5 + Wz / 50 + k / 2;
			double tmp = Pre[k] + normal_distribution(nsample,b,u);//生成一个出力场景
			if(tmp > pref)
			{//对生成的场景添加出力上限
				tmp = pref;
			}
			
			if(tmp < 0)
			{
				tmp = 0;
			}
			output<<tmp<<"   ";//输出
		}
		output<<endl;
	}
	
	input.close();
	output.close();
	return 0;
}
