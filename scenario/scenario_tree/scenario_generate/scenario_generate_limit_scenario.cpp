/*              极限场景法
 *author:朱传林
 *date:2016-01-13
 *function:输入预测值，输出上下极限场景
 *
 */
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <cstdlib>
using namespace std;

#define POS pos_090//定义分位点，本程序采用的是对称分位点

const int NT         = 24;//误差场景时段
const double pref    = 160;//切割点
const double Wz      = 200;//总装机容量

//标准正态分布分位点
const double pos_095 = 1.6449;
const double pos_090 = 1.2816;
const double pos_085 = 1.0364;
const double pos_080 = 0.8416;



double topPos(const double b=0.0475,const double u=0)
{//上极限场景
	return (b * POS + u);
}

double bottomPos(const double b=0.0475,const double u=0)
{//下极限场景
	return (-1 * b * POS + u);
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
	for(int k = 0;k < NT;++k)//上分位点
	{
		b = Pre[k] / 5 + Wz / 50;
		double tmp = Pre[k] + topPos(b,u);//生成一个出力场景
		if(tmp > pref)
		{//对生成的场景添加出力上限
			tmp = pref;
		}
		output<<tmp<<"   ";//输出
	}
	output<<endl;
	
	for(int k = 0;k < NT;++k)//下分位点
	{
		b = Pre[k] / 5 + Wz / 50;
		double tmp = Pre[k] + bottomPos(b,u);//生成一个出力场景
		output<<tmp<<"   ";//输出
	}
	output<<endl;
	
	input.close();
	output.close();
	return 0;
}
