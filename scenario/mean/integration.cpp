/*             	求正态分布积分函数
 *author:朱传林
 *date:2016.01.15以前
 *function:1）本程序采用最原始的积分求解方法：面积逼近法；
 */
#include<fstream>
#include<iostream>
#include<math.h>

using namespace std;

const int u=0;//期望
const int b=1;//方差
const int N=500;//积分分段数,分段越多，结果越精确
const double minvalue=-100;//积分区间
const double maxvalue=100;
const double Pi=3.14159;

double f(double x)//概率密度函数值
{
	double coefficient=1/(sqrt(2*Pi)*b);
	double index=-((x-u)*(x-u)/(2*b*b));
	return coefficient*exp(index);
}

int main()
{
	double distance=(maxvalue-minvalue)/N;//区间大小
	
	double sum(0);
	
	for(int i=0;i<N;++i)
	{
		double length;//矩形长
		double x1,x2;//横坐标
		x1=minvalue+i*distance;
		x2=minvalue+(i+1)*distance;
		length=(f(x1)+f(x2))/2;//计算有重复，可以考虑优化
		sum+=length*distance;
	}
	cout<<"The result of integration:"<<sum<<endl;
	return 0;
}
