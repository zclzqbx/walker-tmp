/*             	非规则分布的期望
 *author:朱传林
 *date:2016.01.15以前
 *function:求单时段的非规则分布的期望
 */
#include<iostream>
#include<math.h>

using namespace std;

const double Wz       = 80;//总装机
const double Pre      = 50;//预测值
const double Poff     = 60;//参考点
const double minvalue = Pre-50;//积分区间
const double maxvalue = Pre+50;
const double Pi       = 3.14159;

//概率密度函数函数
double f(double x,double u,double b)
{
	double coefficient = 1/(sqrt(2*Pi)*b);
	double index       = -((x-u)*(x-u)/(2*b*b));
	return coefficient*exp(index);
}

//求切断点的概率，参数为分段数
double Integration(const double u,const double b,const int N = 300)
{
	double distance = (maxvalue - Poff) / N;//区间大小
	double sum(0);
	
	for(int i = 0;i < N;++i)
	{
		double length;//矩形长
		double x1,x2;//横坐标
		x1 = Poff + i*distance;//横坐标起点是参考点
		x2 = Poff + (i+1)*distance;
		length = (f(x1,u,b)+f(x2,u,b))/2;//计算有重复，可以考虑优化
		sum += length*distance;//求面积
	}
	return sum;
}

//其实求期望和求概率是类似的，只不过是被积函数的改变
double mean(const double u,const double b,const int N = 300)//求（minvalue~Poff）期望，参数为分段数
{
	double distance = (Poff - minvalue)/N;//区间大小
	double sum(0);
	
	for(int i=0;i<N;++i)
	{
		double length;//矩形长
		double x1,x2;//横坐标
		x1 = minvalue + i*distance;//横坐标起点是最小值点
		x2 = minvalue + (i+1)*distance;
		length = (x1*f(x1,u,b)+x2*f(x2,u,b))/2;//计算有重复，可以考虑优化
		sum += length*distance;
	}
	return sum;
}

int main()
{
	double u,b;
	u=Pre;
	b=Pre/5+Wz/50;
	//求0~60之间的积分（实际装机为80，风电出力不可能超过这一值）
	double integ_value(0);
	cout<<Integration(u,b)<<endl;
	integ_value=mean(u,b)+Poff*Integration(u,b);
	cout<<integ_value<<endl;
	return 0;
}