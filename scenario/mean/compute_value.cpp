/*             24时段计算点计算
 *author:朱传林
 *date:2016.01.15以前
 *function:计算24时段的风电出力计算点
 *
 */
#include<fstream>
#include<math.h>

using namespace std;

const int NT      = 24;
const double Wz   = 600;
const double Poff = 450;
const int D       = 500;//横向搜索区间
const double Pi   = 3.14159;

double pre[NT];
double computeResult[NT];

ifstream input("pre_scenarios.txt",ios::in);
ofstream output("compute.txt",ios::ate);

double f(double x,double u,double b)//概率密度函数函数
{
	double coefficient = 1/(sqrt(2*Pi)*b);
	double index = -((x-u)*(x-u)/(2*b*b));
	return coefficient*exp(index);
}

double Integration(const double pre,const double u,const double b,const int N=1000)
{//求概率，参数为分段数
	double distance=(pre+D-Poff)/N;//区间大小
	double sum(0);
	
	for(int i=0;i<N;++i)
	{
		double length;//矩形长
		double x1,x2;//横坐标
		x1=Poff+i*distance;
		x2=Poff+(i+1)*distance;
		length=(f(x1,u,b)+f(x2,u,b))/2;//计算有重复，可以考虑优化
		sum+=length*distance;
	}
	return sum;
}
//其实求期望和求概率是类似的，只不过是被积函数的改变
double mean(const double pre,const double u,const double b,const int N=300)//求（minvalue~Poff）期望，参数为分段数
{
	double distance=(Poff-pre+D)/N;//区间大小
	double sum(0);
	
	for(int i=0;i<N;++i)
	{
		double length;//矩形长
		double x1,x2;//横坐标
		x1 = pre - D + i*distance;
		x2 = pre - D + (i+1)*distance;
		length = (x1 * f(x1,u,b) + x2 * f(x2,u,b)) / 2;//计算有重复，可以考虑优化
		sum += length*distance;
	}
	return sum;
}

double computeValue(double pre)
{
	double u=pre,b;
	b=pre / 5 + Wz / 50;
	return (mean(pre,u,b) + Poff*Integration(pre,u,b));//求的是两段的期望
}
int main()
{
	for(int i=0;i<NT;++i)
	{
		input>>pre[i];
	}
	
	for(int i=0;i<NT;++i)
	{
		computeResult[i] = computeValue(pre[i]);
		output<<computeResult[i]<<"   ";
	}
	
	input.close();
	output.close();
	return 0;
}