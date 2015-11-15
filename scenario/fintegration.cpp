/*             	求任意函数的定积分
 *1）本程序采用最原始的积分求解方法：面积逼近法；
 *2）fun为被积函数，fintegration为积分，积分默认区间为[-1000,1000]，可以自己指定；
 *3）积分区间分段数N关乎精度，N越大，精度越高，但计算量越大；
 */
#include<iostream>
#include<math.h>

using namespace std;
typedef double (*fpoint)(double x);// 函数指针用法

double fun1(double x)//标准正态分布
{
	const int u=0;
	const int b=1;
	const double Pi=3.14159;
	double coefficient=1/(sqrt(2*Pi)*b);
	double index=-((x-u)*(x-u)/(2*b*b));
	return coefficient*exp(index);
}

double fun2(double x)//二次函数
{
	return 3*x*x+4*x;
}

double fun3(double x)//指数函数
{
	return exp(x);
}

double fintegration(fpoint f,const int N=5000,const double minvalue=-1000,const double maxvalue=1000)
{//默认分段数为5000，默认积分区间为[-1000,1000]，被积函数为f，此处f只能是参数为double,返回值为double的一元函数
//如果要让程序可移植性更强，可以采用仿函数
	double distance=(maxvalue-minvalue)/N;//区间大小
	double sum(0);
	for(int i=0;i<N;++i)
	{
		double length;//矩形长
		double x1,x2;//横坐标
		x1=minvalue+i*distance;
		x2=minvalue+(i+1)*distance;
		length=(f(x1)+f(x2))/2;//计算值时有重复，可以考虑优化
		sum+=length*distance;
	}
	return sum;
}

int main()
{	
	cout<<fintegration(fun1)<<endl;
	cout<<fintegration(fun2,1000,0,10)<<endl;
	cout<<fintegration(fun3,1000,0,200)<<endl;
	return 0;
}
