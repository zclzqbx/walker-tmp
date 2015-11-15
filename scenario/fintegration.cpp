/*             	�����⺯���Ķ�����
 *1�������������ԭʼ�Ļ�����ⷽ��������ƽ�����
 *2��funΪ����������fintegrationΪ���֣�����Ĭ������Ϊ[-1000,1000]�������Լ�ָ����
 *3����������ֶ���N�غ����ȣ�NԽ�󣬾���Խ�ߣ���������Խ��
 */
#include<iostream>
#include<math.h>

using namespace std;
typedef double (*fpoint)(double x);// ����ָ���÷�

double fun1(double x)//��׼��̬�ֲ�
{
	const int u=0;
	const int b=1;
	const double Pi=3.14159;
	double coefficient=1/(sqrt(2*Pi)*b);
	double index=-((x-u)*(x-u)/(2*b*b));
	return coefficient*exp(index);
}

double fun2(double x)//���κ���
{
	return 3*x*x+4*x;
}

double fun3(double x)//ָ������
{
	return exp(x);
}

double fintegration(fpoint f,const int N=5000,const double minvalue=-1000,const double maxvalue=1000)
{//Ĭ�Ϸֶ���Ϊ5000��Ĭ�ϻ�������Ϊ[-1000,1000]����������Ϊf���˴�fֻ���ǲ���Ϊdouble,����ֵΪdouble��һԪ����
//���Ҫ�ó������ֲ�Ը�ǿ�����Բ��÷º���
	double distance=(maxvalue-minvalue)/N;//�����С
	double sum(0);
	for(int i=0;i<N;++i)
	{
		double length;//���γ�
		double x1,x2;//������
		x1=minvalue+i*distance;
		x2=minvalue+(i+1)*distance;
		length=(f(x1)+f(x2))/2;//����ֵʱ���ظ������Կ����Ż�
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
