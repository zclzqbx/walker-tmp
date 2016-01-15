/*             	�ǹ���ֲ�������
 *author:�촫��
 *date:2016.01.15��ǰ
 *function:��ʱ�εķǹ���ֲ�������
 */
#include<iostream>
#include<math.h>

using namespace std;

const double Wz       = 80;//��װ��
const double Pre      = 50;//Ԥ��ֵ
const double Poff     = 60;//�ο���
const double minvalue = Pre-50;//��������
const double maxvalue = Pre+50;
const double Pi       = 3.14159;

//�����ܶȺ�������
double f(double x,double u,double b)
{
	double coefficient = 1/(sqrt(2*Pi)*b);
	double index       = -((x-u)*(x-u)/(2*b*b));
	return coefficient*exp(index);
}

//���жϵ�ĸ��ʣ�����Ϊ�ֶ���
double Integration(const double u,const double b,const int N = 300)
{
	double distance = (maxvalue - Poff) / N;//�����С
	double sum(0);
	
	for(int i = 0;i < N;++i)
	{
		double length;//���γ�
		double x1,x2;//������
		x1 = Poff + i*distance;//����������ǲο���
		x2 = Poff + (i+1)*distance;
		length = (f(x1,u,b)+f(x2,u,b))/2;//�������ظ������Կ����Ż�
		sum += length*distance;//�����
	}
	return sum;
}

//��ʵ������������������Ƶģ�ֻ�����Ǳ��������ĸı�
double mean(const double u,const double b,const int N = 300)//��minvalue~Poff������������Ϊ�ֶ���
{
	double distance = (Poff - minvalue)/N;//�����С
	double sum(0);
	
	for(int i=0;i<N;++i)
	{
		double length;//���γ�
		double x1,x2;//������
		x1 = minvalue + i*distance;//�������������Сֵ��
		x2 = minvalue + (i+1)*distance;
		length = (x1*f(x1,u,b)+x2*f(x2,u,b))/2;//�������ظ������Կ����Ż�
		sum += length*distance;
	}
	return sum;
}

int main()
{
	double u,b;
	u=Pre;
	b=Pre/5+Wz/50;
	//��0~60֮��Ļ��֣�ʵ��װ��Ϊ80�������������ܳ�����һֵ��
	double integ_value(0);
	cout<<Integration(u,b)<<endl;
	integ_value=mean(u,b)+Poff*Integration(u,b);
	cout<<integ_value<<endl;
	return 0;
}