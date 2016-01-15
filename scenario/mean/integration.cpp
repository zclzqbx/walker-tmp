/*             	����̬�ֲ����ֺ���
 *author:�촫��
 *date:2016.01.15��ǰ
 *function:1�������������ԭʼ�Ļ�����ⷽ��������ƽ�����
 */
#include<fstream>
#include<iostream>
#include<math.h>

using namespace std;

const int u=0;//����
const int b=1;//����
const int N=500;//���ֶַ���,�ֶ�Խ�࣬���Խ��ȷ
const double minvalue=-100;//��������
const double maxvalue=100;
const double Pi=3.14159;

double f(double x)//�����ܶȺ���ֵ
{
	double coefficient=1/(sqrt(2*Pi)*b);
	double index=-((x-u)*(x-u)/(2*b*b));
	return coefficient*exp(index);
}

int main()
{
	double distance=(maxvalue-minvalue)/N;//�����С
	
	double sum(0);
	
	for(int i=0;i<N;++i)
	{
		double length;//���γ�
		double x1,x2;//������
		x1=minvalue+i*distance;
		x2=minvalue+(i+1)*distance;
		length=(f(x1)+f(x2))/2;//�������ظ������Կ����Ż�
		sum+=length*distance;
	}
	cout<<"The result of integration:"<<sum<<endl;
	return 0;
}
