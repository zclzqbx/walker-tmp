/*             24ʱ�μ�������
 *author:�촫��
 *date:2016.01.15��ǰ
 *function:����24ʱ�εķ����������
 *
 */
#include<fstream>
#include<math.h>

using namespace std;

const int NT      = 24;
const double Wz   = 600;
const double Poff = 450;
const int D       = 500;//������������
const double Pi   = 3.14159;

double pre[NT];
double computeResult[NT];

ifstream input("pre_scenarios.txt",ios::in);
ofstream output("compute.txt",ios::ate);

double f(double x,double u,double b)//�����ܶȺ�������
{
	double coefficient = 1/(sqrt(2*Pi)*b);
	double index = -((x-u)*(x-u)/(2*b*b));
	return coefficient*exp(index);
}

double Integration(const double pre,const double u,const double b,const int N=1000)
{//����ʣ�����Ϊ�ֶ���
	double distance=(pre+D-Poff)/N;//�����С
	double sum(0);
	
	for(int i=0;i<N;++i)
	{
		double length;//���γ�
		double x1,x2;//������
		x1=Poff+i*distance;
		x2=Poff+(i+1)*distance;
		length=(f(x1,u,b)+f(x2,u,b))/2;//�������ظ������Կ����Ż�
		sum+=length*distance;
	}
	return sum;
}
//��ʵ������������������Ƶģ�ֻ�����Ǳ��������ĸı�
double mean(const double pre,const double u,const double b,const int N=300)//��minvalue~Poff������������Ϊ�ֶ���
{
	double distance=(Poff-pre+D)/N;//�����С
	double sum(0);
	
	for(int i=0;i<N;++i)
	{
		double length;//���γ�
		double x1,x2;//������
		x1 = pre - D + i*distance;
		x2 = pre - D + (i+1)*distance;
		length = (x1 * f(x1,u,b) + x2 * f(x2,u,b)) / 2;//�������ظ������Կ����Ż�
		sum += length*distance;
	}
	return sum;
}

double computeValue(double pre)
{
	double u=pre,b;
	b=pre / 5 + Wz / 50;
	return (mean(pre,u,b) + Poff*Integration(pre,u,b));//��������ε�����
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