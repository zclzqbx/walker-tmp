/*              ģ�����ؿ��巽��
 *1)����1����̬�ֲ��������Ҫnsample�����ӣ�-1��1������
 *�ֲ����������nsampleԽ�󣬾���Խ�ߣ���ͬʱҲ��������ӵ�
 *���������ο����ס���̬�ֲ�����������ɷ���_��Ϊ�񡱡�
 *
 */
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <cstdlib>
using namespace std;

const int NT=24;//����ʱ��
const int Num=200;//��������
const int nsample=20;

double uniform_distribution()
{
	return ((double)rand()/RAND_MAX)*2-1;
	//�������ӣ�-1��1��(�ɣ�0��1��-->(-1,1))���ȷֲ��������	
}

long int ran()
{
	static long int a=203;
	a=(a*10001+19)%18307;
	return a;
}

double normal_distribution(const int N=20,const double b=0.0475,const double u=0)
{//����1����̬�ֲ������
	double normal(0);
	srand((unsigned)(time(NULL))+ran());
	for(int i=0;i<N;++i)
	{
		normal+=uniform_distribution();
	}
	return normal*b*sqrt(3.0/N)+u;
}


int main()
{
	double Pre[NT],u=0,b=0.0475;//Ԥ�ⳡ������ֵ������
	double Wz=350;//��װ������
	ifstream input("wind_power.txt",ios::in);
	if(!input)
	{
		cout<<"lack of file wind_power.txt!!!"<<endl;
		return 1;
	}
	ofstream output("mc_scenarios.txt",ios::ate);
	for(int i=0;i<NT;++i)
	{
		input>>Pre[i];
	}	
	for(int j=0;j<Num;++j)
	{
		for(int k=0;k<NT;++k)
		{
			b=Pre[k]/5+Wz/50;
			output<<Pre[k]+normal_distribution(nsample,b,u)<<"   ";
		}
		output<<endl;
	}	
	input.close();
	output.close();
	return 0;
}
