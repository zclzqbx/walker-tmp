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

const int NT      = 24;//����ʱ��
const int Num     = 10;//��������
const int nsample = 20;
const double Poff = 5500;//�и��
const double Wz   = 350;//�����װ������

double uniform_distribution()
{
	//�������ӣ�-1��1��(�ɣ�0��1��-->(-1,1))���ȷֲ��������
	return ((double)rand()/RAND_MAX)*2-1;		
}

long int ran()//���������
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
	double Pre[NT],u=0,b;//Ԥ�ⳡ������ֵ������
	
	ifstream input("forecast_scenarios.txt",ios::in);
	if(!input)
		cout<<"lack of file forecast_scenarios.txt!!!"<<endl;
	ofstream output("mc_scenarios.txt",ios::ate);
	
	for(int i=0;i<NT;++i)
	{
		input>>Pre[i];
	}	

	for(int j=0;j<Num;++j)
	{
		for(int k=0;k<NT;++k)
		{
			b=Pre[k]/10+Wz/50;
			double normal_value(0);
			normal_value=Pre[k]+normal_distribution(nsample,b,u);
			if(normal_value>=Poff)
				normal_value=Poff;
			output<<normal_value<<"   ";
		}
		output<<endl;
	}
	
	input.close();
	output.close();
	return 0;
}
