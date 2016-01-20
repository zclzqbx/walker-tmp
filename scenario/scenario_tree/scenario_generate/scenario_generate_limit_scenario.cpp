/*              ���޳�����
 *author:�촫��
 *date:2016-01-13
 *function:����Ԥ��ֵ��������¼��޳���
 *
 */
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <cstdlib>
using namespace std;

#define POS pos_090//�����λ�㣬��������õ��ǶԳƷ�λ��

const int NT         = 24;//����ʱ��
const double pref    = 160;//�и��
const double Wz      = 200;//��װ������

//��׼��̬�ֲ���λ��
const double pos_095 = 1.6449;
const double pos_090 = 1.2816;
const double pos_085 = 1.0364;
const double pos_080 = 0.8416;



double topPos(const double b=0.0475,const double u=0)
{//�ϼ��޳���
	return (b * POS + u);
}

double bottomPos(const double b=0.0475,const double u=0)
{//�¼��޳���
	return (-1 * b * POS + u);
}

int main()
{
	double Pre[NT],u = 0,b = 0.0475;//Ԥ�ⳡ������ֵ������
	
	ifstream input("wind_power.txt",ios::in);
	if(!input)
	{
		cerr<<"lack of file wind_power.txt!!!"<<endl;
		return -1;
	}
	for(int i=0;i<NT;++i)
	{//��ȡ������Ԥ��ֵ
		input>>Pre[i];
	}
	
	ofstream output("mc_scenarios.txt",ios::ate);
	for(int k = 0;k < NT;++k)//�Ϸ�λ��
	{
		b = Pre[k] / 5 + Wz / 50;
		double tmp = Pre[k] + topPos(b,u);//����һ����������
		if(tmp > pref)
		{//�����ɵĳ�����ӳ�������
			tmp = pref;
		}
		output<<tmp<<"   ";//���
	}
	output<<endl;
	
	for(int k = 0;k < NT;++k)//�·�λ��
	{
		b = Pre[k] / 5 + Wz / 50;
		double tmp = Pre[k] + bottomPos(b,u);//����һ����������
		output<<tmp<<"   ";//���
	}
	output<<endl;
	
	input.close();
	output.close();
	return 0;
}
