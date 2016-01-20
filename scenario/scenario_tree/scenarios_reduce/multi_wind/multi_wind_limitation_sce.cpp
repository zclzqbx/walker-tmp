/*              3��糡����ϳ�������
 *author:�촫��
 *date:2015.10.26
 *function:������ĸ�����糡����ȫ����
 *problem:1���˳�����Ҫȱ������ֻ�ʺ�3��糡������ֲ�Բ
 *
 */
#include <iostream>
#include <fstream>
#include <time.h>

using namespace std;

const int NUM = 2;//ÿ����糡�ĳ�������
const int NT  = 24;//�ܵ�ʱ�θ���
const int NW  = 3;//��糡����

ifstream input_scenario1("scenarios1.txt",ios::in);
ifstream input_scenario2("scenarios2.txt",ios::in);
ifstream input_scenario3("scenarios3.txt",ios::in);
ofstream output("merge_scenarios.txt",ios::ate);

int main()
{
	clock_t start,finish;
	double totaltime;
	start=clock();
	
	double wind_scenario1[NUM][NT],wind_scenario2[NUM][NT],wind_scenario3[NUM][NT];
	//��������
	if(!input_scenario1)
	{
		cerr<<"scenario1.txt is not exist!!!"<<endl;
		return -1;
	}	
	if(!input_scenario2)
	{
		cerr<<"scenario2.txt is not exist!!!"<<endl;
		return -1;
	}
	if(!input_scenario3)
	{
		cerr<<"scenario3.txt is not exist!!!"<<endl;
		return -1;
	}
	
	for(int i = 0;i < NUM;++i)
	{
		for(int t = 0;t < NT;++t)
		{
			input_scenario1>>wind_scenario1[i][t];
			input_scenario2>>wind_scenario2[i][t];
			input_scenario3>>wind_scenario3[i][t];
		}
	}
	input_scenario1.close();
	input_scenario2.close();
	input_scenario3.close();
	
	double temp_storage[NW][NT];
	for(int i=0;i<NUM;++i)//����ѭ��ʵ�����
	{
		for(int t=0;t<NT;++t)	
			temp_storage[0][t]=wind_scenario1[i][t];
		
		for(int j=0;j<NUM;++j)
		{
			for(int t=0;t<NT;++t)	
				temp_storage[1][t]=wind_scenario2[j][t];
			
			for(int k=0;k<NUM;++k)
			{
				for(int t=0;t<NT;++t)	
					temp_storage[2][t]=wind_scenario3[k][t];
				
				for(int w=0;w<NW;++w)//�����������
				{
					for(int t=0;t<NT;++t)
					{
						output<<temp_storage[w][t]<<"   ";
					}
					output<<endl;
				}
				output<<endl;
			}
		}
	}
	
	finish = clock();
	totaltime = (double)(finish-start) / CLOCKS_PER_SEC;
	output<<"totaltime="<<totaltime<<endl;
	
	return 0;
}