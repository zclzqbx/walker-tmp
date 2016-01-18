/*
 *author:�촫��
 *date:2015.11.30
 *function���������糡�ĳ�����ϳ�һ������
 *
 */
#include <iostream>
#include <fstream>

using namespace std;
 
const int WINDNUM = 3;
const int NT      = 24;
const int SCENNUM = 20;

//֮���Էֿ�����ǿ��ǵ�ϵͳ�ڴ������
ifstream input1("scenarios1.txt",ios::in);//��糡1
ifstream input2("scenarios2.txt",ios::in);//��糡2
ifstream input3("scenarios3.txt",ios::in);//��糡3

ofstream output("merge_scenarios.txt",ios::ate);//�ۺϳ���������

int main()
{
	double wind1[NT],wind2[NT],wind3[NT];
	if(!input1)
	{
		cerr<<"scenarios1.txt is not exist!"<<endl;
	}
	if(!input2)
	{
		cerr<<"scenarios2.txt is not exist!"<<endl;
	}
	if(!input3)
	{
		cerr<<"scenarios3.txt is not exist!"<<endl;
	}
	
	for(int sce = 0;sce < SCENNUM;++sce)
	{
		for(int t = 0;t < NT;++t)
		{//��ȡ��糡���ݣ��Ӳ�ͬ���ļ��ж�ȡ
			input1>>wind1[t];
			input2>>wind2[t];
			input3>>wind3[t];
		}
		
		for(int t = 0;t < NT;++t)
		{//�����糡1����
			output<<wind1[t]<<"  ";
		}
		output<<endl;
		
		for(int t = 0;t < NT;++t)
		{//�����糡2����
			output<<wind2[t]<<"  ";
		}
		output<<endl;
		
		for(int t = 0;t < NT;++t)
		{//�����糡3����
			output<<wind3[t]<<"  ";
		}
		output<<endl<<endl;
	}
	return 0;
}

