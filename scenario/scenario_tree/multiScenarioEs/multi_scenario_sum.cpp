/*
 *author:�촫��
 *date:2015.12.3
 *function���������糡�ĳ�����ϳ�һ������(���)
 *
 */
#include <iostream>
#include <fstream>

using namespace std;
 
const int WINDNUM = 3;
const int NT      = 24;
const int SCENNUM = 8;

//֮���Էֿ�����ǿ��ǵ�ϵͳ�ڴ������
ifstream input("scenarios.txt",ios::in);//��糡1

ofstream output("merge_scenarios.txt",ios::ate);//�ۺϳ���������

int main()
{
	double wind1[NT],wind2[NT],wind3[NT];
	if(!input)
	{
		cerr<<"scenarios.txt is not exist!"<<endl;
	}

	
	for(int sce = 0;sce < SCENNUM;++sce)
	{
		for(int t = 0;t < NT;++t)
		{
			input>>wind1[t];
		}
		
		for(int t = 0;t < NT;++t)
		{
			input>>wind2[t];
		}
		
		for(int t = 0;t < NT;++t)
		{
			input>>wind3[t];
		}
		
		for(int t = 0;t < NT;++t)
		{//�����糡1����
			output<<wind1[t]+wind2[t]+wind3[t]<<"  ";
		}
		output<<endl;
	}
	return 0;
}

