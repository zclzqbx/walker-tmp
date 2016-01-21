/*author:朱传林
 *date:2015.12.3
 *function:根据1000个场景产生一个随机场景（平均值）
 *
 */
 
 #include <iostream>
 #include <fstream>
 #include <memory.h>
  
 using namespace std;
 
 const int INPUTSCENUM = 1000;
 const int          NT = 24;
 
 ifstream input("1000scenarios.txt",ios::in);
 ofstream output("random_scenario.txt",ios::out);
 
 int main()
 {
	double resultSce[NT] = {0.0};
	double tmpSce[NT]    = {0.0};
	
	memset(resultSce,0,NT*sizeof(double));
	for(int t = 0;t < NT;++t)
	{
		cout<< resultSce[t]<<" ";
	}
	
	for(int i = 0;i < INPUTSCENUM;++i)
	{
		for(int t = 0;t < NT;++t)
		{
			input>>tmpSce[t];
			resultSce[t] += tmpSce[t];
		}		
	}
	
	for(int t = 0;t < NT;++t)
	{
		resultSce[t] /= 1000;
		output<<resultSce[t]<<"  ";
	}
	output<<endl;
	
	return 0;
 }
 