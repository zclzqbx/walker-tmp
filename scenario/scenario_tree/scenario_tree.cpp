#include<iostream>
#include<fstream>
#include<math.h>

using namespace std;
/***********************************************全局参数及初始化*****************************************************/
const int T=24;//时段数
const int N=200;//原场景总数
int J[T]={0,1,0,0,1,0,2,0,0,1,0,0,2,0,2,0,0,3,0,3,0,4,0,180};//各个时段删除的节点数

double original_scenarios[N][T];//需要读取原场景数据
double rate[200][24];//各个场景的概率

bool flag[200];//场景是否已被删除的标志（初次缩减）
bool flag1[20][24];//标志位（再次缩减）
/*******************************************求场景之间的距离**********************************************************/
double fdistance(const double scenario1,const double scenario2)
{
	double distance;
	distance=sqrt((scenario1-scenario2)*(scenario1-scenario2));
	return distance;
}
/******************************************初步场景缩减，200场景至20场景************************************************/
void delete_scenario(const double original_scenarios[N][T],double reserve_first[20][24],const int num_delete)//第一步删除方法
{
	for(int k=0;k<num_delete;k++)
	{
		int l_delete=0;
		int final_delete=0;
		double final_delete_value=100000000;
		int l1=0;
		int first_cycle=0;
		int first_of_final_delete=0;
		while(1)//直到把所有可能删除的场景试验结束
		{
			if(first_cycle=!0)
			{
			    flag[l_delete]=0;
				l1++;
			}
			for(;l1<200;++l1)
			{
				if(flag[l1]==0)
				{
					l_delete=l1;
					flag[l1]=1;//此for循环有两种结束条件，一是l1==200;一是break;
					break;
				}
			}

			first_cycle=1;
			if(l1>=200)break;//结束循环

			double z_delete[N];//z[i]
			for(int i1=0;i1<N;++i1)
			{
				z_delete[i1]=0;
			}
			int z_num=0;

			for(int l2=0;l2<N;l2++)//求所有J集合的场景与非J集合场景的距离
			{
				if(flag[l2]==1)
				{
					double distance_smallest=100000000;

					for(int j=0;j<N;++j)
					{
						if(flag[j]==0 && fdistance(original_scenarios[l2][23],original_scenarios[j][23])<distance_smallest)
							distance_smallest=fdistance(original_scenarios[l2][23],original_scenarios[j][23]);
					}

					z_delete[z_num++]=rate[l2][23]*distance_smallest;
				}
			}
			double sum_z_delete=0;
			for(int i=0;i<z_num;++i)
			{
				sum_z_delete+=z_delete[i];
			}

			if(first_of_final_delete==0)//判断是否是第一次循环
			{
				final_delete_value=sum_z_delete;
				final_delete=l_delete;
				first_of_final_delete=1;
			}
			else
			{
				if(sum_z_delete<final_delete_value)
				{
					final_delete_value=sum_z_delete;//判断是否要改变最小值，并改变要删除的场景
					final_delete=l_delete;
				}
			}
		}
		flag[final_delete]=1;
	}//此for循环结束时，所有要被删除的场景，其标志位都已变为1；
	
	int reserve_num=0;
	for(int i=0;i<N;++i)//还需要求保留下来场景的概率
	{
		if(flag[i]==0)
		{
			for(int j=0;j<T;++j)
			{
				reserve_first[reserve_num][j]=original_scenarios[i][j];
			}
			reserve_num++;
		}
		
	}
	for(int i=0;i<200;++i)//求概率
	{
		double rate_distance=10000000;
		int rate_final;
		if(flag[i]==1)
		{
			for(int j=0;j<20;++j)
			{
				if(fdistance(original_scenarios[i][23],reserve_first[j][23])<rate_distance)
				{
					rate_distance=fdistance(original_scenarios[i][23],reserve_first[j][23]);
					rate_final=j;
				}
			}
			rate[rate_final][23]+=rate[i][23];
		}
	}

}//第一次缩减，把所有场景缩减至了目标场景

/********************************************再次缩减（20场景至20场景树）****************************************************/
void delete_scenario(double reserve_first[20][24],const int num_delete,const int time)//删除方法
{
	bool flag_new[20];
	for(int i=0;i<20;++i)//初始化工作
		flag_new[i]=0;
	for(int i=0;i<20;++i)
	{
		if(flag1[i][time]==0)
			flag_new[i]=1;//确实没有被删除,相当于判断是否是这一时段剩余场景
	}

    for(int k=0;k<num_delete;++k)//开始删除，直到达到此时段需要删除的场景数
	{
		int l3=0;
		int l_delete=0;
		int first_of_final_delete=0;
        double final_delete_value=0;
		int final_delete=0;

		int first_cycle=0;
		while(1)//每循环一次删除一个场景
		{
			if(first_cycle=!0)
			{
			    flag1[l_delete][time]=0;
				l3++;
			}
			for(;l3<20;++l3)
			{
				if(flag1[l3][time]==0)
				{
					flag1[l3][time]=1;
					l_delete=l3;
					break;
				}
			}
			first_cycle=1;
			if(l3>=20)break;

			double z_delete[20];//while每循环一次，z_delete都会重新赋一次值
		    for(int i=0;i<20;++i)
			     z_delete[i]=0;
			int num_z_delete=0;

			for(int i=0;i<20;++i)
			{
				if(flag1[i][time]==1 && flag_new[i]==1)
				{
					double smallest_distance=10000000;
					for(int j=0;j<20;++j)
					{
						if(flag1[j][time]==0 && fdistance(reserve_first[i][time],reserve_first[j][time])<smallest_distance)
							smallest_distance=fdistance(reserve_first[i][time],reserve_first[j][time]);
					}
					z_delete[num_z_delete++]=rate[i][time]*smallest_distance;
				}
			}
			double sum_z_delete=0;
			for(int i=0;i<num_z_delete;++i)
			{
				sum_z_delete+=z_delete[i];
			}

			if(first_of_final_delete==0)//判断是否是第一次循环
			{
				final_delete_value=sum_z_delete;
				final_delete=l_delete;
				first_of_final_delete=1;
			}
			else
			{
				if(sum_z_delete<final_delete_value)
				{
					final_delete_value=sum_z_delete;//判断是否要改变最小值，并改变要删除的场景
					final_delete=l_delete;
				}
			}
		}

		for(int pre_scenario=0;pre_scenario<=time;++pre_scenario)//把此场景及前面所有场景标志位置1
		{

		    flag1[final_delete][pre_scenario]=1;
		}
	}

	for(int i=0;i<20;++i)//修改场景及其概率
	{
		double rate_distance=1000000;
		int combine_final;
		if(flag1[i][time]==1 && flag_new[i]==1)
		{
			for(int j=0;j<20;++j)
			{
				if(flag1[j][time]==0 && fdistance(reserve_first[i][time],reserve_first[i][time])<rate_distance)
				{
					rate_distance=fdistance(reserve_first[i][time],reserve_first[i][time]);
					combine_final=j;
				}
			}
			for(int p=0;p<time;++p)
			{
				reserve_first[i][p]=reserve_first[combine_final][p];
			}
			rate[combine_final][time]+=rate[i][time];
		}
	}
}
/****************************************以下是主函数部分*******************************************************/

int main()
{
	for(int i=0;i<200;++i)
	{
		rate[i][23]=1/200;//初始概率相等
		flag[i]=0;//原矩阵都未被删除
	}
	
	for(int i=0;i<20;++i)//全部标志位的初始化
		for(int j=0;j<24;++j)
			flag1[i][j]=0;

	double reserve_first[20][24];

	ifstream in("24时段风电场景_200.txt",ios::in);//读取风电场景数据
	for(int i=0;i<N;++i)
		for(int j=0;j<T;++j)
			in>>original_scenarios[i][j];

    delete_scenario(original_scenarios,reserve_first,J[23]);//初次缩减

	for(int i=0;i<20;i++)
	{
		if(flag[i]==0)
		{
			for(int j=22;j>=0;--j)
			{
				rate[i][j]=rate[i][23];
			}
		}
	}

	for(int t=T-2;t>=0;t--)
	{
		if(J[t]!=0)
		{
		    delete_scenario(reserve_first,J[t],t);//再次缩减形成场景树
		}
	}
		
	ofstream out("200to20场景树结果.txt",ios::ate);
	for(int i=0;i<20;++i)
	{
		for(int j=0;j<T;++j)
			out<<reserve_first[i][j]<<"  ";
		out<<endl;
	}

	in.close ();
	out.close();
	return 0;
}
