#include<iostream>
#include<fstream>
#include<math.h>

using namespace std;
/***********************************************ȫ�ֲ�������ʼ��*****************************************************/
const int T=24;//ʱ����
const int N=200;//ԭ��������
int J[T]={0,1,0,0,1,0,2,0,0,1,0,0,2,0,2,0,0,3,0,3,0,4,0,180};//����ʱ��ɾ���Ľڵ���

double original_scenarios[N][T];//��Ҫ��ȡԭ��������
double rate[200][24];//���������ĸ���

bool flag[200];//�����Ƿ��ѱ�ɾ���ı�־������������
bool flag1[20][24];//��־λ���ٴ�������
/*******************************************�󳡾�֮��ľ���**********************************************************/
double fdistance(const double scenario1,const double scenario2)
{
	double distance;
	distance=sqrt((scenario1-scenario2)*(scenario1-scenario2));
	return distance;
}
/******************************************��������������200������20����************************************************/
void delete_scenario(const double original_scenarios[N][T],double reserve_first[20][24],const int num_delete)//��һ��ɾ������
{
	for(int k=0;k<num_delete;k++)
	{
		int l_delete=0;
		int final_delete=0;
		double final_delete_value=100000000;
		int l1=0;
		int first_cycle=0;
		int first_of_final_delete=0;
		while(1)//ֱ�������п���ɾ���ĳ����������
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
					flag[l1]=1;//��forѭ�������ֽ���������һ��l1==200;һ��break;
					break;
				}
			}

			first_cycle=1;
			if(l1>=200)break;//����ѭ��

			double z_delete[N];//z[i]
			for(int i1=0;i1<N;++i1)
			{
				z_delete[i1]=0;
			}
			int z_num=0;

			for(int l2=0;l2<N;l2++)//������J���ϵĳ������J���ϳ����ľ���
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

			if(first_of_final_delete==0)//�ж��Ƿ��ǵ�һ��ѭ��
			{
				final_delete_value=sum_z_delete;
				final_delete=l_delete;
				first_of_final_delete=1;
			}
			else
			{
				if(sum_z_delete<final_delete_value)
				{
					final_delete_value=sum_z_delete;//�ж��Ƿ�Ҫ�ı���Сֵ�����ı�Ҫɾ���ĳ���
					final_delete=l_delete;
				}
			}
		}
		flag[final_delete]=1;
	}//��forѭ������ʱ������Ҫ��ɾ���ĳ��������־λ���ѱ�Ϊ1��
	
	int reserve_num=0;
	for(int i=0;i<N;++i)//����Ҫ�������������ĸ���
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
	for(int i=0;i<200;++i)//�����
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

}//��һ�������������г�����������Ŀ�곡��

/********************************************�ٴ�������20������20��������****************************************************/
void delete_scenario(double reserve_first[20][24],const int num_delete,const int time)//ɾ������
{
	bool flag_new[20];
	for(int i=0;i<20;++i)//��ʼ������
		flag_new[i]=0;
	for(int i=0;i<20;++i)
	{
		if(flag1[i][time]==0)
			flag_new[i]=1;//ȷʵû�б�ɾ��,�൱���ж��Ƿ�����һʱ��ʣ�ೡ��
	}

    for(int k=0;k<num_delete;++k)//��ʼɾ����ֱ���ﵽ��ʱ����Ҫɾ���ĳ�����
	{
		int l3=0;
		int l_delete=0;
		int first_of_final_delete=0;
        double final_delete_value=0;
		int final_delete=0;

		int first_cycle=0;
		while(1)//ÿѭ��һ��ɾ��һ������
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

			double z_delete[20];//whileÿѭ��һ�Σ�z_delete�������¸�һ��ֵ
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

			if(first_of_final_delete==0)//�ж��Ƿ��ǵ�һ��ѭ��
			{
				final_delete_value=sum_z_delete;
				final_delete=l_delete;
				first_of_final_delete=1;
			}
			else
			{
				if(sum_z_delete<final_delete_value)
				{
					final_delete_value=sum_z_delete;//�ж��Ƿ�Ҫ�ı���Сֵ�����ı�Ҫɾ���ĳ���
					final_delete=l_delete;
				}
			}
		}

		for(int pre_scenario=0;pre_scenario<=time;++pre_scenario)//�Ѵ˳�����ǰ�����г�����־λ��1
		{

		    flag1[final_delete][pre_scenario]=1;
		}
	}

	for(int i=0;i<20;++i)//�޸ĳ����������
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
/****************************************����������������*******************************************************/

int main()
{
	for(int i=0;i<200;++i)
	{
		rate[i][23]=1/200;//��ʼ�������
		flag[i]=0;//ԭ����δ��ɾ��
	}
	
	for(int i=0;i<20;++i)//ȫ����־λ�ĳ�ʼ��
		for(int j=0;j<24;++j)
			flag1[i][j]=0;

	double reserve_first[20][24];

	ifstream in("24ʱ�η�糡��_200.txt",ios::in);//��ȡ��糡������
	for(int i=0;i<N;++i)
		for(int j=0;j<T;++j)
			in>>original_scenarios[i][j];

    delete_scenario(original_scenarios,reserve_first,J[23]);//��������

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
		    delete_scenario(reserve_first,J[t],t);//�ٴ������γɳ�����
		}
	}
		
	ofstream out("200to20���������.txt",ios::ate);
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
