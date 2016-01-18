/*                  118���糡�������
 *author:�촫��
 *date:2014.10.9
 *function:1��ֱ�ӽ�ģ���������
 *2����Ӧ���������ͷ�ļ��У�
 */
#include<ilcplex/ilocplex.h>
#include<fstream>
#include<time.h>
#include<iostream>
#include"function.h"

using namespace std;
ILOSTLBEGIN

const IloInt NG     = 54;//�������̨��
const IloInt NT     = 24;//��������ʱ����
const IloInt NW     = 3;//���������
const IloInt Node   = 118;
const IloInt Branch = 186;
const IloInt NL     = 4;//Ŀ�꺯�����Ի��ֶ���
const IloInt Set    = 20;//��糡����

ofstream output("output_data/output_scuc.txt",ios::ate);

IloEnv env;//���л��������н���֮��Ҫ�ر�
IloModel Master_Model(env,"Master_Model");//����ģ��
IloCplex Master_Cplex(Master_Model);//����CPLEX����

IloArray<IloNumVarArray>  P(env,NG);//�����������ά���߱���
IloArray<IloBoolVarArray>  I(env,NG);//������ͣ״̬,��ά���߱���
IloArray< IloArray<IloNumVarArray> > Ps(env,Set);//��Ӧ�����µĻ������

IloArray< IloArray<IloNumVarArray> > Deta(env,NG);//��ά����
IloNumArray2 Z(env,NG),Betaa(env,NG);//�������������Ի�֮����

IloNumArray2 Unit(env,NG),Info_Branch(env,Branch),Pwind(env,NT);
IloArray<IloNumArray2> Pswind(env,Set);
IloNumArray ddtt(env,NG);//deta
IloNumArray Pload(env,NT);//Ԥ�⸺�ɼ������� 
IloNumArray R(env,NT);
///////////////////////////////�밲ȫԼ���йص���///////////////////////////////////////
IloNumArray Sl(env,Node);//ÿ���ڵ�ĸ���״̬,ÿ���ڵ���ռ���ɰٷֱ�
IloIntArray Sw(env,Node);//�����糡�ֲ�

//ͳһ�Խ��Node-1��Ϊƽ��ڵ㣬�ýڵ����Ϊ�㣬��B0Ϊ��Node-1��X��Node-1���ľ���
IloNumArray2 B0(env,Node-1),B0l(env,Node-1);//���ɾ����������
/*******************************************************���ݳ�ʼ��************************************************************/
void define_data(IloEnv env)
{
	Master_Cplex.setParam(IloCplex::EpGap,0.01);
	for(IloInt i=0; i<NG; ++i) //���߱���
	{
		P[i]=IloNumVarArray(env,NT,0.0,IloInfinity,ILOFLOAT);
		I[i]=IloBoolVarArray(env,NT);
	}
	for(IloInt s=0;s<Set;++s)
	{
		Ps[s]=IloArray<IloNumVarArray>(env,NG);
		for(IloInt i=0;i<NG;++i)
			Ps[s][i]=IloNumVarArray(env,NT,0.0,IloInfinity,ILOFLOAT);
	}

	for(IloInt t=0;t<NT;++t)
		Pwind[t]=IloNumArray(env,NW);
		
	for(IloInt s=0;s<Set;++s)
	{
		Pswind[s]=IloNumArray2(env,NT);
		for(IloInt t=0;t<NT;++t)
		{
			Pswind[s][t]=IloNumArray(env,NW);
		}
	}
	
	output<<endl<<"Pwind:"<<endl;
	ifstream windpower("input_data/wind_power.txt",ios::in);//����������
	if(!windpower)
	{	
		output<<"no such file! wind_pwer.txt"<<endl;
	}
	for(IloInt w=0;w<NW;++w)
	{
		for(IloInt t=0; t<NT; ++t)
		{
			windpower>>Pwind[t][w];
			output<<Pwind[t][w]<<"   ";
		}
		output<<endl;
	}
	windpower.close();
		
	output<<endl<<endl<<"Pswind:"<<endl;
	ifstream wind_scenarios("input_data/wind_scenarios.txt",ios::in);//�����糡��
	if(!wind_scenarios)
	{	
		output<<"no such file! wind_scenarios.txt"<<endl;
	}
	for(IloInt s=0;s<Set;++s)
	{	
		for(IloInt w=0;w<NW;++w)
		{
			for(IloInt t=0;t<NT;++t)
			{		
				wind_scenarios>>Pswind[s][t][w];
				output<<Pswind[s][t][w]<<"   ";
			}
			output<<endl;
		}
			output<<endl<<endl;
	}
	wind_scenarios.close();
	
	output<<"Pload:"<<endl;
	ifstream load_input("input_data/load.txt",ios::in);//����
	if(!load_input)
	{
		cerr<<"no such file! load.txt"<<endl;
	}
	for(IloInt t=0; t<NT; ++t)
	{
		load_input>>Pload[t];
		output<<Pload[t]<<"   ";
	}
	output<<endl<<endl;
	load_input.close();
	
	output<<endl<<endl<<"Reserve:"<<endl;
	ifstream reserve("input_data/reserve.txt",ios::in);//����
	if(!reserve)
	{	
		output<<"no such file! reserve.txt"<<endl;
	}
	for(IloInt t=0; t<NT; ++t)
	{
		reserve>>R[t];
		output<<R[t]<<"   ";
	}
	reserve.close();	
	
	output<<endl<<endl<<"Unit:"<<endl;
	ifstream unit_information("input_data/unit_information.txt",ios::in);//�������ֲ�
	if(!unit_information)
	{
		cerr<<"no such file! unit_information.txt"<<endl;
	}
	for(IloInt i=0;i<NG;++i)//���������Ϣ
	{
		Unit[i]=IloNumArray(env,12);//12����Ϣ:��㡢Pmin��Pmax��MinOn��MinOff��RU��RD��Startup Cost��c��b��a��price;
		for(IloInt k=0; k<12; ++k)
		{
			unit_information>>Unit[i][k];
			output<<Unit[i][k]<<"   ";
		}
		output<<endl;
	}
	output<<endl<<endl;
	unit_information.close();
	
	output<<endl<<endl<<"ddtt:"<<endl;
	for(IloInt i=0;i<NG;++i)
	{
		ddtt[i]=Unit[i][5]/2;
		output<<ddtt[i]<<"   ";
	}
	
	output<<endl<<endl<<"Info_Branch:"<<endl;
	ifstream BFile("input_data/Brach_File.txt",ios::in);//��ȡ֧·��Ϣ
	if(!BFile)
	{
		cerr<<"no such file! Brach_File.txt"<<endl;
	}
	for(IloInt k=0; k<Branch; ++k)
	{
		Info_Branch[k]=IloNumArray(env,5);
		for(IloInt h=0; h<5; ++h)
		{
			BFile>>Info_Branch[k][h];
			output<<Info_Branch[k][h]<<"   ";
		}
		output<<endl;
	}
	output<<endl<<endl;
	BFile.close();

	output<<"Sl:"<<endl;
	ifstream load_locate("input_data/load_locate.txt",ios::in);//���ɷֲ�
	if(!load_locate)
	{
		cerr<<"no such file! load_locate.txt"<<endl;
	}
	for(IloInt k=0; k<Node; ++k)
	{
		load_locate>>Sl[k];
		output<<Sl[k]<<"   ";
	}
	output<<endl<<endl;
	load_locate.close();

	output<<"Sw:"<<endl;
	ifstream wind_locate("input_data/wind_locate.txt",ios::in);//��糡�ֲ�
	if(!wind_locate)
	{
		cerr<<"no such file! wind_locate.txt"<<endl;
	}
	for(IloInt k=0; k<Node; ++k)
	{
		wind_locate>>Sw[k];
		output<<Sw[k]<<"   ";
	}
	output<<endl<<endl;
	wind_locate.close();
	
	for(IloInt i=0; i<NG; ++i) //�ֶ����Ի�
	{
		Z[i]=IloNumArray(env,NL+1);
		for(IloInt l=0; l<NL+1; ++l)
		{
			Z[i][l]=Unit[i][1]+(l*(Unit[i][2]-Unit[i][1]))/NL;
		}
	}
	for(IloInt i=0; i<NG; ++i)
	{
		Betaa[i]=IloNumArray(env,NL);
		for(IloInt l=0; l<NL; ++l)
		{
			Betaa[i][l]=( f(env,Z[i][l+1],Unit[i])-f(env,Z[i][l],Unit[i]) )/( Z[i][l+1]-Z[i][l] );
		}
	}
	for(IloInt i=0; i<NG; ++i)
	{
		Deta[i]=IloArray<IloNumVarArray>(env,NL);
		for(IloInt l=0; l<NL; ++l)
		{
			Deta[i][l]=IloNumVarArray(env,NT,0.0,IloInfinity,ILOFLOAT);
		}
	}
	/////������B0��ֵ
	for(IloInt k=0; k<Node-1; ++k)//�ǶԽ���
	{
		B0[k]=IloNumArray(env,Node-1);
		for(IloInt h=0; h<Node-1; ++h)
		{
			if(k==h)
				continue;
			else
			{
				IloInt temp;
				IloInt b=0;
				for(; b<Branch; ++b)
				{
					if((Info_Branch[b][0]-1==k && Info_Branch[b][1]-1==h)
					        || (Info_Branch[b][0]-1==h && Info_Branch[b][1]-1==k))
					{
						temp=b;
						break;
					}
				}
				if(b==Branch)
					B0[k][h]=0;
				else
					B0[k][h]= -1.0/Info_Branch[temp][3];
			}
		}
	}
	for(IloInt k=0; k<Node-1; ++k)//�Խ���
	{
		B0[k][k]=0;
		for(IloInt h=0; h<Node-1; ++h)
		{
			if(k!=h)
			{
				B0[k][k]+=B0[k][h];
			}
		}
		IloNum last=0;
		IloInt b=0;
		IloInt temp;
		for(;b<Branch;++b)
		{
			if((Info_Branch[b][0]-1==k && Info_Branch[b][1]-1==Node-1)
					        || (Info_Branch[b][0]-1==Node-1 && Info_Branch[b][1]-1==k))
			{
				temp=b;
				break;
			}
		}
		if(b==Branch)
			last=0;
		else
			last=-1.0/Info_Branch[temp][3];
		B0[k][k]+=last;
		B0[k][k]=-1*B0[k][k];
	}
	/////��ֵ���	
	for(IloInt k=0; k<Node-1; ++k)
		B0l[k]=IloNumArray(env,Node-1);
		
	IloInvert(env,B0,B0l,Node-1);//�������
}

int main()
{
	clock_t start,finish;
	double totaltime;
	start=clock();
	time_t nowTime=time(0);
	struct tm* nowTimeStruct=localtime(&nowTime);
	output<<"ϵͳ��ǰʱ�䣺"<<1900+nowTimeStruct->tm_year<<"."<<nowTimeStruct->tm_mon<<"."<<
		nowTimeStruct->tm_mday<<"  "<<nowTimeStruct->tm_hour<<":"<<nowTimeStruct->tm_min<<":"<<nowTimeStruct->tm_sec<<endl;

	try
	{
		output<<">>>>>>>>>>>>>>������<<<<<<<<<<<<<<<"<<endl;		
		define_data(env);//���ȳ�ʼ��ȫ�ֱ���
		output<<">>>>>>>>>>>>>>����������<<<<<<<<<<<<<<<"<<endl;	
/*************************************************************Ŀ�꺯��*******************************************************/
		IloNumExpr Cost(env);

		for(IloInt i=0; i<NG; ++i) //���гɱ��������Ի�
			for(IloInt l=0; l<NL; ++l)
				for(IloInt t=0; t<NT; ++t)
				{
					Master_Model.add(Deta[i][l][t]<=Z[i][l+1]-Z[i][l]);
				}
				
		for(IloInt i=0; i<NG; ++i)
			for(IloInt t=0; t<NT; ++t)
			{
				IloNumExpr Deta_Sum(env);
				for(IloInt l=0; l<NL; ++l)
				{
					Deta_Sum+=Deta[i][l][t];
				}
				Master_Model.add(P[i][t]==Deta_Sum+Unit[i][1]*I[i][t]);
				Deta_Sum.end();
			}	
			
		for(IloInt i=0; i<NG; ++i)
		{
			for(IloInt t=0; t<NT; ++t)
			{
				IloNumExpr f1(env);
				f1+=(Unit[i][10]+Unit[i][9]*Unit[i][1]+Unit[i][8]*Unit[i][1]*Unit[i][1])*I[i][t];
				for(IloInt l=0; l<NL; ++l)
					f1+=Betaa[i][l]*Deta[i][l][t];
				Cost+=f1*Unit[i][11];
				f1.end();
			}
		}

		for(IloInt i=0; i<NG; ++i)//�����ɱ�
			for(IloInt t=1; t<NT; ++t)
			{
				Cost+=Unit[i][11]*Unit[i][7]*I[i][t]*(1-I[i][t-1]);
			}

		for(IloInt i=0; i<NG; ++i)//ͣ���ɱ�
			for(IloInt t=1; t<NT; ++t)
			{
				Cost+=Unit[i][11]*Unit[i][7]*I[i][t-1]*(1-I[i][t]);
			}
		Master_Model.add(IloMinimize(env,Cost));
/********************************************************�������������Լ��**************************************************/
		for(IloInt i=0; i<NG; ++i)
			for(IloInt t=0; t<NT; ++t)
			{
				Master_Model.add(P[i][t]<=Unit[i][2]*I[i][t]);//��÷ֿ�д
				Master_Model.add(P[i][t]>=Unit[i][1]*I[i][t]);
			}
/*********************************************************��������Լ��******************************************************/
		for(IloInt i=0; i<NG; ++i)
			for(IloInt t=1; t<NT; ++t)
			{
				Master_Model.add(P[i][t]-P[i][t-1]<=(1-I[i][t]*(1-I[i][t-1]))*Unit[i][5]+I[i][t]*(1-I[i][t-1])*Unit[i][1]);
				Master_Model.add(P[i][t-1]-P[i][t]<=(1-I[i][t-1]*(1-I[i][t]))*Unit[i][6]+I[i][t-1]*(1-I[i][t])*Unit[i][1]);
			}
/*********************************************************���鹦��ƽ��Լ��**************************************************/
		for(IloInt t=0; t<NT; ++t)
		{
			IloNumExpr fire(env);
			IloNum wind(0);
			for(IloInt i=0; i<NG; ++i)
			{
				fire+=P[i][t];
			}
			
			for(IloInt w=0;w<NW;++w)
			{
				wind+=Pwind[t][w];
			}
			
			Master_Model.add(fire==Pload[t]-wind);
			fire.end();
		}
/**********************************************************����Լ��********************************************************/
		for(IloInt t=0; t<NT; ++t) 
		{
			IloNumExpr expr(env);
			for(IloInt i=0; i<NG; ++i)
			{
				expr+=Unit[i][2]*I[i][t];
			}
			Master_Model.add(expr>=R[t]+Pload[t]);
			expr.end();
		}
/**********************************************************��ȫԼ��********************************************************/
		////��ÿ���ڵ��ע�빦��
		for(IloInt t=0; t<NT; ++t)
		{
			IloExprArray Pf(env,Branch);//����
			IloExprArray Psp(env,Node-1);//ע�빦��
			IloExprArray Theta(env,Node);//���
			
			Theta[Node-1]=IloExpr(env);
			
			for(IloInt k=0; k<Node-1; ++k)
			{
				Psp[k]=IloExpr(env);
				IloInt i=0;
				for(;i<NG;++i)
				{
					if(Unit[i][0]-1==k)break;
				}
				if(i<NG)
				{
					Psp[k] +=P[i][t];
				}
				if(Sw[k]>=0)
					Psp[k] += Pwind[t][ Sw[k] ] ;
				Psp[k] -= Sl[k]*Pload[t];
			}
			IloMutiply(env,B0l,Psp,Theta,Node-1);
			
			//���㳱��
			for(IloInt h=0; h<Branch; ++h)
			{
				Pf[h]=IloExpr(env);			
				Pf[h]+=(Theta[(IloInt)Info_Branch[h][0]-1]-Theta[(IloInt)Info_Branch[h][1]-1])/Info_Branch[h][3];
				Master_Model.add(Pf[h]<=Info_Branch[h][4]);
				Master_Model.add(Pf[h]>=-Info_Branch[h][4]);
			}
			Pf.end();
			Psp.end();
			Theta.end();
		}
/************************************************************��С����ʱ��Լ��*****************************************************/
		for(IloInt i=0; i<NG; ++i)
		{
			IloNumExpr expr1(env);
			for(IloInt t=0; t<=Unit[i][3]-1; ++t)
			{
				expr1+=I[i][t];
			}
			Master_Model.add((expr1-I[i][0]*Unit[i][3])>=0);
			expr1.end();

			for(IloInt t=1; t<=NT-Unit[i][3]; ++t)
			{
				IloNumExpr expr2(env);
				for(IloInt k=t; k<=t+Unit[i][3]-1; ++k)
				{
					expr2+=I[i][k];
				}
				Master_Model.add(expr2>=Unit[i][3]*(I[i][t]-I[i][t-1]));
				expr2.end();
			}
			
			for(IloInt t=NT-(IloInt)Unit[i][3]+1; t<=NT-1; ++t)
			{
				IloNumExpr expr2(env);
				for(IloInt h=t; h<=NT-1; ++h)
				{
					expr2+=(I[i][h]-(I[i][t]-I[i][t-1]));
				}
				Master_Model.add(expr2>=0);
				expr2.end();
			}
		}
/*******************************************************��Сͣ��ʱ��Լ��******************************************************/
		for(IloInt i=0; i<NG; ++i)
		{
			IloNumExpr expr1(env);
			for(IloInt t=0; t<=Unit[i][4]-1; ++t)
			{
				expr1+=(1-I[i][t]);
			}
			Master_Model.add((expr1-(1-I[i][0])*Unit[i][4])>=0);
			expr1.end();

			for(IloInt t=1; t<=NT-Unit[i][4]; ++t)
			{
				IloNumExpr expr1(env);
				for(IloInt k=t; k<=t+Unit[i][4]-1; ++k)
				{
					expr1+=(1-I[i][k]);
				}
				Master_Model.add(expr1>=Unit[i][4]*(I[i][t-1]-I[i][t]));
				expr1.end();
			}
			
			for(IloInt t=NT-(IloInt)Unit[i][4]+1; t<=NT-1; ++t)
			{
				IloNumExpr expr2(env);
				for(IloInt h=t; h<=NT-1; ++h)
				{
					expr2+=(1-I[i][h]-(I[i][t-1]-I[i][t]));
				}
				Master_Model.add(expr2>=0);
				expr2.end();
			}
		}

/*******************************************************�����ǳ������Լ��******************************************************/		
		for(IloInt s=0;s<Set;++s)
		{
			for(IloInt t=0;t<NT;++t)//����ƽ��Լ��
			{
				IloNumExpr fire(env);
				for(IloInt i=0; i<NG; ++i)
				{
					fire+=Ps[s][i][t];
				}
				
				IloNum wind(0);
				for(IloInt w=0;w<NW;++w)
					wind+=Pswind[s][t][w];
					
				Master_Model.add(fire==Pload[t]-wind);
				fire.end();
			}
			
			for(IloInt i=0; i<NG; ++i)//�������������Լ��
				for(IloInt t=0; t<NT; ++t)
				{
					Master_Model.add(Ps[s][i][t]<=Unit[i][2]*I[i][t]);//��÷ֿ�д
					Master_Model.add(Ps[s][i][t]>=Unit[i][1]*I[i][t]);
				}
			
			for(IloInt i=0; i<NG; ++i)//������ΪԼ��
				for(IloInt t=0; t<NT; ++t)
				{
					Master_Model.add(Ps[s][i][t]-P[i][t]<=ddtt[i]);
					Master_Model.add(-Ps[s][i][t]+P[i][t]<=ddtt[i]);
				}
			
			for(IloInt t=0; t<NT; ++t)//��ȫԼ��
			{
				IloExprArray Pf(env,Branch);//����
				IloExprArray Psp(env,Node-1);//ע�빦��
				IloExprArray Theta(env,Node);//���
				
				Theta[Node-1]=IloExpr(env);
				for(IloInt k=0; k<Node-1; ++k)
				{
					Psp[k]=IloExpr(env);
					
					IloInt i=0;
					for(;i<NG;++i)
					{
						if(Unit[i][0]-1==k)break;
					}
					if(i<NG)
					{
						Psp[k] +=Ps[s][i][t];
					}
					if(Sw[k]>=0)
						Psp[k] += Pswind[s][t][ Sw[k] ];
					Psp[k] -= Sl[k]*Pload[t];
				}
				IloMutiply(env,B0l,Psp,Theta,Node-1);
				
				//���㳱��
				for(IloInt h=0; h<Branch; ++h)
				{
					Pf[h]=IloExpr(env);			
					Pf[h]+=(Theta[(IloInt)Info_Branch[h][0]-1]-Theta[(IloInt)Info_Branch[h][1]-1])/Info_Branch[h][3];
					Master_Model.add(Pf[h]<=Info_Branch[h][4]);
					Master_Model.add(Pf[h]>=-Info_Branch[h][4]);
				}
				Pf.end();
				Psp.end();
				Theta.end();
			}			
		}		
		Master_Cplex.solve();                        
/************************************************************�����ʾ����**************************************************/
		if (Master_Cplex.getStatus() == IloAlgorithm::Infeasible)//������
			output << "No Solution" << endl;

		output<<"Cost:"<<Master_Cplex.getObjValue()<<endl;
		for(IloInt i=0; i<NG; ++i)
		{
			for(IloInt t=0; t<NT; ++t)
			{
				output<<"���� "<<i+1<<" �� "<<t+1<<" ʱ��״̬:"<<Master_Cplex.getValue(I[i][t])<<"   ����:";
				output<<Master_Cplex.getValue(P[i][t])<<endl;
			}
			output<<endl;
		}
		
		Master_Model.end();
		Master_Cplex.end();
		env.end();
	}
	catch(IloException& ex)//�쳣����
	{
		cerr<<"Error: "<<ex<<endl;
	}
	catch(...)
	{
		cerr << "Error: Unknown exception caught!" << endl;
	}

	finish=clock();
	totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	output<<"totaltime: "<<totaltime<<"s"<<endl<<endl;
	
	output.close();
	return 0;
}