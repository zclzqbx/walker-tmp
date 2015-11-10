/*              ��������ʷֲ����������
 *1���Է������������Ϊ����������ʱ�����Ƿ��Ĳ�ȷ����
 *
 */
#include<fstream>
#include<iostream>
#include<time.h>
#include"function.h"

using namespace std;
ILOSTLBEGIN

const IloInt NG=54;//�������̨��
const IloInt NT=24;//��������ʱ����
const IloInt NL=4;//Ŀ�꺯�����Ի��ֶ���,�����Ϸֶ�Խ��Խ��ȷ���������Ӽ�����
const IloInt NW=1;//��糡����
const IloInt Node=118;//���
const IloInt Branch=186;//֧·

IloEnv env;//���л��������н���֮��Ҫ�ر�
IloModel Master_Model(env,"Master_Model");//����ģ��
IloCplex Master_Cplex(Master_Model);//����CPLEX����

ofstream output("output_fixedwind.txt",ios::ate);

IloArray<IloNumVarArray>  P(env,NG),PS(env,NG);//���鼰��Ӧ�����µĳ���
IloArray<IloBoolVarArray>  I(env,NG);//������ͣ״̬,��ά���߱���
IloArray<IloNumVarArray> PL(env,NT);//NT��ʱ�θ���֧·�ĳ���
IloArray<IloNumVarArray> P_dual(env,NT);
IloNumArray2 P1(env,NG),u(env,NG);//���ڴ洢P��I��ֵ

IloArray< IloArray<IloNumVarArray> > Deta(env,NG);//�ֶ����Ի��������
IloNumArray2 Z(env,NG);
IloNumArray2 Betaa(env,NG);

IloNumArray2 Unit(env,NG),Info_Branch(env,Branch);//���鼰֧·��Ϣ
IloNumArray Pload(env,NT);//���ɼ�������
IloNumArray R(env,NT);//����
IloNumArray2 Pwind(env,NT);//���Ƕ����糡

///////////////////////////////�밲ȫԼ���йص���///////////////////////////////////////
IloNumArray2 B0(env,Node-1),B0l(env,Node-1);//���ɾ���

IloNumArray Sl(env,Node);//ÿ���ڵ�ĸ���״̬,ÿ���ڵ���ռ���ɰٷֱ�
IloIntArray Sw(env,Node);//�����糡�ֲ�

/***********************************************���ݳ�ʼ��*****************************************************/
void define_data(IloEnv env)//���ݳ�ʼ��,��ȫ�ֱ������и�ֵ
{
	Master_Cplex.setParam(IloCplex::EpGap,0.01);//��Сһ��ͻ�out of memory
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>������ʼ��<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	for(IloInt i=0; i<NG; ++i)
	{	
		P[i]=IloNumVarArray(env,NT,0.0,IloInfinity,ILOFLOAT);		
		PS[i]=IloNumVarArray(env,NT,0.0,IloInfinity,ILOFLOAT);
		I[i]=IloBoolVarArray(env,NT);
		P1[i]=IloNumArray(env,NT);
		u[i]=IloNumArray(env,NT);
	}
	
	for(IloInt i=0;i<Node-1;++i)
	{
		B0l[i]=IloNumArray(env,Node-1);
	}
	
	for(IloInt t=0;t<NT;++t)
		Pwind[t]=IloNumArray(env,NW);
		
	for(IloInt t=0;t<NT;++t)
	{
		PL[t]=IloNumVarArray(env,Branch,-IloInfinity,IloInfinity,ILOFLOAT);
		P_dual[t]=IloNumVarArray(env,Node-1,-IloInfinity,IloInfinity,ILOFLOAT);
	}
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>������ʼ������<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>���ݶ�ȡ<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	output<<"Unit:"<<endl;
	ifstream unit_information("unit_information.txt",ios::in);//�������ֲ�
	if(!unit_information)
	{
		output<<"no such file! unit_information.txt"<<endl;
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
	unit_information.close();

	output<<endl<<endl<<"Pload:"<<endl;
	ifstream load("load.txt",ios::in);//����
	if(!load)
	{	
		output<<"no such file! load.txt"<<endl;
	}
	for(IloInt t=0; t<NT; ++t)
	{
		load>>Pload[t];
		output<<Pload[t]<<"   ";
	}
	load.close();
	
	output<<endl<<endl<<"Reserve:"<<endl;
	ifstream reserve("reserve.txt",ios::in);//����
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
	
	output<<endl<<"Pwind:"<<endl;
	ifstream windpower("wind_power.txt",ios::in);//�����㳡��
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
		
	output<<endl<<"Info_Branch:"<<endl;
	ifstream BFile("Brach_File.txt",ios::in);//��ȡ���֧·��Ϣ:��㣬�յ㣬���裬�翹�Լ���������
	if(!BFile)
	{	
		output<<"no such file! Brach_File.txt"<<endl;
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
	BFile.close();
	
	output<<endl<<"Sl:"<<endl;
	ifstream load_locate("load_locate.txt",ios::in);//���ɷֲ�
	if(!load_locate)
	{	
		output<<"no such file! load_locate.txt"<<endl;
	}
	for(IloInt k=0; k<Node; ++k)
	{	
		load_locate>>Sl[k];
		output<<Sl[k]<<"   ";
	}
	load_locate.close();

	output<<endl<<endl<<"Sw:"<<endl;
	ifstream wind_locate("wind_locate.txt",ios::in);//��糡�ֲ�
	if(!wind_locate)
	{
		output<<"no such file! wind_locate.txt"<<endl;
	}
	for(IloInt k=0; k<Node; ++k)
	{
		wind_locate>>Sw[k];
		output<<Sw[k]<<"   ";
	}
	output<<endl<<endl;
	wind_locate.close();
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>���ݶ�ȡ����<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Ŀ�꺯�����Ի�����<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    for(IloInt i=0;i<NG;++i)//�ֶ����Ի�,���Կ���ʹ���������Ի�����������̫��ʱ
	{
		Z[i]=IloNumArray(env,NL+1);
		for(IloInt l=0;l<NL+1;++l)
		{
			Z[i][l]=Unit[i][1]+(l*(Unit[i][2]-Unit[i][1]))/NL;
		}
		
		Betaa[i]=IloNumArray(env,NL);
		for(IloInt l=0;l<NL;++l)
		{
			Betaa[i][l]=( f(env,Z[i][l+1],Unit[i])-f(env,Z[i][l],Unit[i]) )/( Z[i][l+1]-Z[i][l] );
		}

		Deta[i]=IloArray<IloNumVarArray>(env,NL);
		for(IloInt l=0;l<NL;++l)
		{
			Deta[i][l]=IloNumVarArray(env,NT,0.0,IloInfinity,ILOFLOAT);
		}
	}
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Ŀ�꺯�����Ի��������<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>������B0��ֵ<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	for(IloInt k=0; k<Node-1; ++k)//�ǶԽ���,��Ȼ�����һ�������Ϊƽ����
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
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>������B0��ֵ����<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
}
/**************************************************�������****************************************************/
int main()
{
	clock_t start,finish;
	double totaltime;
	start=clock();
	
	try
	{
		define_data(env);//���ȳ�ʼ��ȫ�ֱ���
/*************************************************************������Ŀ�꺯��*******************************************************/
		IloNumExpr Cost(env);
		
		for(IloInt i=0;i<NG;++i)//���гɱ��������Ի�
			for(IloInt l=0;l<NL;++l)
				for(IloInt t=0;t<NT;++t)
				{
					Master_Model.add(Deta[i][l][t]<=Z[i][l+1]-Z[i][l]);
				}
		for(IloInt i=0;i<NG;++i)
			for(IloInt t=0;t<NT;++t)
			{
				IloNumExpr Deta_Sum(env);
				for(IloInt l=0;l<NL;++l)
				{
					Deta_Sum+=Deta[i][l][t];
				}
				Master_Model.add(P[i][t]==Deta_Sum+Unit[i][1]*I[i][t]);
				Deta_Sum.end();
			}
		for(IloInt i=0;i<NG;++i)
		{
			for(IloInt t=0;t<NT;++t)
			{
				IloNumExpr f1(env);
				f1+=(Unit[i][10]+Unit[i][9]*Unit[i][1]+Unit[i][8]*Unit[i][1]*Unit[i][1])*I[i][t];
				for(IloInt l=0;l<NL;++l)
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
		Cost.end();
/********************************************************�������������Լ��**************************************************/
		for(IloInt i=0; i<NG; ++i)
			for(IloInt t=0; t<NT; ++t)
			{
				Master_Model.add(P[i][t]<=Unit[i][2]*I[i][t]);
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
			for(IloInt i=0;i<NG;++i)
			{
				expr+=Unit[i][2]*I[i][t];
			}
			Master_Model.add(expr>=Pload[t]+R[t]);
			expr.end();
		}
/***********************************************************��С����ʱ��Լ��*************************************************/
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
/*******************************************************��Сͣ��ʱ��Լ��**************************************************/
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
/*************************************************��ȫԼ��*********************************************/
		IloInvert(env,B0,B0l,Node-1);//�������
		for(int t=0;t<NT;++t)//У������ʱ�ε�����֧·�ĳ���
		{
			IloNumExprArray Psp(env,Node-1);//ע�빦��
			IloExprArray Theta(env,Node);//���				
			
			Theta[Node-1]=IloExpr(env);
				
			for(IloInt k=0; k<Node-1; ++k)//У�������н��ĳ���
			{
				Psp[k]=IloNumExpr(env);
				
				IloInt i=0;
				for(;i<NG;++i)
				{
					if(Unit[i][0]-1==k)break;
				}
				if(i<NG)
					Master_Model.add(P_dual[t][k]==P[i][t]);	
				else
					Master_Model.add(P_dual[t][k]==0);

				Psp[k] += P_dual[t][k];
				if(Sw[k]>=0)
					Psp[k] += Pwind[t][ Sw[k] ];

				Psp[k] -= Sl[k]*Pload[t];
			}
			
			IloMutiply(env,B0l,Psp,Theta,Node-1);
			//���㳱��
			for(IloInt h=0; h<Branch; ++h)
			{
				Master_Model.add(PL[t][h]==(Theta[(IloInt)Info_Branch[h][0]-1]-Theta[(IloInt)Info_Branch[h][1]-1])/Info_Branch[h][3]);
				Master_Model.add(PL[t][h]<=Info_Branch[h][4]);
				Master_Model.add(-PL[t][h]<=Info_Branch[h][4]);
			}
			Psp.end();
			Theta.end();
		}
		
		Master_Cplex.extract(Master_Model);
		Master_Cplex.solve();
		if (Master_Cplex.getStatus() == IloAlgorithm::Infeasible)//������
		{
			output<<"Master Problem Have No Solution"<<endl;
			goto lable2;
		}
/************************************************************�����ʾ����**************************************************/
		output<<"Cost:"<<Master_Cplex.getObjValue()<<endl;
		for(IloInt t=0;t<NT;++t)
		{
			output<<"ʱ�� "<<t+1<<" �������ͣ(I)��"; 
			for(IloInt i=0;i<NG;++i)
			{
				output<<Master_Cplex.getValue(I[i][t])<<"   ";
			}
			output<<endl;
			output<<"ʱ�� "<<t+1<<" ����ĳ���(P)��";			
			for(IloInt i=0;i<NG;++i)
			{
				output<<Master_Cplex.getValue(P[i][t])<<"   ";
			}
			output<<endl;
			output<<"ʱ�� "<<t+1<<" ��·����(PL)��";			
			for(IloInt b=0;b<Branch;++b)
			{
				output<<Master_Cplex.getValue(PL[t][b])<<"   ";
			}
			output<<endl<<endl;
		}
		Master_Model.end();
		Master_Cplex.end();
		env.end();
	}
	catch(IloException& ex)//�쳣����
	{
		output<<"Error: "<<ex<<endl;
	}
	catch(...)
	{
		output<<"Error: Unknown exception caught!" << endl;
	}
	
	lable2:
	finish=clock();
	totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	output<<"totaltime: "<<totaltime<<"s"<<endl<<endl;

	output.close();	
	return 0;
}
