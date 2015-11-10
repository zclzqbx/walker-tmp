/*              ��������ʷֲ����������
 *1��
 *
 *
 */
#include<ilcplex/ilocplex.h>
#include<fstream>
#include<iostream>
#include<time.h>

using namespace std;
ILOSTLBEGIN

const IloInt NG=54;//�������̨��
const IloInt NT=24;//��������ʱ����
const IloInt NL=4;//Ŀ�꺯�����Ի��ֶ���,�����Ϸֶ�Խ��Խ��ȷ���������Ӽ�����
const IloInt NW=1;//��糡����
const IloInt set=10;//��������
const IloInt Node=118;//���
const IloInt Branch=186;//֧·
const IloNum limit=0.1;//���ظֵ
const IloNum Cost_Modify=1;

IloEnv env;//���л��������н���֮��Ҫ�ر�
IloModel Master_Model(env,"Master_Model");//����ģ��
IloCplex Master_Cplex(Master_Model);//����CPLEX����

ofstream output("output_benders.txt",ios::ate);

IloArray<IloNumVarArray>  P(env,NG),PS(env,NG);//���鼰��Ӧ�����µĳ���
IloArray<IloBoolVarArray>  I(env,NG);//������ͣ״̬,��ά���߱���
IloNumArray2 P1(env,NG),u(env,NG);//���ڴ洢P��I��ֵ

IloArray< IloArray<IloNumVarArray> > Deta(env,NG);//�ֶ����Ի��������
IloNumArray2 Z(env,NG);
IloNumArray2 Betaa(env,NG);

IloNumArray2 Unit(env,NG),Info_Branch(env,Branch);//���鼰֧·��Ϣ
IloNumArray /*Pwind(env,NT),*/Pload(env,NT);//���ɼ�������
IloNumArray2 Pwind(env,NT);
IloArray<IloNumArray2> Pswind(env,set);//��糡��
IloNum detaa(10);//������ΪԼ��,�˱������ɶȱȽϴ�

///////////////////////////////�밲ȫԼ���йص���///////////////////////////////////////
IloNumArray2 B0(env,Node-1),B0l(env,Node-1);//���ɾ���

IloNumArray Sl(env,Node);//ÿ���ڵ�ĸ���״̬,ÿ���ڵ���ռ���ɰٷֱ�
IloIntArray Sw(env,Node);//�����糡�ֲ�

IloNum f(IloNum &Z,IloNumArray &unit)
{
	return unit[10]*Cost_Modify+unit[9]*Z*Cost_Modify+unit[8]*Z*Z*Cost_Modify;
}
/***************************************************��������****************************************************/
void IloInvert(const IloNumArray2 &B0,IloNumArray2 &A,IloInt N)//�������
{
	IloNumArray2 C(env,N);//��չ����

	for(IloInt i=0; i<N; ++i) //���쵥λ��
		for(IloInt j=0; j<N; ++j)
		{
			if(i==j)
				A[i][j]=1;
			else
				A[i][j]=0;
		}

	for(IloInt i=0; i<N; ++i) //���������
	{
		C[i]=IloNumArray(env,2*N); 
		for(IloInt j=0; j<2*N; ++j)
		{
			if(j<N)
				C[i][j]=B0[i][j];
			else
				C[i][j]=A[i][j-N];
		}
	}
		
	for(IloInt i=0; i<N; ++i)
	{
		if(C[i][i]==0)
		{
			IloInt j=i+1;
			for(; j<N; ++j) //�ӱ��е���һ�п�ʼ������Ϊ�����
			{
				if(C[j][i]!=0)//�ҵ�����н���
				{
					IloNum temp2;
					for(IloInt k=i; k<2*N; ++k)
					{
						temp2=C[i][k];
						C[i][k]=C[j][k];
						C[j][k]=temp2;
					}
					break;
				}				
			}
		}
		
		IloNum temp3=C[i][i];//�����׷���ֵ
		
		for(IloInt k=i; k<2*N; ++k)
		{
			C[i][k]/=temp3;
		}
		
		for(IloInt j=0; j<N; ++j) //�ӵ����п�ʼ��
		{
			if(i!=j)
			{
				IloNum temp4=C[j][i];
				for(IloInt k=i; k<2*N; ++k)
				{
					C[j][k]-=(C[i][k]*temp4);
				}
			}
		}
	}
	for(IloInt i=0; i<N; ++i)
		for(IloInt j=N; j<2*N; ++j)
			A[i][j-N]=C[i][j];
}

void IloMutiply(const IloNumArray2 &B0l,IloNumExprArray &Psp,IloExprArray &Theta,const IloInt N)//��˻�����
{
	for(IloInt i=0;i<N;++i)
	{
		Theta[i]=IloExpr(env);
		for(IloInt j=0;j<N;++j)
		{
			Theta[i]+=B0l[i][j]*Psp[j];
		}
	}
}
/***********************************************���ݳ�ʼ��*****************************************************/
void define_data(IloEnv env)//���ݳ�ʼ��
{
	Master_Cplex.setParam(IloCplex::EpGap,0.11);//��Сһ��ͻ�out of memory
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>������ʼ��<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	for(IloInt i=0; i<NG; ++i)
	{	
		P[i]=IloNumVarArray(env,NT,0.0,IloInfinity,ILOFLOAT);		
		PS[i]=IloNumVarArray(env,NT,0.0,IloInfinity,ILOFLOAT);
		P1[i]=IloNumArray(env,NT);
		I[i]=IloBoolVarArray(env,NT);
		u[i]=IloNumArray(env,NT);
	}
	
	for(IloInt i=0;i<Node-1;++i)
	{
		B0l[i]=IloNumArray(env,Node-1);
	}
	
	for(IloInt t=0;t<NT;++t)
		Pwind[t]=IloNumArray(env,NW);
		
	for(IloInt s=0;s<set;++s)
	{
		Pswind[s]=IloNumArray2(env,NT);
		for(IloInt t=0;t<NT;++t)
		{
			Pswind[s][t]=IloNumArray(env,NW);
		}
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
	
	output<<endl<<"Pwind:"<<endl;
	ifstream windpower("wind_power.txt",ios::in);//����������
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
	ifstream wind_scenarios("wind_scenarios.txt",ios::in);//�����糡��
	if(!wind_scenarios)
	{	
		output<<"no such file! wind_scenarios.txt"<<endl;
	}
	for(IloInt s=0;s<set;++s)
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
			Betaa[i][l]=( f(Z[i][l+1],Unit[i])-f(Z[i][l],Unit[i]) )/( Z[i][l+1]-Z[i][l] );
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
		for(IloInt i=0;i<NG;++i)//���гɱ��������Ի����໪�����ģ�
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
		/*for(IloInt t=0; t<NT; ++t) //������ȫ����ʱ��Ԥ����
		{
			IloNumExpr expr1(env),expr2(env);
			for(IloInt i=0; i<NG; ++i)
			{
				expr1+=0.2*P[i][t];
				expr2+=0.2*P[i][t];
			}
			Master_Model.add(expr1>=20);
			Master_Model.add(expr2>=20);
			expr1.end();
			expr2.end();
		}*/
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
/************************************************************��ʼ������������***************************************************************/
		IloInvert(B0,B0l,Node-1);//�������

		IloInt while_interator=0;
		IloBool first_cut_is_add;
		
		while(1)
		{
			first_cut_is_add=IloFalse;
			output<<"���������������"<<while_interator++<<endl;
			Master_Cplex.extract(Master_Model);			
			Master_Cplex.solve();
			
			if (Master_Cplex.getStatus() == IloAlgorithm::Infeasible)//������
				output<<"Master Problem Have No Solution"<<endl;
			else 
				output<<"Master Problem Have Solution, the objective value is: "<<Master_Cplex.getObjValue()<<endl<<endl;
				
			for(IloInt i=0; i<NG; ++i)//���������ͣ״̬
			{
				Master_Cplex.getValues(I[i],u[i]);
				Master_Cplex.getValues(P[i],P1[i]);
			}
/*****************************************************�����ⰲȫԼ��********************************************************/
			for(IloInt t=0;t<NT;++t)//�����ⰲȫԼ����У��
			{
				output<<"Master Problem Security Constraint Check:"<<t+1<<endl;
				IloNumVarArray C1(env,Branch,0,IloInfinity,ILOFLOAT),C2(env,Branch,0,IloInfinity,ILOFLOAT);
				IloNumExprArray PL(env,Branch);//����֧·�ĳ���
				IloNumVarArray P_dual(env,Node-1);
				IloRangeArray range(env,Node-1);
				for(IloInt h=0;h<Branch;++h)
				{
					C1[h]=IloNumVar(env);
					C2[h]=IloNumVar(env);
					PL[h]=IloNumVar(env);
				}				
				IloModel Master_Model_Security(env);
				
				IloNumExpr obj(env);//Ŀ�꺯��
				for(IloInt h=0;h<Branch;++h)
				{
					obj+=(C1[h]+C2[h]);
				}
				Master_Model_Security.add(IloMinimize(env,obj));
				obj.end();
				
				IloNumExprArray Psp(env,Node-1);//ע�빦��
				IloExprArray Theta(env,Node);//���
				
				Theta[Node-1]=IloExpr(env);
					
				for(IloInt k=0; k<Node-1; ++k)//У�������н��ĳ���
				{
					P_dual[k]=IloNumVar(env);
					Psp[k]=IloNumExpr(env);
					
					IloInt i=0;
					for(;i<NG;++i)
					{
						if(Unit[i][0]-1==k)break;
					}
					if(i<NG)
						range[k]=IloRange::IloRange(env,P1[i][t],P_dual[k],P1[i][t]);	
					else
						range[k]=IloRange::IloRange(env,0,P_dual[k],0);

					Master_Model_Security.add(range[k]);
					Psp[k] += P_dual[k];
					if(Sw[k]>=0)
						Psp[k] += Pwind[t][ Sw[k] ];
					Psp[k] -= Sl[k]*Pload[t];
				}				
				IloMutiply(B0l,Psp,Theta,Node-1);
				//���㳱��
				for(IloInt h=0; h<Branch; ++h)
				{
					PL[h]=IloNumExpr(env);
					PL[h]+=(Theta[(IloInt)Info_Branch[h][0]-1]-Theta[(IloInt)Info_Branch[h][1]-1])/Info_Branch[h][3];
					Master_Model_Security.add(PL[h]-C1[h]<=Info_Branch[h][4]);
					Master_Model_Security.add(-PL[h]-C2[h]<=Info_Branch[h][4]);
				}
				Psp.end();
				Theta.end();
				
				IloCplex Master_Security_Cplex(Master_Model_Security);
				Master_Security_Cplex.solve();
					
				if (Master_Security_Cplex.getStatus() == IloAlgorithm::Infeasible)//������
					output << "Master_Security_Problem Have No Solution"<<endl;
				else
					output<<"Master_Security_Problem Have Solution, the objective value is: "<<Master_Security_Cplex.getObjValue()<<endl<<endl;
				
				if(Master_Security_Cplex.getObjValue()>limit)
				{
					IloNumExpr cut(env);
					cut+=Master_Security_Cplex.getObjValue();
					for(IloInt i=0;i<NG;++i)
					{
						if(Unit[i][0]-1>=Node-1)//���һ�����û��
							continue;
						else
							cut+=Master_Security_Cplex.getDual(range[(IloInt)Unit[i][0]-1])*(P[i][t]-P1[i][t]);//+(I[i][t]-u[i][t]));
					}
					Master_Model.add(cut<=0);
					cut.end();

					first_cut_is_add=IloTrue;
				}
				C1.end();
				C2.end();
				PL.end();
				P_dual.end();
				range.end();
				Master_Model_Security.end();
				Master_Security_Cplex.end();
			}
			
			if(first_cut_is_add==IloTrue)
			{
				output<<"�����ⰲȫԼ��Խ�ޣ�����������>>>"<<endl;
				continue;//�������Ӹ�����¿�ʼ
			}
/****************************************************************������У��****************************************************************/			
			IloInt s=0;
			IloBool second_cut_is_add;
			for(; s<set; ++s)//������У��
			{	
				output<<"Num_Senario_Iterator��������:"<<s+1<<endl;
				
				IloNumVarArray S1S(env,NT,0,IloInfinity,ILOFLOAT),S2S(env,NT,0,IloInfinity,ILOFLOAT);//�ڶ��׶��ɳڱ���
								//S3S(env,NT,0,IloInfinity,ILOFLOAT),S4S(env,NT,0,IloInfinity,ILOFLOAT);
				IloArray<IloNumVarArray> S5S(env,NG);
				for(IloInt i=0; i<NG; ++i)
					S5S[i]=IloNumVarArray(env,NT,0.0,IloInfinity,ILOFLOAT);
				
				second_cut_is_add=IloFalse;

				IloModel Sub_Model(env,"Sub_Model");
/********************************������Ŀ�꺯��**************************************/
				IloNumExpr VS(env);
				for(IloInt t=0; t<NT; ++t)
				{
					VS+=(S1S[t]+S2S[t]);					
					for(IloInt i=0; i<NG; ++i)
						VS+=S5S[i][t];
				}
				Sub_Model.add(IloMinimize(env,VS));
				VS.end();
/*****************************������ϵͳƽ��Լ��*********************************/
				for(IloInt t=0; t<NT; ++t)
				{
					IloNumExpr expr1(env);
					for(IloInt i=0; i<NG; ++i)
						expr1+=PS[i][t];
					expr1+=S1S[t]-S2S[t];
					
					IloNum wind(0);
					for(IloInt w=0;w<NW;++w)
						wind+=Pswind[s][t][w];
						
					Sub_Model.add(expr1==Pload[t]-wind);
					expr1.end();
				}
/*******************************�����ⱸ��Լ��**************************************/
				/*for(IloInt t=0; t<NT; ++t)
				{
					IloNumExpr expr2(env),expr3(env);
					for(IloInt i=0; i<NG; ++i)
					{
						expr2+=0.2*PS[i][t];
						expr3+=0.2*PS[i][t];
					}
					Sub_Model.add(expr2+S3S[t]>=20);
					Sub_Model.add(expr3+S4S[t]>=20);
					expr2.end();
					expr3.end();
				}*/
/***************************�����⽨����ΪԼ��*************************************/
				IloArray<IloRangeArray> ramp1(env,NG),ramp2(env,NG); 
				for(IloInt i=0; i<NG; ++i)//����ֵת��Ϊ����Լ��
				{
					ramp1[i]=IloRangeArray(env,NT);
					ramp2[i]=IloRangeArray(env,NT);
					for(IloInt t=0; t<NT; ++t)
					{
						ramp1[i][t]=IloRange::IloRange(env,-IloInfinity,PS[i][t]-S5S[i][t],detaa+P1[i][t],"ramp1");
						ramp2[i][t]=IloRange::IloRange(env,-IloInfinity,-PS[i][t]-S5S[i][t],detaa-P1[i][t],"ramp2");
						
						Sub_Model.add(ramp1[i][t]);
						Sub_Model.add(ramp2[i][t]);
					}
				}
/**************************������������������Լ��*******************************/
				IloArray<IloRangeArray> range_limit1(env,NG),range_limit2(env,NG);
				for(IloInt i=0; i<NG; ++i)
				{
					range_limit1[i]=IloRangeArray(env,NT);
					range_limit2[i]=IloRangeArray(env,NT);
					for(IloInt t=0; t<NT; ++t)
					{
						range_limit1[i][t]=IloRange::IloRange(env,-IloInfinity,PS[i][t],Unit[i][2]*u[i][t],"power_limits1");
						range_limit2[i][t]=IloRange::IloRange(env,-IloInfinity,-PS[i][t],-Unit[i][1]*u[i][t],"power_limits2");
						Sub_Model.add(range_limit1[i][t]);
						Sub_Model.add(range_limit2[i][t]);
					}
				}
/*************************************��Ӹ�************************************/
				IloCplex Sub_Cplex(Sub_Model);
				label1:
				Sub_Cplex.extract(Sub_Model);
				Sub_Cplex.solve();
				
				IloNum objvalue(Sub_Cplex.getObjValue());
				if(Sub_Cplex.getStatus() == IloAlgorithm::Infeasible)//������
				{
					output<< "Senario_Sub_Problem No Solution" << endl;
					output << "���� "<< s+1 <<" �޽�,����ѭ�����������޽⡭��"<<endl;
					goto label2;
				}
				else 
					output<<"Senario_Sub_Problem Have Solution, and objective is: "<<objvalue<<endl<<endl;
				
				if(objvalue>limit)//����������
				{
					IloNumExpr expr(env);
					expr+=objvalue;
					for(IloInt i=0;i<NG;++i)
						for(IloInt t=0;t<NT;++t)
						{
							expr+=0.1*(Sub_Cplex.getDual(range_limit1[i][t])*Unit[i][2]-Sub_Cplex.getDual(range_limit2[i][t])*Unit[i][1])
									*(I[i][t]-u[i][t])+0.1*(Sub_Cplex.getDual(ramp1[i][t])-Sub_Cplex.getDual(ramp2[i][t]))*(P[i][t]-P1[i][t]);
						}
					Master_Model.add(expr<=0);
					range_limit1.end();
					range_limit2.end();
					ramp1.end();
					ramp2.end();
					expr.end();
					output<<"����У�������ⲻ���㣬����������>>>"<<endl;
					break;//��������У��������
				}
				else//����ʼ��ȫУ��
				{
					IloNumArray2 PS1(env,NG);//��ȡ��Ӧ�����µĳ���
					for(IloInt i=0;i<NG;++i)
					{
						PS1[i]=IloNumArray(env,NT);
						Sub_Cplex.getValues(PS[i],PS1[i]);
					}
					
					for(IloInt t=0;t<NT;++t)//������ʱ�ν��а�ȫУ��
					{
						output<<"Sub Problem Security Constraint Check:"<<t+1<<endl;
						IloNumVarArray C1(env,Branch,0,IloInfinity,ILOFLOAT),C2(env,Branch,0,IloInfinity,ILOFLOAT);
						IloNumExprArray PL(env,Branch);
						IloNumVarArray P_dual(env,Node-1);
						IloRangeArray range(env,Node-1);
						for(IloInt l=0;l<Branch;++l)
						{
							C1[l]=IloNumVar(env);
							C2[l]=IloNumVar(env);
							PL[l]=IloNumExpr(env);
						}						
						IloModel SubModel_Security(env);//ֻ�Ǿֲ�����
						
						IloNumExpr obj(env);//Ŀ�꺯��
						for(IloInt l=0;l<Branch;++l)
							obj+=C1[l]+C2[l];
						SubModel_Security.add(IloMinimize(env,obj));
						obj.end();
						
						IloNumExprArray Psp(env,Node-1);//ע�빦��
						IloExprArray Theta(env,Node);//���
						
						Theta[Node-1]=IloExpr(env);
							
						for(IloInt k=0; k<Node-1; ++k)//��Psp
						{
							P_dual[k]=IloNumVar(env);
							Psp[k]=IloNumExpr(env);
							
							IloInt i=0;
							for(;i<NG;++i)
							{
								if(Unit[i][0]-1==k)break;
							}
							if(i<NG)
								range[k]=IloRange::IloRange(env,PS1[i][t],P_dual[k],PS1[i][t]);
							else
								range[k]=IloRange::IloRange(env,0,P_dual[k],0);
								
							SubModel_Security.add(range[k]);
							Psp[k] += P_dual[k];
							// Psp[k] += (Sw[k]+1)*Pswind[s][t];
							if(Sw[k]>=0)
								Psp[k] += Pswind[s][t][ Sw[k] ];
							Psp[k] -= Sl[k]*Pload[t];
						}						
						IloMutiply(B0l,Psp,Theta,Node-1);
						//���㳱��
						for(IloInt h=0; h<Branch; ++h)
						{
							PL[h]=IloNumExpr(env);			
							PL[h]+=(Theta[(IloInt)Info_Branch[h][0]-1]-Theta[(IloInt)Info_Branch[h][1]-1])/Info_Branch[h][3];							
							SubModel_Security.add(PL[h]-C1[h]<=Info_Branch[h][4]);
							SubModel_Security.add(-PL[h]-C2[h]<=Info_Branch[h][4]);
						}
						IloCplex Senario_Security_Cplex(SubModel_Security);
						Senario_Security_Cplex.solve();
						IloNum security_objvalue(Senario_Security_Cplex.getObjValue());
						Psp.end();
						Theta.end();
						PL.end();
						
						if (Senario_Security_Cplex.getStatus() == IloAlgorithm::Infeasible)//������
							output<<"Senario_Security_Problem one No Solution" << endl;
						else 
							output<<"Senario_Security_Problem Have Solution, and objective is: "<<security_objvalue<<endl<<endl;

						if(security_objvalue>limit)
						{
							IloNumExpr cut(env);
							cut+=security_objvalue;
							
							for(IloInt i=0;i<NG;++i)
							{
								if(Unit[i][0]-1>=Node-1)
									continue;
								else
									cut+=Senario_Security_Cplex.getDual(range[(IloInt)Unit[i][0]-1])*(PS[i][t]-PS1[i][t]);//��ȫԼ����
							}
							Sub_Model.add(cut<=0);
							cut.end();							
							second_cut_is_add=IloTrue;
						}
						C1.end();
						C2.end();
						PL.end();
						P_dual.end();
						range.end();
						SubModel_Security.end();
						Senario_Security_Cplex.end();
					}
					
					if(second_cut_is_add==IloTrue)//���ص�ǰ����У��
					{
						output<<"goto label1������У�������ⰲȫԼ��Խ�ޣ����ص�ǰ������>>>"<<endl;
						second_cut_is_add=IloFalse;
						goto label1;
					}
					PS1.end();
				}
				Sub_Model.end();
				Sub_Cplex.end();
			}
			if(s==set)break;
		}
/************************************************************�����ʾ����**************************************************/
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
		label2: ;
	}
	catch(IloException& ex)//�쳣����
	{
		output<<"Error: "<<ex<<endl;
	}
	catch(...)
	{
		output<<"Error: Unknown exception caught!" << endl;
	}
	
	finish=clock();
	totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	output<<"totaltime: "<<totaltime<<"s"<<endl<<endl;
	output.close();	
	return 0;
}
