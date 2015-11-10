/*                ������ͣ�����У��
 *1���ڻ�����ͣ��ʽȷ��������£���Բ�ͬ��糡��
 *У�����Ƿ���ƽ��繦�ʵĲ�����
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
const IloInt NW=1;//��糡����
const IloInt Node=118;//���
const IloInt Branch=186;//֧·

IloEnv env;
IloModel Check_Model(env,"Check_Model");
IloCplex Check_Cplex(Check_Model);
ofstream output("Output_Check_Result.txt",ios::out | ios::app);

/**************���߱���****************/
IloArray<IloNumVarArray> P(env,NG);
IloNumVarArray S1(env,NT,0,IloInfinity,ILOFLOAT),
				S2(env,NT,0,IloInfinity,ILOFLOAT);
IloArray<IloNumVarArray> S3(env,NG);
IloNumVarArray S4(env,Branch,0,IloInfinity,ILOFLOAT),
				S5(env,Branch,0,IloInfinity,ILOFLOAT);


IloNumArray2 Unit(env,NG),Pwind(env,NW);//������Ϣ��Ԥ�ⳡ��
IloNumArray2 I(env,NG),Pre(env,NG);//����״̬��Ԥ�ⳡ���³���
IloNumArray Pload(env,NT),deta(env,NG); 

///////////////////////////////�밲ȫԼ���йص���///////////////////////////////////////
IloNumArray2 Info_Branch(env,Branch);//���֧·��Ϣ,����֧·��㣬�յ㣬���裬�翹�Լ���������
IloNumArray Sl(env,Node);//ÿ���ڵ�ĸ���״̬,ÿ���ڵ���ռ���ɰٷֱ�
IloIntArray Sw(env,Node);//�����糡�ֲ�
IloNumArray2 B0(env,Node-1),B0l(env,Node-1);//���ɾ����������
//ͳһ�Խ��Node-1��Ϊƽ��ڵ㣬�ýڵ����Ϊ�㣬��B0Ϊ��Node-1��X��Node-1���ľ���

void data_define(IloEnv env)
{
	for(IloInt i=0;i<NG;++i)
	{
		P[i]=IloNumVarArray(env,NT,0,IloInfinity,ILOFLOAT);
		S3[i]=IloNumVarArray(env,NT,0,IloInfinity,ILOFLOAT);
	}		
/************************************���ݶ�ȡ******************************************/
	ifstream input_wind("wind_scenarios.txt",ios::in);//�����������ؿ��巽������
	for(IloInt w=0;w<NW;++w)
	{
		Pwind[w]=IloNumArray(env,NT);
		for(IloInt t=0;t<NT;++t)
		{
			input_wind>>Pwind[w][t];
			// output<<Pwind[w][t]<<"   ";
		}
		// output<<endl;
	}
	input_wind.close();

	ifstream input_commitment("uc.txt",ios::in);//������ͣ״̬
	if(!input_commitment)
	{
		cerr<<"uc.txt is not exist!!!"<<endl;
	}
	for(IloInt i=0;i<NG;++i)
	{
		I[i]=IloNumArray(env,NT);
		for(IloInt t=0;t<NT;++t)
		{
			input_commitment>>I[i][t];
		}
	}
	input_commitment.close();
	
	ifstream input_preoutput("Ppower.txt",ios::in);//Ԥ�ⳡ���³���
	if(!input_preoutput)
	{
		cerr<<"Ppower.txt is not exist"<<endl;
	}
	for(IloInt i=0;i<NG;++i)
	{
		Pre[i]=IloNumArray(env,NT);
		for(IloInt t=0;t<NT;++t)
		{
			input_preoutput>>Pre[i][t];
		}
	}
	input_preoutput.close();

	ifstream input_load("load.txt",ios::in);//��ȡ������Ϣ
	if(!input_load)
	{
		cerr<<"load.txt is not exist!!!"<<endl;
	}
	for(IloInt t=0; t<NT; ++t)
	{
		input_load>>Pload[t];
	}
	input_load.close();
			
	ifstream input_BFile("Brach_File.txt",ios::in);//��ȡ���֧·��Ϣ:��㣬�յ㣬���裬�翹�Լ���������
	if(!input_BFile)
	{	
		cerr<<"Brach_File.txt is not exist!!!"<<endl;
	}
	for(IloInt k=0; k<Branch; ++k)
	{	
		Info_Branch[k]=IloNumArray(env,5);
		for(IloInt h=0; h<5; ++h)
		{
			input_BFile>>Info_Branch[k][h];
		}
	}
	input_BFile.close();
	
	ifstream input_unitinfo("unit_information.txt",ios::in);//���������Ϣ
	if(!input_unitinfo)
	{
		cerr<<"unit_information.txt is not exist!!!"<<endl;
	}
	for(IloInt i=0;i<NG;++i)
	{
		Unit[i]=IloNumArray(env,12);//12����Ϣ:��㡢Pmin��Pmax��MinOn��MinOff��RU��RD��Startup Cost��c��b��a��price;
		for(IloInt k=0; k<12; ++k)
		{
			input_unitinfo>>Unit[i][k];
		}
	}
	input_unitinfo.close();
	
	//output<<endl<<endl<<"deta:"<<endl;
	for(IloInt i=0;i<NG;i++)
	{
		deta[i]=Unit[i][5]/10;
		//output<<deta[i]<<"   ";
	}
	
	ifstream input_wind_locate("wind_locate.txt",ios::in);//�糡�ֲ�
	if(!input_wind_locate)
	{
		cerr<<"wind_locate is not exist!!!"<<endl;
	}
	for(IloInt k=0; k<Node; ++k)
	{
		input_wind_locate>>Sw[k];
	}
	input_wind_locate.close();
	
	ifstream input_load_locate("load_locate.txt",ios::in);//���ɷֲ�
	if(!input_load_locate)
	{
		cerr<<"load_locate is not exist!!!"<<endl;
	}
	for(IloInt k=0; k<Node; ++k)
	{
		input_load_locate>>Sl[k];
	}
	input_load_locate.close();
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
}

void IloInvert(const IloNumArray2 &B0,IloNumArray2 &A,IloInt N)//�������
{
	IloNumArray2 C(env,N);//��չ����

	for(IloInt k=0; k<N; ++k)
	{
		C[k]=IloNumArray(env,2*N);
	}

	for(IloInt i=0; i<N; ++i) //���쵥λ��
		for(IloInt j=0; j<N; ++j)
		{
			if(i==j)
				A[i][j]=1;
			else
				A[i][j]=0;
		}

	for(IloInt i=0; i<N; ++i) //���������
		for(IloInt j=0; j<2*N; ++j)
		{
			if(j<N)
				C[i][j]=B0[i][j];
			else
				C[i][j]=A[i][j-N];
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

void IloMutiply(const IloNumArray2 &B0l,const IloExprArray &Psp,IloExprArray &Theta,IloInt N)//��˻�����
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

int main()
{
	clock_t start,finish;
	double totaltime;
	start=clock();	
	
	try{
/************************************������ʼ��*************************************/		
		data_define(env);
		IloInvert(B0,B0l,Node-1);//�������
/************************************ģ��******************************************/		
		IloNumExpr Cost(env);
		for(IloInt t=0;t<NT;++t)
		{
			Cost+=(S1[t]+S2[t]);
			for(IloInt i=0;i<NG;++i)			
			{
				Cost+=S3[i][t];
			}
		}
		for(IloInt h=0;h<Branch;++h)
		{
			Cost+=(S4[h]+S5[h]);
		}
		Check_Model.add(IloMinimize(env,Cost));//Ŀ�꺯��
		
		for(IloInt i=0; i<NG; ++i)//�������������Լ��
			for(IloInt t=0; t<NT; ++t)
			{
				Check_Model.add(P[i][t]<=Unit[i][2]*I[i][t]);
				Check_Model.add(P[i][t]>=Unit[i][1]*I[i][t]);
			}
		
		// for(IloInt i=0; i<NG; ++i)//��������Լ��,������δ�ɳ�
			// for(IloInt t=1; t<NT; ++t)
			// {
				// Check_Model.add(P[i][t]-P[i][t-1]<=(1-I[i][t]*(1-I[i][t-1]))*Unit[i][5]+I[i][t]*(1-I[i][t-1])*Unit[i][1]);
				// Check_Model.add(P[i][t-1]-P[i][t]<=(1-I[i][t-1]*(1-I[i][t]))*Unit[i][6]+I[i][t-1]*(1-I[i][t])*Unit[i][1]);
			// }
		
		for(IloInt t=0; t<NT; ++t)//���鹦��ƽ��Լ��
		{
			IloNumExpr fire(env);
			IloNum wind(0);
			for(IloInt i=0; i<NG; ++i)
			{
				fire+=P[i][t];
			}
			for(IloInt w=0;w<NW;++w)
			{
				wind+=Pwind[w][t];
			}
			Check_Model.add(fire+S1[t]-S2[t]==Pload[t]-wind);
			fire.end();
		}

		for(IloInt i=0;i<NG;++i)//�����������Լ��
			for(IloInt t=0;t<NT;++t)
			{
				Check_Model.add(P[i][t]-Pre[i][t]-S3[i][t]<=deta[i]);
				Check_Model.add(P[i][t]-Pre[i][t]+S3[i][t]>=-deta[i]);
			}
/************************************��ȫԼ��******************************************/		
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
				{
					Psp[k] +=Pwind[ Sw[k] ][t];
				}
				Psp[k] -= Sl[k]*Pload[t];
			}
			IloMutiply(B0l,Psp,Theta,Node-1);
			
			//���㳱��
			for(IloInt h=0; h<Branch; ++h)
			{
				Pf[h]=IloExpr(env);			
				Pf[h]+=(Theta[(IloInt)Info_Branch[h][0]-1]-Theta[(IloInt)Info_Branch[h][1]-1])/Info_Branch[h][3];
				Check_Model.add(Pf[h]-S4[h]<=Info_Branch[h][4]);
				Check_Model.add(Pf[h]+S5[h]>=-Info_Branch[h][4]);
			}
			Pf.end();
			Psp.end();
			Theta.end();
		}
/************************************ģ����⼰���******************************************/
		Check_Cplex.solve();
		if(Check_Cplex.getStatus() == IloAlgorithm::Infeasible)
			output<< "No Solution" << endl;
		
		if(Check_Cplex.getObjValue()>0.2)
			output<<1<<endl;
		else
			output<<0<<endl;
		//output<</*endl<<"Cost:"<<*/Check_Cplex.getObjValue()<<endl;
		// for(IloInt i=0; i<NG; ++i)
		// {
			// for(IloInt t=0; t<NT; ++t)
			// {
				// output<<"���� "<<i+1<<" �� "<<t+1<<" ʱ��״̬:"<<Check_Cplex.getValue(I[i][t])<<"   ����:";
				// output<<Check_Cplex.getValue(P[i][t])<<endl;
			// }
			// output<<endl;
		// }
		Check_Model.end();
		Check_Cplex.end();
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
	
	finish=clock();
	totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	// output<<"totaltime: "<<totaltime<<"s"<<endl<<endl;
	output.close();	
	return 0;
}