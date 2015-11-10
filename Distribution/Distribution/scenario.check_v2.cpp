/*                机组启停满意度校验
 *1）在机组启停方式确定的情况下，针对不同风电场，
 *校验其是否能平衡风功率的波动；
 *
 */
#include<ilcplex/ilocplex.h>
#include<fstream>
#include<iostream>
#include<time.h>

using namespace std;
ILOSTLBEGIN

const IloInt NG=54;//常规机组台数
const IloInt NT=24;//调度周期时段数
const IloInt NW=1;//风电场个数
const IloInt Node=118;//结点
const IloInt Branch=186;//支路

IloEnv env;
IloModel Check_Model(env,"Check_Model");
IloCplex Check_Cplex(Check_Model);
ofstream output("Output_Check_Result.txt",ios::out | ios::app);

/**************决策变量****************/
IloArray<IloNumVarArray> P(env,NG);
IloNumVarArray S1(env,NT,0,IloInfinity,ILOFLOAT),
				S2(env,NT,0,IloInfinity,ILOFLOAT);
IloArray<IloNumVarArray> S3(env,NG);
IloNumVarArray S4(env,Branch,0,IloInfinity,ILOFLOAT),
				S5(env,Branch,0,IloInfinity,ILOFLOAT);


IloNumArray2 Unit(env,NG),Pwind(env,NW);//机组信息及预测场景
IloNumArray2 I(env,NG),Pre(env,NG);//机组状态及预测场景下出力
IloNumArray Pload(env,NT),deta(env,NG); 

///////////////////////////////与安全约束有关的量///////////////////////////////////////
IloNumArray2 Info_Branch(env,Branch);//风电支路信息,包括支路起点，终点，电阻，电抗以及潮流上限
IloNumArray Sl(env,Node);//每个节点的负荷状态,每个节点所占负荷百分比
IloIntArray Sw(env,Node);//常规风电场分布
IloNumArray2 B0(env,Node-1),B0l(env,Node-1);//导纳矩阵及其逆矩阵
//统一以结点Node-1作为平衡节点，该节点相角为零，则B0为（Node-1）X（Node-1）的矩阵

void data_define(IloEnv env)
{
	for(IloInt i=0;i<NG;++i)
	{
		P[i]=IloNumVarArray(env,NT,0,IloInfinity,ILOFLOAT);
		S3[i]=IloNumVarArray(env,NT,0,IloInfinity,ILOFLOAT);
	}		
/************************************数据读取******************************************/
	ifstream input_wind("wind_scenarios.txt",ios::in);//误差场景，由蒙特卡洛方法生成
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

	ifstream input_commitment("uc.txt",ios::in);//机组启停状态
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
	
	ifstream input_preoutput("Ppower.txt",ios::in);//预测场景下出力
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

	ifstream input_load("load.txt",ios::in);//读取负荷信息
	if(!input_load)
	{
		cerr<<"load.txt is not exist!!!"<<endl;
	}
	for(IloInt t=0; t<NT; ++t)
	{
		input_load>>Pload[t];
	}
	input_load.close();
			
	ifstream input_BFile("Brach_File.txt",ios::in);//读取风电支路信息:起点，终点，电阻，电抗以及潮流上限
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
	
	ifstream input_unitinfo("unit_information.txt",ios::in);//读入机组信息
	if(!input_unitinfo)
	{
		cerr<<"unit_information.txt is not exist!!!"<<endl;
	}
	for(IloInt i=0;i<NG;++i)
	{
		Unit[i]=IloNumArray(env,12);//12条信息:结点、Pmin、Pmax、MinOn、MinOff、RU、RD、Startup Cost、c、b、a、price;
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
	
	ifstream input_wind_locate("wind_locate.txt",ios::in);//风场分布
	if(!input_wind_locate)
	{
		cerr<<"wind_locate is not exist!!!"<<endl;
	}
	for(IloInt k=0; k<Node; ++k)
	{
		input_wind_locate>>Sw[k];
	}
	input_wind_locate.close();
	
	ifstream input_load_locate("load_locate.txt",ios::in);//负荷分布
	if(!input_load_locate)
	{
		cerr<<"load_locate is not exist!!!"<<endl;
	}
	for(IloInt k=0; k<Node; ++k)
	{
		input_load_locate>>Sl[k];
	}
	input_load_locate.close();
	/////给矩阵B0赋值
	for(IloInt k=0; k<Node-1; ++k)//非对角线
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
	for(IloInt k=0; k<Node-1; ++k)//对角线
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
	/////赋值完毕	
	for(IloInt k=0; k<Node-1; ++k)
		B0l[k]=IloNumArray(env,Node-1);
}

void IloInvert(const IloNumArray2 &B0,IloNumArray2 &A,IloInt N)//求逆矩阵
{
	IloNumArray2 C(env,N);//扩展矩阵

	for(IloInt k=0; k<N; ++k)
	{
		C[k]=IloNumArray(env,2*N);
	}

	for(IloInt i=0; i<N; ++i) //构造单位阵
		for(IloInt j=0; j<N; ++j)
		{
			if(i==j)
				A[i][j]=1;
			else
				A[i][j]=0;
		}

	for(IloInt i=0; i<N; ++i) //求增广矩阵
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
			for(; j<N; ++j) //从本行的下一行开始搜索不为零的行
			{
				if(C[j][i]!=0)//找到后进行交换
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
		
		IloNum temp3=C[i][i];//保存首非零值
		
		for(IloInt k=i; k<2*N; ++k)
		{
			C[i][k]/=temp3;
		}
		
		for(IloInt j=0; j<N; ++j) //从第零行开始减
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

void IloMutiply(const IloNumArray2 &B0l,const IloExprArray &Psp,IloExprArray &Theta,IloInt N)//求乘积矩阵
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
/************************************变量初始化*************************************/		
		data_define(env);
		IloInvert(B0,B0l,Node-1);//求逆矩阵
/************************************模型******************************************/		
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
		Check_Model.add(IloMinimize(env,Cost));//目标函数
		
		for(IloInt i=0; i<NG; ++i)//机组出力上下限约束
			for(IloInt t=0; t<NT; ++t)
			{
				Check_Model.add(P[i][t]<=Unit[i][2]*I[i][t]);
				Check_Model.add(P[i][t]>=Unit[i][1]*I[i][t]);
			}
		
		// for(IloInt i=0; i<NG; ++i)//机组爬坡约束,此条件未松弛
			// for(IloInt t=1; t<NT; ++t)
			// {
				// Check_Model.add(P[i][t]-P[i][t-1]<=(1-I[i][t]*(1-I[i][t-1]))*Unit[i][5]+I[i][t]*(1-I[i][t-1])*Unit[i][1]);
				// Check_Model.add(P[i][t-1]-P[i][t]<=(1-I[i][t-1]*(1-I[i][t]))*Unit[i][6]+I[i][t-1]*(1-I[i][t])*Unit[i][1]);
			// }
		
		for(IloInt t=0; t<NT; ++t)//机组功率平衡约束
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

		for(IloInt i=0;i<NG;++i)//机组出力爬坡约束
			for(IloInt t=0;t<NT;++t)
			{
				Check_Model.add(P[i][t]-Pre[i][t]-S3[i][t]<=deta[i]);
				Check_Model.add(P[i][t]-Pre[i][t]+S3[i][t]>=-deta[i]);
			}
/************************************安全约束******************************************/		
		for(IloInt t=0; t<NT; ++t)
		{
			IloExprArray Pf(env,Branch);//潮流
			IloExprArray Psp(env,Node-1);//注入功率
			IloExprArray Theta(env,Node);//相角
			
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
			
			//计算潮流
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
/************************************模型求解及输出******************************************/
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
				// output<<"机组 "<<i+1<<" 第 "<<t+1<<" 时段状态:"<<Check_Cplex.getValue(I[i][t])<<"   出力:";
				// output<<Check_Cplex.getValue(P[i][t])<<endl;
			// }
			// output<<endl;
		// }
		Check_Model.end();
		Check_Cplex.end();
		env.end();
	}
	catch(IloException& ex)//异常捕获
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