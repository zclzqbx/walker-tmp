/*              不规则概率分布风电机组组合
 *1）以风电出力计算点作为风电出力，暂时不考虑风电的不确定性
 *
 */
#include<fstream>
#include<iostream>
#include<time.h>
#include"function.h"

using namespace std;
ILOSTLBEGIN

const IloInt NG=54;//常规机组台数
const IloInt NT=24;//调度周期时段数
const IloInt NL=4;//目标函数线性化分段数,理论上分段越多越精确，但会增加计算量
const IloInt NW=1;//风电场个数
const IloInt Node=118;//结点
const IloInt Branch=186;//支路

IloEnv env;//运行环境，运行结束之后要关闭
IloModel Master_Model(env,"Master_Model");//生成模型
IloCplex Master_Cplex(Master_Model);//创建CPLEX环境

ofstream output("output_fixedwind.txt",ios::ate);

IloArray<IloNumVarArray>  P(env,NG),PS(env,NG);//机组及相应场景下的出力
IloArray<IloBoolVarArray>  I(env,NG);//机组启停状态,二维决策变量
IloArray<IloNumVarArray> PL(env,NT);//NT个时段各条支路的潮流
IloArray<IloNumVarArray> P_dual(env,NT);
IloNumArray2 P1(env,NG),u(env,NG);//用于存储P、I的值

IloArray< IloArray<IloNumVarArray> > Deta(env,NG);//分段线性化的相关量
IloNumArray2 Z(env,NG);
IloNumArray2 Betaa(env,NG);

IloNumArray2 Unit(env,NG),Info_Branch(env,Branch);//机组及支路信息
IloNumArray Pload(env,NT);//负荷及风电出力
IloNumArray R(env,NT);//备用
IloNumArray2 Pwind(env,NT);//考虑多个风电场

///////////////////////////////与安全约束有关的量///////////////////////////////////////
IloNumArray2 B0(env,Node-1),B0l(env,Node-1);//导纳矩阵

IloNumArray Sl(env,Node);//每个节点的负荷状态,每个节点所占负荷百分比
IloIntArray Sw(env,Node);//常规风电场分布

/***********************************************数据初始化*****************************************************/
void define_data(IloEnv env)//数据初始化,对全局变量进行赋值
{
	Master_Cplex.setParam(IloCplex::EpGap,0.01);//再小一点就会out of memory
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>变量初始化<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
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
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>变量初始化结束<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>数据读取<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	output<<"Unit:"<<endl;
	ifstream unit_information("unit_information.txt",ios::in);//常规机组分布
	if(!unit_information)
	{
		output<<"no such file! unit_information.txt"<<endl;
	}
	for(IloInt i=0;i<NG;++i)//读入机组信息
	{
		Unit[i]=IloNumArray(env,12);//12条信息:结点、Pmin、Pmax、MinOn、MinOff、RU、RD、Startup Cost、c、b、a、price;
		for(IloInt k=0; k<12; ++k)
		{
			unit_information>>Unit[i][k];
			output<<Unit[i][k]<<"   ";
		}
		output<<endl;
	}
	unit_information.close();

	output<<endl<<endl<<"Pload:"<<endl;
	ifstream load("load.txt",ios::in);//负荷
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
	ifstream reserve("reserve.txt",ios::in);//备用
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
	ifstream windpower("wind_power.txt",ios::in);//风电计算场景
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
	ifstream BFile("Brach_File.txt",ios::in);//读取风电支路信息:起点，终点，电阻，电抗以及潮流上限
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
	ifstream load_locate("load_locate.txt",ios::in);//负荷分布
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
	ifstream wind_locate("wind_locate.txt",ios::in);//风电场分布
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
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>数据读取结束<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>目标函数线性化处理<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
    for(IloInt i=0;i<NG;++i)//分段线性化,可以考虑使用其他线性化方法，这种太耗时
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
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>目标函数线性化处理结束<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>给矩阵B0赋值<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	for(IloInt k=0; k<Node-1; ++k)//非对角线,依然以最后一个结点作为平衡结点
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
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>给矩阵B0赋值结束<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
}
/**************************************************程序入口****************************************************/
int main()
{
	clock_t start,finish;
	double totaltime;
	start=clock();
	
	try
	{
		define_data(env);//首先初始化全局变量
/*************************************************************主问题目标函数*******************************************************/
		IloNumExpr Cost(env);
		
		for(IloInt i=0;i<NG;++i)//运行成本及其线性化
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
		
		for(IloInt i=0; i<NG; ++i)//开机成本
			for(IloInt t=1; t<NT; ++t)
			{
				Cost+=Unit[i][11]*Unit[i][7]*I[i][t]*(1-I[i][t-1]);
			}
			
		for(IloInt i=0; i<NG; ++i)//停机成本
			for(IloInt t=1; t<NT; ++t)
			{
				Cost+=Unit[i][11]*Unit[i][7]*I[i][t-1]*(1-I[i][t]);
			}
		
		Master_Model.add(IloMinimize(env,Cost));
		Cost.end();
/********************************************************机组出力上下限约束**************************************************/
		for(IloInt i=0; i<NG; ++i)
			for(IloInt t=0; t<NT; ++t)
			{
				Master_Model.add(P[i][t]<=Unit[i][2]*I[i][t]);
				Master_Model.add(P[i][t]>=Unit[i][1]*I[i][t]);
			}
/*********************************************************机组爬坡约束******************************************************/
		for(IloInt i=0; i<NG; ++i)
			for(IloInt t=1; t<NT; ++t)
			{
				Master_Model.add(P[i][t]-P[i][t-1]<=(1-I[i][t]*(1-I[i][t-1]))*Unit[i][5]+I[i][t]*(1-I[i][t-1])*Unit[i][1]);
				Master_Model.add(P[i][t-1]-P[i][t]<=(1-I[i][t-1]*(1-I[i][t]))*Unit[i][6]+I[i][t-1]*(1-I[i][t])*Unit[i][1]);
			}
/*********************************************************机组功率平衡约束**************************************************/
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
/**********************************************************备用约束********************************************************/
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
/***********************************************************最小开机时间约束*************************************************/
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
/*******************************************************最小停机时间约束**************************************************/
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
/*************************************************安全约束*********************************************/
		IloInvert(env,B0,B0l,Node-1);//求逆矩阵
		for(int t=0;t<NT;++t)//校验所有时段的所有支路的潮流
		{
			IloNumExprArray Psp(env,Node-1);//注入功率
			IloExprArray Theta(env,Node);//相角				
			
			Theta[Node-1]=IloExpr(env);
				
			for(IloInt k=0; k<Node-1; ++k)//校验了所有结点的潮流
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
			//计算潮流
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
		if (Master_Cplex.getStatus() == IloAlgorithm::Infeasible)//输出结果
		{
			output<<"Master Problem Have No Solution"<<endl;
			goto lable2;
		}
/************************************************************输出显示过程**************************************************/
		output<<"Cost:"<<Master_Cplex.getObjValue()<<endl;
		for(IloInt t=0;t<NT;++t)
		{
			output<<"时段 "<<t+1<<" 机组的启停(I)："; 
			for(IloInt i=0;i<NG;++i)
			{
				output<<Master_Cplex.getValue(I[i][t])<<"   ";
			}
			output<<endl;
			output<<"时段 "<<t+1<<" 机组的出力(P)：";			
			for(IloInt i=0;i<NG;++i)
			{
				output<<Master_Cplex.getValue(P[i][t])<<"   ";
			}
			output<<endl;
			output<<"时段 "<<t+1<<" 线路潮流(PL)：";			
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
	catch(IloException& ex)//异常捕获
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
