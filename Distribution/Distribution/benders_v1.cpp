/*              不规则概率分布风电机组组合
 *1）
 *
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
const IloInt NL=4;//目标函数线性化分段数,理论上分段越多越精确，但会增加计算量
const IloInt NW=1;//风电场个数
const IloInt set=10;//场景个数
const IloInt Node=118;//结点
const IloInt Branch=186;//支路
const IloNum limit=0.1;//返回割阀值
const IloNum Cost_Modify=1;

IloEnv env;//运行环境，运行结束之后要关闭
IloModel Master_Model(env,"Master_Model");//生成模型
IloCplex Master_Cplex(Master_Model);//创建CPLEX环境

ofstream output("output_benders.txt",ios::ate);

IloArray<IloNumVarArray>  P(env,NG),PS(env,NG);//机组及相应场景下的出力
IloArray<IloBoolVarArray>  I(env,NG);//机组启停状态,二维决策变量
IloNumArray2 P1(env,NG),u(env,NG);//用于存储P、I的值

IloArray< IloArray<IloNumVarArray> > Deta(env,NG);//分段线性化的相关量
IloNumArray2 Z(env,NG);
IloNumArray2 Betaa(env,NG);

IloNumArray2 Unit(env,NG),Info_Branch(env,Branch);//机组及支路信息
IloNumArray /*Pwind(env,NT),*/Pload(env,NT);//负荷及风电出力
IloNumArray2 Pwind(env,NT);
IloArray<IloNumArray2> Pswind(env,set);//风电场景
IloNum detaa(10);//建议行为约束,此变量自由度比较大

///////////////////////////////与安全约束有关的量///////////////////////////////////////
IloNumArray2 B0(env,Node-1),B0l(env,Node-1);//导纳矩阵

IloNumArray Sl(env,Node);//每个节点的负荷状态,每个节点所占负荷百分比
IloIntArray Sw(env,Node);//常规风电场分布

IloNum f(IloNum &Z,IloNumArray &unit)
{
	return unit[10]*Cost_Modify+unit[9]*Z*Cost_Modify+unit[8]*Z*Z*Cost_Modify;
}
/***************************************************矩阵求逆****************************************************/
void IloInvert(const IloNumArray2 &B0,IloNumArray2 &A,IloInt N)//求逆矩阵
{
	IloNumArray2 C(env,N);//扩展矩阵

	for(IloInt i=0; i<N; ++i) //构造单位阵
		for(IloInt j=0; j<N; ++j)
		{
			if(i==j)
				A[i][j]=1;
			else
				A[i][j]=0;
		}

	for(IloInt i=0; i<N; ++i) //求增广矩阵
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

void IloMutiply(const IloNumArray2 &B0l,IloNumExprArray &Psp,IloExprArray &Theta,const IloInt N)//求乘积矩阵
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
/***********************************************数据初始化*****************************************************/
void define_data(IloEnv env)//数据初始化
{
	Master_Cplex.setParam(IloCplex::EpGap,0.11);//再小一点就会out of memory
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>变量初始化<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
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
	
	output<<endl<<"Pwind:"<<endl;
	ifstream windpower("wind_power.txt",ios::in);//风电基础场景
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
	ifstream wind_scenarios("wind_scenarios.txt",ios::in);//读入风电场景
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
			Betaa[i][l]=( f(Z[i][l+1],Unit[i])-f(Z[i][l],Unit[i]) )/( Z[i][l+1]-Z[i][l] );
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
		for(IloInt i=0;i<NG;++i)//运行成本及其线性化（青华姐论文）
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
		/*for(IloInt t=0; t<NT; ++t) //条件不全，暂时不预考虑
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
/************************************************************开始整个迭代过程***************************************************************/
		IloInvert(B0,B0l,Node-1);//求逆矩阵

		IloInt while_interator=0;
		IloBool first_cut_is_add;
		
		while(1)
		{
			first_cut_is_add=IloFalse;
			output<<"主问题迭代次数："<<while_interator++<<endl;
			Master_Cplex.extract(Master_Model);			
			Master_Cplex.solve();
			
			if (Master_Cplex.getStatus() == IloAlgorithm::Infeasible)//输出结果
				output<<"Master Problem Have No Solution"<<endl;
			else 
				output<<"Master Problem Have Solution, the objective value is: "<<Master_Cplex.getObjValue()<<endl<<endl;
				
			for(IloInt i=0; i<NG; ++i)//保存机组启停状态
			{
				Master_Cplex.getValues(I[i],u[i]);
				Master_Cplex.getValues(P[i],P1[i]);
			}
/*****************************************************主问题安全约束********************************************************/
			for(IloInt t=0;t<NT;++t)//主问题安全约束的校验
			{
				output<<"Master Problem Security Constraint Check:"<<t+1<<endl;
				IloNumVarArray C1(env,Branch,0,IloInfinity,ILOFLOAT),C2(env,Branch,0,IloInfinity,ILOFLOAT);
				IloNumExprArray PL(env,Branch);//各条支路的潮流
				IloNumVarArray P_dual(env,Node-1);
				IloRangeArray range(env,Node-1);
				for(IloInt h=0;h<Branch;++h)
				{
					C1[h]=IloNumVar(env);
					C2[h]=IloNumVar(env);
					PL[h]=IloNumVar(env);
				}				
				IloModel Master_Model_Security(env);
				
				IloNumExpr obj(env);//目标函数
				for(IloInt h=0;h<Branch;++h)
				{
					obj+=(C1[h]+C2[h]);
				}
				Master_Model_Security.add(IloMinimize(env,obj));
				obj.end();
				
				IloNumExprArray Psp(env,Node-1);//注入功率
				IloExprArray Theta(env,Node);//相角
				
				Theta[Node-1]=IloExpr(env);
					
				for(IloInt k=0; k<Node-1; ++k)//校验了所有结点的场景
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
				//计算潮流
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
					
				if (Master_Security_Cplex.getStatus() == IloAlgorithm::Infeasible)//输出结果
					output << "Master_Security_Problem Have No Solution"<<endl;
				else
					output<<"Master_Security_Problem Have Solution, the objective value is: "<<Master_Security_Cplex.getObjValue()<<endl<<endl;
				
				if(Master_Security_Cplex.getObjValue()>limit)
				{
					IloNumExpr cut(env);
					cut+=Master_Security_Cplex.getObjValue();
					for(IloInt i=0;i<NG;++i)
					{
						if(Unit[i][0]-1>=Node-1)//最后一个结点没算
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
				output<<"主问题安全约束越限，返回主问题>>>"<<endl;
				continue;//如果有添加割，则重新开始
			}
/****************************************************************场景的校验****************************************************************/			
			IloInt s=0;
			IloBool second_cut_is_add;
			for(; s<set; ++s)//场景的校验
			{	
				output<<"Num_Senario_Iterator（场景）:"<<s+1<<endl;
				
				IloNumVarArray S1S(env,NT,0,IloInfinity,ILOFLOAT),S2S(env,NT,0,IloInfinity,ILOFLOAT);//第二阶段松弛变量
								//S3S(env,NT,0,IloInfinity,ILOFLOAT),S4S(env,NT,0,IloInfinity,ILOFLOAT);
				IloArray<IloNumVarArray> S5S(env,NG);
				for(IloInt i=0; i<NG; ++i)
					S5S[i]=IloNumVarArray(env,NT,0.0,IloInfinity,ILOFLOAT);
				
				second_cut_is_add=IloFalse;

				IloModel Sub_Model(env,"Sub_Model");
/********************************子问题目标函数**************************************/
				IloNumExpr VS(env);
				for(IloInt t=0; t<NT; ++t)
				{
					VS+=(S1S[t]+S2S[t]);					
					for(IloInt i=0; i<NG; ++i)
						VS+=S5S[i][t];
				}
				Sub_Model.add(IloMinimize(env,VS));
				VS.end();
/*****************************子问题系统平衡约束*********************************/
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
/*******************************子问题备用约束**************************************/
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
/***************************子问题建议行为约束*************************************/
				IloArray<IloRangeArray> ramp1(env,NG),ramp2(env,NG); 
				for(IloInt i=0; i<NG; ++i)//绝对值转化为两个约束
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
/**************************子问题机组出力上下限约束*******************************/
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
/*************************************添加割************************************/
				IloCplex Sub_Cplex(Sub_Model);
				label1:
				Sub_Cplex.extract(Sub_Model);
				Sub_Cplex.solve();
				
				IloNum objvalue(Sub_Cplex.getObjValue());
				if(Sub_Cplex.getStatus() == IloAlgorithm::Infeasible)//输出结果
				{
					output<< "Senario_Sub_Problem No Solution" << endl;
					output << "场景 "<< s+1 <<" 无解,跳出循环，主问题无解……"<<endl;
					goto label2;
				}
				else 
					output<<"Senario_Sub_Problem Have Solution, and objective is: "<<objvalue<<endl<<endl;
				
				if(objvalue>limit)//返回主问题
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
					output<<"场景校验子问题不满足，返回主问题>>>"<<endl;
					break;//跳出场景校验子问题
				}
				else//否则开始安全校验
				{
					IloNumArray2 PS1(env,NG);//提取相应场景下的出力
					for(IloInt i=0;i<NG;++i)
					{
						PS1[i]=IloNumArray(env,NT);
						Sub_Cplex.getValues(PS[i],PS1[i]);
					}
					
					for(IloInt t=0;t<NT;++t)//对所有时段进行安全校核
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
						IloModel SubModel_Security(env);//只是局部变量
						
						IloNumExpr obj(env);//目标函数
						for(IloInt l=0;l<Branch;++l)
							obj+=C1[l]+C2[l];
						SubModel_Security.add(IloMinimize(env,obj));
						obj.end();
						
						IloNumExprArray Psp(env,Node-1);//注入功率
						IloExprArray Theta(env,Node);//相角
						
						Theta[Node-1]=IloExpr(env);
							
						for(IloInt k=0; k<Node-1; ++k)//求Psp
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
						//计算潮流
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
						
						if (Senario_Security_Cplex.getStatus() == IloAlgorithm::Infeasible)//输出结果
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
									cut+=Senario_Security_Cplex.getDual(range[(IloInt)Unit[i][0]-1])*(PS[i][t]-PS1[i][t]);//安全约束割
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
					
					if(second_cut_is_add==IloTrue)//返回当前场景校验
					{
						output<<"goto label1，场景校验子问题安全约束越限，返回当前子问题>>>"<<endl;
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
/************************************************************输出显示过程**************************************************/
		output<<"Cost:"<<Master_Cplex.getObjValue()<<endl;
		for(IloInt i=0; i<NG; ++i)
		{
			for(IloInt t=0; t<NT; ++t)
			{
				output<<"机组 "<<i+1<<" 第 "<<t+1<<" 时段状态:"<<Master_Cplex.getValue(I[i][t])<<"   出力:";
				output<<Master_Cplex.getValue(P[i][t])<<endl;
			}
			output<<endl;
		}
		Master_Model.end();
		Master_Cplex.end();
		env.end();
		label2: ;
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
	output<<"totaltime: "<<totaltime<<"s"<<endl<<endl;
	output.close();	
	return 0;
}
