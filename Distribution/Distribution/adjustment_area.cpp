/*              不规则概率分布风电机组组合
 *1）此程序用于求解给定机组组合下的风电可接受波动区间问题
 *
 */
#include<ilcplex/ilocplex.h>
#include<fstream>
#include<time.h>
#include"function.h"

ILOSTLBEGIN

const IloInt NG=54;//常规机组台数
const IloInt NW=1;//风电场个数
const IloInt Node=118;//结点
const IloInt Branch=186;//支路


IloEnv env;//运行环境，运行结束之后要关闭
IloModel Master_Model(env,"Master_Model");//生成模型
IloCplex Master_Cplex(Master_Model);//创建CPLEX环境

ofstream output("output_interval.txt",ios::ate);
IloNumVarArray detaP(env,NG,-IloInfinity,IloInfinity,ILOFLOAT),
				detaPw(env,NW,-IloInfinity,IloInfinity,ILOFLOAT);//风电和机组的可调整量

IloNumArray P1(env,NG),u(env,NG);//用于存储P、I的值
IloNumArray PL(env,Branch);
IloNumArray2 Unit(env,NG),Info_Branch(env,Branch);//支路信息

IloNumArray detaa(env,NG);//建议行为约束,其取值应该与机组爬坡率有关

IloNumArray2 B0(env,Node-1),B0l(env,Node-1);//导纳矩阵
IloIntArray Sw(env,Node);//常规风电场分布

/***********************************************数据初始化*****************************************************/
void define_data(IloEnv env)//数据初始化,对全局变量进行赋值
{
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>变量初始化<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	for(IloInt i=0;i<Node-1;++i)
	{
		B0l[i]=IloNumArray(env,Node-1);
	}
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>变量初始化结束<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>数据读取<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	output<<endl<<"Unit:"<<endl;
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
	
	for(IloInt i=0;i<NG;++i)
	{
		detaa[i]=Unit[i][5]/20;
	}
	
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
		if(((k+1)%20)==0)
			output<<endl;
	}
	output<<endl<<endl;
	wind_locate.close();
	
	output<<endl<<"State and Output of Unit:"<<endl;
	ifstream Unit_State("Unit_State.txt",ios::in);//读取机组状态
	output<<"State:"<<endl;
	if(!Unit_State)
	{	
		output<<"no such file! Unit_State.txt"<<endl;
	}
	for(IloInt i=0; i<NG; ++i)//State
	{
		Unit_State>>u[i];
		output<<u[i]<<"   ";
		if(((i+1)%20)==0)
			output<<endl;
	}
	output<<endl;
	
	output<<"Output:"<<endl;
	for(IloInt i=0; i<NG; ++i)//Output
	{
		Unit_State>>P1[i];
		output<<P1[i]<<"   ";
		if(((i+1)%20)==0)
			output<<endl;
	}
	
	output<<endl;
	output<<"PL:"<<endl;
	for(IloInt l=0;l<Branch;++l)
	{
		Unit_State>>PL[l];
		output<<PL[l]<<"   ";
		if(((l+1)%20)==0)
			output<<endl;
	}
	Unit_State.close();

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
	time_t nowTime=time(0);
	struct tm* nowTimeStruct=localtime(&nowTime);
	output<<"系统当前时间："<<1900+nowTimeStruct->tm_year<<"."<<nowTimeStruct->tm_mon<<"."<<
		nowTimeStruct->tm_mday<<"  "<<nowTimeStruct->tm_hour<<":"<<nowTimeStruct->tm_min<<":"<<nowTimeStruct->tm_sec<<endl;
	
	try
	{
		define_data(env);//首先初始化全局变量		
		IloInvert(env,B0,B0l,Node-1);//求逆矩阵
		
		IloNumExpr Cost(env);//目标函数
		for(IloInt w=0;w<NW;++w)
		{
			Cost+=detaPw[w];
		}		
		// Master_Model.add(IloMinimize(env,Cost));
		Master_Model.add(IloMaximize(env,Cost));//目标函数二选一
		Cost.end();

		IloNumExpr expr1(env),expr2(env);//功率平衡约束
		for(IloInt i=0;i<NG;++i)
		{
			expr1+=detaP[i];
		}
		for(IloInt w=0;w<NW;++w)
		{
			expr2+=detaPw[w];
		}
		Master_Model.add(expr1+expr2==0);
		expr1.end();
		expr2.end();
		
		for(IloInt i=0;i<NG;++i)//机组可调节范围
		{
			Master_Model.add(detaP[i]>=Unit[i][1]*u[i]-P1[i]);
			Master_Model.add(detaP[i]<=Unit[i][2]*u[i]-P1[i]);
			
			Master_Model.add(detaP[i]>=-detaa[i]);
			Master_Model.add(detaP[i]<=detaa[i]);
		}
		
		IloNumExprArray detaP_node(env,Node-1),detaPw_node(env,Node-1);//安全约束，实际上安全约束影响不大
		IloNumExprArray Theta(env,Node);
		
		for(IloInt b=0;b<Node-1;++b)
		{
			detaP_node[b]=IloNumExpr(env);
			detaPw_node[b]=IloNumExpr(env);
			IloInt i=0;
			for(;i<NG;++i)
			{
				if(Unit[i][0]==b-1)break;
			}
			
			if(i<NG)
			{
				detaP_node[b]+=detaP[i];
			}	
			
			if(Sw[b]>=0)
			{
				detaPw_node[b]+=detaPw[ Sw[b] ];
			}
		}
		
		for(IloInt b=0;b<Node-1;++b)
		{
			Theta[b]=IloNumExpr(env);
			for(IloInt k=0;k<Node-1;++k)
			{			
				Theta[b]+=B0l[b][k]*(detaP_node[k]+detaPw_node[k]);	
			}
		}
		Theta[Node-1]=IloNumExpr(env);
		
		for(IloInt h=0;h<Branch;++h)
		{
			IloNumExpr exprTheta(env);//莫明其妙的错误
			exprTheta+=(Theta[(IloInt)Info_Branch[h][0]-1]-Theta[(IloInt)Info_Branch[h][1]-1]);
			
			Master_Model.add(exprTheta<=Info_Branch[h][3]*(Info_Branch[h][4]-PL[h]));			
			Master_Model.add(exprTheta>=Info_Branch[h][3]*(-Info_Branch[h][4]-PL[h]));
			exprTheta.end();
			//两个相减的节点顺序没有影响么？
		}
		Theta.end();
		detaP_node.end();
		detaPw_node.end();
		
		Master_Cplex.extract(Master_Model);
		Master_Cplex.solve();
		if (Master_Cplex.getStatus() == IloAlgorithm::Infeasible)//输出结果
		{
			output<<"Master Problem Have No Solution"<<endl;
			goto lable2;
		}
		
/************************************************************输出显示过程**************************************************/
		
		output<<endl<<"Min/Max:"<<Master_Cplex.getObjValue()<<endl;
		
		output<<endl<<"常规机组出力调整量:"<<endl;
		for(IloInt i=0;i<NG;++i)
		{
			output<<Master_Cplex.getValue(detaP[i])<<"  ";
		}
		output<<endl<<"风电机组出力调整量:"<<endl;
		for(IloInt i=0;i<NW;++i)
		{
			output<<Master_Cplex.getValue(detaPw[i])<<"  ";
		}
		output<<endl;
		
		lable2:
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
	
	finish=clock();
	totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	output<<"totaltime: "<<totaltime<<"s"<<endl<<endl;
	output.close();	
	return 0;
}
