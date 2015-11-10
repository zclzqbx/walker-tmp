/*              ��������ʷֲ����������
 *1���˳���������������������µķ��ɽ��ܲ�����������
 *
 */
#include<ilcplex/ilocplex.h>
#include<fstream>
#include<time.h>
#include"function.h"

ILOSTLBEGIN

const IloInt NG=54;//�������̨��
const IloInt NW=1;//��糡����
const IloInt Node=118;//���
const IloInt Branch=186;//֧·


IloEnv env;//���л��������н���֮��Ҫ�ر�
IloModel Master_Model(env,"Master_Model");//����ģ��
IloCplex Master_Cplex(Master_Model);//����CPLEX����

ofstream output("output_interval.txt",ios::ate);
IloNumVarArray detaP(env,NG,-IloInfinity,IloInfinity,ILOFLOAT),
				detaPw(env,NW,-IloInfinity,IloInfinity,ILOFLOAT);//���ͻ���Ŀɵ�����

IloNumArray P1(env,NG),u(env,NG);//���ڴ洢P��I��ֵ
IloNumArray PL(env,Branch);
IloNumArray2 Unit(env,NG),Info_Branch(env,Branch);//֧·��Ϣ

IloNumArray detaa(env,NG);//������ΪԼ��,��ȡֵӦ��������������й�

IloNumArray2 B0(env,Node-1),B0l(env,Node-1);//���ɾ���
IloIntArray Sw(env,Node);//�����糡�ֲ�

/***********************************************���ݳ�ʼ��*****************************************************/
void define_data(IloEnv env)//���ݳ�ʼ��,��ȫ�ֱ������и�ֵ
{
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>������ʼ��<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	for(IloInt i=0;i<Node-1;++i)
	{
		B0l[i]=IloNumArray(env,Node-1);
	}
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>������ʼ������<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	
	/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>���ݶ�ȡ<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
	output<<endl<<"Unit:"<<endl;
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
	
	for(IloInt i=0;i<NG;++i)
	{
		detaa[i]=Unit[i][5]/20;
	}
	
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
		if(((k+1)%20)==0)
			output<<endl;
	}
	output<<endl<<endl;
	wind_locate.close();
	
	output<<endl<<"State and Output of Unit:"<<endl;
	ifstream Unit_State("Unit_State.txt",ios::in);//��ȡ����״̬
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
	time_t nowTime=time(0);
	struct tm* nowTimeStruct=localtime(&nowTime);
	output<<"ϵͳ��ǰʱ�䣺"<<1900+nowTimeStruct->tm_year<<"."<<nowTimeStruct->tm_mon<<"."<<
		nowTimeStruct->tm_mday<<"  "<<nowTimeStruct->tm_hour<<":"<<nowTimeStruct->tm_min<<":"<<nowTimeStruct->tm_sec<<endl;
	
	try
	{
		define_data(env);//���ȳ�ʼ��ȫ�ֱ���		
		IloInvert(env,B0,B0l,Node-1);//�������
		
		IloNumExpr Cost(env);//Ŀ�꺯��
		for(IloInt w=0;w<NW;++w)
		{
			Cost+=detaPw[w];
		}		
		// Master_Model.add(IloMinimize(env,Cost));
		Master_Model.add(IloMaximize(env,Cost));//Ŀ�꺯����ѡһ
		Cost.end();

		IloNumExpr expr1(env),expr2(env);//����ƽ��Լ��
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
		
		for(IloInt i=0;i<NG;++i)//����ɵ��ڷ�Χ
		{
			Master_Model.add(detaP[i]>=Unit[i][1]*u[i]-P1[i]);
			Master_Model.add(detaP[i]<=Unit[i][2]*u[i]-P1[i]);
			
			Master_Model.add(detaP[i]>=-detaa[i]);
			Master_Model.add(detaP[i]<=detaa[i]);
		}
		
		IloNumExprArray detaP_node(env,Node-1),detaPw_node(env,Node-1);//��ȫԼ����ʵ���ϰ�ȫԼ��Ӱ�첻��
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
			IloNumExpr exprTheta(env);//Ī������Ĵ���
			exprTheta+=(Theta[(IloInt)Info_Branch[h][0]-1]-Theta[(IloInt)Info_Branch[h][1]-1]);
			
			Master_Model.add(exprTheta<=Info_Branch[h][3]*(Info_Branch[h][4]-PL[h]));			
			Master_Model.add(exprTheta>=Info_Branch[h][3]*(-Info_Branch[h][4]-PL[h]));
			exprTheta.end();
			//��������Ľڵ�˳��û��Ӱ��ô��
		}
		Theta.end();
		detaP_node.end();
		detaPw_node.end();
		
		Master_Cplex.extract(Master_Model);
		Master_Cplex.solve();
		if (Master_Cplex.getStatus() == IloAlgorithm::Infeasible)//������
		{
			output<<"Master Problem Have No Solution"<<endl;
			goto lable2;
		}
		
/************************************************************�����ʾ����**************************************************/
		
		output<<endl<<"Min/Max:"<<Master_Cplex.getObjValue()<<endl;
		
		output<<endl<<"����������������:"<<endl;
		for(IloInt i=0;i<NG;++i)
		{
			output<<Master_Cplex.getValue(detaP[i])<<"  ";
		}
		output<<endl<<"���������������:"<<endl;
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
