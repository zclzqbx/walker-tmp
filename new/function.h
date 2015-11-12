/*              机组组合相关头文件
 *
 *
 *
 */
#include <ilcplex/ilocplex.h>

static IloNum f(IloEnv env,IloNum &Z,IloNumArray &unit)//成本函数
{
	return unit[10]+unit[9]*Z+unit[8]*Z*Z;
}

static void IloInvert(IloEnv env,const IloNumArray2 &B0,IloNumArray2 &A,IloInt N)//求逆矩阵
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

static void IloMutiply(IloEnv env,const IloNumArray2 &B0l,IloNumExprArray &Psp,
		IloExprArray &Theta,const IloInt N)//求乘积矩阵
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
