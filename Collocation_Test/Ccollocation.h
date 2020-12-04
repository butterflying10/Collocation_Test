#pragma once
#include <conio.h>
#include <stdlib.h>


#include "stdio.h"
#include "math.h"
#include "string.h"

#include <Eigen/Dense>
#include<Eigen/Core>



using namespace Eigen;
class CCollocation
{
public:
	//构造函数
	CCollocation();
	//析构函数
	virtual ~CCollocation();

	int InputData(char* DataFile);


	//构造观测值矩阵
	void Get_L();

	double* GGDD(double A[],int n,double medA);

	//构造观测值权矩阵

	void Get_PL();

	//构造系数阵A
	void Get_A();

	//构造δL阵
	void Get_detaL();









	

	//计算k
	double Get_k(int Pnumber, double* detaLx, double* XYZ);

	//已测信号平差前的协方差
	void Get_QS1();

	//主要是计算高斯拟合函数中的C0，构造信号的协方差阵时也要用到
	void Get_C0();
	//计算两点之间的距离并存到数组中
	void Get_diStance();




	//计算平差后的detax，三个欧拉向量
	void Get_detax();

	//计算平差后的已测信号的  局部区域的形变,并且还会计算未测信号的估值
	void Get_S1_Adjust();

	//计算平差后形变量，倾向参数+  局部区域的形变
	void Get_Lv();

	//计算未测信号估值
	void Get_S2_Adjust();

	//计算未测信号的完全信号
	void Get_Ls2();


	void PrintResult(char * resultfile);


	void Get_unused();

private:
	//系数矩阵
	MatrixXd A;//矩阵A
	MatrixXd AT;//矩阵A的转置矩阵

	//δL -----这个也是观测值，只不过是去掉了倾向参数的观测值(原观测值-中位数)
	
	//X\Y\Z分量
	double* detaLx;
	double* detaLy;
	double* detaLz;

	double C0x;
	double C0y;
	double C0z;



private:
	int n;//观测值总数
	

	char** Pname1;//第一周的点名地址数组
	char** Pname52;//第52周的点名地址数组

	int Pnumber1;//第一周的点数
	int Pnumber52;//第52周的点数

	char** Pname;//参与最小二乘配置的点的点名地址数组
	int Pnumber;//参与最小二乘配置的点数
	double* XYZ;//参与最小二乘配置的点的坐标

	int unusedPnumber;//未参与最小二乘配置的点的个数
	char** unusedPname;
	double* unusedXYZ;



	double* XYZ1;//第一周的点的XYZ坐标数据
	double* XYZ52;//第52周的点的XYZ坐标数据

	double* Qxyz1;//第一周的点的XYZ坐标协因数阵
	double* Qxyz52; //第52周的点的XYZ坐标协因数阵




	
	MatrixXd PL;//观测值权矩阵

	VectorXd L;//L的向量观测值

	VectorXd diStance;//点之间的距离向量

	Vector3d detax;//wx,wy,wz(3个欧拉向量)

	MatrixXd Qdetax;

	MatrixXd QS1;//信号的协方差阵（已测信号）

	MatrixXd QS12;//

	MatrixXd QS1_Adjust;//精度评定

	VectorXd S1_Adjust;//平差后的已测信号 局部区域的形变量

	VectorXd S2_Adjust;//未测信号的估值，推估

	VectorXd Lv;//平差后的形变量   倾向参数+局部区域形变
	VectorXd Ls2;//未测点的完全信号




	

};
