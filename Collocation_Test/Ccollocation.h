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
	//���캯��
	CCollocation();
	//��������
	virtual ~CCollocation();

	int InputData(char* DataFile);


	//����۲�ֵ����
	void Get_L();

	double* GGDD(double A[],int n,double medA);

	//����۲�ֵȨ����

	void Get_PL();

	//����ϵ����A
	void Get_A();

	//�����L��
	void Get_detaL();









	

	//����k
	double Get_k(int Pnumber, double* detaLx, double* XYZ);

	//�Ѳ��ź�ƽ��ǰ��Э����
	void Get_QS1();

	//��Ҫ�Ǽ����˹��Ϻ����е�C0�������źŵ�Э������ʱҲҪ�õ�
	void Get_C0();
	//��������֮��ľ��벢�浽������
	void Get_diStance();




	//����ƽ����detax������ŷ������
	void Get_detax();

	//����ƽ�����Ѳ��źŵ�  �ֲ�������α�,���һ������δ���źŵĹ�ֵ
	void Get_S1_Adjust();

	//����ƽ����α������������+  �ֲ�������α�
	void Get_Lv();

	//����δ���źŹ�ֵ
	void Get_S2_Adjust();

	//����δ���źŵ���ȫ�ź�
	void Get_Ls2();


	void PrintResult(char * resultfile);


	void Get_unused();

private:
	//ϵ������
	MatrixXd A;//����A
	MatrixXd AT;//����A��ת�þ���

	//��L -----���Ҳ�ǹ۲�ֵ��ֻ������ȥ������������Ĺ۲�ֵ(ԭ�۲�ֵ-��λ��)
	
	//X\Y\Z����
	double* detaLx;
	double* detaLy;
	double* detaLz;

	double C0x;
	double C0y;
	double C0z;



private:
	int n;//�۲�ֵ����
	

	char** Pname1;//��һ�ܵĵ�����ַ����
	char** Pname52;//��52�ܵĵ�����ַ����

	int Pnumber1;//��һ�ܵĵ���
	int Pnumber52;//��52�ܵĵ���

	char** Pname;//������С�������õĵ�ĵ�����ַ����
	int Pnumber;//������С�������õĵ���
	double* XYZ;//������С�������õĵ������

	int unusedPnumber;//δ������С�������õĵ�ĸ���
	char** unusedPname;
	double* unusedXYZ;



	double* XYZ1;//��һ�ܵĵ��XYZ��������
	double* XYZ52;//��52�ܵĵ��XYZ��������

	double* Qxyz1;//��һ�ܵĵ��XYZ����Э������
	double* Qxyz52; //��52�ܵĵ��XYZ����Э������




	
	MatrixXd PL;//�۲�ֵȨ����

	VectorXd L;//L�������۲�ֵ

	VectorXd diStance;//��֮��ľ�������

	Vector3d detax;//wx,wy,wz(3��ŷ������)

	MatrixXd Qdetax;

	MatrixXd QS1;//�źŵ�Э�������Ѳ��źţ�

	MatrixXd QS12;//

	MatrixXd QS1_Adjust;//��������

	VectorXd S1_Adjust;//ƽ�����Ѳ��ź� �ֲ�������α���

	VectorXd S2_Adjust;//δ���źŵĹ�ֵ���ƹ�

	VectorXd Lv;//ƽ�����α���   �������+�ֲ������α�
	VectorXd Ls2;//δ������ȫ�ź�




	

};
