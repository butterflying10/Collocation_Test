#pragma once
#include<fstream>
#include <iostream>

#include <string>

#include < iomanip >//保留小数点后几位的

#include<cmath>

#include <Eigen/Dense>
#include<Eigen/Core>
#include <Eigen/Sparse>//稀疏矩阵
#include "CCollocation.h"
#include "public.h" //里面有计算中位数的方法


using namespace Eigen;
using namespace std;

//构造函数
CCollocation::CCollocation()
{

}

//析构函数
CCollocation::~CCollocation()
{

}

/*输入原始数据(第1周和第52周的各点坐标和其协因数)
* 存储到D1、D52中
* 
*/
int CCollocation::InputData(char* DataFile)
{
    ifstream infile(DataFile, ios::in | ios::_Nocreate);


    if (!infile)
    {
        cerr << "打开原始数据失败" << endl;
        return 0;
    }

    string note;//备注
    getline(infile, note);
    

    infile >> note;//读取“点数”
    

    
    infile >> Pnumber1;

    //cout << Pnumber1 << endl;

    Pname1 = new char* [Pnumber1];//第一周的点名地址数据
    for (int i = 0; i < Pnumber1; i++)Pname1[i] = NULL;
    XYZ1 = new double[Pnumber1 * 3];//第一周 点的XYZ坐标
    Qxyz1 = new double[Pnumber1 * 3];//第一周 点的XYZ坐标的方差


    getline(infile, note);
    getline(infile, note);
    getline(infile, note);
    getline(infile, note);
    getline(infile, note);

    /*
    读第一周的数据点 XYZ坐标  方差
    */

    for (int i = 0; i < Pnumber1; i++)
    {
        int note;
        infile >> note;

        //cout << note << endl;

        char pname[20];
        infile >> pname;
        int len = strlen(pname);
        Pname1[i] = new char[len + 1];//加入字符串结束的标志‘/0’
        strcpy_s(Pname1[i], strlen(pname) + 1, pname);
    
      

        //存储第一周各点的坐标
        infile >> XYZ1[3 * i];
        infile >> XYZ1[3 * i + 1];
        infile >> XYZ1[3 * i + 2];

        //cout << XYZ1[3 * i] << setw(20) << XYZ1[3 * i + 1] << setw(20) << XYZ1[3 * i + 2] << endl;

        //存储第一周各点的方差
        infile >> Qxyz1[3 * i];
        infile >> Qxyz1[3 * i + 1];
        infile >> Qxyz1[3 * i + 2];

        //cout << Qxyz1[3 * i] << setw(20) << Qxyz1[3 * i + 1] << setw(20) << Qxyz1[3 * i + 2] << endl;


    }


    getline(infile, note);
    getline(infile, note);
    getline(infile, note);

    infile >> note;//读取“点数”

    //int Pnumber52;
    infile >> Pnumber52;//第五十二周点数

    //cout << Pnumber52 << endl;

    Pname52 = new char* [Pnumber52];//第52周的点名地址数据
    for (int i = 0; i < Pnumber52; i++)Pname52[i] = NULL;
    XYZ52 = new double[Pnumber52 * 3];//第52周 点的XYZ坐标
    Qxyz52 = new double[Pnumber52 * 3];//第52周 点的XYZ坐标的方差

    getline(infile, note);
    getline(infile, note);
    getline(infile, note);
    getline(infile, note);
    getline(infile, note);


    /*
    读第52周的数据点 XYZ坐标  方差
    */

    for (int i = 0; i < Pnumber52; i++)
    {
        int note;
        infile >> note;

        //cout << note << endl;

        char pname[20];
        infile >> pname;
        int len = strlen(pname);
        Pname52[i] = new char[len + 1];
        strcpy_s(Pname52[i], strlen(pname) + 1, pname);
        


        //存储第一周各点的坐标
        infile >> XYZ52[3 * i];
        infile >> XYZ52[3 * i + 1];
        infile >> XYZ52[3 * i + 2];

        //cout << XYZ52[3 * i] << setw(20) << XYZ52[3 * i + 1] << setw(20) << XYZ52[3 * i + 2] << endl;

        //存储第一周各点的方差
        infile >> Qxyz52[3 * i];
        infile >> Qxyz52[3 * i + 1];
        infile >> Qxyz52[3 * i + 2];

        //cout << Qxyz52[3 * i] << setw(20) << Qxyz52[3 * i + 1] << setw(20) << Qxyz52[3 * i + 2] << endl;

    }


    infile.close();

}

void CCollocation::Get_unused()
{
    //获取未测点信号的坐标和未测点信号的个数
    unusedPnumber = Pnumber1 + Pnumber52 - Pnumber * 2;

    unusedPname = new char* [unusedPnumber];
    for (int i = 0; i < unusedPnumber; i++)unusedPname[i] = NULL;
    unusedXYZ = new double[unusedPnumber * 3];

    int k = 0;
    for (int i = 0; i < Pnumber1; i++)
    {
        int flag = 0;
        for (int j = 0; j < Pnumber; j++)
        {
            if (strcmp(Pname1[i], Pname[j]) == 0)

            {
                flag = 1;
                break;

            }
            if (flag == 0 && j == Pnumber - 1)
            {
                unusedPname[k] = Pname1[i];
                unusedXYZ[3 * k] = XYZ1[3 * i];
                unusedXYZ[3 * k + 1] = XYZ1[3 * i + 1];
                unusedXYZ[3 * k + 2] = XYZ1[3 * i + 2];

                cout << unusedPname[k] << setw(20) << unusedXYZ[3 * k] << setw(20) << unusedXYZ[3 * k + 1] << setw(20) << unusedXYZ[3 * k + 2] << endl;

                k++;

            }

        }

    }
    for (int i = 0; i < Pnumber52; i++)
    {

        for (int j = 0; j < Pnumber; j++)
        {
            if (strcmp(Pname52[i], Pname[j]) == 0)
            {

                break;

            }

            if (strcmp(Pname52[i], Pname[j]) != 0 && j == Pnumber - 1)
            {
                unusedPname[k] = Pname52[i];
                unusedXYZ[3 * k] = XYZ52[3 * i];
                unusedXYZ[3 * k + 1] = XYZ52[3 * i + 1];
                unusedXYZ[3 * k + 2] = XYZ52[3 * i + 2];

                cout << unusedPname[k] << setw(20) << unusedXYZ[3 * k] << setw(20) << unusedXYZ[3 * k + 1] << setw(20) << unusedXYZ[3 * k + 2] << endl;

                k++;
            }

        }

    }

}

/*
构造观测值矩阵和观测值协因数阵
*/

void CCollocation::Get_L()
{

    /*
    参与最小二乘配置的点数n
    */
    Pnumber = 0;
    unusedPnumber = 0;

    
    for (int i = 0; i < Pnumber1; i++)
    {
        for (int j = 0; j < Pnumber52; j++)
        {
            if (strcmp(Pname1[i], Pname52[j]) ==0)
            { 
                Pnumber++;
                break;
            }
        }  
    }
    



    Pname = new char* [Pnumber];
    XYZ = new double[Pnumber * 3];



    L.resize(Pnumber * 3);

    int k = 0;
    for (int i = 0; i < Pnumber1; i++)
    {
        for (int j = 0; j < Pnumber52; j++)
        {
            if (strcmp(Pname1[i], Pname52[j]) == 0)
            {
                Pname[k] = Pname1[i];
                //cout << Pname[k++] << endl; 
                L(3 * k) = XYZ52[3 * j] - XYZ1[3 * i];
                L(3 * k + 1) = XYZ52[3 * j + 1] - XYZ1[3 * i + 1];
                L(3 * k + 2) = XYZ52[3 * j + 2] - XYZ1[3 * i + 2];

                //cout << L[3 * k] << setw(20) << L[3 * k + 1] << setw(20) << L[3 * k + 2] << endl;

                /*v(3*k)= Qxyz1[3 * i]/10000.0 + Qxyz1[3 * j]/10000.0;
                v(3 * k + 1) = Qxyz1[3 * i + 1]/10000.0 + Qxyz1[3 * j + 1]/10000.0;
                v(3 * k + 2) = Qxyz1[3 * i + 2]/10000.0 + Qxyz1[3 * j + 2]/10000.0;*/
                XYZ[3 * k] = (XYZ52[3 * j] + XYZ1[3 * i]) / 2.0;
                XYZ[3 * k + 1] = (XYZ52[3 * j + 1] + XYZ1[3 * i + 1]) / 2.0;
                XYZ[3 * k + 2] = (XYZ52[3 * j + 2] + XYZ1[3 * i + 2]) / 2.0;

                //cout<< XYZ[3 * k] << setw(20) << XYZ[3 * k + 1] << setw(20) << XYZ[3 * k + 2] << endl;     

                k++;
            }

        }
    }


    //QL要换一种方法定义，而不是读取秩亏自由网平差的结果，这个结果不是一个简单的对角阵
    ////cout << v << endl;
    //QL = v.asDiagonal();
    //PL = QL.inverse();
        ////归一化   归一化后的矩阵的每一行或者列的元素的平方和为1
    //PL.normalize();
    //cout << PL << endl;



    

}
/*构造观测值权矩阵*/

void CCollocation::Get_PL()
{
    

    double* Lx = new double[Pnumber];//存储原观测值X分量的数据
    double* Ly = new double[Pnumber];//存储原观测值Y分量的数据
    double* Lz = new double[Pnumber];//存储原观测值Z分量的数据




    detaLx = new double[Pnumber];//存储detaL   X分量的数据
    detaLy = new double[Pnumber];//存储detaL   Y分量的数据
    detaLz = new double[Pnumber];//存储detaL   Z分量的数据



   /*
   观测值
   */
    for (int i = 0; i < Pnumber; i++)
    {
        Lx[i] = L[3 * i];
        Ly[i] = L[3 * i + 1];
        Lz[i] = L[3 * i + 2];

        
    }

    //求取Lx数组的中位数
    double medLx = Median(Lx, Pnumber, false);
    double medLy = Median(Ly, Pnumber, false);
    double medLz = Median(Lz, Pnumber, false);

    //这个在计算C0和k的时候用到了
    for (int i = 0; i < Pnumber; i++)
    {
        detaLx[i] = (Lx[i] - medLx);
        detaLy[i] = (Ly[i] - medLy);
        detaLz[i] = (Lz[i] - medLz);
    }



    double* Px = new double[Pnumber];
    Px = GGDD(detaLx, Pnumber, medLx);
    cout << "=================================" << endl;

    double* Py = new double[Pnumber];
    Py = GGDD(detaLy, Pnumber, medLy);
    cout << "=================================" << endl;
    double* Pz = new double[Pnumber];
    Pz = GGDD(detaLz, Pnumber, medLz);

    //for(int i=0;i<Pnumber;i++) cout <<Px[i]<<setw(20)<<Py[i]<<setw(20)<< Pz[i] << endl;


    //开始构造PL
    PL.resize(Pnumber * 3, Pnumber * 3);

    VectorXd pp(Pnumber * 3);

    for (int i = 0; i < Pnumber; i++)
    {
        pp(3*i) = Px[i];
        pp(3*i + 1) = Py[i];
        pp(3*i + 2) = Pz[i];
    
    }
    PL = pp.asDiagonal();

}

/*
A:::观测值某一分量的观测值
n:::长度
medA:::中位数
*/

double* CCollocation::GGDD(double A[], int n, double medA)
{
    //构造PL的初始矩阵  单位阵，但是为了方便操作，只存放对角线元素，成为一个向量

    double* P = new double[n];
    for (int j = 0; j < n; j++) P[j] = 1;
    double* detaA = new double[n];

    int KK = 2.0;
    for (int ii = 0; ii < 100; ii++)
    {
        for (int i = 0; i < n; i++)
        {
            detaA[i] = A[i] - medA;
        }
        //计算初始方差因子
        double Sigmax0 = Median(detaA, n, true) / 0.6745;

        //cout << "Sigmax0：：：" << Sigmax0 << endl;


        for (int i = 0; i < n; i++)
        {
            if (abs(detaA[i] / Sigmax0) <= KK)
            {
                P[i] = P[i];
            }
            else
            {
                P[i] = P[i] * KK / abs(detaA[i] / Sigmax0);
            }

        }
        double newmedA = 0.0;
        double sumP = 0.0;
       
        //计算加权平均值
        for (int i = 0; i < n; i++)
        {
            newmedA += P[i] * A[i];
            sumP += P[i];

        }
        newmedA = newmedA / sumP;

       // cout << "GGDD:::迭代次数:::" << ii << "  medA:::" << medA << "  newmedA:::" << newmedA << endl;

        if (abs(newmedA - medA) < 0.0000001)
        {
           
            return P;
            delete[]P;
            break;
            
        }
        else if(ii==99)
        {
            return P;
            delete[]P;
            break;
        }
        else
        {
            medA = newmedA;
        }


    }


}

/*
构造系数阵A
*/
void CCollocation::Get_A()
{
    //A = Eigen::Matrix<double, Dynamic, Dynamic>();
    A.resize(Pnumber * 3, 3);
    AT.resize(3, Pnumber * 3);

    for (int i = 0; i < Pnumber; i++)
    {
        /*
        Eigen中矩阵的赋值不能通过循环一个一个给值
        */
        /*A << 0.0, XYZ[3 * i + 2], -XYZ[3 * i + 1],
            -XYZ[3 * i + 2], 0.0, XYZ[3 * i],
            XYZ[3 * i + 1], -XYZ[3 * i], 0.0;*/

        //分块初始化赋值
        A.block(3 * i, 0, 3, 3)
            << 0.0, XYZ[3 * i + 2], -XYZ[3 * i + 1],
            -XYZ[3 * i + 2], 0.0, XYZ[3 * i],
            XYZ[3 * i + 1], -XYZ[3 * i], 0.0;
    }
    

    //转置
    AT = A.transpose();

    //cout << A << endl;
    //cout << endl<<endl;
}

void CCollocation::Get_detax()
{
    MatrixXd ATPLA_N(3, 3);

    ATPLA_N = (AT * PL * A).inverse();

    Qdetax.resize(3, 3);
    Qdetax = ATPLA_N;

    MatrixXd ATPLL(3, 1);

    ATPLL = AT * PL * L;

    detax = ATPLA_N * ATPLL;

    cout <<"detax:::"<< detax << endl;
}

void CCollocation::Get_C0()
{
    //计算信号的初始方差因子
    C0x = 0.0;
    for (int i = 0; i < Pnumber; i++) C0x = C0x + detaLx[i] * detaLx[i];
    C0x = C0x / Pnumber;

    C0y = 0.0;
    for (int i = 0; i < Pnumber; i++) C0y = C0y + detaLy[i] * detaLy[i];
    C0y = C0y / Pnumber;

    C0z = 0.0;
    for (int i = 0; i < Pnumber; i++) C0z = C0z + detaLz[i] * detaLz[i];
    C0z = C0z / Pnumber;
}

//void CCollocation::Get_diStance()
//{
//
//    diStance.resize(Pnumber * (Pnumber - 1) / 2);
//    int k = 0;
//    //第i个点
//    for (int i = 0; i < Pnumber; i++)
//    {
//        //第j个点
//        for (int j = i + 1; j < Pnumber; j++)
//        {
//            //算出来毫无意义啊啊啊啊啊，又不是真的距离,也是有意义的，
//
//            double Sx = (XYZ[3 * i] - XYZ[3 * j]) * (XYZ[3 * i] - XYZ[3 * j]);
//            double Sy = (XYZ[3 * i + 1] - XYZ[3 * j + 1]) * (XYZ[3 * i + 1] - XYZ[3 * j + 1]);
//            double Sz = (XYZ[3 * i + 2] - XYZ[3 * j + 2]) * (XYZ[3 * i + 2] - XYZ[3 * j + 2]);
//
//            // cout << Sx << setw(20) << Sy << setw(20) << Sz << endl;
//
//             //把S换成千米级，公里
//            diStance(k)= (Sx + Sy + Sz) / 1000000.0;
//         
//            k++;
//        }
//    }
//
//    //cout << diStance << endl << endl << endl;
//
//
//
//}


/*计算k
*/
double  CCollocation::Get_k(int Pnumber,double * detaLx,double *XYZ)
{

    //计算信号的初始方差因子
    double c=0.0;
    for (int i = 0; i < Pnumber; i++) c = c+ detaLx[i]*detaLx[i];
    c = c / Pnumber;

    

    //计算信号的初始协方差
     double* CSx = new  double[Pnumber * (Pnumber - 1) / 2];
     double* S = new  double[Pnumber * (Pnumber - 1) / 2];

    int k = 0;
    double* Kx = new double[Pnumber * (Pnumber - 1) / 2];

    //先把所有的值弄成负数
    for (int i = 0; i < Pnumber * (Pnumber - 1) / 2; i++) Kx[i] = -1.0;

    //看Kx这个数组中有多少个正数
    int positiveNumber = 0;

    //第i个点
    for (int i = 0; i < Pnumber; i++)
    {   
        //第j个点
        for (int j = i+1; j < Pnumber; j++)
        {
            CSx[k] = detaLx[i] * detaLx[j];

            double S = ((XYZ[3 * j] - XYZ[3 * i]) * (XYZ[3 * j] - XYZ[3 * i]) + (XYZ[3 * j + 1] - XYZ[3 * i + 1]) * (XYZ[3 * j + 1] - XYZ[3 * i + 1]) + (XYZ[3 * j + 2] - XYZ[3 * i + 2]) * (XYZ[3 * j + 2] - XYZ[3 * i + 2]))/1000000.0;
            //计算x分量下的k1,k2,k3,k4,k5.........

            if (CSx[k] > 0) Kx[k] = (log(c) - log(CSx[k])) / S;

            if (Kx[k] > 0) positiveNumber++;
            
            
            k++;
        }   
    }
    //cout <<"正数的个数为："<< Number << endl;
    //去掉kx中的负值并获取中位数
    //构造只含正数的Kx
    double *positiveKx=new double[positiveNumber];

    int j = 0;
    for (int i = 0; i < Pnumber * (Pnumber - 1) / 2; i++)
    {
        if (Kx[i] > 0)
        {
            positiveKx[j] = Kx[i];
            j++;
        }
    
    }


    
    double medKx=Median(positiveKx,positiveNumber,false);


    //接下来我们需要迭代得到medkx
     //构造降权矩阵
    //它是一个对角阵，所以可以用向量或者数组来表示
     VectorXd P=VectorXd::Constant(positiveNumber,1);
     
    // cout << P << endl;

     for (int ii = 0; ii < 50; ii++)
     {
         double* detaKx = new double[positiveNumber];
         for (int i = 0; i < positiveNumber; i++)
         {

             detaKx[i] = abs(positiveKx[i] - medKx);
         }


         double detaK0 = Median(detaKx, positiveNumber,false);

         

         

         //构造降权矩阵
         //它是一个对角阵，所以可以用向量或者数组来表示
        


         double KK = 2.0;
         for (int i = 0; i < positiveNumber; i++)
         {
             if (detaKx[i] / detaK0 <= KK) P(i) = P(i);
             else
             {
                 P(i) = P(i) * KK / abs(detaKx[i] / detaK0);
             }

         }
         //cout << P << endl;

         //计算属于k的加权平均值

         double newmedKx = 0.0;
         for (int i = 0; i < positiveNumber; i++)
         {
             newmedKx += P(i) * positiveKx[i];

         }
         newmedKx = newmedKx / P.sum();


        // cout <<"迭代次数:::"<<ii<< "   medKx:::" << medKx << "     " << "newmedKx" << newmedKx << endl;
         


         if (abs(medKx - newmedKx) < 0.000000001)
         {
             break;
         }
         else
         {
             medKx = newmedKx;
         }
     }

    return medKx;

}
/*构造已测信号的方差协方差阵*/
void CCollocation::Get_QS1()
{

    Get_C0();
    //Get_diStance();
    //计算k
    double kx=Get_k(Pnumber, detaLx, XYZ);


    cout << endl << endl << "--------------------------------" << endl;

    double ky= Get_k(Pnumber, detaLy, XYZ);

    cout << endl << endl << "--------------------------------" << endl;
    double kz= Get_k(Pnumber, detaLz, XYZ);

    cout <<"k::::::"<< kx << setw(20) << ky << setw(20) << kz << endl;

    //QS
    
    QS1 = MatrixXd::Zero(Pnumber*3, Pnumber*3);
    int k = 0;
    for (int i = 0; i < Pnumber; i++)
    {
        for (int j = i; j < Pnumber; j++)
        {
            if (j == i)
            {
                
                QS1.block(3 * i, 3 * j, 3, 3) << 
                    C0x / 2.0, 0, 0,
                    0, C0y / 2.0, 0,
                    0, 0, C0z / 2.0;
            }
            if (j > i)
            {
                double S = ((XYZ[3 * j] - XYZ[3 * i]) * (XYZ[3 * j] - XYZ[3 * i]) + (XYZ[3 * j + 1] - XYZ[3 * i + 1]) * (XYZ[3 * j + 1] - XYZ[3 * i + 1]) + (XYZ[3 * j + 2] - XYZ[3 * i + 2]) * (XYZ[3 * j + 2] - XYZ[3 * i + 2]))/1000000.0;

                QS1.block(3 * i, 3 * j, 3, 3) <<
                    C0x * exp(-kx * S), 0, 0,
                    0, C0y* exp(-ky * S), 0,
                    0, 0, C0z* exp(-kz * S);
                k++;
            
            }
        }
   }
    MatrixXd QST = QS1.transpose();
    QS1 = QS1 + QST;

    //验证对不对
    //cout << QS.bottomRightCorner(30, 30) << endl;


    QS12 = MatrixXd::Zero(unusedPnumber * 3, Pnumber * 3);

    for (int i = 0; i < unusedPnumber; i++)
    {
        for (int j = 0; j < Pnumber; j++)
        {
            double S = ((XYZ[3 * j] - unusedXYZ[3 * i]) * (XYZ[3 * j] - unusedXYZ[3 * i])+ (XYZ[3 * j+1] - unusedXYZ[3 * i+1]) * (XYZ[3 * j+1] - unusedXYZ[3 * i+1])+ (XYZ[3 * j+2] - unusedXYZ[3 * i+2]) * (XYZ[3 * j+2] - unusedXYZ[3 * i+2]))/1000000.0;
            //cout << "S::" << S << endl;

            QS12.block(3 * i, 3 * j, 3, 3) <<
                C0x * exp(-kx *S ), 0, 0,
                0, C0y* exp(-ky * S), 0,
                0, 0, C0z* exp(-kz * S);
        
        
        }
    
    }
    //cout << "QS12" << endl;
    //cout << QS12 << endl;


}

void CCollocation::Get_S1_Adjust()
{
    S1_Adjust.resize(Pnumber * 3);

    S1_Adjust = QS1 * PL * (L - A * detax);

    S2_Adjust.resize(unusedPnumber * 3);

    S2_Adjust = QS12 * PL * (L - A * detax);


    //QS = MatrixXd::Zero(Pnumber * 3, Pnumber * 3);
    QS1_Adjust.resize(Pnumber * 3, Pnumber * 3);

    MatrixXd G(Pnumber * 3, Pnumber * 3);
    G = PL - PL * A * Qdetax * AT * PL;
    QS1_Adjust = QS1 - QS1 * G * QS1;

    //cout << QS1_Adjust << endl<<endl;

}

void CCollocation::Get_Lv()
{
    Lv.resize(Pnumber * 3);
    Lv = A * detax + S1_Adjust;
    cout << "点名                   已测信号S                     L平差前             L平差后" << endl;
    for (int i = 0; i < Pnumber; i++)
    {
        cout << endl <<Pname[i]<<setw(20)<< S1_Adjust(3*i)<<setw(20)<< L(3*i) <<setw(20)<<Lv(3*i)<< endl;

        cout << endl <<  setw(20) << S1_Adjust(3 * i+1) << setw(20) << L(3 * i+1) << setw(20) << Lv(3 * i+1) << endl;

        cout << endl << setw(20) << S1_Adjust(3 * i + 2) << setw(20) << L(3 * i + 2) << setw(20) << Lv(3 * i + 2) << endl;

    }
    //cout << endl << Lv << endl;
}

void CCollocation::PrintResult(char* resultfile)
{
    ofstream outfile(resultfile, ios::out);
    outfile << "-------------------------平差结果--------------------------\n";
    outfile << "----------------------各站一年内XYZ方向运动速度------------------\n";
    

    outfile << "\n站点       方向          平差前速度            平差后速度------(单位:mm/year)\n";

    for (int i = 0; i < Pnumber; i++)
    {
       

        outfile <<"  "<<  setw(20) << "  X  " << setw(20) << L(3 * i) * 1000.0 << setw(20) << Lv(3 * i) * 1000.0 << endl;

        outfile  << Pname[i] << setw(15) << "  Y  " << setw(20) << L(3 * i+1) * 1000.0 << setw(20) << Lv(3 * i+1) * 1000.0 << endl;

        outfile  <<"  "<< setw(20) << "  Z  " << setw(20) << L(3 * i+2) * 1000.0 << setw(20) << Lv(3 * i+2) * 1000.0 << endl;

        outfile << "-------------------------------------------------------------------------" << endl;


    }

    outfile << "\n\n-----------------------欧拉矢量---------------------\n";

    outfile << "wx: " << detax(0) << setw(20) << "wy:  " << detax(1) << setw(20) << "wz:  " << detax(2) << endl;

    outfile << "\n\n--------------------已测点信号估值----------------------\n";

    outfile << "点名       方向        已测点信号估值------(单位:mm/year) \n";

    for (int i = 0; i < Pnumber; i++)
    {
        outfile <<" " << setw(20) <<" X "<<setw(20)<< S1_Adjust(3*i)*1000.0 << endl;
        outfile << Pname[i] << setw(15) <<" Y "<<setw(20)<< S1_Adjust(3 * i+1)*1000.0 << endl;
        outfile <<" " << setw(20) <<" Z "<<setw(20)<< S1_Adjust(3 * i+2)*1000.0 << endl;
        outfile << "-----------------------------------------------" << endl;
    
    }


    outfile << "\n\n--------------------未测点信号估值----------------------\n";

    outfile << "点名       方向        未测点信号估值------(单位:mm/year) \n";

    for (int i = 0; i < unusedPnumber; i++)
    {
        outfile << " " << setw(20) << " X " << setw(20) << S2_Adjust(3 * i) * 1000.0 << endl;
        outfile << unusedPname[i] << setw(15) << " Y " << setw(20) << S2_Adjust(3 * i + 1) * 1000.0 << endl;
        outfile << " " << setw(20) << " Z " << setw(20) << S2_Adjust(3 * i + 2) * 1000.0 << endl;
        outfile << "-----------------------------------------------" << endl;

    }
}