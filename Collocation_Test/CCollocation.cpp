#pragma once
#include<fstream>
#include <iostream>

#include <string>

#include < iomanip >//����С�����λ��

#include<cmath>

#include <Eigen/Dense>
#include<Eigen/Core>
#include <Eigen/Sparse>//ϡ�����
#include "CCollocation.h"
#include "public.h" //�����м�����λ���ķ���


using namespace Eigen;
using namespace std;

//���캯��
CCollocation::CCollocation()
{

}

//��������
CCollocation::~CCollocation()
{

}

/*����ԭʼ����(��1�ܺ͵�52�ܵĸ����������Э����)
* �洢��D1��D52��
* 
*/
int CCollocation::InputData(char* DataFile)
{
    ifstream infile(DataFile, ios::in | ios::_Nocreate);


    if (!infile)
    {
        cerr << "��ԭʼ����ʧ��" << endl;
        return 0;
    }

    string note;//��ע
    getline(infile, note);
    

    infile >> note;//��ȡ��������
    

    
    infile >> Pnumber1;

    //cout << Pnumber1 << endl;

    Pname1 = new char* [Pnumber1];//��һ�ܵĵ�����ַ����
    for (int i = 0; i < Pnumber1; i++)Pname1[i] = NULL;
    XYZ1 = new double[Pnumber1 * 3];//��һ�� ���XYZ����
    Qxyz1 = new double[Pnumber1 * 3];//��һ�� ���XYZ����ķ���


    getline(infile, note);
    getline(infile, note);
    getline(infile, note);
    getline(infile, note);
    getline(infile, note);

    /*
    ����һ�ܵ����ݵ� XYZ����  ����
    */

    for (int i = 0; i < Pnumber1; i++)
    {
        int note;
        infile >> note;

        //cout << note << endl;

        char pname[20];
        infile >> pname;
        int len = strlen(pname);
        Pname1[i] = new char[len + 1];//�����ַ��������ı�־��/0��
        strcpy_s(Pname1[i], strlen(pname) + 1, pname);
    
      

        //�洢��һ�ܸ��������
        infile >> XYZ1[3 * i];
        infile >> XYZ1[3 * i + 1];
        infile >> XYZ1[3 * i + 2];

        //cout << XYZ1[3 * i] << setw(20) << XYZ1[3 * i + 1] << setw(20) << XYZ1[3 * i + 2] << endl;

        //�洢��һ�ܸ���ķ���
        infile >> Qxyz1[3 * i];
        infile >> Qxyz1[3 * i + 1];
        infile >> Qxyz1[3 * i + 2];

        //cout << Qxyz1[3 * i] << setw(20) << Qxyz1[3 * i + 1] << setw(20) << Qxyz1[3 * i + 2] << endl;


    }


    getline(infile, note);
    getline(infile, note);
    getline(infile, note);

    infile >> note;//��ȡ��������

    //int Pnumber52;
    infile >> Pnumber52;//����ʮ���ܵ���

    //cout << Pnumber52 << endl;

    Pname52 = new char* [Pnumber52];//��52�ܵĵ�����ַ����
    for (int i = 0; i < Pnumber52; i++)Pname52[i] = NULL;
    XYZ52 = new double[Pnumber52 * 3];//��52�� ���XYZ����
    Qxyz52 = new double[Pnumber52 * 3];//��52�� ���XYZ����ķ���

    getline(infile, note);
    getline(infile, note);
    getline(infile, note);
    getline(infile, note);
    getline(infile, note);


    /*
    ����52�ܵ����ݵ� XYZ����  ����
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
        


        //�洢��һ�ܸ��������
        infile >> XYZ52[3 * i];
        infile >> XYZ52[3 * i + 1];
        infile >> XYZ52[3 * i + 2];

        //cout << XYZ52[3 * i] << setw(20) << XYZ52[3 * i + 1] << setw(20) << XYZ52[3 * i + 2] << endl;

        //�洢��һ�ܸ���ķ���
        infile >> Qxyz52[3 * i];
        infile >> Qxyz52[3 * i + 1];
        infile >> Qxyz52[3 * i + 2];

        //cout << Qxyz52[3 * i] << setw(20) << Qxyz52[3 * i + 1] << setw(20) << Qxyz52[3 * i + 2] << endl;

    }


    infile.close();

}

void CCollocation::Get_unused()
{
    //��ȡδ����źŵ������δ����źŵĸ���
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
����۲�ֵ����͹۲�ֵЭ������
*/

void CCollocation::Get_L()
{

    /*
    ������С�������õĵ���n
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


    //QLҪ��һ�ַ������壬�����Ƕ�ȡ�ȿ�������ƽ��Ľ��������������һ���򵥵ĶԽ���
    ////cout << v << endl;
    //QL = v.asDiagonal();
    //PL = QL.inverse();
        ////��һ��   ��һ����ľ����ÿһ�л����е�Ԫ�ص�ƽ����Ϊ1
    //PL.normalize();
    //cout << PL << endl;



    

}
/*����۲�ֵȨ����*/

void CCollocation::Get_PL()
{
    

    double* Lx = new double[Pnumber];//�洢ԭ�۲�ֵX����������
    double* Ly = new double[Pnumber];//�洢ԭ�۲�ֵY����������
    double* Lz = new double[Pnumber];//�洢ԭ�۲�ֵZ����������




    detaLx = new double[Pnumber];//�洢detaL   X����������
    detaLy = new double[Pnumber];//�洢detaL   Y����������
    detaLz = new double[Pnumber];//�洢detaL   Z����������



   /*
   �۲�ֵ
   */
    for (int i = 0; i < Pnumber; i++)
    {
        Lx[i] = L[3 * i];
        Ly[i] = L[3 * i + 1];
        Lz[i] = L[3 * i + 2];

        
    }

    //��ȡLx�������λ��
    double medLx = Median(Lx, Pnumber, false);
    double medLy = Median(Ly, Pnumber, false);
    double medLz = Median(Lz, Pnumber, false);

    //����ڼ���C0��k��ʱ���õ���
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


    //��ʼ����PL
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
A:::�۲�ֵĳһ�����Ĺ۲�ֵ
n:::����
medA:::��λ��
*/

double* CCollocation::GGDD(double A[], int n, double medA)
{
    //����PL�ĳ�ʼ����  ��λ�󣬵���Ϊ�˷��������ֻ��ŶԽ���Ԫ�أ���Ϊһ������

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
        //�����ʼ��������
        double Sigmax0 = Median(detaA, n, true) / 0.6745;

        //cout << "Sigmax0������" << Sigmax0 << endl;


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
       
        //�����Ȩƽ��ֵ
        for (int i = 0; i < n; i++)
        {
            newmedA += P[i] * A[i];
            sumP += P[i];

        }
        newmedA = newmedA / sumP;

       // cout << "GGDD:::��������:::" << ii << "  medA:::" << medA << "  newmedA:::" << newmedA << endl;

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
����ϵ����A
*/
void CCollocation::Get_A()
{
    //A = Eigen::Matrix<double, Dynamic, Dynamic>();
    A.resize(Pnumber * 3, 3);
    AT.resize(3, Pnumber * 3);

    for (int i = 0; i < Pnumber; i++)
    {
        /*
        Eigen�о���ĸ�ֵ����ͨ��ѭ��һ��һ����ֵ
        */
        /*A << 0.0, XYZ[3 * i + 2], -XYZ[3 * i + 1],
            -XYZ[3 * i + 2], 0.0, XYZ[3 * i],
            XYZ[3 * i + 1], -XYZ[3 * i], 0.0;*/

        //�ֿ��ʼ����ֵ
        A.block(3 * i, 0, 3, 3)
            << 0.0, XYZ[3 * i + 2], -XYZ[3 * i + 1],
            -XYZ[3 * i + 2], 0.0, XYZ[3 * i],
            XYZ[3 * i + 1], -XYZ[3 * i], 0.0;
    }
    

    //ת��
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
    //�����źŵĳ�ʼ��������
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
//    //��i����
//    for (int i = 0; i < Pnumber; i++)
//    {
//        //��j����
//        for (int j = i + 1; j < Pnumber; j++)
//        {
//            //������������尡�����������ֲ�����ľ���,Ҳ��������ģ�
//
//            double Sx = (XYZ[3 * i] - XYZ[3 * j]) * (XYZ[3 * i] - XYZ[3 * j]);
//            double Sy = (XYZ[3 * i + 1] - XYZ[3 * j + 1]) * (XYZ[3 * i + 1] - XYZ[3 * j + 1]);
//            double Sz = (XYZ[3 * i + 2] - XYZ[3 * j + 2]) * (XYZ[3 * i + 2] - XYZ[3 * j + 2]);
//
//            // cout << Sx << setw(20) << Sy << setw(20) << Sz << endl;
//
//             //��S����ǧ�׼�������
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


/*����k
*/
double  CCollocation::Get_k(int Pnumber,double * detaLx,double *XYZ)
{

    //�����źŵĳ�ʼ��������
    double c=0.0;
    for (int i = 0; i < Pnumber; i++) c = c+ detaLx[i]*detaLx[i];
    c = c / Pnumber;

    

    //�����źŵĳ�ʼЭ����
     double* CSx = new  double[Pnumber * (Pnumber - 1) / 2];
     double* S = new  double[Pnumber * (Pnumber - 1) / 2];

    int k = 0;
    double* Kx = new double[Pnumber * (Pnumber - 1) / 2];

    //�Ȱ����е�ֵŪ�ɸ���
    for (int i = 0; i < Pnumber * (Pnumber - 1) / 2; i++) Kx[i] = -1.0;

    //��Kx����������ж��ٸ�����
    int positiveNumber = 0;

    //��i����
    for (int i = 0; i < Pnumber; i++)
    {   
        //��j����
        for (int j = i+1; j < Pnumber; j++)
        {
            CSx[k] = detaLx[i] * detaLx[j];

            double S = ((XYZ[3 * j] - XYZ[3 * i]) * (XYZ[3 * j] - XYZ[3 * i]) + (XYZ[3 * j + 1] - XYZ[3 * i + 1]) * (XYZ[3 * j + 1] - XYZ[3 * i + 1]) + (XYZ[3 * j + 2] - XYZ[3 * i + 2]) * (XYZ[3 * j + 2] - XYZ[3 * i + 2]))/1000000.0;
            //����x�����µ�k1,k2,k3,k4,k5.........

            if (CSx[k] > 0) Kx[k] = (log(c) - log(CSx[k])) / S;

            if (Kx[k] > 0) positiveNumber++;
            
            
            k++;
        }   
    }
    //cout <<"�����ĸ���Ϊ��"<< Number << endl;
    //ȥ��kx�еĸ�ֵ����ȡ��λ��
    //����ֻ��������Kx
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


    //������������Ҫ�����õ�medkx
     //���콵Ȩ����
    //����һ���Խ������Կ���������������������ʾ
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

         

         

         //���콵Ȩ����
         //����һ���Խ������Կ���������������������ʾ
        


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

         //��������k�ļ�Ȩƽ��ֵ

         double newmedKx = 0.0;
         for (int i = 0; i < positiveNumber; i++)
         {
             newmedKx += P(i) * positiveKx[i];

         }
         newmedKx = newmedKx / P.sum();


        // cout <<"��������:::"<<ii<< "   medKx:::" << medKx << "     " << "newmedKx" << newmedKx << endl;
         


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
/*�����Ѳ��źŵķ���Э������*/
void CCollocation::Get_QS1()
{

    Get_C0();
    //Get_diStance();
    //����k
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

    //��֤�Բ���
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
    cout << "����                   �Ѳ��ź�S                     Lƽ��ǰ             Lƽ���" << endl;
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
    outfile << "-------------------------ƽ����--------------------------\n";
    outfile << "----------------------��վһ����XYZ�����˶��ٶ�------------------\n";
    

    outfile << "\nվ��       ����          ƽ��ǰ�ٶ�            ƽ����ٶ�------(��λ:mm/year)\n";

    for (int i = 0; i < Pnumber; i++)
    {
       

        outfile <<"  "<<  setw(20) << "  X  " << setw(20) << L(3 * i) * 1000.0 << setw(20) << Lv(3 * i) * 1000.0 << endl;

        outfile  << Pname[i] << setw(15) << "  Y  " << setw(20) << L(3 * i+1) * 1000.0 << setw(20) << Lv(3 * i+1) * 1000.0 << endl;

        outfile  <<"  "<< setw(20) << "  Z  " << setw(20) << L(3 * i+2) * 1000.0 << setw(20) << Lv(3 * i+2) * 1000.0 << endl;

        outfile << "-------------------------------------------------------------------------" << endl;


    }

    outfile << "\n\n-----------------------ŷ��ʸ��---------------------\n";

    outfile << "wx: " << detax(0) << setw(20) << "wy:  " << detax(1) << setw(20) << "wz:  " << detax(2) << endl;

    outfile << "\n\n--------------------�Ѳ���źŹ�ֵ----------------------\n";

    outfile << "����       ����        �Ѳ���źŹ�ֵ------(��λ:mm/year) \n";

    for (int i = 0; i < Pnumber; i++)
    {
        outfile <<" " << setw(20) <<" X "<<setw(20)<< S1_Adjust(3*i)*1000.0 << endl;
        outfile << Pname[i] << setw(15) <<" Y "<<setw(20)<< S1_Adjust(3 * i+1)*1000.0 << endl;
        outfile <<" " << setw(20) <<" Z "<<setw(20)<< S1_Adjust(3 * i+2)*1000.0 << endl;
        outfile << "-----------------------------------------------" << endl;
    
    }


    outfile << "\n\n--------------------δ����źŹ�ֵ----------------------\n";

    outfile << "����       ����        δ����źŹ�ֵ------(��λ:mm/year) \n";

    for (int i = 0; i < unusedPnumber; i++)
    {
        outfile << " " << setw(20) << " X " << setw(20) << S2_Adjust(3 * i) * 1000.0 << endl;
        outfile << unusedPname[i] << setw(15) << " Y " << setw(20) << S2_Adjust(3 * i + 1) * 1000.0 << endl;
        outfile << " " << setw(20) << " Z " << setw(20) << S2_Adjust(3 * i + 2) * 1000.0 << endl;
        outfile << "-----------------------------------------------" << endl;

    }
}