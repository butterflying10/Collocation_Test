
#include <iostream>
#include <Eigen/Dense>
#include<Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include "conio.h"
#include "CCollocation.h"

#include "public.h"
using namespace std;
using namespace Eigen;

/*最小二乘配置*/
int main()
{
    CCollocation ccollocation;

    char* datafile = "C:\\Butterflying\\grade_first\\adjust\\Collocation_Test\\collocationdata.txt";

    ccollocation.InputData(datafile);
    

    ccollocation.Get_L();
    ccollocation.Get_unused();
    ccollocation.Get_PL();
    ccollocation.Get_A();
    ccollocation.Get_detax();
    
    
    
    ccollocation.Get_QS1();
    ccollocation.Get_S1_Adjust();
    ccollocation.Get_Lv();

    ccollocation.PrintResult("C:\\Butterflying\\grade_first\\adjust\\Collocation_Test\\result.txt");
    


}

