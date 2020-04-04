#include <iostream> 
#include <Eigen/Dense>
#include <complex>
#include <Eigen/Eigenvalues>

#define pai 3.1415926535898
using namespace Eigen;
using namespace std;

class typeofrlt_plcf{
       public:
        MatrixXcd Damp,Freq,Eigenvalue;
    } result_plcf;//用来存储返回矩阵的对象

typeofrlt_plcf PLCF(const MatrixXcd & FRF,const int & Fs,const int & ModelOrder,const int & DOF_Impulse,const int & Delay_est){
    int DOF_measure=FRF.rows();
    int Num_samples=FRF.cols();
    int Order_d=ModelOrder;
    int Order_n=ModelOrder+1+Delay_est;
    MatrixXcd Y(DOF_measure*Num_samples,DOF_Impulse);
    MatrixXcd HH(DOF_measure*Num_samples,DOF_measure*Order_n+DOF_Impulse*Order_d);
    MatrixXcd Para_est(DOF_measure*Order_n+DOF_Impulse*Order_d,DOF_Impulse);//Para_est.col=Y.col;Para_est.row=HH.col
    MatrixXcd tmp(DOF_Impulse,DOF_Impulse*Order_d);
    MatrixXcd CoefficientMatrix_NormalForm(DOF_Impulse*Order_d,DOF_Impulse*Order_d);
    
    MatrixXcd H(DOF_measure,DOF_Impulse*Order_d),I(DOF_measure,DOF_measure*Order_n);

    Y.setZero(DOF_measure*Num_samples,DOF_Impulse);
    HH.setZero(DOF_measure*Num_samples,DOF_measure*Order_n+DOF_Impulse*Order_d);

    
    for(int i = 0;i < Num_samples;i++){
        H.setZero(DOF_measure,DOF_Impulse*Order_d);
        I.setZero(DOF_measure,DOF_measure*Order_n);

        for(int j = 0;j < Order_d;j++){
            H.block(0,j*DOF_Impulse,DOF_measure,DOF_Impulse)=exp(1i*2.0*pai*double(i*j/Num_samples))*FRF.block(0,i*DOF_Impulse,DOF_measure,DOF_Impulse);
            //按列扩展矩阵
        }
        for(int k = 0;k < Order_n;k++){
            I.block(0,k*DOF_measure,DOF_measure,DOF_measure)=-exp(1i*2.0*pai*double(i*k/Num_samples))*MatrixXcd::Identity(DOF_measure,DOF_measure);
            //按列扩展矩阵
        }
        //HH.block(i*DOF_measure,0,DOF_measure,DOF_measure*Order_n+DOF_Impulse*Order_d) << H,I;
        HH.block(i*DOF_measure,0,DOF_measure,DOF_Impulse*Order_d) = H;
        HH.block(i*DOF_measure,DOF_Impulse*Order_d,DOF_measure,DOF_measure*Order_n)=I;
        //按行扩展矩阵
        Y.block(i*DOF_measure,0,DOF_measure,DOF_Impulse)=-exp(1i*2.0*pai*double(i*Order_d/Num_samples))*FRF.block(0,i*DOF_Impulse,DOF_measure,DOF_Impulse);
        //按行扩展矩阵
    }
    //得到ARMA的参数估计
    Para_est=(HH.adjoint()*HH).inverse()*(HH.adjoint()*Y);//
    //Para_est=HH.bdcSvd(ComputeThinU | ComputeThinV).solve(Y);
    //Para_est=(HH.adjoint() * HH).ldlt().solve(HH.adjoint() * Y);//Normal equation 
    //Para_est=HH.colPivHouseholderQr().solve(Y);//QR
    //cout<<Para_est<<endl;
    tmp.setZero(DOF_Impulse,Order_d*DOF_Impulse);
    //组成标准型
    for(int m = 0;m < Order_d-1;m++){
        tmp.block(0,(Order_d-m-1)*DOF_Impulse,DOF_Impulse,DOF_Impulse)=-Para_est.block(m*DOF_Impulse,0,DOF_Impulse,DOF_Impulse);
        //tmp.block(0,(Order_d-m-1)*DOF_Impulse,DOF_Impulse,DOF_Impulse)=
    }
    CoefficientMatrix_NormalForm << tmp,MatrixXd::Identity(DOF_Impulse*(Order_d-1),DOF_Impulse*(Order_d-1)),MatrixXd::Zero(DOF_Impulse*(Order_d-1),DOF_Impulse);
    ComplexEigenSolver<MatrixXcd> es(CoefficientMatrix_NormalForm);
    MatrixXcd ee=(es.eigenvalues().array().log())*double(Fs);
    result_plcf.Eigenvalue=ee;
    result_plcf.Freq=ee.imag()/2./pai;
    result_plcf.Damp=-ee.real();
    
    return result_plcf;
}
