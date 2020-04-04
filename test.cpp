#include <iostream> 
#include <Eigen/Dense>
#include <Eigen/LU>
#include <complex>
#include <Eigen/Eigenvalues>

#include"ModeShape.h"

#define pai 3.1415926535898
using namespace Eigen;
using namespace std;


int main()  
{  
    double f=100;
    double f2=300;
    double f3=600;
    double sigma=10;
    double sigma2=20;
    double sigma3=30;
    double fs=5000;
    int N=5000;   
    VectorXd t(N);
    MatrixXd x(1,N);
    MatrixXd x_tdomin(3,N);
    VectorXcd ModalParameters(6);
    MatrixXcd frf(3,N),x_in_f(1,N),din(3,N);
    t.setLinSpaced(N,0,1-1/N);   
    //cout<<t<<endl;
    x=((-sigma*t).array().exp()).cwiseProduct((2.*pai*f*t).array().sin())+((-sigma2*t).array().exp()).cwiseProduct((2.*pai*f2*t).array().sin())+((-sigma3*t).array().exp()).cwiseProduct((2.*pai*f3*t).array().sin());
    
    x_tdomin.row(0)=x.transpose();
    x_tdomin.row(1)=2*x.transpose();
    x_tdomin.row(2)=3*x.transpose();
    MatrixXd x_in(1,N);
    x_in.setZero(1,N);
    x_in(0)=1;
    ModalParameters << -sigma+1i*2.*pai*f , -sigma2+1i*2.*pai*f2 , -sigma3+1i*2.*pai*f3 , -sigma-1i*2.*pai*f , -sigma2-1i*2.*pai*f2 , -sigma3-1i*2.*pai*f3;
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            for(int k=0;k<3;k++){
                frf(k,i)=frf(k,i)+exp(-1i*2.*pai*double(i*j/N))*x_tdomin(k,j);
            } 
            x_in_f(i)=x_in_f(i)+exp(-1i*2.*pai*double(i*j/N))*x_in(j);       
        }
    }
    din.row(0)=frf.row(0);
    din.row(1)=frf.row(1);
    din.row(2)=frf.row(2);
    
    int Delay=0;
    typeofrlt_ms rr=Mode_shape(din,ModalParameters,fs,Delay);
    
    //typeofrlt_plcf poles=PLCF(din,fs,20,1,Delay);
    
    cout << "ModalParametersmeters: \n"<<ModalParameters << endl;
    cout << "FRF: \n"<< frf << endl;
    cout << "Modeshape: \n"<<rr.MS << endl;
    cout << "Residue: \n"<<rr.Residue << endl;
    //cout << "Eigenvalue: \n"<<poles.Eigenvalue << endl;
    //cout << "Damp Coefficients: \n"<<poles.Damp << endl;
    //cout << "Natuarl Frequency: \n"<<poles.Freq << endl;
    
    
}  