#include <iostream> 
#include <Eigen/Dense>
#include <complex>
#include <Eigen/Eigenvalues>


#define pai 3.1415926535898
using namespace Eigen;
using namespace std;

class typeofrlt_ms{
       public:
        MatrixXcd MS,Residue;
    } result;//用来存储返回矩阵的对象

typeofrlt_ms Mode_shape( const MatrixXcd & FRF,const VectorXcd & Mode_ModalParameters,const int & Fre_sample,const int & Delay_est){
    
    int sample_num=FRF.cols();//
    int ModalParameters_num=Mode_ModalParameters.rows();
    int DOF_measure=FRF.rows();//degree of freedom 
    MatrixXcd ms(DOF_measure,ModalParameters_num),residue(DOF_measure,ModalParameters_num),tmp_res(DOF_measure,ModalParameters_num+Delay_est),base_func(ModalParameters_num+Delay_est,sample_num);
    
    for(int i = 0;i < ModalParameters_num+Delay_est;i++){
        //基函数组成矩阵的行数=模态参数个数+空间抖动裕度
        for(int j = 0;j < sample_num;j++){
            if (i < ModalParameters_num){
                base_func(i,j)=1.0/(1.0-exp(-1i*2.0*pai/double(sample_num)*double(j))*exp(Mode_ModalParameters(i)/double(Fre_sample)));
            }
            else{
                base_func(i,j)=exp(-1i*2.0*pai/double(sample_num)*double(j*(i-ModalParameters_num)));//
            }
        } 
    }
    tmp_res=FRF*base_func.adjoint()*(base_func*base_func.adjoint()).inverse();//
    
    residue=tmp_res.leftCols(ModalParameters_num);//对tmp_res切片前ModalParameters_num列,
    for(int i= 0 ;i < ModalParameters_num;i++){
        ms.col(i)=residue.col(i)/residue.col(i).norm();//get the modeshape by nomalizing residue,norm对于vertor来说是2范数，对矩阵来说是Frobinus范数
    }
    
    result.MS=ms;
    result.Residue=residue;
    return result;
}