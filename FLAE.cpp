//Author: Jin Wu, Zebo Zhou et al.
//e-mail: jin_wu_uestc@hotmail.com;	klinsmann.zhou@gmail.com
//The C++ codes of FLAE, a Fast Linear Attitude Estimator for Wahba’s problem.

#include "StartUP.h"
#include "FLAE.h"

static double Power(double a,int b);

FLAE::FLAE():
      N(3)
{
}

void FLAE::setParams(
           int n,
           double R[],
           double B[],
           double w[])
{ 
  N=n;
  for(int i=0;i<N;++i)
  {
    ref[i](0)=   R[i*3+0];        ref[i](1)=   R[i*3+1];            ref[i](2)=   R[i*3+2];
    body[i](0)=  B[i*3+0];        body[i](1)=  B[i*3+1];            body[i](2)=  B[i*3+2];
    weights[i]=  w[i];
  }
}

void FLAE::SolverNewton(void)
{
  MM.setZero();
  
  for(int i=0;i<N;++i)
  {
    MM=MM+weights[i]*body[i]*(ref[i].transpose());
  }
  
  double Hx1,Hx2,Hx3,Hy1,Hy2,Hy3,Hz1,Hz2,Hz3;
  Hx1=MM(0,0);    Hx2=MM(0,1);    Hx3=MM(0,2);
  Hy1=MM(1,0);    Hy2=MM(1,1);    Hy3=MM(1,2);
  Hz1=MM(2,0);    Hz2=MM(2,1);    Hz3=MM(2,2);
  
  W(0,0)=Hx1+Hy2+Hz3;   W(0,1)=-Hy3+Hz2;        W(0,2)=-Hz1+Hx3;        W(0,3)=-Hx2+Hy1;
  W(1,0)=-Hy3+Hz2;      W(1,1)= Hx1-Hy2-Hz3;    W(1,2)=Hx2+Hy1;         W(1,3)=Hx3+Hz1;
  W(2,0)=-Hz1+Hx3;      W(2,1)=Hx2+Hy1;         W(2,2)=Hy2-Hx1-Hz3;     W(2,3)=Hy3+Hz2;
  W(3,0)=-Hx2+Hy1;      W(3,1)=Hx3+Hz1;         W(3,2)=Hy3+Hz2;         W(3,3)=Hz3-Hy2-Hx1;
  
  double Hx1_4,Hx1_2,Hx2_4,Hx2_2,Hx3_4,Hx3_2;
  double Hy1_4,Hy1_2,Hy2_4,Hy2_2,Hy3_4,Hy3_2;
  double Hz1_4,Hz1_2,Hz2_4,Hz2_2,Hz3_4,Hz3_2;
  
  Hx1_2=Hx1*Hx1;        Hx1_4=Hx1_2*Hx1_2;
  Hx2_2=Hx2*Hx2;        Hx2_4=Hx2_2*Hx2_2;
  Hx3_2=Hx3*Hx3;        Hx3_4=Hx3_2*Hx3_2;
  
  Hy1_2=Hy1*Hy1;        Hy1_4=Hy1_2*Hy1_2;
  Hy2_2=Hy2*Hy2;        Hy2_4=Hy2_2*Hy2_2;
  Hy3_2=Hy3*Hy3;        Hy3_4=Hy3_2*Hy3_2;
  
  Hz1_2=Hz1*Hz1;        Hz1_4=Hz1_2*Hz1_2;
  Hz2_2=Hz2*Hz2;        Hz2_4=Hz2_2*Hz2_2;
  Hz3_2=Hz3*Hz3;        Hz3_4=Hz3_2*Hz3_2;
  
  detW=
   Hx1_4 + Hx2_4 + Hx3_4 + Hy1_4 + Hy2_4 + Hy3_4 + Hz1_4 + Hz2_4 + Hz3_4 +
   2.0*(Hx1_2*Hx2_2 + 
   Hx1_2*Hx3_2 + Hx2_2*Hx3_2 + 
   Hx1_2*Hy1_2 - Hx2_2*Hy1_2 - 
   Hx3_2*Hy1_2 + 4.0*Hx1*Hx2*Hy1*Hy2 - 
   Hx1_2*Hy2_2 + Hx2_2*Hy2_2 - 
   Hx3_2*Hy2_2 + Hy1_2*Hy2_2 + 
   4.0*Hx1*Hx3*Hy1*Hy3 + 4.0*Hx2*Hx3*Hy2*Hy3 - Hx1_2*Hy3_2 - 
   Hx2_2*Hy3_2 + Hx3_2*Hy3_2 + 
   Hy1_2*Hy3_2 + Hy2_2*Hy3_2 + 
   Hx1_2*Hz1_2 - Hx2_2*Hz1_2 - 
   Hx3_2*Hz1_2 + Hy1_2*Hz1_2 - 
   Hy2_2*Hz1_2 - Hy3_2*Hz1_2 + 
   4.0*Hx1*Hx2*Hz1*Hz2 + 4.0*Hy1*Hy2*Hz1*Hz2 - Hx1_2*Hz2_2 + 
   Hx2_2*Hz2_2 - Hx3_2*Hz2_2 - 
   Hy1_2*Hz2_2 + Hy2_2*Hz2_2 - 
   Hy3_2*Hz2_2 + Hz1_2*Hz2_2 + 
   4.0*Hx1*Hx3*Hz1*Hz3 + 4.0*Hy1*Hy3*Hz1*Hz3 + 4.0*Hx2*Hx3*Hz2*Hz3 + 4.0*Hy2*Hy3*Hz2*Hz3 - 
   Hx1_2*Hz3_2 - Hx2_2*Hz3_2 + 
   Hx3_2*Hz3_2 - Hy1_2*Hz3_2 - 
   Hy2_2*Hz3_2 + Hy3_2*Hz3_2 + 
   Hz1_2*Hz3_2 + Hz2_2*Hz3_2 );
  
  c=detW;
  b=8.0*(Hx3*Hy2*Hz1 - Hx2*Hy3*Hz1 - Hx3*Hy1*Hz2 + Hx1*Hy3*Hz2 + Hx2*Hy1*Hz3 - Hx1*Hy2*Hz3);
  a= -2.0*(Hx1_2+Hx2_2+Hx3_2+
                  Hy1_2+Hy2_2+Hy3_2+
                  Hz1_2+Hz2_2+Hz3_2);
  
  lambda=1.0;
  double old_lambda=0.0;

  while(fabsf(old_lambda-lambda)>1e-5)
  {
     old_lambda=lambda;
     lambda=lambda-((Power(lambda,4)+a*lambda*lambda+b*lambda+c)/(4.0*Power(lambda,3)+2.0*a*lambda+b));
  }
  
  for(int i=0;i<4;++i)
    for(int j=0;j<4;++j)
    {
      if(i==j)
        G(i,j)=(double)(W(i,j)-lambda);
      else
        G(i,j)=(double)W(i,j);
    }
  
  pivot = G(0,0);  
  G(0,0) /= pivot;              G(0,1) /= pivot;                G(0,2) /= pivot;                G(0,3) /= pivot;
  G(1,0) -= G(1,0)*G(0,0);      G(1,1) -= G(1,0)*G(0,1);        G(1,2) -= G(1,0)*G(0,2);        G(1,3) -= G(1,0)*G(0,3);
  G(2,0) -= G(2,0)*G(0,0);      G(2,1) -= G(2,0)*G(0,1);        G(2,2) -= G(2,0)*G(0,2);        G(2,3) -= G(2,0)*G(0,3);
  G(3,0) -= G(3,0)*G(0,0);      G(3,1) -= G(3,0)*G(0,1);        G(3,2) -= G(3,0)*G(0,2);        G(3,3) -= G(3,0)*G(0,3);

  pivot = G(1,1);
  G(1,0) /= pivot;              G(1,1) /= pivot;                G(1,2) /= pivot;                G(1,3) /= pivot;
  G(0,0) -= G(0,1)*G(1,0);      G(0,1) -= G(0,1)*G(1,1);        G(0,2) -= G(0,1)*G(1,2);        G(0,3) -= G(0,1)*G(1,3);
  G(2,0) -= G(2,1)*G(1,0);      G(2,1) -= G(2,1)*G(1,1);        G(2,2) -= G(2,1)*G(1,2);        G(2,3) -= G(2,1)*G(1,3);
  G(3,0) -= G(3,1)*G(1,0);      G(3,1) -= G(3,1)*G(1,1);        G(3,2) -= G(3,1)*G(1,2);        G(3,3) -= G(3,1)*G(1,3);
   
  pivot = G(2,2);
  G(2,0) /= pivot;              G(2,1) /= pivot;                G(2,2) /= pivot;                G(2,3) /= pivot;
  G(0,0) -= G(0,2)*G(2,0);      G(0,1) -= G(0,2)*G(2,1);        G(0,2) -= G(0,2)*G(2,2);        G(0,3) -= G(0,2)*G(2,3);
  G(1,0) -= G(1,2)*G(2,0);      G(1,1) -= G(1,2)*G(2,1);        G(1,2) -= G(1,2)*G(2,2);        G(1,3) -= G(1,2)*G(2,3);
  G(3,0) -= G(3,2)*G(2,0);      G(3,1) -= G(3,2)*G(2,1);        G(3,2) -= G(3,2)*G(2,2);        G(3,3) -= G(3,2)*G(2,3);

  quaternion(0)=G(0,3);         quaternion(1)=G(1,3);           quaternion(2)=G(2,3);           quaternion(3)=-1.0;
  quaternion.normalize();
}

matrix::Quaternion<double>     FLAE::getQ(void)
{
  return quaternion;
}

static double Power(double a,int b)
{
  double res=1.0f;
  for(int i=0;i<b;++i)
    res*=a;
  return (res);
}

