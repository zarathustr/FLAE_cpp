#ifndef __FLAE_H
#define __FLAE_H

#include "MATH_LIB.h"

#ifdef __cplusplus
extern "C"{
#endif

class FLAE
{
public:
  FLAE();
  ~FLAE(){};
  
  //void SolverSymbolic(void);
  void setParams(
       int n,
       double R[],
       double B[],
       double w[]);
  void SolverNewton(void);
  matrix::Quaternion<double>   getQ(void);
  
private:
  int                           N;             //Numbers of Vector Observation Pairs
  matrix::Vector3<double>              ref[4];           //Reference Vectors
  matrix::Vector3<double>              body[4];          //Body observations
  double                         weights[4];       //Weights
  
  double                                 lambda;
  matrix::Quaternion<double>             quaternion;    //Quaternion solution 
  matrix::SquareMatrix<double,4>         W;
  matrix::SquareMatrix<double,4>         G;
  matrix::SquareMatrix<double,3>         MM;
  double pivot;
  double a;
  double b;
  double c;
  double detW;
};

#ifdef __cplusplus
}
#endif

#endif