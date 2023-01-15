//------------------------------------------------------------------------------
// GaussJacobiQuadrature.hpp: Definition of the class, which contains the
//                            variables and functions to create the integration
//                            points for Gauss Jacobi quadrature.
//------------------------------------------------------------------------------

#pragma once

// Include files needed.

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>

// Definition of the class CGaussJacobiQuadrature.
class CGaussJacobiQuadrature {
public:

  /*!
   * \brief Constructor of the class, nothing to be done.
   */
  CGaussJacobiQuadrature(){}

  /*!
   * \brief Destructor of the class, nothing to be done.
   */
  ~CGaussJacobiQuadrature(){}

  /*!
   * \brief Function, which serves as the API to compute the integration points
            and weights.
   * \param[in]     alpha     Parameter in the weighting function (b-x)^alpha*(x-a)^beta
                              in the Gauss Jacobi rule.
   * \param[in]     beta      Parameter in the weighting function (b-x)^alpha*(x-a)^beta
                              in the Gauss Jacobi rule.
   * \param[in]     a         Lower bound of the integration interval, usually -1.0.
   * \param[in]     b         Upper bound of the integration interval, usually  1.0.
   * \param[in,out] GJPoints  Location of the Gauss-Jacobi integration points.
   * \param[in,out] GJWeights Weights of the Gauss-Jacobi integration points.
   */
  void GetQuadraturePoints(const double   alpha,          const double   beta,
                           const double   a,              const double   b,
                           std::vector<double> &GJPoints, std::vector<double> &GJWeights);
private:
  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void cdgqf(int nt, int kind, double alpha, double beta, double t[],
             double wts[]);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void cgqf(int nt, int kind, double alpha, double beta, double a,
            double b, double t[], double wts[]);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  double class_matrix(int kind, int m, double alpha, double beta,
                         double aj[], double bj[]);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void imtqlx(int n, double d[], double e[], double z[]);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void parchk(int kind, int m, double alpha, double beta);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  double r8_epsilon();

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  double r8_sign(double x);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void scqf(int nt, double t[], int mlt[], double wts[], int nwts, int ndx[],
            double swts[], double st[], int kind, double alpha,
            double beta, double a, double b);

  /*!
   * \brief Function in the original implementation of John Burkardt to compute
            the integration points of the Gauss-Jacobi quadrature rule.
   */
  void sgqf(int nt, double aj[], double bj[], double zemu, double t[],
            double wts[]);
};
