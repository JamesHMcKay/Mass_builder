// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

// File generated at Sat 2 Sep 2017 19:00:27

#include "MSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

namespace {

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(double n, const Eigen::MatrixBase<Derived>& m)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m - Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(double n, const Eigen::MatrixBase<Derived>& m)
{
   return - m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

} // anonymous namespace

/**
 * Calculates the one-loop beta function of MassG.
 *
 * @return one-loop beta function
 */
double MSSM_soft_parameters::calc_beta_MassG_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MassG;

   beta_MassG = Re(-6*MassG*oneOver16PiSqr*Sqr(g3));


   return beta_MassG;
}

/**
 * Calculates the two-loop beta function of MassG.
 *
 * @return two-loop beta function
 */
double MSSM_soft_parameters::calc_beta_MassG_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;


   double beta_MassG;

   beta_MassG = Re(0.4*twoLoop*Sqr(g3)*(11*(MassB + MassG)*Sqr(g1) + 5*(9
      *(MassG + MassWB)*Sqr(g2) + 4*(traceAdjYdTYd + traceAdjYuTYu - MassG*
      traceYdAdjYd - MassG*traceYuAdjYu + 7*MassG*Sqr(g3)))));


   return beta_MassG;
}

/**
 * Calculates the three-loop beta function of MassG.
 *
 * @return three-loop beta function
 */
double MSSM_soft_parameters::calc_beta_MassG_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   const double traceAdjYdYd = TRACE_STRUCT.traceAdjYdYd;
   const double traceAdjYuYu = TRACE_STRUCT.traceAdjYuYu;
   const double traceTYdAdjYd = TRACE_STRUCT.traceTYdAdjYd;
   const double traceAdjYeYe = TRACE_STRUCT.traceAdjYeYe;
   const double traceTYeAdjYe = TRACE_STRUCT.traceTYeAdjYe;
   const double traceTYuAdjYu = TRACE_STRUCT.traceTYuAdjYu;
   const double traceAdjYdYdAdjYdYd = TRACE_STRUCT.traceAdjYdYdAdjYdYd;
   const double traceAdjYdTYdAdjYdYd = TRACE_STRUCT.traceAdjYdTYdAdjYdYd;
   const double traceAdjYuYuAdjYdYd = TRACE_STRUCT.traceAdjYuYuAdjYdYd;
   const double traceAdjYuYuAdjYuYu = TRACE_STRUCT.traceAdjYuYuAdjYuYu;
   const double traceAdjYuTYuAdjYdYd = TRACE_STRUCT.traceAdjYuTYuAdjYdYd;
   const double traceAdjYuTYuAdjYuYu = TRACE_STRUCT.traceAdjYuTYuAdjYuYu;
   const double traceTYdAdjYuYuAdjYd = TRACE_STRUCT.traceTYdAdjYuYuAdjYd;


   double beta_MassG;

   beta_MassG = Re(0.02666666666666667*threeLoop*Sqr(g3)*(-1702*Power(g1,
      4)*(2*MassB + MassG) - 5*Sqr(g1)*(4*(8*MassB*traceAdjYdYd + 8*MassG*
      traceAdjYdYd + 11*MassB*traceAdjYuYu + 11*MassG*traceAdjYuYu - 8*
      traceTYdAdjYd - 11*traceTYuAdjYu) + 9*(MassB + MassG + MassWB)*Sqr(g2) -
      22*(MassB + 2*MassG)*Sqr(g3)) - 25*(-1041*Power(g3,4)*MassG + 81*Power(g2
      ,4)*(MassG + 2*MassWB) + 104*(2*MassG*(traceAdjYdYd + traceAdjYuYu) -
      traceTYdAdjYd - traceTYuAdjYu)*Sqr(g3) - 18*Sqr(g2)*(-2*(MassG*(
      traceAdjYdYd + traceAdjYuYu) + MassWB*(traceAdjYdYd + traceAdjYuYu) -
      traceTYdAdjYd - traceTYuAdjYu) + (2*MassG + MassWB)*Sqr(g3)) + 6*(12*
      traceAdjYdTYdAdjYdYd + 4*traceAdjYuTYuAdjYdYd + 12*traceAdjYuTYuAdjYuYu +
      18*traceAdjYdYd*traceTYdAdjYd + 3*traceAdjYeYe*traceTYdAdjYd + 4*
      traceTYdAdjYuYuAdjYd + 3*traceAdjYdYd*traceTYeAdjYe + 18*traceAdjYuYu*
      traceTYuAdjYu - MassG*(6*traceAdjYdYdAdjYdYd + 3*traceAdjYdYd*
      traceAdjYeYe + 4*traceAdjYuYuAdjYdYd + 6*traceAdjYuYuAdjYuYu + 9*Sqr(
      traceAdjYdYd) + 9*Sqr(traceAdjYuYu))))));


   return beta_MassG;
}

} // namespace flexiblesusy
