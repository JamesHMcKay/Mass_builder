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

// File generated at Sat 2 Sep 2017 18:59:14

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
 * Calculates the one-loop beta function of TYe.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MSSM_soft_parameters::calc_beta_TYe_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (oneOver16PiSqr*(0.4*Ye*(9*MassB*Sqr(g1) + 5*(3*
      traceAdjYdTYd + traceAdjYeTYe + 3*MassWB*Sqr(g2))) + (3*traceYdAdjYd +
      traceYeAdjYe - 1.8*Sqr(g1) - 3*Sqr(g2))*TYe + 4*(Ye*Ye.adjoint()*TYe) + 5
      *(TYe*Ye.adjoint()*Ye))).real();


   return beta_TYe;
}

/**
 * Calculates the two-loop beta function of TYe.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MSSM_soft_parameters::calc_beta_TYe_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYe;

   beta_TYe = (twoLoop*(0.1*(-4*Ye*(135*Power(g1,4)*MassB + Sqr(g1)*(2*(
      traceAdjYdTYd - 3*traceAdjYeTYe - MassB*traceYdAdjYd + 3*MassB*
      traceYeAdjYe) + 9*(MassB + MassWB)*Sqr(g2)) + 5*(15*Power(g2,4)*MassWB +
      3*(6*traceYdAdjYdTYdAdjYd + traceYdAdjYuTYuAdjYd + 2*traceYeAdjYeTYeAdjYe
      + traceYuAdjYdTYdAdjYu) - 16*(traceAdjYdTYd - MassG*traceYdAdjYd)*Sqr(g3
      ))) + (135*Power(g1,4) + 2*Sqr(g1)*(-2*traceYdAdjYd + 6*traceYeAdjYe + 9*
      Sqr(g2)) + 5*(15*Power(g2,4) - 6*(3*traceYdAdjYdYdAdjYd +
      traceYdAdjYuYuAdjYd + traceYeAdjYeYeAdjYe) + 32*traceYdAdjYd*Sqr(g3)))*
      TYe) - 6*(3*traceAdjYdTYd + traceAdjYeTYe + 2*MassWB*Sqr(g2))*(Ye*
      Ye.adjoint()*Ye) + (-4*(3*traceYdAdjYd + traceYeAdjYe) + 1.2*Sqr(g1) + 6*
      Sqr(g2))*(Ye*Ye.adjoint()*TYe) + (-5*(3*traceYdAdjYd + traceYeAdjYe) -
      1.2*Sqr(g1) + 12*Sqr(g2))*(TYe*Ye.adjoint()*Ye) - 6*(Ye*Ye.adjoint()*Ye*
      Ye.adjoint()*TYe) - 8*(Ye*Ye.adjoint()*TYe*Ye.adjoint()*Ye) - 6*(TYe*
      Ye.adjoint()*Ye*Ye.adjoint()*Ye))).real();


   return beta_TYe;
}

/**
 * Calculates the three-loop beta function of TYe.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MSSM_soft_parameters::calc_beta_TYe_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   const double traceAdjYdYd = TRACE_STRUCT.traceAdjYdYd;
   const double traceAdjYeYe = TRACE_STRUCT.traceAdjYeYe;
   const double traceAdjYuYu = TRACE_STRUCT.traceAdjYuYu;
   const double traceTYdAdjYd = TRACE_STRUCT.traceTYdAdjYd;
   const double traceTYeAdjYe = TRACE_STRUCT.traceTYeAdjYe;
   const double traceTYuAdjYu = TRACE_STRUCT.traceTYuAdjYu;
   const double traceAdjYdYdAdjYdYd = TRACE_STRUCT.traceAdjYdYdAdjYdYd;
   const double traceAdjYdTYdAdjYdYd = TRACE_STRUCT.traceAdjYdTYdAdjYdYd;
   const double traceAdjYeYeAdjYeYe = TRACE_STRUCT.traceAdjYeYeAdjYeYe;
   const double traceAdjYeTYeAdjYeYe = TRACE_STRUCT.traceAdjYeTYeAdjYeYe;
   const double traceAdjYuYuAdjYdYd = TRACE_STRUCT.traceAdjYuYuAdjYdYd;
   const double traceAdjYuTYuAdjYdYd = TRACE_STRUCT.traceAdjYuTYuAdjYdYd;
   const double traceTYdAdjYuYuAdjYd = TRACE_STRUCT.traceTYdAdjYuYuAdjYd;
   const double traceAdjYdYdAdjYdYdAdjYdYd =
      TRACE_STRUCT.traceAdjYdYdAdjYdYdAdjYdYd;
   const double traceAdjYdYdAdjYdTYdAdjYdYd =
      TRACE_STRUCT.traceAdjYdYdAdjYdTYdAdjYdYd;
   const double traceAdjYeYeAdjYeYeAdjYeYe =
      TRACE_STRUCT.traceAdjYeYeAdjYeYeAdjYeYe;
   const double traceAdjYeYeAdjYeTYeAdjYeYe =
      TRACE_STRUCT.traceAdjYeYeAdjYeTYeAdjYeYe;
   const double traceAdjYuYuAdjYuYuAdjYdYd =
      TRACE_STRUCT.traceAdjYuYuAdjYuYuAdjYdYd;
   const double traceAdjYuYuAdjYuTYuAdjYdYd =
      TRACE_STRUCT.traceAdjYuYuAdjYuTYuAdjYdYd;
   const double traceAdjYuTYuAdjYuYuAdjYdYd =
      TRACE_STRUCT.traceAdjYuTYuAdjYuYuAdjYdYd;
   const double traceTYdAdjYuYuAdjYuYuAdjYd =
      TRACE_STRUCT.traceTYdAdjYuYuAdjYuYuAdjYd;


   Eigen::Matrix<double,3,3> beta_TYe;

   const Eigen::Matrix<double,3,3> beta_TYe_1 = ((-0.0013333333333333333*
      threeLoop*Ye*(449874*Power(g1,6)*MassB + 50*Sqr(g1)*(135*Power(g2,4)*(
      MassB + 2*MassWB) - 180*traceAdjYdTYdAdjYdYd + 90*MassB*
      traceAdjYdYdAdjYdYd - 540*traceAdjYeTYeAdjYeYe + 270*MassB*
      traceAdjYeYeAdjYeYe + 72*traceAdjYuTYuAdjYdYd - 72*MassB*
      traceAdjYuYuAdjYdYd + 72*traceTYdAdjYuYuAdjYd - 9*(MassWB*traceAdjYdYd +
      27*MassWB*traceAdjYeYe + MassB*(traceAdjYdYd + 27*traceAdjYeYe) -
      traceTYdAdjYd - 27*traceTYeAdjYe)*Sqr(g2) - 568*MassB*traceAdjYdYd*Sqr(g3
      ) - 568*MassG*traceAdjYdYd*Sqr(g3) + 568*traceTYdAdjYd*Sqr(g3)) + 5*Power
      (g1,4)*(5022*(2*MassB + MassWB)*Sqr(g2) + 5*(-2*MassB*(1505*traceAdjYdYd
      + 2619*traceAdjYeYe + 1404*traceAdjYuYu) + 1505*traceTYdAdjYd + 2619*
      traceTYeAdjYe + 1404*traceTYuAdjYu + 4752*(2*MassB + MassG)*Sqr(g3))) +
      125*(6210*Power(g2,6)*MassWB + 45*Power(g2,4)*(-2*MassWB*(21*traceAdjYdYd
      + 7*traceAdjYeYe + 12*traceAdjYuYu) + 21*traceTYdAdjYd + 7*traceTYeAdjYe
      + 12*traceTYuAdjYu + 48*(MassG + 2*MassWB)*Sqr(g3)) - 36*Sqr(g2)*(6*
      traceAdjYdTYdAdjYdYd - 3*MassWB*traceAdjYdYdAdjYdYd + 2*
      traceAdjYeTYeAdjYeYe - MassWB*traceAdjYeYeAdjYeYe + 6*
      traceAdjYuTYuAdjYdYd - 6*MassWB*traceAdjYuYuAdjYdYd + 6*
      traceTYdAdjYuYuAdjYd + 44*(MassG*traceAdjYdYd + MassWB*traceAdjYdYd -
      traceTYdAdjYd)*Sqr(g3)) - 4*(160*Power(g3,4)*(2*MassG*traceAdjYdYd -
      traceTYdAdjYd) + 9*(3*traceAdjYdYdAdjYdTYdAdjYdYd + 12*traceAdjYdYd*
      traceAdjYeTYeAdjYeYe + 4*traceAdjYeTYeAdjYeYe*traceAdjYeYe + 12*
      traceAdjYdTYdAdjYdYd*(3*traceAdjYdYd + traceAdjYeYe) +
      traceAdjYeYeAdjYeTYeAdjYeYe + 3*traceAdjYuTYuAdjYuYuAdjYdYd + 6*
      traceAdjYuTYuAdjYdYd*traceAdjYuYu + 3*traceAdjYuYuAdjYuTYuAdjYdYd + 18*
      traceAdjYdYdAdjYdYd*traceTYdAdjYd + 6*traceAdjYeYeAdjYeYe*traceTYdAdjYd +
      6*traceAdjYuYu*traceTYdAdjYuYuAdjYd + 3*traceTYdAdjYuYuAdjYuYuAdjYd + 6*
      traceAdjYdYdAdjYdYd*traceTYeAdjYe + 2*traceAdjYeYeAdjYeYe*traceTYeAdjYe +
      6*traceAdjYuYuAdjYdYd*traceTYuAdjYu) + 72*(6*traceAdjYdTYdAdjYdYd +
      traceAdjYuTYuAdjYdYd - MassG*(3*traceAdjYdYdAdjYdYd + traceAdjYuYuAdjYdYd
      ) + traceTYdAdjYuYuAdjYd)*Sqr(g3)))) + 0.016*threeLoop*(16119*Power(g1,6)
      *MassB + 5*Power(g1,4)*(MassB*(77*traceAdjYdYd - 81*traceAdjYeYe) + 729*(
      2*MassB + MassWB)*Sqr(g2) + 2376*(2*MassB + MassG)*Sqr(g3)) + 25*Sqr(g1)*
      (81*Power(g2,4)*(MassB + 2*MassWB) + 45*(MassB + MassWB)*traceAdjYdYd*Sqr
      (g2) - 2*(-54*traceAdjYdTYdAdjYdYd + 27*(MassB*traceAdjYdYdAdjYdYd + 2*
      traceAdjYeTYeAdjYeYe) + 56*(MassB + MassG)*traceAdjYdYd*Sqr(g3))) - 125*(
      945*Power(g2,6)*MassWB - 27*Power(g2,4)*(7*MassWB*traceAdjYdYd + 8*(MassG
      + 2*MassWB)*Sqr(g3)) - 18*Sqr(g2)*(6*traceAdjYdTYdAdjYdYd - 3*MassWB*
      traceAdjYdYdAdjYdYd + 2*traceAdjYeTYeAdjYeYe - 8*(MassG + MassWB)*
      traceAdjYdYd*Sqr(g3)) - 2*(16*Power(g3,4)*MassG*traceAdjYdYd + 27*
      traceAdjYdYdAdjYdTYdAdjYdYd - 72*(2*traceAdjYdTYdAdjYdYd - MassG*
      traceAdjYdYdAdjYdYd)*Sqr(g3))))*(Ye*1.2020569031595942))*UNITMATRIX(3))
      .real();
   const Eigen::Matrix<double,3,3> beta_TYe_2 = ((-0.04*threeLoop*(Power(
      g1,4)*(77*traceTYdAdjYd - 81*traceTYeAdjYe) - 25*(63*Power(g2,4)*(2*
      MassWB*traceAdjYeYe - 3*traceTYdAdjYd - traceTYeAdjYe) - 36*Sqr(g2)*(
      MassWB*traceAdjYeYeAdjYeYe - 8*traceTYdAdjYd*Sqr(g3)) + 4*(9*
      traceAdjYeYeAdjYeTYeAdjYeYe - 8*Sqr(g3)*(3*traceAdjYuTYuAdjYdYd - 3*MassG
      *traceAdjYuYuAdjYdYd + 3*traceTYdAdjYuYuAdjYd + traceTYdAdjYd*Sqr(g3))))
      + 10*Sqr(g1)*(9*(9*MassB*traceAdjYeYe + 9*MassWB*traceAdjYeYe + 5*
      traceTYdAdjYd - 9*traceTYeAdjYe)*Sqr(g2) - 2*(3*MassB*(9*
      traceAdjYeYeAdjYeYe - 7*traceAdjYuYuAdjYdYd) + 7*(3*traceAdjYuTYuAdjYdYd
      + 3*traceTYdAdjYuYuAdjYd + 8*traceTYdAdjYd*Sqr(g3)))))*(Ye*
      1.2020569031595942) - 0.004*threeLoop*(10746*Power(g1,6) + 5*Power(g1,4)*
      (77*traceAdjYdYd - 81*traceAdjYeYe + 1458*Sqr(g2) + 4752*Sqr(g3)) + 50*
      Sqr(g1)*(81*Power(g2,4) + 9*(5*traceAdjYdYd - 9*traceAdjYeYe)*Sqr(g2) - 2
      *(27*traceAdjYdYdAdjYdYd - 27*traceAdjYeYeAdjYeYe + 21*
      traceAdjYuYuAdjYdYd + 56*traceAdjYdYd*Sqr(g3))) - 125*(630*Power(g2,6) -
      9*Power(g2,4)*(7*(3*traceAdjYdYd + traceAdjYeYe) + 48*Sqr(g3)) + 36*Sqr(
      g2)*(3*traceAdjYdYdAdjYdYd + traceAdjYeYeAdjYeYe + 8*traceAdjYdYd*Sqr(g3)
      ) - 4*(8*Power(g3,4)*traceAdjYdYd - 3*(3*traceAdjYdYdAdjYdYdAdjYdYd +
      traceAdjYeYeAdjYeYeAdjYeYe) + 24*(3*traceAdjYdYdAdjYdYd +
      traceAdjYuYuAdjYdYd)*Sqr(g3))))*(TYe*1.2020569031595942) + 0.12*threeLoop
      *(1917*Power(g1,4)*MassB + 5*Sqr(g1)*(-98*MassB*traceAdjYdYd - 6*MassB*
      traceAdjYeYe + 98*traceTYdAdjYd + 6*traceTYeAdjYe + 117*(MassB + MassWB)*
      Sqr(g2)) + 25*(73*Power(g2,4)*MassWB - 10*(MassWB*(3*traceAdjYdYd +
      traceAdjYeYe) - 3*traceTYdAdjYd - traceTYeAdjYe)*Sqr(g2) + 4*(18*
      traceAdjYdTYdAdjYdYd + 6*traceAdjYeTYeAdjYeYe + 3*traceAdjYuTYuAdjYdYd -
      9*traceAdjYdYd*traceTYdAdjYd - 3*traceAdjYeYe*traceTYdAdjYd + 3*
      traceTYdAdjYuYuAdjYd - 3*traceAdjYdYd*traceTYeAdjYe - traceAdjYeYe*
      traceTYeAdjYe + 16*(MassG*traceAdjYdYd - traceTYdAdjYd)*Sqr(g3))))*(Ye*
      Ye.adjoint()*Ye) - 0.04*threeLoop*(2124*Power(g1,4) + 5*Sqr(g1)*(-187*
      traceAdjYdYd - 9*traceAdjYeYe + 216*Sqr(g2)) + 25*(66*Power(g2,4) - 21*(3
      *traceAdjYdYd + traceAdjYeYe)*Sqr(g2) + 4*(-18*traceAdjYdYdAdjYdYd + 6*
      traceAdjYdYd*traceAdjYeYe - 6*traceAdjYeYeAdjYeYe - 6*traceAdjYuYuAdjYdYd
      + 32*traceAdjYdYd*Sqr(g3) + 9*Sqr(traceAdjYdYd) + Sqr(traceAdjYeYe))))*(
      Ye*Ye.adjoint()*TYe) - 0.01*threeLoop*(8757*Power(g1,4) + 10*Sqr(g1)*(
      -508*traceAdjYdYd - 36*traceAdjYeYe + 621*Sqr(g2)) + 25*(393*Power(g2,4)
      - 96*(3*traceAdjYdYd + traceAdjYeYe)*Sqr(g2) + 20*(-18*
      traceAdjYdYdAdjYdYd + 6*traceAdjYdYd*traceAdjYeYe - 6*traceAdjYeYeAdjYeYe
      - 6*traceAdjYuYuAdjYdYd + 32*traceAdjYdYd*Sqr(g3) + 9*Sqr(traceAdjYdYd)
      + Sqr(traceAdjYeYe))))*(TYe*Ye.adjoint()*Ye) + 0.72*threeLoop*(81*Power(
      g1,4)*MassB + 225*Power(g2,4)*MassWB - 5*Sqr(g1)*(-2*MassB*traceAdjYdYd +
      27*(MassB + MassWB)*Sqr(g2)))*(Ye*Ye.adjoint()*Ye*1.2020569031595942))*
      UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYe_3 = ((0.0006666666666666666*
      threeLoop*(149958*Power(g1,6) + 5*Power(g1,4)*(-7525*traceAdjYdYd - 13095
      *traceAdjYeYe - 7020*traceAdjYuYu + 5022*Sqr(g2) + 23760*Sqr(g3)) + 50*
      Sqr(g1)*(135*Power(g2,4) + 90*traceAdjYdYdAdjYdYd + 270*
      traceAdjYeYeAdjYeYe - 72*traceAdjYuYuAdjYdYd - 9*(traceAdjYdYd + 27*
      traceAdjYeYe)*Sqr(g2) - 568*traceAdjYdYd*Sqr(g3)) + 125*(2070*Power(g2,6)
      + 45*Power(g2,4)*(-21*traceAdjYdYd - 7*traceAdjYeYe - 12*traceAdjYuYu +
      48*Sqr(g3)) - 36*Sqr(g2)*(-3*traceAdjYdYdAdjYdYd - traceAdjYeYeAdjYeYe -
      6*traceAdjYuYuAdjYdYd + 44*traceAdjYdYd*Sqr(g3)) + 4*(-160*Power(g3,4)*
      traceAdjYdYd + 3*(3*traceAdjYdYdAdjYdYdAdjYdYd + 18*traceAdjYdYdAdjYdYd*
      traceAdjYeYe + 6*traceAdjYeYe*traceAdjYeYeAdjYeYe + 18*traceAdjYdYd*(3*
      traceAdjYdYdAdjYdYd + traceAdjYeYeAdjYeYe) + traceAdjYeYeAdjYeYeAdjYeYe +
      18*traceAdjYuYu*traceAdjYuYuAdjYdYd + 9*traceAdjYuYuAdjYuYuAdjYdYd) + 72
      *(3*traceAdjYdYdAdjYdYd + traceAdjYuYuAdjYdYd)*Sqr(g3))))*TYe - 7.2*
      threeLoop*((3*MassB*traceAdjYeYe + traceTYdAdjYd - 3*traceTYeAdjYe)*Sqr(
      g1) - 5*(3*MassWB*traceAdjYdYd + MassWB*traceAdjYeYe - 3*traceTYdAdjYd -
      traceTYeAdjYe)*Sqr(g2) + 40*(MassG*traceAdjYdYd - traceTYdAdjYd)*Sqr(g3))
      *(Ye*Ye.adjoint()*Ye*1.2020569031595942) - 0.24*threeLoop*(54*Power(g1,4)
      - 5*Sqr(g1)*(-13*traceAdjYdYd + 9*traceAdjYeYe + 54*Sqr(g2)) + 25*(12*
      Power(g2,4) + 3*(3*traceAdjYdYd + traceAdjYeYe)*Sqr(g2) - 32*traceAdjYdYd
      *Sqr(g3)))*(Ye*Ye.adjoint()*TYe*1.2020569031595942) - 0.06*threeLoop*(513
      *Power(g1,4) - 10*Sqr(g1)*(8*traceAdjYdYd + 36*traceAdjYeYe + 135*Sqr(g2)
      ) + 25*(33*Power(g2,4) + 24*(3*traceAdjYdYd + traceAdjYeYe)*Sqr(g2) - 160
      *traceAdjYdYd*Sqr(g3)))*(TYe*Ye.adjoint()*Ye*1.2020569031595942) +
      threeLoop*(8*(3*traceTYdAdjYd + traceTYeAdjYe) - 21.6*MassB*Sqr(g1) - 12*
      MassWB*Sqr(g2))*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*Ye) + threeLoop*(19.8*
      Sqr(g1) + 3*(6*traceAdjYdYd + 2*traceAdjYeYe + Sqr(g2)))*(Ye*Ye.adjoint()
      *Ye*Ye.adjoint()*TYe) + threeLoop*(24*traceAdjYdYd + 8*traceAdjYeYe +
      21.6*Sqr(g1) + 12*Sqr(g2))*(Ye*Ye.adjoint()*TYe*Ye.adjoint()*Ye) +
      threeLoop*(18*traceAdjYdYd + 6*traceAdjYeYe + 12.6*Sqr(g1) + 15*Sqr(g2))*
      (TYe*Ye.adjoint()*Ye*Ye.adjoint()*Ye) + threeLoop*(-10.8*Sqr(g1) + 18*Sqr
      (g2))*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*TYe*1.2020569031595942) + 3.6*
      threeLoop*(3*Sqr(g1) - 5*Sqr(g2))*(TYe*Ye.adjoint()*Ye*Ye.adjoint()*Ye*
      1.2020569031595942) + 6*threeLoop*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*Ye*
      Ye.adjoint()*TYe) + 12*threeLoop*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*TYe*
      Ye.adjoint()*Ye) + 12*threeLoop*(Ye*Ye.adjoint()*TYe*Ye.adjoint()*Ye*
      Ye.adjoint()*Ye) + 12*threeLoop*(TYe*Ye.adjoint()*Ye*Ye.adjoint()*Ye*
      Ye.adjoint()*Ye) + 24*threeLoop*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*Ye*
      Ye.adjoint()*TYe*1.2020569031595942) + 36*threeLoop*(Ye*Ye.adjoint()*Ye*
      Ye.adjoint()*TYe*Ye.adjoint()*Ye*1.2020569031595942) + 36*threeLoop*(Ye*
      Ye.adjoint()*TYe*Ye.adjoint()*Ye*Ye.adjoint()*Ye*1.2020569031595942) + 30
      *threeLoop*(TYe*Ye.adjoint()*Ye*Ye.adjoint()*Ye*Ye.adjoint()*Ye*
      1.2020569031595942))*UNITMATRIX(3)).real();

   beta_TYe = beta_TYe_1 + beta_TYe_2 + beta_TYe_3;


   return beta_TYe;
}

} // namespace flexiblesusy
