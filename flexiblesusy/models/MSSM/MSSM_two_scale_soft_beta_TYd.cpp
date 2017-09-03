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

// File generated at Sat 2 Sep 2017 18:59:01

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
 * Calculates the one-loop beta function of TYd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MSSM_soft_parameters::calc_beta_TYd_one_loop(const Soft_traces& soft_traces) const
{
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (oneOver16PiSqr*(0.13333333333333333*Yd*(7*MassB*Sqr(g1) +
      5*(9*traceAdjYdTYd + 3*traceAdjYeTYe + 9*MassWB*Sqr(g2) + 16*MassG*Sqr(g3
      ))) + (3*traceYdAdjYd + traceYeAdjYe - 0.4666666666666667*Sqr(g1) - 3*Sqr
      (g2) - 5.333333333333333*Sqr(g3))*TYd + 4*(Yd*Yd.adjoint()*TYd) + 2*(Yd*
      Yu.adjoint()*TYu) + 5*(TYd*Yd.adjoint()*Yd) + TYd*Yu.adjoint()*Yu)).real(
      );


   return beta_TYd;
}

/**
 * Calculates the two-loop beta function of TYd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MSSM_soft_parameters::calc_beta_TYd_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceYdAdjYdTYdAdjYd = TRACE_STRUCT.traceYdAdjYdTYdAdjYd;
   const double traceYdAdjYuTYuAdjYd = TRACE_STRUCT.traceYdAdjYuTYuAdjYd;
   const double traceYeAdjYeTYeAdjYe = TRACE_STRUCT.traceYeAdjYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjYu = TRACE_STRUCT.traceYuAdjYdTYdAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_TYd;

   beta_TYd = (twoLoop*(-0.044444444444444446*Yd*(287*Power(g1,4)*MassB +
      5*(-32*Power(g3,4)*MassG + 135*Power(g2,4)*MassWB + 27*(6*
      traceYdAdjYdTYdAdjYd + traceYdAdjYuTYuAdjYd + 2*traceYeAdjYeTYeAdjYe +
      traceYuAdjYdTYdAdjYu) - 144*(traceAdjYdTYd - MassG*traceYdAdjYd)*Sqr(g3)
      + 72*(MassG + MassWB)*Sqr(g2)*Sqr(g3)) + Sqr(g1)*(45*(MassB + MassWB)*Sqr
      (g2) + 2*(9*(traceAdjYdTYd - 3*traceAdjYeTYe - MassB*traceYdAdjYd + 3*
      MassB*traceYeAdjYe) + 20*(MassB + MassG)*Sqr(g3)))) + (3.188888888888889*
      Power(g1,4) + 7.5*Power(g2,4) - 1.7777777777777777*Power(g3,4) - 9*
      traceYdAdjYdYdAdjYd - 3*traceYdAdjYuYuAdjYd - 3*traceYeAdjYeYeAdjYe + Sqr
      (g1)*(-0.4*traceYdAdjYd + 1.2*traceYeAdjYe + Sqr(g2) + 0.8888888888888888
      *Sqr(g3)) + 16*traceYdAdjYd*Sqr(g3) + 8*Sqr(g2)*Sqr(g3))*TYd - 0.4*(4*
      MassB*Sqr(g1) + 15*(3*traceAdjYdTYd + traceAdjYeTYe + 2*MassWB*Sqr(g2)))*
      (Yd*Yd.adjoint()*Yd) + (-4*(3*traceYdAdjYd + traceYeAdjYe) + 1.2*Sqr(g1)
      + 6*Sqr(g2))*(Yd*Yd.adjoint()*TYd) + (-6*traceAdjYuTYu - 1.6*MassB*Sqr(g1
      ))*(Yd*Yu.adjoint()*Yu) + (-6*traceYuAdjYu + 1.6*Sqr(g1))*(Yd*Yu.adjoint(
      )*TYu) + (-5*(3*traceYdAdjYd + traceYeAdjYe) + 1.2*Sqr(g1) + 12*Sqr(g2))*
      (TYd*Yd.adjoint()*Yd) + (-3*traceYuAdjYu + 0.8*Sqr(g1))*(TYd*Yu.adjoint()
      *Yu) - 6*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd) - 8*(Yd*Yd.adjoint()*TYd*
      Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*TYd) - 4*(Yd*
      Yu.adjoint()*Yu*Yu.adjoint()*TYu) - 4*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*
      Yd) - 4*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*Yu) - 6*(TYd*Yd.adjoint()*Yd*
      Yd.adjoint()*Yd) - 4*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*(TYd*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_TYd;
}

/**
 * Calculates the three-loop beta function of TYd.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MSSM_soft_parameters::calc_beta_TYd_three_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjYuYuAdjYuYu = TRACE_STRUCT.traceAdjYuYuAdjYuYu;
   const double traceAdjYuTYuAdjYdYd = TRACE_STRUCT.traceAdjYuTYuAdjYdYd;
   const double traceAdjYuTYuAdjYuYu = TRACE_STRUCT.traceAdjYuTYuAdjYuYu;
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


   Eigen::Matrix<double,3,3> beta_TYd;

   const Eigen::Matrix<double,3,3> beta_TYd_1 = ((-0.00044444444444444447
      *threeLoop*Yd*(389302*Power(g1,6)*MassB + 5*Power(g1,4)*(-15*(1890*MassB*
      traceAdjYdYd + 3798*MassB*traceAdjYeYe + 728*MassB*traceAdjYuYu - 945*
      traceTYdAdjYd - 1899*traceTYeAdjYe - 364*traceTYuAdjYu) + 1962*(2*MassB +
      MassWB)*Sqr(g2) + 15568*(2*MassB + MassG)*Sqr(g3)) + 50*Sqr(g1)*(765*
      Power(g2,4)*(MassB + 2*MassWB) - 9*Sqr(g2)*(3*(MassB*traceAdjYdYd +
      MassWB*traceAdjYdYd + 27*MassB*traceAdjYeYe + 27*MassWB*traceAdjYeYe -
      traceTYdAdjYd - 27*traceTYeAdjYe) + 16*(MassB + MassG + MassWB)*Sqr(g3))
      + 2*(1060*Power(g3,4)*(MassB + 2*MassG) + 27*(-10*traceAdjYdTYdAdjYdYd +
      5*MassB*traceAdjYdYdAdjYdYd - 30*traceAdjYeTYeAdjYeYe + 15*MassB*
      traceAdjYeYeAdjYeYe + 4*traceAdjYuTYuAdjYdYd - 4*MassB*
      traceAdjYuYuAdjYdYd + 4*traceTYdAdjYuYuAdjYd) - 852*(MassB*traceAdjYdYd +
      MassG*traceAdjYdYd - traceTYdAdjYd)*Sqr(g3))) + 125*(18630*Power(g2,6)*
      MassWB + 45*Power(g2,4)*(-6*MassWB*(21*traceAdjYdYd + 7*traceAdjYeYe + 12
      *traceAdjYuYu) + 63*traceTYdAdjYd + 21*traceTYeAdjYe + 36*traceTYuAdjYu +
      112*(MassG + 2*MassWB)*Sqr(g3)) + 36*Sqr(g2)*(68*Power(g3,4)*(2*MassG +
      MassWB) + 3*(-6*traceAdjYdTYdAdjYdYd + MassWB*(3*traceAdjYdYdAdjYdYd +
      traceAdjYeYeAdjYeYe + 6*traceAdjYuYuAdjYdYd) - 2*(traceAdjYeTYeAdjYeYe +
      3*(traceAdjYuTYuAdjYdYd + traceTYdAdjYuYuAdjYd))) - 132*(MassG*
      traceAdjYdYd + MassWB*traceAdjYdYd - traceTYdAdjYd)*Sqr(g3)) + 4*(5440*
      Power(g3,6)*MassG - 480*Power(g3,4)*(2*MassG*(2*traceAdjYdYd +
      traceAdjYuYu) - 2*traceTYdAdjYd - traceTYuAdjYu) - 27*(3*
      traceAdjYdYdAdjYdTYdAdjYdYd + 12*traceAdjYdYd*traceAdjYeTYeAdjYeYe + 4*
      traceAdjYeTYeAdjYeYe*traceAdjYeYe + 12*traceAdjYdTYdAdjYdYd*(3*
      traceAdjYdYd + traceAdjYeYe) + traceAdjYeYeAdjYeTYeAdjYeYe + 3*
      traceAdjYuTYuAdjYuYuAdjYdYd + 6*traceAdjYuTYuAdjYdYd*traceAdjYuYu + 3*
      traceAdjYuYuAdjYuTYuAdjYdYd + 18*traceAdjYdYdAdjYdYd*traceTYdAdjYd + 6*
      traceAdjYeYeAdjYeYe*traceTYdAdjYd + 6*traceAdjYuYu*traceTYdAdjYuYuAdjYd +
      3*traceTYdAdjYuYuAdjYuYuAdjYd + 6*traceAdjYdYdAdjYdYd*traceTYeAdjYe + 2*
      traceAdjYeYeAdjYeYe*traceTYeAdjYe + 6*traceAdjYuYuAdjYdYd*traceTYuAdjYu)
      - 216*(6*traceAdjYdTYdAdjYdYd + traceAdjYuTYuAdjYdYd - MassG*(3*
      traceAdjYdYdAdjYdYd + traceAdjYuYuAdjYdYd) + traceTYdAdjYuYuAdjYd)*Sqr(g3
      )))) + 0.016*threeLoop*(4179*Power(g1,6)*MassB + 25*(81*Power(g2,4)*(
      MassB + 2*MassWB) + 4*(44*Power(g3,4)*(MassB + 2*MassG) + 27*
      traceAdjYdTYdAdjYdYd))*Sqr(g1) - 375*(640*Power(g3,6)*MassG + 315*Power(
      g2,6)*MassWB - 12*(4*Power(g3,4)*(2*MassG + MassWB) + 3*
      traceAdjYdTYdAdjYdYd)*Sqr(g2) - 72*Power(g2,4)*(MassG + 2*MassWB)*Sqr(g3)
      + 96*traceAdjYdTYdAdjYdYd*Sqr(g3)) + 35*Power(g1,4)*(27*(2*MassB +
      MassWB)*Sqr(g2) + 11*(MassB*traceAdjYdYd + 8*(2*MassB + MassG)*Sqr(g3))))
      *(Yd*1.2020569031595942))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_2 = ((-0.04*threeLoop*(Power(
      g1,4)*(162*MassB*traceAdjYeYe + 77*traceTYdAdjYd - 81*traceTYeAdjYe) - 10
      *Sqr(g1)*(9*(5*MassB*traceAdjYdYd + 5*MassWB*traceAdjYdYd - 9*MassB*
      traceAdjYeYe - 9*MassWB*traceAdjYeYe - 5*traceTYdAdjYd + 9*traceTYeAdjYe)
      *Sqr(g2) - 2*(3*(18*traceAdjYeTYeAdjYeYe + MassB*(9*traceAdjYdYdAdjYdYd -
      9*traceAdjYeYeAdjYeYe + 7*traceAdjYuYuAdjYdYd) - 7*(traceAdjYuTYuAdjYdYd
      + traceTYdAdjYuYuAdjYd)) + 56*(MassB*traceAdjYdYd + MassG*traceAdjYdYd -
      traceTYdAdjYd)*Sqr(g3))) - 25*(63*Power(g2,4)*(2*MassWB*(3*traceAdjYdYd
      + traceAdjYeYe) - 3*traceTYdAdjYd - traceTYeAdjYe) - 36*Sqr(g2)*(3*MassWB
      *traceAdjYdYdAdjYdYd - 2*traceAdjYeTYeAdjYeYe + MassWB*
      traceAdjYeYeAdjYeYe + 8*(MassG*traceAdjYdYd + MassWB*traceAdjYdYd -
      traceTYdAdjYd)*Sqr(g3)) + 4*(9*(3*traceAdjYdYdAdjYdTYdAdjYdYd +
      traceAdjYeYeAdjYeTYeAdjYeYe) + 8*Power(g3,4)*(2*MassG*traceAdjYdYd -
      traceTYdAdjYd) - 24*(traceAdjYuTYuAdjYdYd - MassG*(3*traceAdjYdYdAdjYdYd
      + traceAdjYuYuAdjYdYd) + traceTYdAdjYuYuAdjYd)*Sqr(g3))))*(Yd*
      1.2020569031595942) - 0.004*threeLoop*(2786*Power(g1,6) + 5*Power(g1,4)*(
      77*traceAdjYdYd - 81*traceAdjYeYe + 378*Sqr(g2) + 1232*Sqr(g3)) + 50*Sqr(
      g1)*(81*Power(g2,4) + 9*(5*traceAdjYdYd - 9*traceAdjYeYe)*Sqr(g2) + 2*(88
      *Power(g3,4) - 27*traceAdjYdYdAdjYdYd + 27*traceAdjYeYeAdjYeYe - 21*
      traceAdjYuYuAdjYdYd - 56*traceAdjYdYd*Sqr(g3))) - 125*(630*Power(g2,6) -
      9*Power(g2,4)*(7*(3*traceAdjYdYd + traceAdjYeYe) + 48*Sqr(g3)) - 36*Sqr(
      g2)*(8*Power(g3,4) - 3*traceAdjYdYdAdjYdYd - traceAdjYeYeAdjYeYe - 8*
      traceAdjYdYd*Sqr(g3)) + 4*(320*Power(g3,6) - 8*Power(g3,4)*traceAdjYdYd +
      3*(3*traceAdjYdYdAdjYdYdAdjYdYd + traceAdjYeYeAdjYeYeAdjYeYe) - 24*(3*
      traceAdjYdYdAdjYdYd + traceAdjYuYuAdjYdYd)*Sqr(g3))))*(TYd*
      1.2020569031595942) + 0.013333333333333334*threeLoop*(5269*Power(g1,4)*
      MassB + 5*Sqr(g1)*(381*(MassB + MassWB)*Sqr(g2) + 2*(-153*MassB*
      traceAdjYdYd + 69*MassB*traceAdjYeYe + 153*traceTYdAdjYd - 69*
      traceTYeAdjYe + 148*(MassB + MassG)*Sqr(g3))) + 75*(219*Power(g2,4)*
      MassWB + 2*Sqr(g2)*(-15*(3*MassWB*traceAdjYdYd + MassWB*traceAdjYeYe - 3*
      traceTYdAdjYd - traceTYeAdjYe) + 92*(MassG + MassWB)*Sqr(g3)) - 4*(8*
      Power(g3,4)*MassG + 3*(-18*traceAdjYdTYdAdjYdYd - 6*traceAdjYeTYeAdjYeYe
      - 3*traceAdjYuTYuAdjYdYd + 9*traceAdjYdYd*traceTYdAdjYd + 3*traceAdjYeYe*
      traceTYdAdjYd - 3*traceTYdAdjYuYuAdjYd + 3*traceAdjYdYd*traceTYeAdjYe +
      traceAdjYeYe*traceTYeAdjYe) - 12*(MassG*traceAdjYdYd - MassG*traceAdjYeYe
      - traceTYdAdjYd + traceTYeAdjYe)*Sqr(g3))))*(Yd*Yd.adjoint()*Yd) -
      0.013333333333333334*threeLoop*(1792*Power(g1,4) + 5*Sqr(g1)*(-201*
      traceAdjYdYd + 252*Sqr(g2) + 224*Sqr(g3)) + 25*(198*Power(g2,4) - 32*
      Power(g3,4) + 9*Sqr(g2)*(-21*traceAdjYdYd + 32*Sqr(g3))))*(Yd*Yd.adjoint(
      )*TYd))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_3 = ((-0.2*threeLoop*(-360*
      traceAdjYdYdAdjYdYd + 120*traceAdjYdYd*traceAdjYeYe - 120*
      traceAdjYeYeAdjYeYe - 120*traceAdjYuYuAdjYdYd + 31*traceAdjYeYe*Sqr(g1) -
      105*traceAdjYeYe*Sqr(g2) + 160*(traceAdjYdYd - traceAdjYeYe)*Sqr(g3) +
      180*Sqr(traceAdjYdYd) + 20*Sqr(traceAdjYeYe))*(Yd*Yd.adjoint()*TYd) +
      0.013333333333333334*threeLoop*(3767*Power(g1,4)*MassB + 15*Sqr(g1)*(59*(
      MassB + MassWB)*Sqr(g2) + 4*(-5*MassB*traceAdjYuYu + 5*traceTYuAdjYu + 34
      *(MassB + MassG)*Sqr(g3))) + 25*(135*Power(g2,4)*MassWB + 12*Sqr(g2)*(-9*
      MassWB*traceAdjYuYu + 9*traceTYuAdjYu + 2*(MassG + MassWB)*Sqr(g3)) - 4*(
      8*Power(g3,4)*MassG - 9*(traceAdjYuTYuAdjYdYd + 6*traceAdjYuTYuAdjYuYu +
      traceTYdAdjYuYuAdjYd - 3*traceAdjYuYu*traceTYuAdjYu) - 12*(MassG*
      traceAdjYuYu - traceTYuAdjYu)*Sqr(g3))))*(Yd*Yu.adjoint()*Yu) -
      0.006666666666666667*threeLoop*(3767*Power(g1,4) + 30*Sqr(g1)*(-20*
      traceAdjYuYu + 59*Sqr(g2) + 136*Sqr(g3)) + 25*(135*Power(g2,4) + 24*Sqr(
      g2)*(-9*traceAdjYuYu + 2*Sqr(g3)) - 4*(8*Power(g3,4) + 18*
      traceAdjYuYuAdjYdYd + 54*traceAdjYuYuAdjYuYu - 24*traceAdjYuYu*Sqr(g3) -
      27*Sqr(traceAdjYuYu))))*(Yd*Yu.adjoint()*TYu) + 0.0033333333333333335*
      threeLoop*(-8639*Power(g1,4) - 29475*Power(g2,4) + 3600*Sqr(g2)*(6*
      traceAdjYdYd + 2*traceAdjYeYe - 15*Sqr(g3)) - 10*Sqr(g1)*(-516*
      traceAdjYdYd + 228*traceAdjYeYe + 639*Sqr(g2) + 440*Sqr(g3)) + 500*(8*
      Power(g3,4) - 24*(traceAdjYdYd - traceAdjYeYe)*Sqr(g3) - 3*(-18*
      traceAdjYdYdAdjYdYd + 6*traceAdjYdYd*traceAdjYeYe - 6*traceAdjYeYeAdjYeYe
      - 6*traceAdjYuYuAdjYdYd + 9*Sqr(traceAdjYdYd) + Sqr(traceAdjYeYe))))*(
      TYd*Yd.adjoint()*Yd) + (-12.556666666666667*Power(g1,4)*threeLoop - 11.25
      *Power(g2,4)*threeLoop + 2.6666666666666665*Power(g3,4)*threeLoop + 6*
      threeLoop*traceAdjYuYuAdjYdYd + 18*threeLoop*traceAdjYuYuAdjYuYu - 8*
      threeLoop*traceAdjYuYu*Sqr(g3) - 0.1*threeLoop*Sqr(g1)*(-20*traceAdjYuYu
      + 59*Sqr(g2) + 136*Sqr(g3)) + Sqr(g2)*(18*threeLoop*traceAdjYuYu - 4*
      threeLoop*Sqr(g3)) - 9*threeLoop*Sqr(traceAdjYuYu))*(TYd*Yu.adjoint()*Yu)
      - 0.08*threeLoop*(7*Power(g1,4)*MassB + 15*Sqr(g1)*(17*(MassB + MassWB)*
      Sqr(g2) + 2*(-9*MassB*traceAdjYdYd + 7*MassB*traceAdjYeYe + 9*
      traceTYdAdjYd - 7*traceTYeAdjYe + 16*(MassB + MassG)*Sqr(g3))) - 25*(81*
      Power(g2,4)*MassWB + 16*Sqr(g3)*(-9*MassG*traceAdjYdYd + 9*traceTYdAdjYd
      + 34*MassG*Sqr(g3)) - 6*Sqr(g2)*(-9*MassWB*traceAdjYdYd - 3*MassWB*
      traceAdjYeYe + 9*traceTYdAdjYd + 3*traceTYeAdjYe + 16*(MassG + MassWB)*
      Sqr(g3))))*(Yd*Yd.adjoint()*Yd*1.2020569031595942) + 0.02666666666666667*
      threeLoop*(14*Power(g1,4) + 5*Sqr(g1)*(-117*traceAdjYdYd + 81*
      traceAdjYeYe + 90*Sqr(g2) + 208*Sqr(g3)) - 25*(108*Power(g2,4) - 9*Sqr(g2
      )*(-9*traceAdjYdYd + 16*Sqr(g3)) + 32*Sqr(g3)*(-9*traceAdjYdYd + 17*Sqr(
      g3))))*(Yd*Yd.adjoint()*TYd*1.2020569031595942))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_TYd_4 = ((-18*threeLoop*
      traceAdjYeYe*Sqr(g2)*(Yd*Yd.adjoint()*TYd*1.2020569031595942) -
      0.02666666666666667*threeLoop*(143*Power(g1,4)*MassB - 25*(544*Power(g3,4
      )*MassG + 189*Power(g2,4)*MassWB - 144*(MassG*traceAdjYuYu -
      traceTYuAdjYu)*Sqr(g3)) + 5*Sqr(g1)*(135*(MassB + MassWB)*Sqr(g2) + 8*(-9
      *MassB*traceAdjYuYu + 9*traceTYuAdjYu + 16*(MassB + MassG)*Sqr(g3))))*(Yd
      *Yu.adjoint()*Yu*1.2020569031595942) + 0.013333333333333334*threeLoop*(
      143*Power(g1,4) + 10*Sqr(g1)*(-72*traceAdjYuYu + 135*Sqr(g2) + 128*Sqr(g3
      )) - 25*(189*Power(g2,4) + 544*Power(g3,4) - 288*traceAdjYuYu*Sqr(g3)))*(
      Yd*Yu.adjoint()*TYu*1.2020569031595942) + 0.006666666666666667*threeLoop*
      (7*Power(g1,4) + 10*Sqr(g1)*(-252*traceAdjYdYd + 216*traceAdjYeYe + 279*
      Sqr(g2) + 448*Sqr(g3)) - 25*(297*Power(g2,4) - 72*Sqr(g2)*(-9*
      traceAdjYdYd - 3*traceAdjYeYe + 16*Sqr(g3)) + 160*(17*Power(g3,4) - 9*
      traceAdjYdYd*Sqr(g3))))*(TYd*Yd.adjoint()*Yd*1.2020569031595942) +
      0.006666666666666667*threeLoop*(143*Power(g1,4) + 10*Sqr(g1)*(-72*
      traceAdjYuYu + 135*Sqr(g2) + 128*Sqr(g3)) - 25*(189*Power(g2,4) + 544*
      Power(g3,4) - 288*traceAdjYuYu*Sqr(g3)))*(TYd*Yu.adjoint()*Yu*
      1.2020569031595942) - 0.26666666666666666*threeLoop*(MassB*Sqr(g1) + 5*(
      -18*traceTYdAdjYd - 6*traceTYeAdjYe + 9*MassWB*Sqr(g2) + 64*MassG*Sqr(g3)
      ))*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) + (0.6*threeLoop*Sqr(g1) +
      threeLoop*(18*traceAdjYdYd + 6*traceAdjYeYe + 3*Sqr(g2) + 64*Sqr(g3)))*(
      Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd) + 0.26666666666666666*threeLoop*(Sqr
      (g1) + 5*(18*traceAdjYdYd + 6*traceAdjYeYe + 9*Sqr(g2) + 64*Sqr(g3)))*(Yd
      *Yd.adjoint()*TYd*Yd.adjoint()*Yd) + 0.13333333333333333*threeLoop*(29*
      MassB*Sqr(g1) - 5*(18*traceTYdAdjYd + 6*traceTYeAdjYe - 36*traceTYuAdjYu
      + 27*MassWB*Sqr(g2) + 64*MassG*Sqr(g3)))*(Yd*Yu.adjoint()*Yu*Yd.adjoint()
      *Yd) + (-6*threeLoop*traceAdjYdYd - 2*threeLoop*traceAdjYeYe + 12*
      threeLoop*traceAdjYuYu - 1.9333333333333333*threeLoop*Sqr(g1) + 9*
      threeLoop*Sqr(g2) + 21.333333333333332*threeLoop*Sqr(g3))*(Yd*Yu.adjoint(
      )*Yu*Yd.adjoint()*TYd) - 0.6666666666666666*threeLoop*(11*MassB*Sqr(g1) -
      9*(2*traceTYuAdjYu + MassWB*Sqr(g2)) + 64*MassG*Sqr(g3))*(Yd*Yu.adjoint(
      )*Yu*Yu.adjoint()*Yu) + 0.6666666666666666*threeLoop*(18*traceAdjYuYu +
      11*Sqr(g1) - 9*Sqr(g2) + 64*Sqr(g3))*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu
      ) - 0.13333333333333333*threeLoop*(29*Sqr(g1) - 5*(-6*(3*traceAdjYdYd +
      traceAdjYeYe - 6*traceAdjYuYu) + 27*Sqr(g2) + 64*Sqr(g3)))*(Yd*Yu.adjoint
      ()*TYu*Yd.adjoint()*Yd) + 0.6666666666666666*threeLoop*(18*traceAdjYuYu +
      11*Sqr(g1) - 9*Sqr(g2) + 64*Sqr(g3))*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*
      Yu) + (-0.2*threeLoop*Sqr(g1) + threeLoop*(18*traceAdjYdYd + 6*
      traceAdjYeYe + 15*Sqr(g2) + 64*Sqr(g3)))*(TYd*Yd.adjoint()*Yd*Yd.adjoint(
      )*Yd) - 0.13333333333333333*threeLoop*(29*Sqr(g1) - 5*(-6*(3*traceAdjYdYd
      + traceAdjYeYe - 6*traceAdjYuYu) + 27*Sqr(g2) + 64*Sqr(g3)))*(TYd*
      Yu.adjoint()*Yu*Yd.adjoint()*Yd) + 0.3333333333333333*threeLoop*(18*
      traceAdjYuYu + 11*Sqr(g1) - 9*Sqr(g2) + 64*Sqr(g3))*(TYd*Yu.adjoint()*Yu*
      Yu.adjoint()*Yu) - 1.2*threeLoop*(Sqr(g1) - 15*Sqr(g2))*(Yd*Yd.adjoint()*
      Yd*Yd.adjoint()*TYd*1.2020569031595942) - 7.2*MassB*threeLoop*Sqr(g1)*(Yd
      *Yu.adjoint()*Yu*Yd.adjoint()*Yd*1.2020569031595942))*UNITMATRIX(3)).real
      ();
   const Eigen::Matrix<double,3,3> beta_TYd_5 = ((0.00007407407407407407*
      threeLoop*(389302*Power(g1,6) + 15*Power(g1,4)*(-15*(945*traceAdjYdYd +
      1899*traceAdjYeYe + 364*traceAdjYuYu) + 1962*Sqr(g2) + 15568*Sqr(g3)) +
      150*Sqr(g1)*(765*Power(g2,4) + 2120*Power(g3,4) + 54*(5*
      traceAdjYdYdAdjYdYd + 15*traceAdjYeYeAdjYeYe - 4*traceAdjYuYuAdjYdYd) -
      1704*traceAdjYdYd*Sqr(g3) - 9*Sqr(g2)*(3*traceAdjYdYd + 81*traceAdjYeYe +
      16*Sqr(g3))) + 125*(18630*Power(g2,6) + 135*Power(g2,4)*(-3*(21*
      traceAdjYdYd + 7*traceAdjYeYe + 12*traceAdjYuYu) + 112*Sqr(g3)) + 108*Sqr
      (g2)*(68*Power(g3,4) + 3*(3*traceAdjYdYdAdjYdYd + traceAdjYeYeAdjYeYe + 6
      *traceAdjYuYuAdjYdYd) - 132*traceAdjYdYd*Sqr(g3)) + 4*(5440*Power(g3,6) -
      1440*Power(g3,4)*(2*traceAdjYdYd + traceAdjYuYu) + 27*(3*
      traceAdjYdYdAdjYdYdAdjYdYd + 18*traceAdjYdYdAdjYdYd*traceAdjYeYe + 6*
      traceAdjYeYe*traceAdjYeYeAdjYeYe + 18*traceAdjYdYd*(3*traceAdjYdYdAdjYdYd
      + traceAdjYeYeAdjYeYe) + traceAdjYeYeAdjYeYeAdjYeYe + 18*traceAdjYuYu*
      traceAdjYuYuAdjYdYd + 9*traceAdjYuYuAdjYuYuAdjYdYd) + 648*(3*
      traceAdjYdYdAdjYdYd + traceAdjYuYuAdjYdYd)*Sqr(g3))))*TYd + 36*MassWB*
      threeLoop*Sqr(g2)*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd*1.2020569031595942)
      + 3.6*threeLoop*(Sqr(g1) - 5*Sqr(g2))*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*
      TYd*1.2020569031595942) + 12*threeLoop*(MassB*Sqr(g1) - 3*MassWB*Sqr(g2))
      *(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu*1.2020569031595942) - 12*threeLoop*(
      Sqr(g1) - 3*Sqr(g2))*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu*
      1.2020569031595942) + 7.2*threeLoop*(Sqr(g1) - 5*Sqr(g2))*(Yd*Yu.adjoint(
      )*TYu*Yd.adjoint()*Yd*1.2020569031595942) - 12*threeLoop*(Sqr(g1) - 3*Sqr
      (g2))*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*Yu*1.2020569031595942) + 1.2*
      threeLoop*(Sqr(g1) - 15*Sqr(g2))*(TYd*Yd.adjoint()*Yd*Yd.adjoint()*Yd*
      1.2020569031595942) + 7.2*threeLoop*(Sqr(g1) - 5*Sqr(g2))*(TYd*Yu.adjoint
      ()*Yu*Yd.adjoint()*Yd*1.2020569031595942) - 6*threeLoop*(Sqr(g1) - 3*Sqr(
      g2))*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu*1.2020569031595942) + 6*
      threeLoop*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd) + 12*
      threeLoop*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd) + 6*
      threeLoop*(Yd*Yd.adjoint()*Yd*Yu.adjoint()*Yu*Yd.adjoint()*TYd) + 4*
      threeLoop*(Yd*Yd.adjoint()*Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd) + 12*
      threeLoop*(Yd*Yd.adjoint()*TYd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) + 4*
      threeLoop*(Yd*Yd.adjoint()*TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*
      threeLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd*Yd.adjoint()*TYd) + 8*
      threeLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd*Yu.adjoint()*TYu) - 4*
      threeLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*TYd*Yd.adjoint()*Yd) + 8*
      threeLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*TYd*Yu.adjoint()*Yu) + 6*
      threeLoop*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu*Yd.adjoint()*TYd) + 12*
      threeLoop*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu*Yd.adjoint()*Yd) - 4*
      threeLoop*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) + 8*
      threeLoop*(Yd*Yu.adjoint()*TYu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) + 12*
      threeLoop*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*Yu*Yd.adjoint()*Yd) + 12*
      threeLoop*(TYd*Yd.adjoint()*Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*
      threeLoop*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) + 4*
      threeLoop*(TYd*Yu.adjoint()*Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) + 12*
      threeLoop*(TYd*Yu.adjoint()*Yu*Yu.adjoint()*Yu*Yd.adjoint()*Yd) + 24*
      threeLoop*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd*
      1.2020569031595942) + 36*threeLoop*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*TYd*
      Yd.adjoint()*Yd*1.2020569031595942) + 36*threeLoop*(Yd*Yd.adjoint()*TYd*
      Yd.adjoint()*Yd*Yd.adjoint()*Yd*1.2020569031595942) + 12*threeLoop*(Yd*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu*Yu.adjoint()*TYu*1.2020569031595942) + 12
      *threeLoop*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*TYu*Yu.adjoint()*Yu*
      1.2020569031595942) + 12*threeLoop*(Yd*Yu.adjoint()*TYu*Yu.adjoint()*Yu*
      Yu.adjoint()*Yu*1.2020569031595942) + 30*threeLoop*(TYd*Yd.adjoint()*Yd*
      Yd.adjoint()*Yd*Yd.adjoint()*Yd*1.2020569031595942) + 6*threeLoop*(TYd*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu*1.2020569031595942))*
      UNITMATRIX(3)).real();

   beta_TYd = beta_TYd_1 + beta_TYd_2 + beta_TYd_3 + beta_TYd_4 +
      beta_TYd_5;


   return beta_TYd;
}

} // namespace flexiblesusy
