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

// File generated at Tue 7 Nov 2017 11:38:51

#include "MSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MSSM_susy_parameters::calc_beta_Yd_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (oneOver16PiSqr*(Yd*(3*traceYdAdjYd + traceYeAdjYe -
      0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3)) + 3*(
      Yd*Yd.adjoint()*Yd) + Yd*Yu.adjoint()*Yu)).real();


   return beta_Yd;
}

/**
 * Calculates the two-loop beta function of Yd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MSSM_susy_parameters::calc_beta_Yd_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (twoLoop*(0.011111111111111112*Yd*(287*Power(g1,4) + 2*Sqr(
      g1)*(-18*traceYdAdjYd + 54*traceYeAdjYe + 45*Sqr(g2) + 40*Sqr(g3)) + 5*(
      135*Power(g2,4) + 144*Sqr(g2)*Sqr(g3) - 2*(16*Power(g3,4) + 27*(3*
      traceYdAdjYdYdAdjYd + traceYdAdjYuYuAdjYd + traceYeAdjYeYeAdjYe) - 144*
      traceYdAdjYd*Sqr(g3)))) + (-3*(3*traceYdAdjYd + traceYeAdjYe) + 0.8*Sqr(
      g1) + 6*Sqr(g2))*(Yd*Yd.adjoint()*Yd) + (-3*traceYuAdjYu + 0.8*Sqr(g1))*(
      Yd*Yu.adjoint()*Yu) - 4*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yd*
      Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu)
      )).real();


   return beta_Yd;
}

/**
 * Calculates the three-loop beta function of Yd.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MSSM_susy_parameters::calc_beta_Yd_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   const double traceAdjYdYd = TRACE_STRUCT.traceAdjYdYd;
   const double traceAdjYeYe = TRACE_STRUCT.traceAdjYeYe;
   const double traceAdjYuYu = TRACE_STRUCT.traceAdjYuYu;
   const double traceAdjYdYdAdjYdYd = TRACE_STRUCT.traceAdjYdYdAdjYdYd;
   const double traceAdjYeYeAdjYeYe = TRACE_STRUCT.traceAdjYeYeAdjYeYe;
   const double traceAdjYuYuAdjYdYd = TRACE_STRUCT.traceAdjYuYuAdjYdYd;
   const double traceAdjYuYuAdjYuYu = TRACE_STRUCT.traceAdjYuYuAdjYuYu;
   const double traceAdjYdYdAdjYdYdAdjYdYd =
      TRACE_STRUCT.traceAdjYdYdAdjYdYdAdjYdYd;
   const double traceAdjYeYeAdjYeYeAdjYeYe =
      TRACE_STRUCT.traceAdjYeYeAdjYeYeAdjYeYe;
   const double traceAdjYuYuAdjYuYuAdjYdYd =
      TRACE_STRUCT.traceAdjYuYuAdjYuYuAdjYdYd;


   Eigen::Matrix<double,3,3> beta_Yd;

   const Eigen::Matrix<double,3,3> beta_Yd_1 = ((0.00007407407407407407*
      threeLoop*Yd*(389302*Power(g1,6) + 15*Power(g1,4)*(-15*(945*traceAdjYdYd
      + 1899*traceAdjYeYe + 364*traceAdjYuYu) + 1962*Sqr(g2) + 15568*Sqr(g3)) +
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
      traceAdjYdYdAdjYdYd + traceAdjYuYuAdjYdYd)*Sqr(g3)))) - 0.004*threeLoop*(
      2786*Power(g1,6) + 5*Power(g1,4)*(77*traceAdjYdYd - 81*traceAdjYeYe + 378
      *Sqr(g2) + 1232*Sqr(g3)) + 50*Sqr(g1)*(81*Power(g2,4) + 9*(5*traceAdjYdYd
      - 9*traceAdjYeYe)*Sqr(g2) + 2*(88*Power(g3,4) - 27*traceAdjYdYdAdjYdYd +
      27*traceAdjYeYeAdjYeYe - 21*traceAdjYuYuAdjYdYd - 56*traceAdjYdYd*Sqr(g3
      ))) - 125*(630*Power(g2,6) - 9*Power(g2,4)*(7*(3*traceAdjYdYd +
      traceAdjYeYe) + 48*Sqr(g3)) - 36*Sqr(g2)*(8*Power(g3,4) - 3*
      traceAdjYdYdAdjYdYd - traceAdjYeYeAdjYeYe - 8*traceAdjYdYd*Sqr(g3)) + 4*(
      320*Power(g3,6) - 8*Power(g3,4)*traceAdjYdYd + 3*(3*
      traceAdjYdYdAdjYdYdAdjYdYd + traceAdjYeYeAdjYeYeAdjYeYe) - 24*(3*
      traceAdjYdYdAdjYdYd + traceAdjYuYuAdjYdYd)*Sqr(g3))))*(Yd*
      1.2020569031595942) - 0.0033333333333333335*threeLoop*(5269*Power(g1,4) +
      10*Sqr(g1)*(-306*traceAdjYdYd + 138*traceAdjYeYe + 381*Sqr(g2) + 296*Sqr
      (g3)) + 75*(219*Power(g2,4) + 4*Sqr(g2)*(-15*(3*traceAdjYdYd +
      traceAdjYeYe) + 92*Sqr(g3)) - 4*(8*Power(g3,4) - 24*(traceAdjYdYd -
      traceAdjYeYe)*Sqr(g3) - 3*(-18*traceAdjYdYdAdjYdYd + 6*traceAdjYdYd*
      traceAdjYeYe - 6*traceAdjYeYeAdjYeYe - 6*traceAdjYuYuAdjYdYd + 9*Sqr(
      traceAdjYdYd) + Sqr(traceAdjYeYe)))))*(Yd*Yd.adjoint()*Yd) + (
      -12.556666666666667*Power(g1,4)*threeLoop - 11.25*Power(g2,4)*threeLoop +
      2.6666666666666665*Power(g3,4)*threeLoop + 6*threeLoop*
      traceAdjYuYuAdjYdYd + 18*threeLoop*traceAdjYuYuAdjYuYu - 8*threeLoop*
      traceAdjYuYu*Sqr(g3) - 0.1*threeLoop*Sqr(g1)*(-20*traceAdjYuYu + 59*Sqr(
      g2) + 136*Sqr(g3)) + Sqr(g2)*(18*threeLoop*traceAdjYuYu - 4*threeLoop*Sqr
      (g3)) - 9*threeLoop*Sqr(traceAdjYuYu))*(Yd*Yu.adjoint()*Yu) + 0.02*
      threeLoop*(7*Power(g1,4) + 30*Sqr(g1)*(17*Sqr(g2) + 32*Sqr(g3)) + 75*Sqr(
      g2)*(-27*Sqr(g2) + 64*Sqr(g3)))*(Yd*Yd.adjoint()*Yd*1.2020569031595942))*
      UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_Yd_2 = ((-0.4*threeLoop*(680*
      Power(g3,4) + 3*((9*traceAdjYdYd - 7*traceAdjYeYe)*Sqr(g1) + 15*(3*
      traceAdjYdYd + traceAdjYeYe)*Sqr(g2)) - 360*traceAdjYdYd*Sqr(g3))*(Yd*
      Yd.adjoint()*Yd*1.2020569031595942) + 0.006666666666666667*threeLoop*(143
      *Power(g1,4) + 10*Sqr(g1)*(-72*traceAdjYuYu + 135*Sqr(g2) + 128*Sqr(g3))
      - 25*(189*Power(g2,4) + 544*Power(g3,4) - 288*traceAdjYuYu*Sqr(g3)))*(Yd*
      Yu.adjoint()*Yu*1.2020569031595942) + 0.13333333333333333*threeLoop*(Sqr(
      g1) + 5*(18*traceAdjYdYd + 6*traceAdjYeYe + 9*Sqr(g2) + 64*Sqr(g3)))*(Yd*
      Yd.adjoint()*Yd*Yd.adjoint()*Yd) + 0.06666666666666667*threeLoop*(-29*Sqr
      (g1) + 5*(-6*(3*traceAdjYdYd + traceAdjYeYe - 6*traceAdjYuYu) + 27*Sqr(g2
      ) + 64*Sqr(g3)))*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) +
      0.3333333333333333*threeLoop*(18*traceAdjYuYu + 11*Sqr(g1) - 9*Sqr(g2) +
      64*Sqr(g3))*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu) + 3.6*threeLoop*(Sqr(g1)
      - 5*Sqr(g2))*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd*1.2020569031595942) - 6
      *threeLoop*(Sqr(g1) - 3*Sqr(g2))*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu*
      1.2020569031595942) + 6*threeLoop*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd*
      Yd.adjoint()*Yd) + 2*threeLoop*(Yd*Yd.adjoint()*Yd*Yu.adjoint()*Yu*
      Yd.adjoint()*Yd) - 2*threeLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd*
      Yd.adjoint()*Yd) + 4*threeLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd*
      Yu.adjoint()*Yu) + 6*threeLoop*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu*
      Yd.adjoint()*Yd) + 18*threeLoop*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd*
      Yd.adjoint()*Yd*1.2020569031595942) + 6*threeLoop*(Yd*Yu.adjoint()*Yu*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu*1.2020569031595942))*UNITMATRIX(3)).real(
      );

   beta_Yd = beta_Yd_1 + beta_Yd_2;


   return beta_Yd;
}

} // namespace flexiblesusy
