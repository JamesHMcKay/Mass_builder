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

// File generated at Tue 7 Nov 2017 11:38:55

#include "MSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MSSM_susy_parameters::calc_beta_Yu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(-0.06666666666666667*Yu*(-45*traceYuAdjYu +
      13*Sqr(g1) + 45*Sqr(g2) + 80*Sqr(g3)) + Yu*Yd.adjoint()*Yd + 3*(Yu*
      Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the two-loop beta function of Yu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MSSM_susy_parameters::calc_beta_Yu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(Yu*(6.095555555555555*Power(g1,4) + 7.5*Power(g2,4
      ) - 1.7777777777777777*Power(g3,4) - 3*traceYdAdjYuYuAdjYd - 9*
      traceYuAdjYuYuAdjYu + 16*traceYuAdjYu*Sqr(g3) + 8*Sqr(g2)*Sqr(g3) + Sqr(
      g1)*(Sqr(g2) + 0.08888888888888889*(9*traceYuAdjYu + 34*Sqr(g3)))) + (-3*
      traceYdAdjYd - traceYeAdjYe + 0.4*Sqr(g1))*(Yu*Yd.adjoint()*Yd) + (-9*
      traceYuAdjYu + 0.4*Sqr(g1) + 6*Sqr(g2))*(Yu*Yu.adjoint()*Yu) - 2*(Yu*
      Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu)
      - 4*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the three-loop beta function of Yu.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MSSM_susy_parameters::calc_beta_Yu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   const double traceAdjYdYd = TRACE_STRUCT.traceAdjYdYd;
   const double traceAdjYeYe = TRACE_STRUCT.traceAdjYeYe;
   const double traceAdjYuYu = TRACE_STRUCT.traceAdjYuYu;
   const double traceAdjYuYuAdjYdYd = TRACE_STRUCT.traceAdjYuYuAdjYdYd;
   const double traceAdjYdYdAdjYdYd = TRACE_STRUCT.traceAdjYdYdAdjYdYd;
   const double traceAdjYeYeAdjYeYe = TRACE_STRUCT.traceAdjYeYeAdjYeYe;
   const double traceAdjYuYuAdjYuYu = TRACE_STRUCT.traceAdjYuYuAdjYuYu;
   const double traceAdjYdYdAdjYuYuAdjYdYd =
      TRACE_STRUCT.traceAdjYdYdAdjYuYuAdjYdYd;
   const double traceAdjYuYuAdjYuYuAdjYuYu =
      TRACE_STRUCT.traceAdjYuYuAdjYuYuAdjYuYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   const Eigen::Matrix<double,3,3> beta_Yu_1 = ((0.00007407407407407407*
      threeLoop*Yu*(704194*Power(g1,6) + 15*Power(g1,4)*(-15*(364*traceAdjYdYd
      + 468*traceAdjYeYe + 2565*traceAdjYuYu) + 6822*Sqr(g2) + 21232*Sqr(g3)) +
      150*Sqr(g1)*(765*Power(g2,4) - 9*Sqr(g2)*(57*traceAdjYuYu + 16*Sqr(g3))
      + 2*(436*Power(g3,4) + 54*traceAdjYuYuAdjYdYd + 513*traceAdjYuYuAdjYuYu -
      1860*traceAdjYuYu*Sqr(g3))) + 125*(18630*Power(g2,6) + 135*Power(g2,4)*(
      -3*(12*traceAdjYdYd + 4*traceAdjYeYe + 21*traceAdjYuYu) + 112*Sqr(g3)) +
      108*Sqr(g2)*(68*Power(g3,4) + 9*(2*traceAdjYuYuAdjYdYd +
      traceAdjYuYuAdjYuYu) - 132*traceAdjYuYu*Sqr(g3)) + 4*(5440*Power(g3,6) -
      1440*Power(g3,4)*(traceAdjYdYd + 2*traceAdjYuYu) + 81*(3*
      traceAdjYdYdAdjYuYuAdjYdYd + 6*traceAdjYdYd*traceAdjYuYuAdjYdYd + 2*
      traceAdjYeYe*traceAdjYuYuAdjYdYd + 18*traceAdjYuYu*traceAdjYuYuAdjYuYu +
      traceAdjYuYuAdjYuYuAdjYuYu) + 648*(traceAdjYuYuAdjYdYd + 3*
      traceAdjYuYuAdjYuYu)*Sqr(g3)))) - 0.004*threeLoop*(5174*Power(g1,6) + 65*
      Power(g1,4)*(5*traceAdjYuYu + 54*Sqr(g2) + 176*Sqr(g3)) + 50*Sqr(g1)*(81*
      Power(g2,4) + 176*Power(g3,4) - 6*traceAdjYuYuAdjYdYd + 18*
      traceAdjYuYuAdjYuYu - 63*traceAdjYuYu*Sqr(g2) - 208*traceAdjYuYu*Sqr(g3))
      - 125*(630*Power(g2,6) - 27*Power(g2,4)*(7*traceAdjYuYu + 16*Sqr(g3)) -
      36*Sqr(g2)*(8*Power(g3,4) - 3*traceAdjYuYuAdjYuYu - 8*traceAdjYuYu*Sqr(g3
      )) + 4*(320*Power(g3,6) - 8*Power(g3,4)*traceAdjYuYu + 9*
      traceAdjYuYuAdjYuYuAdjYuYu - 24*(traceAdjYuYuAdjYdYd + 3*
      traceAdjYuYuAdjYuYu)*Sqr(g3))))*(Yu*1.2020569031595942) -
      0.0033333333333333335*threeLoop*(1899*Power(g1,4) + 10*Sqr(g1)*(-96*
      traceAdjYdYd + 48*traceAdjYeYe + 123*Sqr(g2) + 152*Sqr(g3)) + 25*(135*
      Power(g2,4) + 24*Sqr(g2)*(-3*(3*traceAdjYdYd + traceAdjYeYe) + 2*Sqr(g3))
      - 4*(8*Power(g3,4) - 24*(traceAdjYdYd - traceAdjYeYe)*Sqr(g3) - 3*(-18*
      traceAdjYdYdAdjYdYd + 6*traceAdjYdYd*traceAdjYeYe - 6*traceAdjYeYeAdjYeYe
      - 6*traceAdjYuYuAdjYdYd + 9*Sqr(traceAdjYdYd) + Sqr(traceAdjYeYe)))))*(
      Yu*Yd.adjoint()*Yd) + (-28.536666666666665*Power(g1,4)*threeLoop - 54.75*
      Power(g2,4)*threeLoop - 0.03333333333333333*threeLoop*Sqr(g1)*(-270*
      traceAdjYuYu + 579*Sqr(g2) + 424*Sqr(g3)) + Sqr(g2)*(45*threeLoop*
      traceAdjYuYu - 92*threeLoop*Sqr(g3)) + threeLoop*(8*Power(g3,4) + 18*
      traceAdjYuYuAdjYdYd + 54*traceAdjYuYuAdjYuYu - 24*traceAdjYuYu*Sqr(g3) -
      27*Sqr(traceAdjYuYu)))*(Yu*Yu.adjoint()*Yu) + 0.03333333333333333*
      threeLoop*(7*Power(g1,4) + 2*Sqr(g1)*(-72*traceAdjYdYd + 36*traceAdjYeYe
      + 27*Sqr(g2) + 128*Sqr(g3)) - 5*(189*Power(g2,4) + 544*Power(g3,4) - 288*
      traceAdjYdYd*Sqr(g3)))*(Yu*Yd.adjoint()*Yd*1.2020569031595942) - 0.02*
      threeLoop*(117*Power(g1,4) - 10*Sqr(g1)*(18*traceAdjYuYu + 123*Sqr(g2) -
      32*Sqr(g3)) + 25*(81*Power(g2,4) - 12*Sqr(g2)*(-9*traceAdjYuYu + 16*Sqr(
      g3)) + 32*Sqr(g3)*(-9*traceAdjYuYu + 17*Sqr(g3))))*(Yu*Yu.adjoint()*Yu*
      1.2020569031595942))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_Yu_2 = ((0.06666666666666667*
      threeLoop*(7*Sqr(g1) + 5*(18*traceAdjYdYd + 6*traceAdjYeYe - 9*Sqr(g2) +
      64*Sqr(g3)))*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) + (1.2666666666666666*
      threeLoop*Sqr(g1) + 9*threeLoop*Sqr(g2) + 0.6666666666666666*threeLoop*(
      18*traceAdjYdYd + 6*traceAdjYeYe - 9*traceAdjYuYu + 32*Sqr(g3)))*(Yu*
      Yd.adjoint()*Yd*Yu.adjoint()*Yu) + 0.6666666666666666*threeLoop*(18*
      traceAdjYuYu + 5*Sqr(g1) + 9*Sqr(g2) + 64*Sqr(g3))*(Yu*Yu.adjoint()*Yu*
      Yu.adjoint()*Yu) - 1.2*threeLoop*(Sqr(g1) - 15*Sqr(g2))*(Yu*Yd.adjoint()*
      Yd*Yd.adjoint()*Yd*1.2020569031595942) + 3.6*threeLoop*(Sqr(g1) - 5*Sqr(
      g2))*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu*1.2020569031595942) + 6*
      threeLoop*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd*Yu.adjoint()*Yu) + 4*
      threeLoop*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*
      threeLoop*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu) + 2*
      threeLoop*(Yu*Yu.adjoint()*Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) + 6*
      threeLoop*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu) + 6*
      threeLoop*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd*
      1.2020569031595942) + 18*threeLoop*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu*
      Yu.adjoint()*Yu*1.2020569031595942))*UNITMATRIX(3)).real();

   beta_Yu = beta_Yu_1 + beta_Yu_2;


   return beta_Yu;
}

} // namespace flexiblesusy
