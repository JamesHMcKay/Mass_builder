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

// File generated at Wed 15 Nov 2017 19:56:11

#include "EW_triplet_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> EW_triplet_susy_parameters::calc_beta_Yu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(Yu*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu - 0.85*Sqr(g1) - 2.25*Sqr(g2) - 8*Sqr(g3)) - 1.5*(Yu*
      Yd.adjoint()*Yd) + 1.5*(Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the two-loop beta function of Yu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> EW_triplet_susy_parameters::calc_beta_Yu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(0.0016666666666666668*Yu*(1187*Power(g1,4) + 5*Sqr
      (g1)*(75*traceYdAdjYd + 225*traceYeAdjYe + 255*traceYuAdjYu - 54*Sqr(g2)
      + 152*Sqr(g3)) - 75*(38*Power(g2,4) - 3*Sqr(g2)*(5*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu) + 24*Sqr(g3)) + 2*(432*Power(g3,4) - 80*(
      traceYdAdjYd + traceYuAdjYu)*Sqr(g3) + 3*(9*traceYdAdjYdYdAdjYd - 2*
      traceYdAdjYuYuAdjYd + 3*traceYeAdjYeYeAdjYe + 9*traceYuAdjYuYuAdjYu - 8*
      Sqr(LamH))))) + 0.0125*(-43*Sqr(g1) + 5*(20*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu) + 9*Sqr(g2) - 256*Sqr(g3)))*(Yu*Yd.adjoint
      ()*Yd) + 0.0125*(223*Sqr(g1) + 5*(135*Sqr(g2) + 4*(-3*(16*LamH + 9*
      traceYdAdjYd + 3*traceYeAdjYe + 9*traceYuAdjYu) + 64*Sqr(g3))))*(Yu*
      Yu.adjoint()*Yu) + 2.75*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 0.25*(Yu*
      Yd.adjoint()*Yd*Yu.adjoint()*Yu) - Yu*Yu.adjoint()*Yu*Yd.adjoint()*Yd +
      1.5*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the three-loop beta function of Yu.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> EW_triplet_susy_parameters::calc_beta_Yu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

} // namespace flexiblesusy
