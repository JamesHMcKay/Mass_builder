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

// File generated at Sat 2 Sep 2017 18:58:37

#include "MSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Mu.
 *
 * @return one-loop beta function
 */
double MSSM_susy_parameters::calc_beta_Mu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Mu;

   beta_Mu = Re(oneOver16PiSqr*(-0.6*Mu*Sqr(g1) + Mu*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu - 3*Sqr(g2))));


   return beta_Mu;
}

/**
 * Calculates the two-loop beta function of Mu.
 *
 * @return two-loop beta function
 */
double MSSM_susy_parameters::calc_beta_Mu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Mu;

   beta_Mu = Re(0.02*twoLoop*Mu*(207*Power(g1,4) + 10*Sqr(g1)*(-2*
      traceYdAdjYd + 6*traceYeAdjYe + 4*traceYuAdjYu + 9*Sqr(g2)) + 25*(15*
      Power(g2,4) + 2*(-3*(3*traceYdAdjYdYdAdjYd + 2*traceYdAdjYuYuAdjYd +
      traceYeAdjYeYeAdjYe + 3*traceYuAdjYuYuAdjYu) + 16*(traceYdAdjYd +
      traceYuAdjYu)*Sqr(g3)))));


   return beta_Mu;
}

/**
 * Calculates the three-loop beta function of Mu.
 *
 * @return three-loop beta function
 */
double MSSM_susy_parameters::calc_beta_Mu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   const double traceAdjYuYu = TRACE_STRUCT.traceAdjYuYu;
   const double traceAdjYdYdAdjYdYd = TRACE_STRUCT.traceAdjYdYdAdjYdYd;
   const double traceAdjYeYeAdjYeYe = TRACE_STRUCT.traceAdjYeYeAdjYeYe;
   const double traceAdjYuYuAdjYdYd = TRACE_STRUCT.traceAdjYuYuAdjYdYd;
   const double traceAdjYeYe = TRACE_STRUCT.traceAdjYeYe;
   const double traceAdjYdYd = TRACE_STRUCT.traceAdjYdYd;
   const double traceAdjYuYuAdjYuYu = TRACE_STRUCT.traceAdjYuYuAdjYuYu;
   const double traceAdjYdYdAdjYdYdAdjYdYd =
      TRACE_STRUCT.traceAdjYdYdAdjYdYdAdjYdYd;
   const double traceAdjYdYdAdjYuYuAdjYdYd =
      TRACE_STRUCT.traceAdjYdYdAdjYuYuAdjYdYd;
   const double traceAdjYeYeAdjYeYeAdjYeYe =
      TRACE_STRUCT.traceAdjYeYeAdjYeYeAdjYeYe;
   const double traceAdjYuYuAdjYuYuAdjYdYd =
      TRACE_STRUCT.traceAdjYuYuAdjYuYuAdjYdYd;
   const double traceAdjYuYuAdjYuYuAdjYuYu =
      TRACE_STRUCT.traceAdjYuYuAdjYuYuAdjYuYu;


   double beta_Mu;

   beta_Mu = Re(threeLoop*Mu*(19.556928691529336*Power(g1,6) +
      551.1479244952723*Power(g2,6) - 72.56624378388685*Power(g3,4)*
      traceAdjYdYd + 54*traceAdjYdYd*traceAdjYdYdAdjYdYd + 24.637024256872696*
      traceAdjYdYdAdjYdYdAdjYdYd + 9*traceAdjYdYdAdjYuYuAdjYdYd + 18*
      traceAdjYdYdAdjYdYd*traceAdjYeYe + 18*traceAdjYdYd*traceAdjYeYeAdjYeYe +
      6*traceAdjYeYe*traceAdjYeYeAdjYeYe + 8.212341418957566*
      traceAdjYeYeAdjYeYeAdjYeYe - 72.56624378388685*Power(g3,4)*traceAdjYuYu +
      18*traceAdjYdYd*traceAdjYuYuAdjYdYd + 6*traceAdjYeYe*traceAdjYuYuAdjYdYd
      + 18*traceAdjYuYu*traceAdjYuYuAdjYdYd + 54*traceAdjYuYu*
      traceAdjYuYuAdjYuYu + 9*traceAdjYuYuAdjYuYuAdjYdYd + 24.637024256872696*
      traceAdjYuYuAdjYuYuAdjYuYu + Power(g2,4)*(-192.34437734858165*
      traceAdjYdYd - 64.11479244952722*traceAdjYeYe - 192.34437734858165*
      traceAdjYuYu - 79.64429108247236*Sqr(g3)) + Power(g1,4)*(
      -18.534500964199108*traceAdjYdYd - 30.902667816881458*traceAdjYeYe -
      40.84600730744081*traceAdjYuYu - 11.143993098711256*Sqr(g2) -
      11.681162692095945*Sqr(g3)) - 101.09619405498157*traceAdjYdYdAdjYdYd*Sqr(
      g3) - 67.39746270332105*traceAdjYuYuAdjYdYd*Sqr(g3) - 101.09619405498157*
      traceAdjYuYuAdjYuYu*Sqr(g3) + Sqr(g1)*(-14.973321831185427*Power(g2,4) +
      15.982214554123619*traceAdjYdYdAdjYdYd - 3.982214554123619*
      traceAdjYeYeAdjYeYe + 10.339746270332105*traceAdjYuYuAdjYdYd +
      7.072595148625461*traceAdjYuYuAdjYuYu + (-11.118512128436349*traceAdjYdYd
      + 11.373321831185427*traceAdjYeYe + 9.445916979810885*traceAdjYuYu)*Sqr(
      g2) + (7.992741297441576*traceAdjYdYd + 8.67223383810579*traceAdjYuYu)*
      Sqr(g3)) + Sqr(g2)*(73.91107277061809*traceAdjYdYdAdjYdYd +
      24.637024256872696*traceAdjYeYeAdjYeYe + 36*traceAdjYuYuAdjYdYd +
      73.91107277061809*traceAdjYuYuAdjYuYu + (41.09619405498157*traceAdjYdYd +
      41.09619405498157*traceAdjYuYu)*Sqr(g3))));


   return beta_Mu;
}

} // namespace flexiblesusy
