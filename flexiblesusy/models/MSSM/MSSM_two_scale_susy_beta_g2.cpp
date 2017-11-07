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

// File generated at Tue 7 Nov 2017 11:38:56

#include "MSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of g2.
 *
 * @return one-loop beta function
 */
double MSSM_susy_parameters::calc_beta_g2_one_loop(const Susy_traces& susy_traces) const
{


   double beta_g2;

   beta_g2 = Re(Power(g2,3)*oneOver16PiSqr);


   return beta_g2;
}

/**
 * Calculates the two-loop beta function of g2.
 *
 * @return two-loop beta function
 */
double MSSM_susy_parameters::calc_beta_g2_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2;

   beta_g2 = Re(0.2*Power(g2,3)*twoLoop*(9*Sqr(g1) + 5*(-2*(3*
      traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu) + 25*Sqr(g2) + 24*Sqr(g3)))
      );


   return beta_g2;
}

/**
 * Calculates the three-loop beta function of g2.
 *
 * @return three-loop beta function
 */
double MSSM_susy_parameters::calc_beta_g2_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   const double traceAdjYdYd = TRACE_STRUCT.traceAdjYdYd;
   const double traceAdjYeYe = TRACE_STRUCT.traceAdjYeYe;
   const double traceAdjYuYu = TRACE_STRUCT.traceAdjYuYu;
   const double traceAdjYdYdAdjYdYd = TRACE_STRUCT.traceAdjYdYdAdjYdYd;
   const double traceAdjYeYeAdjYeYe = TRACE_STRUCT.traceAdjYeYeAdjYeYe;
   const double traceAdjYuYuAdjYdYd = TRACE_STRUCT.traceAdjYuYuAdjYdYd;
   const double traceAdjYuYuAdjYuYu = TRACE_STRUCT.traceAdjYuYuAdjYuYu;


   double beta_g2;

   beta_g2 = Re(0.04*Power(g2,3)*threeLoop*(-457*Power(g1,4) + 5*Sqr(g1)*
      (-11*traceAdjYdYd - 21*traceAdjYeYe - 29*traceAdjYuYu + 9*Sqr(g2) - 8*Sqr
      (g3)) + 25*(35*Power(g2,4) + Sqr(g2)*(-11*(3*traceAdjYdYd + traceAdjYeYe
      + 3*traceAdjYuYu) + 24*Sqr(g3)) + 2*(22*Power(g3,4) + 12*
      traceAdjYdYdAdjYdYd + 6*traceAdjYdYd*traceAdjYeYe + 4*traceAdjYeYeAdjYeYe
      + 6*traceAdjYuYuAdjYdYd + 12*traceAdjYuYuAdjYuYu - 16*(traceAdjYdYd +
      traceAdjYuYu)*Sqr(g3) + 9*Sqr(traceAdjYdYd) + Sqr(traceAdjYeYe) + 9*Sqr(
      traceAdjYuYu)))));


   return beta_g2;
}

} // namespace flexiblesusy
