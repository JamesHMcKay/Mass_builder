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

// File generated at Sat 2 Sep 2017 18:58:38

#include "MSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of g3.
 *
 * @return one-loop beta function
 */
double MSSM_susy_parameters::calc_beta_g3_one_loop(const Susy_traces& susy_traces) const
{


   double beta_g3;

   beta_g3 = Re(-3*Power(g3,3)*oneOver16PiSqr);


   return beta_g3;
}

/**
 * Calculates the two-loop beta function of g3.
 *
 * @return two-loop beta function
 */
double MSSM_susy_parameters::calc_beta_g3_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g3;

   beta_g3 = Re(0.2*Power(g3,3)*twoLoop*(11*Sqr(g1) + 5*(-4*traceYdAdjYd
      - 4*traceYuAdjYu + 9*Sqr(g2) + 14*Sqr(g3))));


   return beta_g3;
}

/**
 * Calculates the three-loop beta function of g3.
 *
 * @return three-loop beta function
 */
double MSSM_susy_parameters::calc_beta_g3_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   const double traceAdjYdYd = TRACE_STRUCT.traceAdjYdYd;
   const double traceAdjYeYe = TRACE_STRUCT.traceAdjYeYe;
   const double traceAdjYuYu = TRACE_STRUCT.traceAdjYuYu;
   const double traceAdjYdYdAdjYdYd = TRACE_STRUCT.traceAdjYdYdAdjYdYd;
   const double traceAdjYuYuAdjYdYd = TRACE_STRUCT.traceAdjYuYuAdjYdYd;
   const double traceAdjYuYuAdjYuYu = TRACE_STRUCT.traceAdjYuYuAdjYuYu;


   double beta_g3;

   beta_g3 = Re(0.013333333333333334*Power(g3,3)*threeLoop*(-1702*Power(
      g1,4) - 5*Sqr(g1)*(32*traceAdjYdYd + 44*traceAdjYuYu + 9*Sqr(g2) - 22*Sqr
      (g3)) - 25*(81*Power(g2,4) - 347*Power(g3,4) + 104*(traceAdjYdYd +
      traceAdjYuYu)*Sqr(g3) - 18*Sqr(g2)*(-2*(traceAdjYdYd + traceAdjYuYu) +
      Sqr(g3)) - 6*(6*traceAdjYdYdAdjYdYd + 3*traceAdjYdYd*traceAdjYeYe + 4*
      traceAdjYuYuAdjYdYd + 6*traceAdjYuYuAdjYuYu + 9*Sqr(traceAdjYdYd) + 9*Sqr
      (traceAdjYuYu)))));


   return beta_g3;
}

} // namespace flexiblesusy
