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
 * Calculates the one-loop beta function of g1.
 *
 * @return one-loop beta function
 */
double MSSM_susy_parameters::calc_beta_g1_one_loop(const Susy_traces& susy_traces) const
{


   double beta_g1;

   beta_g1 = Re(6.6*Power(g1,3)*oneOver16PiSqr);


   return beta_g1;
}

/**
 * Calculates the two-loop beta function of g1.
 *
 * @return two-loop beta function
 */
double MSSM_susy_parameters::calc_beta_g1_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g1;

   beta_g1 = Re(0.04*Power(g1,3)*twoLoop*(199*Sqr(g1) + 5*(-14*
      traceYdAdjYd - 18*traceYeAdjYe - 26*traceYuAdjYu + 27*Sqr(g2) + 88*Sqr(g3
      ))));


   return beta_g1;
}

/**
 * Calculates the three-loop beta function of g1.
 *
 * @return three-loop beta function
 */
double MSSM_susy_parameters::calc_beta_g1_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   const double traceAdjYdYd = TRACE_STRUCT.traceAdjYdYd;
   const double traceAdjYeYe = TRACE_STRUCT.traceAdjYeYe;
   const double traceAdjYuYu = TRACE_STRUCT.traceAdjYuYu;
   const double traceAdjYdYdAdjYdYd = TRACE_STRUCT.traceAdjYdYdAdjYdYd;
   const double traceAdjYeYeAdjYeYe = TRACE_STRUCT.traceAdjYeYeAdjYeYe;
   const double traceAdjYuYuAdjYdYd = TRACE_STRUCT.traceAdjYuYuAdjYdYd;
   const double traceAdjYuYuAdjYuYu = TRACE_STRUCT.traceAdjYuYuAdjYuYu;


   double beta_g1;

   beta_g1 = Re(-0.0026666666666666666*Power(g1,3)*threeLoop*(32117*Power
      (g1,4) + 5*Sqr(g1)*(49*traceAdjYdYd + 243*traceAdjYeYe + 169*traceAdjYuYu
      + 207*Sqr(g2) + 1096*Sqr(g3)) + 25*(243*Power(g2,4) + 9*Sqr(g2)*(11*
      traceAdjYdYd + 21*traceAdjYeYe + 29*traceAdjYuYu + 8*Sqr(g3)) - 2*(242*
      Power(g3,4) - 16*(8*traceAdjYdYd + 11*traceAdjYuYu)*Sqr(g3) + 3*(27*
      traceAdjYdYdAdjYdYd + 42*traceAdjYdYd*traceAdjYeYe + 27*
      traceAdjYeYeAdjYeYe + 29*traceAdjYuYuAdjYdYd + 42*traceAdjYuYuAdjYuYu +
      18*Sqr(traceAdjYdYd) + 12*Sqr(traceAdjYeYe) + 45*Sqr(traceAdjYuYu))))));


   return beta_g1;
}

} // namespace flexiblesusy
