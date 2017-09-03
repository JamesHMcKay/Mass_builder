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

// File generated at Thu 24 Aug 2017 11:18:21

#include "EW_triplet_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces
#define TRACE_STRUCT_TYPE Soft_traces
#define CALCULATE_TRACES() calc_soft_traces(TRACE_STRUCT);

const int EW_triplet_soft_parameters::numberOfParameters;

EW_triplet_soft_parameters::EW_triplet_soft_parameters(const EW_triplet_input_parameters& input_)
   : EW_triplet_susy_parameters(input_)
   , Yc(0), mu2(0), v(0)

{
   set_number_of_parameters(numberOfParameters);
}

EW_triplet_soft_parameters::EW_triplet_soft_parameters(
   const EW_triplet_susy_parameters& susy_model
   , double Yc_, double mu2_, double v_

)
   : EW_triplet_susy_parameters(susy_model)
   , Yc(Yc_), mu2(mu2_), v(v_)

{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd EW_triplet_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

EW_triplet_soft_parameters EW_triplet_soft_parameters::calc_beta() const
{
   double beta_Yc = 0.;
   double beta_mu2 = 0.;
   double beta_v = 0.;

   if (get_loops() > 0) {
      TRACE_STRUCT_TYPE TRACE_STRUCT;
      CALCULATE_TRACES();

      beta_Yc += calc_beta_Yc_one_loop(TRACE_STRUCT);
      beta_mu2 += calc_beta_mu2_one_loop(TRACE_STRUCT);
      beta_v += calc_beta_v_one_loop(TRACE_STRUCT);

      if (get_loops() > 1) {
         beta_Yc += calc_beta_Yc_two_loop(TRACE_STRUCT);
         beta_mu2 += calc_beta_mu2_two_loop(TRACE_STRUCT);
         beta_v += calc_beta_v_two_loop(TRACE_STRUCT);

         if (get_loops() > 2) {

         }
      }
   }


   const EW_triplet_susy_parameters susy_betas(EW_triplet_susy_parameters::calc_beta());

   return EW_triplet_soft_parameters(susy_betas, beta_Yc, beta_mu2, beta_v);
}

EW_triplet_soft_parameters EW_triplet_soft_parameters::calc_beta(unsigned loops) const
{
   EW_triplet_soft_parameters p(*this);
   p.set_loops(loops);

   return p.calc_beta();
}

void EW_triplet_soft_parameters::clear()
{
   EW_triplet_susy_parameters::clear();

   Yc = 0.;
   mu2 = 0.;
   v = 0.;

}

Eigen::ArrayXd EW_triplet_soft_parameters::get() const
{
   Eigen::ArrayXd pars(EW_triplet_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(31) = Yc;
   pars(32) = mu2;
   pars(33) = v;


   return pars;
}

void EW_triplet_soft_parameters::print(std::ostream& ostr) const
{
   EW_triplet_susy_parameters::print(ostr);
   ostr << "soft parameters at Q = " << get_scale() << ":\n";
   ostr << "Yc = " << Yc << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "v = " << v << '\n';

}

void EW_triplet_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   EW_triplet_susy_parameters::set(pars);

   Yc = pars(31);
   mu2 = pars(32);
   v = pars(33);

}

void EW_triplet_soft_parameters::calc_soft_traces(Soft_traces& soft_traces) const
{
   if (get_loops() > 0) {
      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());

   }

   if (get_loops() > 1) {
      TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint())
         .trace());
      TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint())
         .trace());
      TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint())
         .trace());
      TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint())
         .trace());

   }

   if (get_loops() > 2) {

   }
}

std::ostream& operator<<(std::ostream& ostr, const EW_triplet_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
