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

// File generated at Wed 15 Nov 2017 19:56:13

#ifndef EW_triplet_TWO_SCALE_soft_parameters_H
#define EW_triplet_TWO_SCALE_soft_parameters_H

#include "EW_triplet_two_scale_susy_parameters.hpp"

#include <iosfwd>

namespace flexiblesusy {

#ifdef TRACE_STRUCT_TYPE
   #undef TRACE_STRUCT_TYPE
#endif
#define TRACE_STRUCT_TYPE Soft_traces

class EW_triplet_soft_parameters : public EW_triplet_susy_parameters {
public:
   explicit EW_triplet_soft_parameters(const EW_triplet_input_parameters& input_ = EW_triplet_input_parameters());
   EW_triplet_soft_parameters(const EW_triplet_susy_parameters& , double Yc_, double mu2_, double v_
);
   virtual ~EW_triplet_soft_parameters() {}
   virtual Eigen::ArrayXd beta() const;
   virtual Eigen::ArrayXd get() const;
   virtual void print(std::ostream&) const;
   virtual void set(const Eigen::ArrayXd&);

   EW_triplet_soft_parameters calc_beta() const;
   EW_triplet_soft_parameters calc_beta(unsigned) const;
   virtual void clear();

   void set_Yc(double Yc_) { Yc = Yc_; }
   void set_mu2(double mu2_) { mu2 = mu2_; }
   void set_v(double v_) { v = v_; }

   double get_Yc() const { return Yc; }
   double get_mu2() const { return mu2; }
   double get_v() const { return v; }


protected:
   double Yc;
   double mu2;
   double v;


private:
   static const int numberOfParameters = 34;

   struct Soft_traces {
      double traceYdAdjYd;
      double traceYeAdjYe;
      double traceYuAdjYu;
      double traceYdAdjYdYdAdjYd;
      double traceYdAdjYuYuAdjYd;
      double traceYeAdjYeYeAdjYe;
      double traceYuAdjYuYuAdjYu;

   };
   void calc_soft_traces(Soft_traces&) const;

   double calc_beta_Yc_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Yc_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_Yc_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_mu2_three_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_one_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_two_loop(const TRACE_STRUCT_TYPE&) const;
   double calc_beta_v_three_loop(const TRACE_STRUCT_TYPE&) const;

};

std::ostream& operator<<(std::ostream&, const EW_triplet_soft_parameters&);

} // namespace flexiblesusy

#undef TRACE_STRUCT_TYPE

#endif
