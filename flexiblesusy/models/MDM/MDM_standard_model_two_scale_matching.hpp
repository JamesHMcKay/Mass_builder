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

#ifndef MDM_TWO_SCALE_MATCHING_H
#define MDM_TWO_SCALE_MATCHING_H

#include "two_scale_matching.hpp"

namespace flexiblesusy {

class Two_scale;
class Two_scale_model;

template <class T> class Constraint;
template <class T> class MDM;
template <class T> class MDM_standard_model_Matching;

namespace standard_model {
template <class T> class StandardModel;
}

template<>
class MDM_standard_model_Matching<Two_scale> : public Matching<Two_scale>
{
public:
   MDM_standard_model_Matching();

   MDM_standard_model_Matching(standard_model::StandardModel<Two_scale>*,
                                       MDM<Two_scale>*,
                                       Constraint<Two_scale>*,
                                       unsigned, unsigned, unsigned);

   virtual ~MDM_standard_model_Matching();

   virtual double get_scale() const;
   virtual void set_models(Two_scale_model*, Two_scale_model*);
   virtual void match_high_to_low_scale_model();
   virtual void match_low_to_high_scale_model();

   unsigned get_higgs_index() const;
   unsigned get_loop_order_up() const;
   unsigned get_loop_order_down() const;
   void set_higgs_index(unsigned);
   void set_loop_order_up(unsigned);
   void set_loop_order_down(unsigned);
   void set_scale(double);
   void set_constraint(Constraint<Two_scale>*);
   void match_high_to_low_scale_model_tree_level();
   void match_low_to_high_scale_model_tree_level();

private:
   MDM<Two_scale>* model;
   standard_model::StandardModel<Two_scale>* eft;
   Constraint<Two_scale>* constraint;
   double scale;
   unsigned loop_order_up, loop_order_down;
   unsigned higgs_idx; ///< Higgs index
};

} // namespace flexiblesusy

#endif
