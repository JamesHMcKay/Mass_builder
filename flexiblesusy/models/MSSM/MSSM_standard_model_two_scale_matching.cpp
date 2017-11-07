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

// File generated at Tue 7 Nov 2017 11:41:20

#include "MSSM_standard_model_two_scale_matching.hpp"
#include "MSSM_standard_model_matching.hpp"
#include "MSSM_two_scale_model.hpp"
#include "standard_model_two_scale_model.hpp"
#include "error.hpp"
#include "two_scale_constraint.hpp"

#define CLASSNAME MSSM_standard_model_Matching<Two_scale>

namespace flexiblesusy {

CLASSNAME::MSSM_standard_model_Matching()
   : model(0), eft(0), constraint(nullptr), scale(0.)
   , loop_order_up(0), loop_order_down(0)
   , higgs_idx(0)
{}

CLASSNAME::MSSM_standard_model_Matching(
   standard_model::StandardModel<Two_scale>* low_,
   MSSM<Two_scale>* high_,
   Constraint<Two_scale>* constraint_,
   unsigned loop_order_up_, unsigned loop_order_down_,
   unsigned higgs_idx_)
   : model(high_)
   , eft(low_)
   , constraint(constraint_)
   , scale(0.)
   , higgs_idx(higgs_idx_)
{
   set_loop_order_up(loop_order_up_);
   set_loop_order_down(loop_order_down_);
}

CLASSNAME::~MSSM_standard_model_Matching()
{}

double CLASSNAME::get_scale() const
{
   if (scale != 0.)
      return scale;

   if (!constraint)
      throw SetupError("Constraint pointer in matching class is NULL!");

   return constraint->get_scale();
}

void CLASSNAME::set_models(Two_scale_model* low, Two_scale_model* high)
{
   eft = cast_model<standard_model::StandardModel<Two_scale>*>(low);
   model = cast_model<MSSM<Two_scale>*>(high);
}

void CLASSNAME::match_high_to_low_scale_model()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   if (model->get_thresholds() && loop_order_down)
      MSSM_standard_model_matching::match_high_to_low_scale_model(*eft, *model, loop_order_down, higgs_idx);
   else
      MSSM_standard_model_matching::match_high_to_low_scale_model_tree_level(*eft, *model, higgs_idx);
}

void CLASSNAME::match_low_to_high_scale_model()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   if (model->get_thresholds() && loop_order_up)
      MSSM_standard_model_matching::match_low_to_high_scale_model(*model, *eft, loop_order_up);
   else
      MSSM_standard_model_matching::match_low_to_high_scale_model_tree_level(*model, *eft);
}

void CLASSNAME::match_high_to_low_scale_model_tree_level()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   MSSM_standard_model_matching::match_high_to_low_scale_model_tree_level(*eft, *model, higgs_idx);
}

void CLASSNAME::match_low_to_high_scale_model_tree_level()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   MSSM_standard_model_matching::match_low_to_high_scale_model_tree_level(*model, *eft);
}

unsigned CLASSNAME::get_higgs_index() const
{
   return higgs_idx;
}

unsigned CLASSNAME::get_loop_order_up() const
{
   return loop_order_up;
}

unsigned CLASSNAME::get_loop_order_down() const
{
   return loop_order_down;
}

void CLASSNAME::set_higgs_index(unsigned idx)
{
   higgs_idx = idx;
}

void CLASSNAME::set_loop_order_up(unsigned loop_order_up_)
{
   loop_order_up = loop_order_up_;
}

void CLASSNAME::set_loop_order_down(unsigned loop_order_down_)
{
   if (loop_order_down_ > 1) {
      WARNING("Matching loop order " << loop_order_down_
              << " for downwards matching currently not"
              " supported!  I'm using 1-loop matching.");
      loop_order_down_ = 1;
   }

   loop_order_down = loop_order_down_;
}

void CLASSNAME::set_constraint(Constraint<Two_scale>* constraint_)
{
   constraint = constraint_;
}

void CLASSNAME::set_scale(double scale_)
{
   scale = scale_;
}

} // namespace flexiblesusy
