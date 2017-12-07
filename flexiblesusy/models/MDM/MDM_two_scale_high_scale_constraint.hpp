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

// File generated at Mon 27 Nov 2017 21:10:29

#ifndef MDM_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H
#define MDM_TWO_SCALE_HIGH_SCALE_CONSTRAINT_H

#include "MDM_high_scale_constraint.hpp"
#include "MDM_input_parameters.hpp"
#include "two_scale_constraint.hpp"

namespace flexiblesusy {

template <class T>
class MDM;

class Two_scale;

template<>
class MDM_high_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   MDM_high_scale_constraint();
   MDM_high_scale_constraint(MDM<Two_scale>*);
   virtual ~MDM_high_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   double get_initial_scale_guess() const;
   const MDM_input_parameters& get_input_parameters() const;
   MDM<Two_scale>* get_model() const;
   void initialize();
   void set_scale(double); ///< fix unification scale (0 = unfixed)

protected:
   void update_scale();
   bool check_non_perturbative();

private:
   double scale;
   double initial_scale_guess;
   MDM<Two_scale>* model;
};

} // namespace flexiblesusy

#endif
