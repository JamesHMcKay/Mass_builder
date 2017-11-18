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

// File generated at Wed 15 Nov 2017 15:15:50

/**
 * @file MDM_two_scale_model.hpp
 * @brief contains class for model with routines needed to solve boundary
 *        value problem using the two_scale solver by solvingt EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Wed 15 Nov 2017 15:15:50 with FlexibleSUSY
 * 1.7.4 (git commit: unknown) and SARAH 4.12.2 .
 */

#ifndef MDM_TWO_SCALE_H
#define MDM_TWO_SCALE_H

#include "MDM_model.hpp"
#include "MDM_mass_eigenstates.hpp"
#include "two_scale_model.hpp"

namespace flexiblesusy {

class Two_scale;
/**
 * @class MDM<Two_scale>
 * @brief model class with routines for determing masses and mixinga and EWSB
 */
template<>
class MDM<Two_scale> : public Two_scale_model, public MDM_mass_eigenstates {
public:
   explicit MDM(const MDM_input_parameters& input_ = MDM_input_parameters());
   virtual ~MDM();

   // interface functions
   virtual void calculate_spectrum();
   virtual void clear_problems();
   virtual std::string name() const;
   virtual void run_to(double scale, double eps = -1.0);
   virtual void print(std::ostream& out = std::cout) const;
   virtual void set_precision(double);
};

std::ostream& operator<<(std::ostream&, const MDM<Two_scale>&);

} // namespace flexiblesusy

#endif
