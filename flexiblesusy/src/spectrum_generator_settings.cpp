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

#include "spectrum_generator_settings.hpp"

#include <iostream>

namespace flexiblesusy {

/**
 * Default constructor
 *
 * Calls reset() to initialize all spectrum generator settings to
 * their default values.
 */
Spectrum_generator_settings::Spectrum_generator_settings()
{
   reset();
}

double Spectrum_generator_settings::get(Settings o) const
{
   return values.at(o);
}

void Spectrum_generator_settings::set(Settings o, double value)
{
   values.at(o) = value;
}

/**
 * Resets all spectrum generator settings to their defaults.
 *
 * | enum                             | possible values              | default value   |
 * |----------------------------------|------------------------------|-----------------|
 * | precision                        | any positive double          | 1.0e-4          |
 * | max_iterations                   | any positive double          | 0 (= automatic) |
 * | algorithm                        | 0 (two-scale) or 1 (lattice) | 0 (= two-scale) |
 * | calculate_sm_masses              | 0 (no) or 1 (yes)            | 0 (= no)        |
 * | pole_mass_loop_order             | 0, 1, 2                      | 2 (= 2-loop)    |
 * | ewsb_loop_order                  | 0, 1, 2                      | 2 (= 2-loop)    |
 * | beta_loop_order                  | 0, 1, 2, 3                   | 2 (= 2-loop)    |
 * | threshold_corrections_loop_order | 0, 1, 2                      | 2 (= 2-loop)    |
 * | higgs_2loop_correction_at_as     | 0, 1                         | 1 (= enabled)   |
 * | higgs_2loop_correction_ab_as     | 0, 1                         | 1 (= enabled)   |
 * | higgs_2loop_correction_at_at     | 0, 1                         | 1 (= enabled)   |
 * | higgs_2loop_correction_atau_atau | 0, 1                         | 1 (= enabled)   |
 * | force_output                     | 0 (no) or 1 (yes)            | 0 (= no)        |
 * | top_pole_qcd_corrections         | 0 (1L), 1 (2L), 2 (3L)       | 1 (= 2L QCD)    |
 * | beta_zero_threshold              | any positive double          | 1.0e-11         |
 * | calculate_observables            | 0 (no) or 1 (yes)            | 0 (= no)        |
 * | force_positive_masses            | 0 (no) or 1 (yes)            | 0 (= no)        |
 * | pole_mass_scale                  | any positive double          | 0 (= SUSY scale)|
 * | eft_pole_mass_scale              | any positive double          | 0 (= minimum of {Mt, SUSY scale})|
 * | eft_matching_scale               | any positive double          | 0 (= SUSY scale)|
 * | eft_matching_loop_order_up       | 0, 1, 2                      | 2 (= 2-loop)    |
 * | eft_matching_loop_order_down     | 0, 1                         | 1 (= 1-loop)    |
 * | eft_higgs_index                  | any integer >= 0             | 0 (= lightest)  |
 * | calculate_bsm_masses             | 0 (no) or 1 (yes)            | 1 (= yes)       |
 */
void Spectrum_generator_settings::reset()
{
   values[precision]             = 1.0e-4;
   values[max_iterations]        = 0.; // 0 = automatic
   values[algorithm]             = 0.; // 0 = two-scale
   values[calculate_sm_masses]   = 0.; // 0 = false
   values[pole_mass_loop_order]  = 2.;
   values[ewsb_loop_order]       = 2.;
   values[beta_loop_order]       = 2.;
   values[threshold_corrections_loop_order] = 2.;
   values[higgs_2loop_correction_at_as]     = 1.;
   values[higgs_2loop_correction_ab_as]     = 1.;
   values[higgs_2loop_correction_at_at]     = 1.;
   values[higgs_2loop_correction_atau_atau] = 1.;
   values[calculate_sm_masses]   = 0.; // 0 = false
   values[top_pole_qcd_corrections]         = 1.;
   values[beta_zero_threshold]              = 1.0e-11;
   values[calculate_observables]            = 0;
   values[force_positive_masses]            = 0;
   values[pole_mass_scale]                  = 0;
   values[eft_pole_mass_scale]              = 0;
   values[eft_matching_scale]               = 0;
   values[eft_matching_loop_order_up]       = 2.;
   values[eft_matching_loop_order_down]     = 1.;
   values[eft_higgs_index]                  = 0;
   values[calculate_bsm_masses]             = 1.;
   values[force_output]                     = 0;
}

Two_loop_corrections Spectrum_generator_settings::get_two_loop_corrections() const
{
   Two_loop_corrections two_loop_corrections;
   two_loop_corrections.higgs_at_as     = get(higgs_2loop_correction_at_as);
   two_loop_corrections.higgs_ab_as     = get(higgs_2loop_correction_ab_as);
   two_loop_corrections.higgs_at_at     = get(higgs_2loop_correction_at_at);
   two_loop_corrections.higgs_atau_atau = get(higgs_2loop_correction_atau_atau);
   two_loop_corrections.top_qcd         = get(top_pole_qcd_corrections);

   return two_loop_corrections;
}

void Spectrum_generator_settings::set_two_loop_corrections(
   const Two_loop_corrections& two_loop_corrections)
{
   set(higgs_2loop_correction_at_as, two_loop_corrections.higgs_at_as);
   set(higgs_2loop_correction_ab_as, two_loop_corrections.higgs_ab_as);
   set(higgs_2loop_correction_at_at, two_loop_corrections.higgs_at_at);
   set(higgs_2loop_correction_atau_atau, two_loop_corrections.higgs_atau_atau);
   set(top_pole_qcd_corrections, two_loop_corrections.top_qcd);
}

std::ostream& operator<<(std::ostream& ostr, const Spectrum_generator_settings& sgs)
{
   ostr << "(";

   for (unsigned i = 0; i < Spectrum_generator_settings::NUMBER_OF_OPTIONS; i++) {
      ostr << sgs.get(static_cast<Spectrum_generator_settings::Settings>(i));
      if (i < Spectrum_generator_settings::NUMBER_OF_OPTIONS - 1)
         ostr << ", ";
   }

   ostr << ")";

   return ostr;
}

} // namespace flexiblesusy
