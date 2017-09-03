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


#ifndef STANDARD_MODEL_EFFECTIVE_COUPLINGS_H
#define STANDARD_MODEL_EFFECTIVE_COUPLINGS_H

#include "standard_model.hpp"
#include "lowe.h"
#include "physical_input.hpp"

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

namespace standard_model {

class Standard_model_effective_couplings {
public:
   Standard_model_effective_couplings(const Standard_model&,
                                   const softsusy::QedQcd&,
                                   const Physical_input&);
   ~Standard_model_effective_couplings();

   void do_run_couplings(bool flag) { rg_improve = flag; }
   bool do_run_couplings() const { return rg_improve; }
   void do_include_qcd_corrections(bool flag) { include_qcd_corrections = flag; }
   bool do_include_qcd_corrections() const { return include_qcd_corrections; }
   void set_physical_inputs(const Physical_input& inputs_) { physical_input = inputs_; }
   void set_low_energy_data(const softsusy::QedQcd& qedqcd_) { qedqcd = qedqcd_; }
   void set_model(const Standard_model& model_);

   double get_hhVPVP_partial_width() const;
   double get_hhVGVG_partial_width() const;
   double get_AhVPVP_partial_width() const;
   double get_AhVGVG_partial_width() const;
   std::complex<double> get_eff_CphhVPVP() const { return eff_CphhVPVP; }
   std::complex<double> get_eff_CphhVGVG() const { return eff_CphhVGVG; }
   std::complex<double> get_eff_CpAhVPVP() const { return eff_CpAhVPVP; }
   std::complex<double> get_eff_CpAhVGVG() const { return eff_CpAhVGVG; }

   void calculate_effective_couplings();

   std::complex<double> CpFdhhbarFdPL(unsigned gt1, unsigned gt3) const;
   std::complex<double> CpFuhhbarFuPL(unsigned gt1, unsigned gt3) const;
   std::complex<double> CpFehhbarFePL(unsigned gt1, unsigned gt3) const;
   double CphhVWpconjVWp() const;
   std::complex<double> CpAhFdbarFdPL(unsigned gt2, unsigned gt3) const;
   std::complex<double> CpAhFubarFuPL(unsigned gt2, unsigned gt3) const;
   std::complex<double> CpAhFebarFePL(unsigned gt2, unsigned gt3) const;
   void calculate_eff_CphhVPVP();
   void calculate_eff_CphhVGVG();
   void calculate_eff_CpAhVPVP();
   void calculate_eff_CpAhVGVG();

private:
   Standard_model model;
   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   bool rg_improve;
   bool include_qcd_corrections;

   void copy_mixing_matrices_from_model();

   standard_model::Standard_model initialise_SM() const;
   void run_SM_strong_coupling_to(standard_model::Standard_model, double m);

   // higher order corrections to the amplitudes for
   // effective coupling to photons
   std::complex<double> scalar_scalar_qcd_factor(double, double) const;
   std::complex<double> scalar_fermion_qcd_factor(double, double) const;
   std::complex<double> pseudoscalar_fermion_qcd_factor(double, double) const;

   // higher order corrections to the leading order
   // effective couplings to gluons
   double number_of_active_flavours(double) const;
   double scalar_scaling_factor(double) const;
   double pseudoscalar_scaling_factor(double) const;

   Eigen::Matrix<std::complex<double>,3,3> Vd;
   Eigen::Matrix<std::complex<double>,3,3> Ud;
   Eigen::Matrix<std::complex<double>,3,3> Vu;
   Eigen::Matrix<std::complex<double>,3,3> Uu;
   Eigen::Matrix<std::complex<double>,3,3> Ve;
   Eigen::Matrix<std::complex<double>,3,3> Ue;
   Eigen::Matrix<double,2,2> ZZ;

   std::complex<double> eff_CphhVPVP;
   std::complex<double> eff_CphhVGVG;
   std::complex<double> eff_CpAhVPVP;
   std::complex<double> eff_CpAhVGVG;

};

} // namespace standard_model
} // namespace flexiblesusy

#endif
