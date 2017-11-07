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

// File generated at Tue 7 Nov 2017 11:40:45

#ifndef MSSM_INPUT_PARAMETERS_H
#define MSSM_INPUT_PARAMETERS_H

#include <complex>
#include <Eigen/Core>

namespace flexiblesusy {

struct MSSM_input_parameters {
   double TanBeta;
   int SignMu;
   double Qin;
   double QSUSY;
   double mHd2IN;
   double mHu2IN;
   Eigen::Matrix<double,3,3> Aeij;
   Eigen::Matrix<double,3,3> Adij;
   Eigen::Matrix<double,3,3> Auij;
   Eigen::Matrix<double,3,3> mq2Input;
   Eigen::Matrix<double,3,3> ml2Input;
   Eigen::Matrix<double,3,3> md2Input;
   Eigen::Matrix<double,3,3> mu2Input;
   Eigen::Matrix<double,3,3> me2Input;
   double MassBInput;
   double MassWBInput;
   double MassGInput;

   MSSM_input_parameters()
      : TanBeta(0), SignMu(1), Qin(0), QSUSY(0), mHd2IN(0), mHu2IN(0), Aeij(
   Eigen::Matrix<double,3,3>::Zero()), Adij(Eigen::Matrix<double,3,3>::Zero()),
   Auij(Eigen::Matrix<double,3,3>::Zero()), mq2Input(Eigen::Matrix<double,3,3>
   ::Zero()), ml2Input(Eigen::Matrix<double,3,3>::Zero()), md2Input(
   Eigen::Matrix<double,3,3>::Zero()), mu2Input(Eigen::Matrix<double,3,3>::Zero
   ()), me2Input(Eigen::Matrix<double,3,3>::Zero()), MassBInput(0), MassWBInput
   (0), MassGInput(0)

   {}

   Eigen::ArrayXd get() const;
   void set(const Eigen::ArrayXd&);
};

std::ostream& operator<<(std::ostream&, const MSSM_input_parameters&);

} // namespace flexiblesusy

#endif
