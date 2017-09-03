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

// File generated at Thu 24 Aug 2017 11:18:34

#include "EW_triplet_physical.hpp"
#include "slha_io.hpp"

#include <iostream>

#define LOCALPHYSICAL(p) p

namespace flexiblesusy {

EW_triplet_physical::EW_triplet_physical()
   :
    MVG(0), MHp(0), MFv(Eigen::Array<double,3,1>::Zero()), MFc(0), MFn(0), MAh
       (0), Mhh(0), MFd(Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<double
       ,3,1>::Zero()), MFe(Eigen::Array<double,3,1>::Zero()), MVWp(0), MVP(0), MVZ
       (0)

   , Vd(Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ud(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), Vu(Eigen::Matrix<std::complex<double>,3,
      3>::Zero()), Uu(Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ve(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), Ue(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZZ(Eigen::Matrix<double,2,2>::Zero())

{
}

void EW_triplet_physical::clear()
{
   MVG = 0.;
   MHp = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MFc = 0.;
   MFn = 0.;
   MAh = 0.;
   Mhh = 0.;
   MFd = Eigen::Matrix<double,3,1>::Zero();
   Vd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ud = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   Vu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Uu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   Ve = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ue = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MVWp = 0.;
   MVP = 0.;
   MVZ = 0.;

}

/**
 * Convert masses and mixing matrices to Haber-Kane convention:
 * Fermion masses are always positive and mixing matrices are allowed
 * to be complex.
 */
void EW_triplet_physical::convert_to_hk()
{

}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 */
void EW_triplet_physical::convert_to_slha()
{

}

Eigen::ArrayXd EW_triplet_physical::get() const
{
   Eigen::ArrayXd pars(get_masses());

   pars.conservativeResize(133);

   pars(21) = Re(Vd(0,0));
   pars(22) = Im(Vd(0,0));
   pars(23) = Re(Vd(0,1));
   pars(24) = Im(Vd(0,1));
   pars(25) = Re(Vd(0,2));
   pars(26) = Im(Vd(0,2));
   pars(27) = Re(Vd(1,0));
   pars(28) = Im(Vd(1,0));
   pars(29) = Re(Vd(1,1));
   pars(30) = Im(Vd(1,1));
   pars(31) = Re(Vd(1,2));
   pars(32) = Im(Vd(1,2));
   pars(33) = Re(Vd(2,0));
   pars(34) = Im(Vd(2,0));
   pars(35) = Re(Vd(2,1));
   pars(36) = Im(Vd(2,1));
   pars(37) = Re(Vd(2,2));
   pars(38) = Im(Vd(2,2));
   pars(39) = Re(Ud(0,0));
   pars(40) = Im(Ud(0,0));
   pars(41) = Re(Ud(0,1));
   pars(42) = Im(Ud(0,1));
   pars(43) = Re(Ud(0,2));
   pars(44) = Im(Ud(0,2));
   pars(45) = Re(Ud(1,0));
   pars(46) = Im(Ud(1,0));
   pars(47) = Re(Ud(1,1));
   pars(48) = Im(Ud(1,1));
   pars(49) = Re(Ud(1,2));
   pars(50) = Im(Ud(1,2));
   pars(51) = Re(Ud(2,0));
   pars(52) = Im(Ud(2,0));
   pars(53) = Re(Ud(2,1));
   pars(54) = Im(Ud(2,1));
   pars(55) = Re(Ud(2,2));
   pars(56) = Im(Ud(2,2));
   pars(57) = Re(Vu(0,0));
   pars(58) = Im(Vu(0,0));
   pars(59) = Re(Vu(0,1));
   pars(60) = Im(Vu(0,1));
   pars(61) = Re(Vu(0,2));
   pars(62) = Im(Vu(0,2));
   pars(63) = Re(Vu(1,0));
   pars(64) = Im(Vu(1,0));
   pars(65) = Re(Vu(1,1));
   pars(66) = Im(Vu(1,1));
   pars(67) = Re(Vu(1,2));
   pars(68) = Im(Vu(1,2));
   pars(69) = Re(Vu(2,0));
   pars(70) = Im(Vu(2,0));
   pars(71) = Re(Vu(2,1));
   pars(72) = Im(Vu(2,1));
   pars(73) = Re(Vu(2,2));
   pars(74) = Im(Vu(2,2));
   pars(75) = Re(Uu(0,0));
   pars(76) = Im(Uu(0,0));
   pars(77) = Re(Uu(0,1));
   pars(78) = Im(Uu(0,1));
   pars(79) = Re(Uu(0,2));
   pars(80) = Im(Uu(0,2));
   pars(81) = Re(Uu(1,0));
   pars(82) = Im(Uu(1,0));
   pars(83) = Re(Uu(1,1));
   pars(84) = Im(Uu(1,1));
   pars(85) = Re(Uu(1,2));
   pars(86) = Im(Uu(1,2));
   pars(87) = Re(Uu(2,0));
   pars(88) = Im(Uu(2,0));
   pars(89) = Re(Uu(2,1));
   pars(90) = Im(Uu(2,1));
   pars(91) = Re(Uu(2,2));
   pars(92) = Im(Uu(2,2));
   pars(93) = Re(Ve(0,0));
   pars(94) = Im(Ve(0,0));
   pars(95) = Re(Ve(0,1));
   pars(96) = Im(Ve(0,1));
   pars(97) = Re(Ve(0,2));
   pars(98) = Im(Ve(0,2));
   pars(99) = Re(Ve(1,0));
   pars(100) = Im(Ve(1,0));
   pars(101) = Re(Ve(1,1));
   pars(102) = Im(Ve(1,1));
   pars(103) = Re(Ve(1,2));
   pars(104) = Im(Ve(1,2));
   pars(105) = Re(Ve(2,0));
   pars(106) = Im(Ve(2,0));
   pars(107) = Re(Ve(2,1));
   pars(108) = Im(Ve(2,1));
   pars(109) = Re(Ve(2,2));
   pars(110) = Im(Ve(2,2));
   pars(111) = Re(Ue(0,0));
   pars(112) = Im(Ue(0,0));
   pars(113) = Re(Ue(0,1));
   pars(114) = Im(Ue(0,1));
   pars(115) = Re(Ue(0,2));
   pars(116) = Im(Ue(0,2));
   pars(117) = Re(Ue(1,0));
   pars(118) = Im(Ue(1,0));
   pars(119) = Re(Ue(1,1));
   pars(120) = Im(Ue(1,1));
   pars(121) = Re(Ue(1,2));
   pars(122) = Im(Ue(1,2));
   pars(123) = Re(Ue(2,0));
   pars(124) = Im(Ue(2,0));
   pars(125) = Re(Ue(2,1));
   pars(126) = Im(Ue(2,1));
   pars(127) = Re(Ue(2,2));
   pars(128) = Im(Ue(2,2));
   pars(129) = ZZ(0,0);
   pars(130) = ZZ(0,1);
   pars(131) = ZZ(1,0);
   pars(132) = ZZ(1,1);


   return pars;
}

void EW_triplet_physical::set(const Eigen::ArrayXd& pars)
{
   set_masses(pars);

   Vd(0,0) = std::complex<double>(pars(21), pars(22));
   Vd(0,1) = std::complex<double>(pars(23), pars(24));
   Vd(0,2) = std::complex<double>(pars(25), pars(26));
   Vd(1,0) = std::complex<double>(pars(27), pars(28));
   Vd(1,1) = std::complex<double>(pars(29), pars(30));
   Vd(1,2) = std::complex<double>(pars(31), pars(32));
   Vd(2,0) = std::complex<double>(pars(33), pars(34));
   Vd(2,1) = std::complex<double>(pars(35), pars(36));
   Vd(2,2) = std::complex<double>(pars(37), pars(38));
   Ud(0,0) = std::complex<double>(pars(39), pars(40));
   Ud(0,1) = std::complex<double>(pars(41), pars(42));
   Ud(0,2) = std::complex<double>(pars(43), pars(44));
   Ud(1,0) = std::complex<double>(pars(45), pars(46));
   Ud(1,1) = std::complex<double>(pars(47), pars(48));
   Ud(1,2) = std::complex<double>(pars(49), pars(50));
   Ud(2,0) = std::complex<double>(pars(51), pars(52));
   Ud(2,1) = std::complex<double>(pars(53), pars(54));
   Ud(2,2) = std::complex<double>(pars(55), pars(56));
   Vu(0,0) = std::complex<double>(pars(57), pars(58));
   Vu(0,1) = std::complex<double>(pars(59), pars(60));
   Vu(0,2) = std::complex<double>(pars(61), pars(62));
   Vu(1,0) = std::complex<double>(pars(63), pars(64));
   Vu(1,1) = std::complex<double>(pars(65), pars(66));
   Vu(1,2) = std::complex<double>(pars(67), pars(68));
   Vu(2,0) = std::complex<double>(pars(69), pars(70));
   Vu(2,1) = std::complex<double>(pars(71), pars(72));
   Vu(2,2) = std::complex<double>(pars(73), pars(74));
   Uu(0,0) = std::complex<double>(pars(75), pars(76));
   Uu(0,1) = std::complex<double>(pars(77), pars(78));
   Uu(0,2) = std::complex<double>(pars(79), pars(80));
   Uu(1,0) = std::complex<double>(pars(81), pars(82));
   Uu(1,1) = std::complex<double>(pars(83), pars(84));
   Uu(1,2) = std::complex<double>(pars(85), pars(86));
   Uu(2,0) = std::complex<double>(pars(87), pars(88));
   Uu(2,1) = std::complex<double>(pars(89), pars(90));
   Uu(2,2) = std::complex<double>(pars(91), pars(92));
   Ve(0,0) = std::complex<double>(pars(93), pars(94));
   Ve(0,1) = std::complex<double>(pars(95), pars(96));
   Ve(0,2) = std::complex<double>(pars(97), pars(98));
   Ve(1,0) = std::complex<double>(pars(99), pars(100));
   Ve(1,1) = std::complex<double>(pars(101), pars(102));
   Ve(1,2) = std::complex<double>(pars(103), pars(104));
   Ve(2,0) = std::complex<double>(pars(105), pars(106));
   Ve(2,1) = std::complex<double>(pars(107), pars(108));
   Ve(2,2) = std::complex<double>(pars(109), pars(110));
   Ue(0,0) = std::complex<double>(pars(111), pars(112));
   Ue(0,1) = std::complex<double>(pars(113), pars(114));
   Ue(0,2) = std::complex<double>(pars(115), pars(116));
   Ue(1,0) = std::complex<double>(pars(117), pars(118));
   Ue(1,1) = std::complex<double>(pars(119), pars(120));
   Ue(1,2) = std::complex<double>(pars(121), pars(122));
   Ue(2,0) = std::complex<double>(pars(123), pars(124));
   Ue(2,1) = std::complex<double>(pars(125), pars(126));
   Ue(2,2) = std::complex<double>(pars(127), pars(128));
   ZZ(0,0) = pars(129);
   ZZ(0,1) = pars(130);
   ZZ(1,0) = pars(131);
   ZZ(1,1) = pars(132);

}

Eigen::ArrayXd EW_triplet_physical::get_masses() const
{
   Eigen::ArrayXd pars(21);

   pars(0) = MVG;
   pars(1) = MHp;
   pars(2) = MFv(0);
   pars(3) = MFv(1);
   pars(4) = MFv(2);
   pars(5) = MFc;
   pars(6) = MFn;
   pars(7) = MAh;
   pars(8) = Mhh;
   pars(9) = MFd(0);
   pars(10) = MFd(1);
   pars(11) = MFd(2);
   pars(12) = MFu(0);
   pars(13) = MFu(1);
   pars(14) = MFu(2);
   pars(15) = MFe(0);
   pars(16) = MFe(1);
   pars(17) = MFe(2);
   pars(18) = MVWp;
   pars(19) = MVP;
   pars(20) = MVZ;

   return pars;
}

void EW_triplet_physical::set_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MHp = pars(1);
   MFv(0) = pars(2);
   MFv(1) = pars(3);
   MFv(2) = pars(4);
   MFc = pars(5);
   MFn = pars(6);
   MAh = pars(7);
   Mhh = pars(8);
   MFd(0) = pars(9);
   MFd(1) = pars(10);
   MFd(2) = pars(11);
   MFu(0) = pars(12);
   MFu(1) = pars(13);
   MFu(2) = pars(14);
   MFe(0) = pars(15);
   MFe(1) = pars(16);
   MFe(2) = pars(17);
   MVWp = pars(18);
   MVP = pars(19);
   MVZ = pars(20);

}

void EW_triplet_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "pole masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MHp = " << MHp << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MFc = " << MFc << '\n';
   ostr << "MFn = " << MFn << '\n';
   ostr << "MAh = " << MAh << '\n';
   ostr << "Mhh = " << Mhh << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MVWp = " << MVWp << '\n';

   ostr << "----------------------------------------\n"
           "pole mass mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "Vd = " << Vd << '\n';
   ostr << "Ud = " << Ud << '\n';
   ostr << "Vu = " << Vu << '\n';
   ostr << "Uu = " << Uu << '\n';
   ostr << "Ve = " << Ve << '\n';
   ostr << "Ue = " << Ue << '\n';
   ostr << "ZZ = " << ZZ << '\n';

}

std::ostream& operator<<(std::ostream& ostr, const EW_triplet_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy