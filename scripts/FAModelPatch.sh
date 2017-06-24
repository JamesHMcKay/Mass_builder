#!/bin/sh

originalfile=/Users/jamesmckay/Documents/Programs/Mass_builder/models/QED/QED_backup.mod
file=/Users/jamesmckay/Documents/Programs/Mass_builder/models/QED/QED.mod

cp $originalfile $file

sed -i '' -e  "s/FourVector/FAFourVector/g"  $file
sed -i '' -e  "s/DiracSpinor/FADiracSpinor/g"  $file
sed -i '' -e  "s/DiracSlash/FADiracSlash/g"  $file
sed -i '' -e  "s/DiracMatrix/FADiracMatrix/g"  $file
sed -i '' -e  "s/DiracTrace/FADiracTrace/g"  $file
sed -i '' -e  "s/MetricTensor/FAMetricTensor/g"  $file
sed -i '' -e  "s/ScalarProduct/FAScalarProduct/g"  $file
sed -i '' -e  "s/ChiralityProjector/FAChiralityProjector/g"  $file
sed -i '' -e  "s/GS/FAGS/g"  $file
sed -i '' -e  "s/SUNT/FASUNT/g"  $file
sed -i '' -e  "s/SUNF/FASUNF/g"  $file
sed -i '' -e  "s/LeviCivita/FALeviCivita/g"  $file
sed -i '' -e  "s/Loop/FALoop/g"  $file
sed -i '' -e  "s/NonCommutative/FANonCommutative/g"  $file
sed -i '' -e  "s/PolarizationVector/FAPolarizationVector/g"  $file
sed -i '' -e  "s/FeynAmp/FAFeynAmp/g"  $file
sed -i '' -e  "s/PropagatorDenominator/FAPropagatorDenominator/g"  $file
sed -i '' -e  "s/GaugeXi/FAGaugeXi/g"  $file
sed -i '' -e  "s/FeynAmpDenominator/FAFeynAmpDenominator/g"  $file
#sed -i '' -e  "s/FeynAmpList/FAFeynAmpList/g"  $file
#sed -i '' -e  "s/ME/FCGV[\"ME\"]/g"  $file
#sed -i '' -e  "s/MM/FCGV[\"MM\"]/g"  $file
#sed -i '' -e  "s/ML/FCGV[\"ML\"]/g"  $file
#sed -i '' -e  "s/MU/FCGV[\"MU\"]/g"  $file
#sed -i '' -e  "s/MD/FCGV[\"MD\"]/g"  $file
#sed -i '' -e  "s/MC/FCGV[\"MC\"]/g"  $file
#sed -i '' -e  "s/MS/FCGV[\"MS\"]/g"  $file
#sed -i '' -e  "s/MT/FCGV[\"MT\"]/g"  $file
#sed -i '' -e  "s/MB/FCGV[\"MB\"]/g"  $file
#sed -i '' -e  "s/MH/FCGV[\"MH\"]/g"  $file
#sed -i '' -e  "s/MZ/FCGV[\"MZ\"]/g"  $file
#sed -i '' -e  "s/MW/FCGV[\"MW\"]/g"  $file
#sed -i '' -e  "s/EL/FCGV[\"EL\"]/g"  $file
#sed -i '' -e  "s/CW/FCGV[\"CW\"]/g"  $file
#sed -i '' -e  "s/SW/FCGV[\"SW\"]/g"  $file

originalfile=/Users/jamesmckay/Documents/Programs/Mass_builder/models/QED/QED_backup.gen
file=/Users/jamesmckay/Documents/Programs/Mass_builder/models/QED/QED.gen

cp $originalfile $file

sed -i '' -e  "s/FourVector/FAFourVector/g"  $file
sed -i '' -e  "s/DiracSpinor/FADiracSpinor/g"  $file
sed -i '' -e  "s/DiracSlash/FADiracSlash/g"  $file
sed -i '' -e  "s/DiracMatrix/FADiracMatrix/g"  $file
sed -i '' -e  "s/DiracTrace/FADiracTrace/g"  $file
sed -i '' -e  "s/MetricTensor/FAMetricTensor/g"  $file
sed -i '' -e  "s/ScalarProduct/FAScalarProduct/g"  $file
sed -i '' -e  "s/ChiralityProjector/FAChiralityProjector/g"  $file
sed -i '' -e  "s/GS/FAGS/g"  $file
sed -i '' -e  "s/SUNT/FASUNT/g"  $file
sed -i '' -e  "s/SUNF/FASUNF/g"  $file
sed -i '' -e  "s/LeviCivita/FALeviCivita/g"  $file
sed -i '' -e  "s/Loop/FALoop/g"  $file
sed -i '' -e  "s/NonCommutative/FANonCommutative/g"  $file
sed -i '' -e  "s/PolarizationVector/FAPolarizationVector/g"  $file
sed -i '' -e  "s/FeynAmp/FAFeynAmp/g"  $file
sed -i '' -e  "s/PropagatorDenominator/FAPropagatorDenominator/g"  $file
sed -i '' -e  "s/GaugeXi/FAGaugeXi/g"  $file
sed -i '' -e  "s/FeynAmpDenominator/FAFeynAmpDenominator/g"  $file
#sed -i '' -e  "s/FeynAmpList/FAFeynAmpList/g"  $file
#sed -i '' -e  "s/ME/FCGV[\"ME\"]/g"  $file
#sed -i '' -e  "s/MM/FCGV[\"MM\"]/g"  $file
#sed -i '' -e  "s/ML/FCGV[\"ML\"]/g"  $file
#sed -i '' -e  "s/MU/FCGV[\"MU\"]/g"  $file
#sed -i '' -e  "s/MD/FCGV[\"MD\"]/g"  $file
#sed -i '' -e  "s/MC/FCGV[\"MC\"]/g"  $file
#sed -i '' -e  "s/MS/FCGV[\"MS\"]/g"  $file
#sed -i '' -e  "s/MT/FCGV[\"MT\"]/g"  $file
#sed -i '' -e  "s/MB/FCGV[\"MB\"]/g"  $file
#sed -i '' -e  "s/MH/FCGV[\"MH\"]/g"  $file
#sed -i '' -e  "s/MZ/FCGV[\"MZ\"]/g"  $file
#sed -i '' -e  "s/MW/FCGV[\"MW\"]/g"  $file
#sed -i '' -e  "s/EL/FCGV[\"EL\"]/g"  $file
#sed -i '' -e  "s/CW/FCGV[\"CW\"]/g"  $file
#sed -i '' -e  "s/SW/FCGV[\"SW\"]/g"  $file
