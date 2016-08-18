#!/bin/sh
cd /Users/jamesmckay/Documents/Programs/Mass_builder/Generator


tag=$(head tag.txt)

cp names_updated.txt output/integrals.txt
cp output_tidy.txt output/coeff_integrals_"$tag".txt




rm names_updated.txt
rm names_updated_temp.txt

sed -n '/= 0 ;/!p' output.txt > output_tidy2.txt

cut -d' ' -f3 output_tidy2.txt >> names_updated_temp.txt

cut -c2-20 names_updated_temp.txt >> names_updated.txt


cp names_updated.txt output/products.txt