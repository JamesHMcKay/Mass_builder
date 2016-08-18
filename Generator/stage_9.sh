#!/bin/sh
cd /Users/jamesmckay/Documents/Programs/Mass_builder/Generator

rm names_updated.txt
rm names_updated_temp.txt
rm output_products_2_tidy.txt

tag=$(head tag.txt)

sed -n '/= 0 ;/!p' output_products.txt > output_products_tidy.txt

cut -d' ' -f3 output_products_tidy.txt >> names_updated_temp.txt

cut -c2-20 names_updated_temp.txt >> names_updated.txt

cp output_products_tidy.txt output/coeff_products_"$tag".txt
cp names_updated.txt output/products.txt

sed -n '/= 0 ;/!p' output_products_2.txt > output_products_2_tidy.txt
rm names_updated.txt

cut -c2-20 output_products_2_tidy.txt >> names_updated.txt