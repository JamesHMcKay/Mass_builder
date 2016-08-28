#!/bin/sh
#cd /Users/jamesmckay/Documents/Programs/Mass_builder/build

rm output/names_updated.txt
rm output/names_updated_temp.txt

tag=$(head output/tag.txt)
model=$(head output/model.txt)

sed -n '/= 0 ;/!p' output/output_products.txt > output/output_products_tidy.txt

cut -d' ' -f3 output/output_products_tidy.txt >> output/names_updated_temp.txt

cut -c2-20 output/names_updated_temp.txt >> output/names_updated.txt

cp output/output_products_tidy.txt models/"$model"/output/coeff_products_"$tag".txt
cp output/names_updated.txt output/products.txt

sed -n '/= 0 ;/!p' output/output_products_2.txt > output/output_products_3_tidy.txt
rm output/names_updated.txt

cut -c2-20 output/output_products_3_tidy.txt >> output/names_updated.txt