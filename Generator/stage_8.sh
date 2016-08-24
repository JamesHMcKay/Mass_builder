#!/bin/sh
cd /Users/jamesmckay/Documents/Programs/Mass_builder/Generator


tag=$(head tag.txt)
model=$(head model.txt)
if [ ! -e "models/"$model"/output/" ]; then
  mkdir models/"$model"/output
fi

cp output_tidy.txt models/"$model"/output/coeff_integrals_"$tag".txt




rm names_updated.txt
rm names_updated_temp.txt

sed -n '/= 0 ;/!p' output.txt > output_tidy2.txt

cut -d' ' -f3 output_tidy2.txt >> names_updated_temp.txt

cut -c2-20 names_updated_temp.txt >> names_updated.txt