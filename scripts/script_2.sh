#!/bin/sh
#cd /Users/jamesmckay/Documents/Programs/Mass_builder/build


tag=$(head output/tag.txt)
model=$(head output/model.txt)
if [ ! -e "models/"$model"/output/" ]; then
  mkdir models/"$model"/output
fi

cp output/output_tidy_new.txt models/"$model"/output/coeff_integrals_"$tag".txt




rm output/names_updated.txt
rm output/names_updated_temp.txt

sed -n '/= 0 ;/!p' output/output.txt > output/output_tidy2.txt

cut -d' ' -f3 output/output_tidy2.txt >> output/names_updated_temp.txt

cut -c2-20 output/names_updated_temp.txt >> output/names_updated.txt