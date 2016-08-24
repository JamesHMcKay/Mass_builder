#!/bin/sh
cd /Users/jamesmckay/Documents/Programs/Mass_builder/Generator
rm names_updated.txt
rm names_updated_temp.txt

sed -n '/= 0 ;/!p' output.txt > output_tidy.txt

cut -d' ' -f3 output_tidy.txt >> names_updated_temp.txt

cut -c2-20 names_updated_temp.txt >> names_updated.txt


FAoutput="output_tidy.txt"

sed -i '' -e 's/TBI(4,Power(p,2),List(List(1,mc),List(1,ma)))/1.0L*Bac/g' $FAoutput

sed -i '' -e 's/TBI(4,Power(p,2),List(List(1,mz),List(1,mc)))/1.0L*Bcz/g' $FAoutput

sed -i '' -e 's/TBI(4,Power(p,2),List(List(1,mw),List(1,mc)))/1.0L*Bcw/g' $FAoutput

sed -i '' -e 's/TBI(4,Power(p,2),List(List(1,ms),List(1,ms)))/1.0L*Bss/g' $FAoutput

sed -i '' -e 's/TAI(4,0,List(List(1,mw)))/1.0L*Aw/g' $FAoutput

sed -i '' -e 's/TAI(4,0,List(List(1,mz)))/1.0L*Az/g' $FAoutput

sed -i '' -e 's/TAI(4,0,List(List(1,mc)))/1.0L*Ac/g' $FAoutput

sed -i '' -e 's/TAI(4,0,List(List(1,ma)))/1.0L*Aa/g' $FAoutput

sed -i '' -e 's/TAI(4,0,List(List(1,ms)))/1.0L*As/g' $FAoutput






