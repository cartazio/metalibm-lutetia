#/bin/bash

unique="$2-$$"
beforedir=`pwd`
cd $1
echo "#include \"expansion.h\"" > "impl-$unique.c"
echo " " >> "impl-$unique.c"
cat "$beforedir/$2" | sed -e 's/void/static inline void/g;' | sed -e "s/static inline void $4(/void $4(/g;" >> "impl-$unique.c"
gcc -Winline -finline-limit=1200 -std=c99 -O2 -D$3 -DPOLYNOMIALNAME=$4 $6 -fPIC  -c "impl-$unique.c" 
gcc -Winline -finline-limit=1200 -std=c99 -D$3 -DPOLYNOMIALNAME=$4 $6 -fPIC  -c expansion.c 
gcc -std=c99 -shared -o "$5" "`basename "impl-$unique.c" .c`.o"  expansion.o -lgmp -lmpfr 
rm "`basename "impl-$unique.c" .c`.o"
rm "impl-$unique.c"
rm expansion.o
cd $beforedir


