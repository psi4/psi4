#!/bin/sh

header=CInt.h
file=cint_def.h
olddef=CINT_DEF_H
newdef=__CINT_H__
include=("")
struct=(ERD OED BasisSet)

cp -r ${file} _temp
sed -i '/'"$olddef"'/d' _temp
sed -i '/#include/d' _temp

echo -n "#ifndef " > ${header}
echo ${newdef} >> ${header}
echo -n "#define " >> ${header}
echo ${newdef} >> ${header}
echo "" >> ${header}
echo "" >> ${header}

for str in "${include[@]}"
do
    echo ${include} >> ${header}
done
echo "" >> ${header}
echo "" >> ${header}

for str in "${struct[@]}"
do
    echo -n "struct " >> ${header}
    echo -n ${str} >> ${header}
    echo ";" >> ${header}
done

cat _temp >> ${header}

echo -n "#endif /* " >> ${header}
echo -n ${newdef} >> ${header}
echo " */" >> ${header}

rm -r _temp
cp -r ${header} ../include/
rm -r ${header}
