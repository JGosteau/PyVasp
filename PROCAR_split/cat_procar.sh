file1=PROCAR_1
file2=PROCAR_2
filef=PROCAR

file1=$1
file2=$2
filef=$3

l_kpt=`grep "k-point" $file1 | awk '{val = $2}END{print val}'`

awk 'NR>2 {print $0}' $file1 > $file1"_tmp"

nkpt=`awk -v l_kpt=$l_kpt -v file=$file2"_tmp" 'BEGIN{i=l_kpt+1; print "" > file} NR>2 && $1=="k-point"{print " "$1" "i" "$3"    "$4" "$5" "$6"     "$7" "$8" "$9 > file; i++}NR>3 && $1 != "k-point"{print $0 > file}END{print i-1}' $file2 `

head -2 $file1 | awk -v nkpt=$nkpt 'NR==1{print $0}NR==2{print $1" "$2" "$3"  "nkpt"       "$5" "$6" "$7"   "$8"       "$9" "$10" "$11"    "$12}' > $filef"_ini"

cat $filef"_ini" $file1"_tmp" $file2"_tmp" > $filef
rm $file1"_tmp" $file2"_tmp" $filef"_ini"
