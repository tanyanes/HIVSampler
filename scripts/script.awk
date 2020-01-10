echo "3 RBOX C"
numline=$(wc -l $1)
set -- $numline
echo $1
awk '{ print "     " $1 " " $2 " " $3 }' $2
