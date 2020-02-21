for i in *data_*
do
    mv "$i" "`echo $i | sed 's/data_//'`"
done
