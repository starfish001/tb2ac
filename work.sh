#less data1.txt_backup |perl -p -i -e "s/^\s+//g;s/-//g" |egrep "chr|genome" -v|perl -p -i -e "s/^/1\t/g;s/ +/\t/g" >data1.txt
./tb2ac -x 1000 -f data1.txt -m 2
