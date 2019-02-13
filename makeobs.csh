#!/bin/csh 

#set rad_name = ( KDDC KFSD KVNX KICT KINX KTWX KGLD KOAX KARX KDVN KSGF KDMX KLNX KUEX KEAX )
#set rad_name = ( KDDC )

set start = "2018,5,29,18"
set stop  = "2018,5,30,03"
set newse_radar_file = "radars.20180529.csh"
set dir = "29May2018"

prep_nexrad.py --start $start --end $stop --newse $newse_radar_file >& out_newse

# Do one radar for a test
#prep_nexrad.py --start $start --end $stop --radar KDDC 

# move stuff into its own dir...

mkdir $dir
mv K* $dir
mv log_* $dir
mv obs_seq* $dir
mv out_* $dir
mv roth* $dir
