#!/bin/csh 

set rad_name = ( KDDC KFSD KVNX KICT KINX KTWX KGLD KOAX KARX KDVN KSGF KDMX KLNX KUEX KEAX )

set start = "2018,5,01,18"
set stop  = "2018,5,02,03"

foreach radar ($rad_name)
  echo $radar
  get_nexrad.py --start $start --end $stop -r $radar >& out_$radar &
end
