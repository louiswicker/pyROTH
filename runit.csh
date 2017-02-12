#!/bin/csh 

set rad_name = ( KDDC KTLX KVNX KICT KLZK KINX KAMA KSRX KTWX KFDR KGLD KOAX KSGF KDMX KLNX KUEX KEAX )
set start = "2016,5,23,18"
set stop  = "2016,5,24,05"

foreach radar ($rad_name)
  echo $radar
  getnexrad.py --start $start --end $stop -r $radar >& out_$radar &
end
