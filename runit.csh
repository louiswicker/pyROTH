#!/bin/csh 

set rad_name = ( KDDC )
set start = "2016,5,23,18"
set stop  = "2016,5,23,35"

foreach radar ($rad_name)
  echo $radar
  get_nexrad.py --start $start --end $stop -r $radar >& out_$radar &
end
