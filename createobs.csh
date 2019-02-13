#!/bin/csh 

set rad_name = ( KDDC KFSD KVNX KICT KINX KTWX KGLD KOAX KARX KDVN KSGF KDMX KLNX KUEX KEAX )

set start = "2018,5,01,19,57,00"
set stop  = "2018,5,01,20,05,00"

foreach radar ($rad_name)
  echo $radar
  pyDart.py -f "obs_seq_2018_05_01_"$radar"_VR.h5" --start $start --end $end --search --hdf2ascii 
end
