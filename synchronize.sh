#!/bin/bash
# Read Password
echo -n "Password for itsshift:" 
read -s password

while true
do
	echo -n "$password" | kinit itsshift
	rsync -avruh Plots/ShiftReport24h_THR_* itsshift@lxplus.cern.ch:/eos/user/i/itsshift/Reports24h/THR
	rsync -avruh Plots/ShiftReport24h_FHR_* itsshift@lxplus.cern.ch:/eos/user/i/itsshift/Reports24h/FHR
	echo "Sleeping 1h"
	sleep 1h
done

