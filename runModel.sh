#!/bin/bash
#chmod -R 777 /HYDRA/mount
find /HYDRA/mount -type f | xargs dos2unix

arg1=${1}
arg2=${2}
arg3=${3}

/HYDRA/hydra_sim -ind /HYDRA/mount/$arg1 -ainp /HYDRA/mount/$arg2 -sim $arg3