#!/bin/bash
find /HYDRA/mount -type f | xargs dos2unix
./hydra_sim -ind /HYDRA/mount/hydra_sim.dat -ainp /HYDRA/mount/hydra_sim.pin