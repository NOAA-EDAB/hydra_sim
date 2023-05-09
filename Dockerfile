FROM ubuntu:18.04

RUN apt-get update && apt-get install -yq build-essential autoconf libnetcdf-dev libxml2-dev libproj-dev valgrind wget unzip git nano dos2unix

# pulls ADBM from github and unzips in folder ADMBcode
RUN mkdir /ADMBcode 
RUN wget https://github.com/admb-project/admb/releases/download/admb-12.2/admb-12.2-linux.zip
RUN mv admb-12.2-linux.zip /ADMBcode
RUN unzip ADMBcode/admb-12.2-linux.zip -d /ADMBcode

# pulls hydra repo from github into folder HYDRA
RUN mkdir /HYDRA
RUN git clone https://github.com/NOAA-EDAB/hydra_sim.git /HYDRA

WORKDIR /HYDRA

# in case of windows stuff
RUN chmod -R 777 /HYDRA

# compile the model
RUN /ADMBcode/admb-12.2/admb hydra_sim.tpl

# create a directory in which all output will be created.
This will be the dst directory for mounted files  
RUN mkdir /HYDRA/mount
RUN mv /HYDRA/hydra_sim.dat /HYDRA/mount/hydra_sim.dat
RUN mv /HYDRA/hydra_sim.pin /HYDRA/mount/hydra_sim.pin
RUN mv /HYDRA/runModel.sh /HYDRA/mount/runModel.sh

WORKDIR /HYDRA/mount

ENTRYPOINT ["sh"]
CMD ["runModel.sh"]
