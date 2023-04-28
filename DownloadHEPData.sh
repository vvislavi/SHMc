#First, fetch published Pb--Pb
#D0, D*, Dch:
mkdir PbPb_Published && cd PbPb_Published
wget -c https://www.hepdata.net/download/submission/ins1946131/1/root -O - | tar -xv
#Ds:
wget -c https://www.hepdata.net/download/submission/ins1946931/2/root -O - | tar -xv
#Lc
wget -c https://www.hepdata.net/download/submission/ins1990765/2/root -O - | tar -xv
#Then, get pp:
cd .. && mkdir pp_Published && cd pp_Published
#d0,dch,ds
wget -c https://www.hepdata.net/download/submission/ins1848990/1/root -O - | tar -xv
#d*
wget -c https://www.hepdata.net/download/submission/ins1716440/1/root -O - | tar -xv
#Lc
wget -c https://www.hepdata.net/download/submission/ins1829739/1/root -O - | tar -xv
cd ..
