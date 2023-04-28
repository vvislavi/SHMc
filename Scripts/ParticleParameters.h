#ifndef __PARTICLE_PARAMETERS__
#define __PARTICLE_PARAMETERS__
#include "TString.h"
struct Particle {
  TString ppHEPFile; //Link to HEP data file
  TString kernelName; //Particle name in kernel
  TString pLabel; //Particle label
  double xLow; //Lower pT of the spectrum
  double xHigh; //Higher pT of the spectrum
  double Sigma; //Integrated x-section
  double dSigma;//Error on x-section
  double degen; //Degeneracy of particle
  double scaling; //Scaling factor. Relevant e.g. when adding new states for Lc
  double thermalFD; //Feeddown on thermal spectra -- needed to scale gamma for dif. FO surface
  double ratAM; //Ratio for integrated yields by AM
  double ratAA; //Ratio for integrated yields by AA
  int    prefCol=kBlack; //Preferred color when drawing (where available)
  bool operator==(Particle& rhs) { if(kernelName.EqualTo(rhs.kernelName)) return true; return false; };
  TString paperRef;
  vector<TString> spectraPath;
  vector<TString> raaPath;
  vector<TString> ratioPath;
};
/*Template for a particle for copy-paste:
Particle someName = {.ppHEPFile="",
               // .ppHEPTabNo=,
               .kernelName="",
               .pLabel="",
               .xLow=,
               .xHigh=,
               .Sigma=,
               .dSigma=,
               .degen=,
               .scaling=,
               .thermalFD=,
               .ratAM=,
               .ratAA=
};
*/
Particle D0 = {
               .ppHEPFile="pp_Published/HEPData-ins1848990-v1-root.root/Table 4/Graph1D_y1",
               .kernelName="Dc1865zer",
               .pLabel="D^{0}",
               .xLow=0,
               .xHigh=16,
               .Sigma=440,
               .dSigma=18,
               .degen=1,
               .scaling=1,
               .thermalFD=4.115,
               .ratAM=4.115,
               .ratAA=4.088,
               .prefCol=kViolet,
               .paperRef="JHEP 01 (2022) 174",
               .spectraPath={"PbPb_Published/HEPData-ins1946131-v1-root.root/Table 1/Graph1D_y1",
                             "PbPb_Published/HEPData-ins1946131-v1-root.root/Table 2/Graph1D_y1"},
               .raaPath={"PbPb_Published/HEPData-ins1946131-v1-root.root/Table 11/Graph1D_y1",
                             "PbPb_Published/HEPData-ins1946131-v1-root.root/Table 14/Graph1D_y1"},
               .ratioPath={}
};
Particle DStar = {
               .ppHEPFile="pp_Published/HEPData-ins1716440-v1-root.root/Table 3/Graph1D_y1",
               .kernelName="Dc2010plu",
               .pLabel="D^{*+}",
               .xLow=0,
               .xHigh=16,
               .Sigma=178,
               .dSigma=18,
               .degen=3,
               .scaling=1,
               .thermalFD=1.233,
               .ratAM=1.233,
               .ratAA=1.220,
               .prefCol=kGreen+3,
               .paperRef="JHEP 01 (2022) 174",
               .spectraPath={"PbPb_Published/HEPData-ins1946131-v1-root.root/Table 5/Graph1D_y1",
                             "PbPb_Published/HEPData-ins1946131-v1-root.root/Table 6/Graph1D_y1"},
               .raaPath={"PbPb_Published/HEPData-ins1946131-v1-root.root/Table 13/Graph1D_y1",
                             "PbPb_Published/HEPData-ins1946131-v1-root.root/Table 16/Graph1D_y1"},
               .ratioPath={"PbPb_Published/HEPData-ins1946131-v1-root.root/Table 9/Graph1D_y1",
                             "PbPb_Published/HEPData-ins1946131-v1-root.root/Table 10/Graph1D_y1"}

};
Particle Dch = {
               .ppHEPFile="pp_Published/HEPData-ins1848990-v1-root.root/Table 5/Graph1D_y1",
               .kernelName="Dc1869plu",
               .pLabel="D^{+}",
               .xLow=0,
               .xHigh=16,
               .Sigma=195,
               .dSigma=18,
               .degen=1,
               .scaling=1,
               .thermalFD=1.831,
               .ratAM=1.831,
               .ratAA=1.860,
               .prefCol=kBlue+1,
               .paperRef="JHEP 01 (2022) 174",
               .spectraPath={"PbPb_Published/HEPData-ins1946131-v1-root.root/Table 3/Graph1D_y1",
                             "PbPb_Published/HEPData-ins1946131-v1-root.root/Table 4/Graph1D_y1"},
               .raaPath={"PbPb_Published/HEPData-ins1946131-v1-root.root/Table 12/Graph1D_y1",
                             "PbPb_Published/HEPData-ins1946131-v1-root.root/Table 15/Graph1D_y1"},
               .ratioPath={"PbPb_Published/HEPData-ins1946131-v1-root.root/Table 7/Graph1D_y1",
                             "PbPb_Published/HEPData-ins1946131-v1-root.root/Table 8/Graph1D_y1"}

};
Particle Ds = {
               .ppHEPFile="pp_Published/HEPData-ins1848990-v1-root.root/Table 6/Graph1D_y1",
               .kernelName="Ds1968plu",
               .pLabel="D_{s}^{+}",
               .xLow=0,
               .xHigh=16,
               .Sigma=82,
               .dSigma=18,
               .degen=1,
               .scaling=1,
               .thermalFD=2.628,
               .ratAM=2.628,
               .ratAA=2.627,
               .prefCol=kRed+1,
               .paperRef="Phys. Lett. B 827 (2022) 136986",
               .spectraPath={"PbPb_Published/HEPData-ins1946931-v2-root.root/Table 1/Graph1D_y1",
                             "PbPb_Published/HEPData-ins1946931-v2-root.root/Table 2/Graph1D_y1"},
               .raaPath={"PbPb_Published/HEPData-ins1946931-v2-root.root/Table 3/Graph1D_y1",
                             "PbPb_Published/HEPData-ins1946931-v2-root.root/Table 4/Graph1D_y1"},
               .ratioPath={"PbPb_Published/HEPData-ins1946931-v2-root.root/Table 5/Graph1D_y1",
                             "PbPb_Published/HEPData-ins1946931-v2-root.root/Table 6/Graph1D_y1"}
};
Particle Lc = {
               .ppHEPFile="pp_Published/HEPData-ins1829739-v1-root.root/Table 1/Graph1D_y1",
               .kernelName="Lc2286plu",
               .pLabel="#Lambda_{c}",
               .xLow=0,
               .xHigh=12,
               .Sigma=278,
               .dSigma=18,
               .degen=1,
               .scaling=1, //For hidden states, this would be 2.146
               .thermalFD=4.938,
               .ratAM=4.935,
               .ratAA=4.939,
               .prefCol=kBlack,
               .paperRef="JHEP 01 (2022) 174",
               .spectraPath={"PbPb_Published/HEPData-ins1990765-v2-root.root/Table 1a/Graph1D_y1",
                             "PbPb_Published/HEPData-ins1990765-v2-root.root/Table 1a/Graph1D_y2"},
               .raaPath={"PbPb_Published/HEPData-ins1990765-v2-root.root/Table 2/Graph1D_y1",
                             "PbPb_Published/HEPData-ins1990765-v2-root.root/Table 2/Graph1D_y2"},
               .ratioPath={"PbPb_Published/HEPData-ins1990765-v2-root.root/Table 1b/Graph1D_y1",
                             "PbPb_Published/HEPData-ins1990765-v2-root.root/Table 1b/Graph1D_y2"}

};
#endif
