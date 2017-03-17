//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.2.2, 2014-11-06
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================


#ifndef CPPProcess_DY_gg_to_ddbar_cc
#define CPPProcess_DY_gg_to_ddbar_cc

#include "MEM/MEMAlgo/interface/CPPProcess_DY_gg_to_ddbar.h"
#include "MEM/MEMAlgo/interface/HelAmps_sm.h"

using namespace MG5_sm; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > d d~ ta- ta+ WEIGHTED=6
// Process: g g > s s~ ta- ta+ WEIGHTED=6

//--------------------------------------------------------------------------
// Initialize process.

void CPPProcess_DY_gg_to_ddbar::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_sm::getInstance(); 
  SLHAReader slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings(); 
  pars->printIndependentParameters(); 
  pars->printIndependentCouplings(); 
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->mdl_MTA); 
  mME.push_back(pars->mdl_MTA); 
  jamp2[0] = new double[2]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void CPPProcess_DY_gg_to_ddbar::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 
  static bool firsttime = true; 
  if (firsttime)
  {
    pars->printDependentParameters(); 
    pars->printDependentCouplings(); 
    firsttime = false; 
  }

  // Reset color flows
  for(int i = 0; i < 2; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 64; 
  static bool goodhel[64] = {64 * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[64]; 
  static int jhel; 
  //std::complex<double> * * wfs; 
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[64][6] = {{-1, -1, -1, -1, -1, -1},
      {-1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, 1, -1}, {-1, -1, -1, -1, 1, 1},
      {-1, -1, -1, 1, -1, -1}, {-1, -1, -1, 1, -1, 1}, {-1, -1, -1, 1, 1, -1},
      {-1, -1, -1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, -1, 1},
      {-1, -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1, 1}, {-1, -1, 1, 1, -1, -1},
      {-1, -1, 1, 1, -1, 1}, {-1, -1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1}, {-1,
      1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1}, {-1,
      1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, -1}, {-1, 1, -1, 1, -1, 1}, {-1, 1,
      -1, 1, 1, -1}, {-1, 1, -1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1}, {-1, 1, 1,
      -1, -1, 1}, {-1, 1, 1, -1, 1, -1}, {-1, 1, 1, -1, 1, 1}, {-1, 1, 1, 1,
      -1, -1}, {-1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, -1}, {-1, 1, 1, 1, 1, 1},
      {1, -1, -1, -1, -1, -1}, {1, -1, -1, -1, -1, 1}, {1, -1, -1, -1, 1, -1},
      {1, -1, -1, -1, 1, 1}, {1, -1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1}, {1,
      -1, -1, 1, 1, -1}, {1, -1, -1, 1, 1, 1}, {1, -1, 1, -1, -1, -1}, {1, -1,
      1, -1, -1, 1}, {1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, 1}, {1, -1, 1, 1,
      -1, -1}, {1, -1, 1, 1, -1, 1}, {1, -1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1},
      {1, 1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, 1}, {1, 1, -1, -1, 1, -1}, {1,
      1, -1, -1, 1, 1}, {1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1}, {1, 1, -1,
      1, 1, -1}, {1, 1, -1, 1, 1, 1}, {1, 1, 1, -1, -1, -1}, {1, 1, 1, -1, -1,
      1}, {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, 1}, {1, 1, 1, 1, -1, -1}, {1,
      1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {256}; 

  ntry = ntry + 1; 

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i; 
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]); 
        t[0] = matrix_gg_ddxtamtap(); 

        double tsum = 0; 
        for(int iproc = 0; iproc < nprocesses; iproc++ )
        {
          matrix_element[iproc] += t[iproc]; 
          tsum += t[iproc]; 
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel])
        {
          goodhel[ihel] = true; 
          ngood++; 
          igood[ngood] = ihel; 
        }
      }
    }
    jhel = 0; 
    sum_hel = min(sum_hel, ngood); 
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0; j < sum_hel; j++ )
    {
      jhel++; 
      if (jhel >= ngood)
        jhel = 0; 
      double hwgt = double(ngood)/double(sum_hel); 
      int ihel = igood[jhel]; 
      calculate_wavefunctions(perm, helicities[ihel]); 
      t[0] = matrix_gg_ddxtamtap(); 

      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i]; 



}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double CPPProcess_DY_gg_to_ddbar::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 21 && id2 == 21)
  {
    // Add matrix elements for processes with beams (21, 21)
    return matrix_element[0] * 2; 
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.; 
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void CPPProcess_DY_gg_to_ddbar::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  //int i, j; 

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]); 
  vxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
  oxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  ixxxxx(p[perm[5]], mME[5], hel[5], -1, w[5]); 
  VVV1P0_1(w[0], w[1], pars->GC_10, pars->ZERO, pars->ZERO, w[6]); 
  FFV1P0_3(w[5], w[4], pars->GC_3, pars->ZERO, pars->ZERO, w[7]); 
  FFV1_1(w[2], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[8]); 
  FFV1_2(w[3], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[9]); 
  FFV2_4_3(w[5], w[4], pars->GC_50, pars->GC_59, pars->mdl_MZ, pars->mdl_WZ,
      w[10]);
  FFV1_1(w[2], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[11]); 
  FFV1_2(w[3], w[1], pars->GC_11, pars->ZERO, pars->ZERO, w[12]); 
  FFV1_1(w[11], w[1], pars->GC_11, pars->ZERO, pars->ZERO, w[13]); 
  FFV1_2(w[3], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[14]); 
  FFV1_1(w[2], w[1], pars->GC_11, pars->ZERO, pars->ZERO, w[15]); 
  FFV1_2(w[14], w[1], pars->GC_11, pars->ZERO, pars->ZERO, w[16]); 
  FFV1_1(w[15], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[17]); 
  FFV1_2(w[12], w[0], pars->GC_11, pars->ZERO, pars->ZERO, w[18]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[3], w[8], w[7], pars->GC_1, amp[0]); 
  FFV1_0(w[9], w[2], w[7], pars->GC_1, amp[1]); 
  FFV2_3_0(w[3], w[8], w[10], pars->GC_50, pars->GC_58, amp[2]); 
  FFV2_3_0(w[9], w[2], w[10], pars->GC_50, pars->GC_58, amp[3]); 
  FFV1_0(w[12], w[11], w[7], pars->GC_1, amp[4]); 
  FFV2_3_0(w[12], w[11], w[10], pars->GC_50, pars->GC_58, amp[5]); 
  FFV1_0(w[3], w[13], w[7], pars->GC_1, amp[6]); 
  FFV2_3_0(w[3], w[13], w[10], pars->GC_50, pars->GC_58, amp[7]); 
  FFV1_0(w[14], w[15], w[7], pars->GC_1, amp[8]); 
  FFV2_3_0(w[14], w[15], w[10], pars->GC_50, pars->GC_58, amp[9]); 
  FFV1_0(w[16], w[2], w[7], pars->GC_1, amp[10]); 
  FFV2_3_0(w[16], w[2], w[10], pars->GC_50, pars->GC_58, amp[11]); 
  FFV1_0(w[3], w[17], w[7], pars->GC_1, amp[12]); 
  FFV2_3_0(w[3], w[17], w[10], pars->GC_50, pars->GC_58, amp[13]); 
  FFV1_0(w[18], w[2], w[7], pars->GC_1, amp[14]); 
  FFV2_3_0(w[18], w[2], w[10], pars->GC_50, pars->GC_58, amp[15]); 

}
double CPPProcess_DY_gg_to_ddbar::matrix_gg_ddxtamtap() 
{
  int i, j; 
  // Local variables
  //const int ngraphs = 16; 
  const int ncolor = 2; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[2] = {3, 3}; 
  static const double cf[2][2] = {{16, -2}, {-2, 16}}; 

  // Calculate color flows
  jamp[0] = -std::complex<double> (0, 1) * amp[0] - std::complex<double> (0, 1)
      * amp[1] - std::complex<double> (0, 1) * amp[2] - std::complex<double>
      (0, 1) * amp[3] + amp[4] + amp[5] + amp[6] + amp[7] + amp[14] + amp[15];
  jamp[1] = +std::complex<double> (0, 1) * amp[0] + std::complex<double> (0, 1)
      * amp[1] + std::complex<double> (0, 1) * amp[2] + std::complex<double>
      (0, 1) * amp[3] + amp[8] + amp[9] + amp[10] + amp[11] + amp[12] + amp[13];

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return matrix; 
}



#endif
