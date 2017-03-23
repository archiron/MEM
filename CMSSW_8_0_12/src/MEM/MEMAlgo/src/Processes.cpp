//# include "Config/PyConfig.h"

# include "MEM/MEMAlgo/interface/LHAPDF.h"

# include "MEM/MEMAlgo/interface/Processes.h"

#include "MEM/MEMAlgo/interface/CPPProcess_ttH_gg.h"
#include "MEM/MEMAlgo/interface/CPPProcess_ttH_qqbar.h"

#include "MEM/MEMAlgo/interface/CPPProcess_ttZ_gg.h"
#include "MEM/MEMAlgo/interface/CPPProcess_ttZ_ddbar.h"
#include "MEM/MEMAlgo/interface/CPPProcess_ttZ_uubar.h"

#include "MEM/MEMAlgo/interface/CPPProcess_ttZ_Zonly_gg.h"
#include "MEM/MEMAlgo/interface/CPPProcess_ttZ_Zonly_uubar.h"

#include "MEM/MEMAlgo/interface/CPPProcess_ttZ_Zonly_Zll_gg.h"
#include "MEM/MEMAlgo/interface/CPPProcess_ttZ_Zonly_Zll_uubar.h"

#include "MEM/MEMAlgo/interface/CPPProcess_dubar_to_ttbarWm.h"

#include "MEM/MEMAlgo/interface/CPPProcess_gd_to_ttbaruWm.h"
#include "MEM/MEMAlgo/interface/CPPProcess_gubar_to_ttbardbarWm.h"
#include "MEM/MEMAlgo/interface/CPPProcess_dubar_to_ttbargWm.h"

#include "MEM/MEMAlgo/interface/CPPProcess_gg_to_ttbarg.h"
#include "MEM/MEMAlgo/interface/CPPProcess_gu_to_ttbaru.h"
#include "MEM/MEMAlgo/interface/CPPProcess_uubar_to_ttbarg.h"

#include "MEM/MEMAlgo/interface/CPPProcess_gg_to_ttbar.h"
#include "MEM/MEMAlgo/interface/CPPProcess_uubar_to_ttbar.h"

#include "MEM/MEMAlgo/interface/Constants.h"


//static double pdfQ_;

// Main ttH processes
static CPPProcess_ttH_gg process_ttH_gg;
static CPPProcess_ttH_qqbar process_ttH_qqbar;

// Main ttZ processes
static CPPProcess_ttZ_gg process_ttZ_gg;
static CPPProcess_ttZ_ddbar process_ttZ_ddbar;
static CPPProcess_ttZ_uubar process_ttZ_uubar;

static CPPProcess_ttZ_Zonly_gg process_ttZ_Zonly_gg;
static CPPProcess_ttZ_Zonly_uubar process_ttZ_Zonly_uubar;

static CPPProcess_ttZ_Zonly_gg process_ttZ_Zonly_Zll_gg;
static CPPProcess_ttZ_Zonly_uubar process_ttZ_Zonly_Zll_uubar;

// Main ttW processes
static CPPProcess_dubar_to_ttbarWm process_ttWm_dubar;

// Main ttWjets processes
static CPPProcess_gd_to_ttbaruWm process_ttWmJets_gd;
static CPPProcess_gubar_to_ttbardbarWm process_ttWmJets_gubar;
static CPPProcess_dubar_to_ttbargWm process_ttWmJets_dubar;

// Main ttjets processes
static CPPProcess_gg_to_ttbarg process_gg_ttg;
static CPPProcess_gu_to_ttbaru process_gq_ttq;
static CPPProcess_uubar_to_ttbarg process_qqbar_ttg;

// Main ttbar processes
static CPPProcess_gg_to_ttbar process_gg_ttbar;
static CPPProcess_uubar_to_ttbar process_qqbar_ttbar;


// Init with python configuration file
//void init_ttH_HTauTauMEProcesses( const char *configname ) {
void init_ttH_HTauTauMEProcesses( RunConfig runConfig ) {

//  std::cout << "configname : "  << runConfig.configName_ << std::endl;
//
//  std::cout << "MGParamCardFile : " << runConfig.MGParamCardFileName_ << std::endl;
//  std::cout << "Q_ : " << runConfig.Q_ << std::endl;

  //init_ttH_HTauTauMEProcesses( paramFileName.c_str(), _Q );
  //std::cout << "init_ttH_HTauTauMEProcesses( RunConfig runConfig ) : appel init_ttH_HTauTauMEProcesses" << std::endl; // TEMPORAIRE
  init_ttH_HTauTauMEProcesses( runConfig.MGParamCardFileName_, runConfig.Q_); 
  
/*  Py_XDECREF( pVariable );
  closePyConfig( pModule );*/
}

//void init_ttH_HTauTauMEProcesses( const char *_fileName, double pdfQ ) {
void init_ttH_HTauTauMEProcesses( string fileName, double pdfQ ) {
    //std::cout << "init_ttH_HTauTauMEProcesses" << std::endl; // TEMPORAIRE
  
/*  std::cout << "_fileName : "  << _fileName << std::endl;
  pdfQ_ = pdfQ;
  string fileName( _fileName );*/
  
//  std::cout << "process_ttH_gg.initProc(" << _fileName  << ") " << std::endl;
  std::cout << "process_ttH_gg.initProc(" << fileName  << ") " << std::endl;
  // ttH
  process_ttH_gg.initProc(fileName);
  process_ttH_qqbar.initProc(fileName);

  // ttZ
  process_ttZ_gg.initProc(fileName);
  process_ttZ_ddbar.initProc(fileName);
  process_ttZ_uubar.initProc(fileName);

  process_ttZ_Zonly_gg.initProc(fileName);
  process_ttZ_Zonly_uubar.initProc(fileName);

  process_ttZ_Zonly_Zll_gg.initProc(fileName);
  process_ttZ_Zonly_Zll_uubar.initProc(fileName);
  
  
  //ttW
  process_ttWm_dubar.initProc(fileName);

  // ttWjets
  process_ttWmJets_gd.initProc(fileName);
  process_ttWmJets_gubar.initProc(fileName);
  process_ttWmJets_dubar.initProc(fileName);

  // ttjets
  process_gg_ttg.initProc(fileName);
  process_gq_ttq.initProc(fileName);
  process_qqbar_ttg.initProc(fileName);

  // ttbar
  process_gg_ttbar.initProc(fileName);
  process_qqbar_ttbar.initProc(fileName);

}


// ttH
double getME_ttH_gg( vector<double*> &p ) {
  process_ttH_gg.setMomenta( p );
  process_ttH_gg.sigmaKin();
  return
  process_ttH_gg.getMatrixElements()[0];
}

double getME_ttH_qqbar( vector<double*> &p ) {
  process_ttH_qqbar.setMomenta( p );
  process_ttH_qqbar.sigmaKin();
  return
  process_ttH_qqbar.getMatrixElements()[0];
}

// ttZ
double getME_ttZ_gg( vector<double*> &p ) {
  process_ttZ_gg.setMomenta( p );
  process_ttZ_gg.sigmaKin();
  return
  process_ttZ_gg.getMatrixElements()[0];
}

double getME_ttZ_ddbar( vector<double*> &p ) {
  process_ttZ_ddbar.setMomenta( p );
  process_ttZ_ddbar.sigmaKin();
  return
  process_ttZ_ddbar.getMatrixElements()[0];
}

double getME_ttZ_uubar( vector<double*> &p ) {
  process_ttZ_uubar.setMomenta( p );
  process_ttZ_uubar.sigmaKin();
  return
  process_ttZ_uubar.getMatrixElements()[0];
}

double getME_ttZ_Zonly_gg( vector<double*> &p ) {
  process_ttZ_Zonly_gg.setMomenta( p );
  process_ttZ_Zonly_gg.sigmaKin();
  return
  process_ttZ_Zonly_gg.getMatrixElements()[0];
}

double getME_ttZ_Zonly_uubar( vector<double*> &p ) {
  process_ttZ_Zonly_uubar.setMomenta( p );
  process_ttZ_Zonly_uubar.sigmaKin();
  return
  process_ttZ_Zonly_uubar.getMatrixElements()[0];
}

double getME_ttZ_Zonly_Zll_gg( vector<double*> &p ) {
  process_ttZ_Zonly_Zll_gg.setMomenta( p );
  process_ttZ_Zonly_Zll_gg.sigmaKin();
  return
  process_ttZ_Zonly_Zll_gg.getMatrixElements()[0];
}

double getME_ttZ_Zonly_Zll_uubar( vector<double*> &p ) {
  process_ttZ_Zonly_Zll_uubar.setMomenta( p );
  process_ttZ_Zonly_Zll_uubar.sigmaKin();
  return
  process_ttZ_Zonly_Zll_uubar.getMatrixElements()[0];
}


// ttW
double getME_ttWm_dubar( vector<double*> &p ) {
  process_ttWm_dubar.setMomenta( p );
  process_ttWm_dubar.sigmaKin();
  return
  process_ttWm_dubar.getMatrixElements()[0];
}


// ttWjets
double getME_ttWmJets_gd( vector<double*> &p ) {
  process_ttWmJets_gd.setMomenta( p );
  process_ttWmJets_gd.sigmaKin();
  return
  process_ttWmJets_gd.getMatrixElements()[0];
}

double getME_ttWmJets_gubar( vector<double*> &p ) {
  process_ttWmJets_gubar.setMomenta( p );
  process_ttWmJets_gubar.sigmaKin();
  return
  process_ttWmJets_gubar.getMatrixElements()[0];
}

double getME_ttWmJets_dubar( vector<double*> &p ) {
  process_ttWmJets_dubar.setMomenta( p );
  process_ttWmJets_dubar.sigmaKin();
  return
  process_ttWmJets_dubar.getMatrixElements()[0];
}

// ttjets

double getME_gg_ttg( vector<double*> &p ) {
  process_gg_ttg.setMomenta( p );
  process_gg_ttg.sigmaKin();
  return
  process_gg_ttg.getMatrixElements()[0];
}

double getME_gq_ttq( vector<double*> &p ) {
  process_gq_ttq.setMomenta( p );
  process_gq_ttq.sigmaKin();
  return
  process_gq_ttq.getMatrixElements()[0];
}

double getME_qqbar_ttg( vector<double*> &p ) {
  process_qqbar_ttg.setMomenta( p );
  process_qqbar_ttg.sigmaKin();
  return
  process_qqbar_ttg.getMatrixElements()[0];
}


// ttbar

double getME_gg_ttbar( vector<double*> &p ) {
  process_gg_ttbar.setMomenta( p );
  process_gg_ttbar.sigmaKin();
  return
  process_gg_ttbar.getMatrixElements()[0];
}

double getME_qqbar_ttbar( vector<double*> &p ) {
  process_qqbar_ttbar.setMomenta( p );
  process_qqbar_ttbar.sigmaKin();
  return
  process_qqbar_ttbar.getMatrixElements()[0];
}


double get_ttHWeightedME( double x0, double x1, const vector<double*> &p) {
  
  double pdfQ_ = Physics::mtop + 0.5*Physics::mHiggs; 
  double current_ME2, pdf_sum;
  double ME2 = 0;
  
  // The 2 permutations
  vector<double*> p0=p;
  vector<double*> p1=p;  
  p1[0]=p[1];
  p1[1]=p[0];

  
  // Parton codes:
  //    t  b  c  s  u  d  g  d  u  s  c  b  t
  //   -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
  double fg0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 );
  double fu0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 );
  double fd0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 );
  double fc0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 );
  double fs0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 );
  double fubar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 );
  double fdbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 );
  double fcbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 );
  double fsbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 );

  double fg1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  0 );
  double fu1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  2 );
  double fd1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  1 );
  double fc1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  4 );
  double fs1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  3 );
  double fubar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -2 );
  double fdbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -1 );
  double fcbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -4 );
  double fsbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -3 );
  
  //gg->ttH
  current_ME2 = getME_ttH_gg( p0 );
  pdf_sum     = fg0*fg1;
  ME2        += pdf_sum * current_ME2;

  //qqbar->ttH
  current_ME2 = getME_ttH_qqbar( p0 );
  pdf_sum     = fu0*fubar1+fd0*fdbar1+fc0*fcbar1+fs0*fsbar1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_ttH_qqbar( p1 );
  pdf_sum     = fu1*fubar0+fd1*fdbar0+fc1*fcbar0+fs1*fsbar0;
  ME2        += pdf_sum * current_ME2;

  return ME2 / (x0 * x1);

}




double get_ttZWeightedME( double x0, double x1, const vector<double*> &p) {
  
  double pdfQ_ = Physics::mtop + 0.5*Physics::mZ; 
  double current_ME2, pdf_sum;
  double ME2 = 0;
  
  // The 2 permutations
  vector<double*> p0=p;
  vector<double*> p1=p;  
  p1[0]=p[1];
  p1[1]=p[0];

  
  // Parton codes:
  //    t  b  c  s  u  d  g  d  u  s  c  b  t
  //   -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
  double fg0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 );
  double fu0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 );
  double fd0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 );
  double fc0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 );
  double fs0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 );
  double fubar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 );
  double fdbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 );
  double fcbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 );
  double fsbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 );

  double fg1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  0 );
  double fu1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  2 );
  double fd1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  1 );
  double fc1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  4 );
  double fs1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  3 );
  double fubar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -2 );
  double fdbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -1 );
  double fcbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -4 );
  double fsbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -3 );

  //gg->ttZ
  current_ME2 = getME_ttZ_gg( p0 );
  pdf_sum     = fg0*fg1;
  ME2        += pdf_sum * current_ME2;

  //uubar->ttZ
  current_ME2 = getME_ttZ_uubar( p0 );
  pdf_sum     = fu0*fubar1+fc0*fcbar1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_ttZ_uubar( p1 );
  pdf_sum     = fu1*fubar0+fc1*fcbar0;
  ME2        += pdf_sum * current_ME2; 

  //ddbar->ttZ
  current_ME2 = getME_ttZ_ddbar( p0 );
  pdf_sum     = fd0*fdbar1+fs0*fsbar1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_ttZ_ddbar( p1 );
  pdf_sum     = fd1*fdbar0+fs1*fsbar0;
  ME2        += pdf_sum * current_ME2; 
 

  return ME2 / (x0 * x1);

}






double get_ttZZonlyWeightedME( double x0, double x1, const vector<double*> &p) {
  
  double pdfQ_ = Physics::mtop + 0.5*Physics::mZ; 
  double current_ME2, pdf_sum;
  double ME2 = 0;
  
  // The 2 permutations
  vector<double*> p0=p;
  vector<double*> p1=p;  
  p1[0]=p[1];
  p1[1]=p[0];

  
  // Parton codes:
  //    t  b  c  s  u  d  g  d  u  s  c  b  t
  //   -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
  double fg0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 );
  double fu0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 );
  double fd0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 );
  double fc0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 );
  double fs0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 );
  double fubar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 );
  double fdbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 );
  double fcbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 );
  double fsbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 );

  double fg1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  0 );
  double fu1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  2 );
  double fd1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  1 );
  double fc1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  4 );
  double fs1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  3 );
  double fubar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -2 );
  double fdbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -1 );
  double fcbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -4 );
  double fsbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -3 );

  //gg->ttZ
  current_ME2 = getME_ttZ_Zonly_gg( p0 );
  pdf_sum     = fg0*fg1;
  ME2        += pdf_sum * current_ME2;

  //qqbar->ttZ
  current_ME2 = getME_ttZ_Zonly_uubar( p0 );
  pdf_sum     = fu0*fubar1+fc0*fcbar1+fd0*fdbar1+fs0*fsbar1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_ttZ_Zonly_uubar( p1 );
  pdf_sum     = fu1*fubar0+fc1*fcbar0+fd1*fdbar0+fs1*fsbar0;
  ME2        += pdf_sum * current_ME2; 

  return ME2 / (x0 * x1);

}






double get_ttZZonlyZllWeightedME( double x0, double x1, const vector<double*> &p) {
  
  double pdfQ_ = Physics::mtop + 0.5*Physics::mZ; 
  double current_ME2, pdf_sum;
  double ME2 = 0;
  
  // The 2 permutations
  vector<double*> p0=p;
  vector<double*> p1=p;  
  p1[0]=p[1];
  p1[1]=p[0];

  
  // Parton codes:
  //    t  b  c  s  u  d  g  d  u  s  c  b  t
  //   -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
  double fg0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 );
  double fu0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 );
  double fd0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 );
  double fc0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 );
  double fs0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 );
  double fubar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 );
  double fdbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 );
  double fcbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 );
  double fsbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 );

  double fg1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  0 );
  double fu1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  2 );
  double fd1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  1 );
  double fc1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  4 );
  double fs1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  3 );
  double fubar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -2 );
  double fdbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -1 );
  double fcbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -4 );
  double fsbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -3 );

  std::cout << "PDFSet=" << PDFSet << ", x0=" << x0 << ", pdfQ_=" << pdfQ_ << std::endl;
  std::cout << "fg0 :" << LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 ) << std::endl;
  std::cout << "fu0 :" << LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 ) << std::endl;
  std::cout << "fd0 :" << LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 ) << std::endl;
  std::cout << "fc0 :" << LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 ) << std::endl;
  std::cout << "fs0 :" << LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 ) << std::endl;
  std::cout << "fubar0 :" << LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 ) << std::endl;
  std::cout << "fdbar0 :" << LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 ) << std::endl;
  std::cout << "fcbar0 :" << LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 ) << std::endl;
  std::cout << "fsbar0 :" << LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 ) << std::endl;

  //gg->ttZ
  current_ME2 = getME_ttZ_Zonly_Zll_gg( p0 );
  pdf_sum     = fg0*fg1;
  ME2        += pdf_sum * current_ME2;

  //qqbar->ttZ
  current_ME2 = getME_ttZ_Zonly_Zll_uubar( p0 );
  pdf_sum     = fu0*fubar1+fc0*fcbar1+fd0*fdbar1+fs0*fsbar1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_ttZ_Zonly_Zll_uubar( p1 );
  pdf_sum     = fu1*fubar0+fc1*fcbar0+fd1*fdbar0+fs1*fsbar0;
  ME2        += pdf_sum * current_ME2; 

  return ME2 / (x0 * x1);

}






double get_ttWmWeightedME( double x0, double x1, const vector<double*> &p) {
  
  double pdfQ_ = Physics::mtop + 0.5*Physics::mW; 
  double current_ME2, pdf_sum;
  double ME2 = 0;
  
  // The 2 permutations
  vector<double*> p0=p;
  vector<double*> p1=p;  
  p1[0]=p[1];
  p1[1]=p[0];

  
  // Parton codes:
  //    t  b  c  s  u  d  g  d  u  s  c  b  t
  //   -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
//  double fg0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 );
//  double fu0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 );
  double fd0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 );
//  double fc0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 );
  double fs0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 );
  double fubar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 );
//  double fdbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 );
  double fcbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 );
//  double fsbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 );

//  double fg1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  0 );
//  double fu1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  2 );
  double fd1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  1 );
//  double fc1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  4 );
  double fs1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  3 );
  double fubar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -2 );
//  double fdbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -1 );
  double fcbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -4 );
//  double fsbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -3 );


  //dubar->ttWm, scbar->ttWm
  current_ME2 = getME_ttWm_dubar( p0 );
  pdf_sum     = fd0*fubar1+fs0*fcbar1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_ttWm_dubar( p1 );
  pdf_sum     = fd1*fubar0+fs1*fcbar0;
  ME2        += pdf_sum * current_ME2; 
 

  return ME2 / (x0 * x1);

}





double get_ttWpWeightedME( double x0, double x1, const vector<double*> &p) {
  
  double pdfQ_ = Physics::mtop + 0.5*Physics::mW; 
  double current_ME2, pdf_sum;
  double ME2 = 0;
  
  // The 2 permutations
  vector<double*> p0=p;
  vector<double*> p1=p;  
  p1[0]=p[1];
  p1[1]=p[0];

  
  // Parton codes:
  //    t  b  c  s  u  d  g  d  u  s  c  b  t
  //   -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
//  double fg0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 );
  double fu0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 );
//  double fd0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 );
  double fc0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 );
//  double fs0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 );
//  double fubar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 );
  double fdbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 );
//  double fcbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 );
  double fsbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 );

//  double fg1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  0 );
  double fu1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  2 );
//  double fd1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  1 );
  double fc1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  4 );
//  double fs1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  3 );
//  double fubar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -2 );
  double fdbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -1 );
//  double fcbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -4 );
  double fsbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -3 );


  //dbaru->ttWp, sbarc->ttWp
  current_ME2 = getME_ttWm_dubar( p0 );
  pdf_sum     = fdbar0*fu1+fsbar0*fc1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_ttWm_dubar( p1 );
  pdf_sum     = fdbar1*fu0+fsbar1*fc0;
  ME2        += pdf_sum * current_ME2; 
 

  return ME2 / (x0 * x1);

}




double get_ttWmJetsWeightedME( double x0, double x1, const vector<double*> &p) {
  
  double pdfQ_ = Physics::mtop + 0.5*Physics::mW; 
  double current_ME2, pdf_sum;
  double ME2 = 0;
  
  // The 2 permutations
  vector<double*> p0=p;
  vector<double*> p1=p;  
  p1[0]=p[1];
  p1[1]=p[0];

  
  // Parton codes:
  //    t  b  c  s  u  d  g  d  u  s  c  b  t
  //   -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
  double fg0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 );
//  double fu0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 );
  double fd0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 );
//  double fc0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 );
  double fs0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 );
  double fubar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 );
//  double fdbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 );
  double fcbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 );
//  double fsbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 );

  double fg1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  0 );
//  double fu1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  2 );
  double fd1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  1 );
//  double fc1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  4 );
  double fs1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  3 );
  double fubar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -2 );
//  double fdbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -1 );
  double fcbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -4 );
//  double fsbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -3 );

  //gd->ttuWm, gs->ttcWm
  current_ME2 = getME_ttWmJets_gd( p0 );
  pdf_sum     = fg0*fd1 + fg0*fs1;
  ME2        += pdf_sum * current_ME2;

  current_ME2 = getME_ttWmJets_gd( p1 );
  pdf_sum     = fg1*fd0 + fg1*fs0;
  ME2        += pdf_sum * current_ME2;

  //gubar->ttdbarWm, gcbar->ttsbarWm
  current_ME2 = getME_ttWmJets_gubar( p0 );
  pdf_sum     = fg0*fubar1+fg0*fcbar1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_ttWmJets_gubar( p1 );
  pdf_sum     = fg1*fubar0+fg1*fcbar0;
  ME2        += pdf_sum * current_ME2; 

  //dubar->ttgWm, scbar->ttgWm
  current_ME2 = getME_ttWmJets_dubar( p0 );
  pdf_sum     = fd0*fubar1+fs0*fcbar1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_ttWmJets_dubar( p1 );
  pdf_sum     = fd1*fubar0+fs1*fcbar0;
  ME2        += pdf_sum * current_ME2; 
 

  return ME2 / (x0 * x1);

}





double get_ttWpJetsWeightedME( double x0, double x1, const vector<double*> &p) {
  
  double pdfQ_ = Physics::mtop + 0.5*Physics::mW; 
  double current_ME2, pdf_sum;
  double ME2 = 0;
  
  // The 2 permutations
  vector<double*> p0=p;
  vector<double*> p1=p;  
  p1[0]=p[1];
  p1[1]=p[0];

  
  // Parton codes:
  //    t  b  c  s  u  d  g  d  u  s  c  b  t
  //   -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
  double fg0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 );
  double fu0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 );
//  double fd0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 );
  double fc0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 );
//  double fs0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 );
//  double fubar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 );
  double fdbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 );
//  double fcbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 );
  double fsbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 );

  double fg1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  0 );
  double fu1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  2 );
//  double fd1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  1 );
  double fc1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  4 );
//  double fs1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  3 );
//  double fubar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -2 );
  double fdbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -1 );
//  double fcbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -4 );
  double fsbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -3 );

  //gdbar->ttubarWp, gsbar->ttcbarWp
  current_ME2 = getME_ttWmJets_gd( p0 );
  pdf_sum     = fg0*fdbar1 + fg0*fsbar1;
  ME2        += pdf_sum * current_ME2;

  current_ME2 = getME_ttWmJets_gd( p1 );
  pdf_sum     = fg1*fdbar0 + fg1*fsbar0;
  ME2        += pdf_sum * current_ME2;

  //gu->ttdWm, gc->ttsWp
  current_ME2 = getME_ttWmJets_gubar( p0 );
  pdf_sum     = fg0*fu1+fg0*fc1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_ttWmJets_gubar( p1 );
  pdf_sum     = fg1*fu0+fg1*fc0;
  ME2        += pdf_sum * current_ME2; 

  //dbaru->ttgWp, sbarc->ttgWp
  current_ME2 = getME_ttWmJets_dubar( p0 );
  pdf_sum     = fdbar0*fu1+fsbar0*fc1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_ttWmJets_dubar( p1 );
  pdf_sum     = fdbar1*fu0+fsbar1*fc0;
  ME2        += pdf_sum * current_ME2; 
 

  return ME2 / (x0 * x1);

}






double get_ttjetsWeightedME( double x0, double x1, const vector<double*> &p) {
  
  double pdfQ_ = Physics::mtop; 
  double current_ME2, pdf_sum;
  double ME2 = 0;
  
  // The 2 permutations
  vector<double*> p0=p;
  vector<double*> p1=p;  
  p1[0]=p[1];
  p1[1]=p[0];

  
  // Parton codes:
  //    t  b  c  s  u  d  g  d  u  s  c  b  t
  //   -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
  double fg0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 );
  double fu0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 );
  double fd0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 );
  double fc0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 );
  double fs0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 );
  double fubar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 );
  double fdbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 );
  double fcbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 );
  double fsbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 );

  double fg1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  0 );
  double fu1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  2 );
  double fd1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  1 );
  double fc1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  4 );
  double fs1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  3 );
  double fubar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -2 );
  double fdbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -1 );
  double fcbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -4 );
  double fsbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -3 );

  //gg->ttg
  current_ME2 = getME_gg_ttg( p0 );
  pdf_sum     = fg0*fg1;
  ME2        += pdf_sum * current_ME2;

  //gq->ttq
  current_ME2 = getME_gq_ttq( p0 );
  pdf_sum     = fg0*(fu1+fd1+fs1+fc1+fubar1+fdbar1+fsbar1+fcbar1);
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_gq_ttq( p1 );
  pdf_sum     = fg1*(fu0+fd0+fs0+fc0+fubar0+fdbar0+fsbar0+fcbar0);
  ME2        += pdf_sum * current_ME2;


  //qqbar->ttg
  current_ME2 = getME_qqbar_ttg( p0 );
  pdf_sum     = fu0*fubar1+fd0*fdbar1+fc0*fcbar1+fs0*fsbar1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_qqbar_ttg( p1 );
  pdf_sum     = fu1*fubar0+fd1*fdbar0+fc1*fcbar0+fs1*fsbar0;
  ME2        += pdf_sum * current_ME2;
 
  return ME2 / (x0 * x1);

}








double get_ttbarWeightedME( double x0, double x1, const vector<double*> &p) {
  
  double pdfQ_ = Physics::mtop; 
  double current_ME2, pdf_sum;
  double ME2 = 0;
  
  // The 2 permutations
  vector<double*> p0=p;
  vector<double*> p1=p;  
  p1[0]=p[1];
  p1[1]=p[0];

  
  // Parton codes:
  //    t  b  c  s  u  d  g  d  u  s  c  b  t
  //   -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6
  double fg0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  0 );
  double fu0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  2 );
  double fd0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  1 );
  double fc0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  4 );
  double fs0    = LHAPDF::xfx( PDFSet, x0, pdfQ_,  3 );
  double fubar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -2 );
  double fdbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -1 );
  double fcbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -4 );
  double fsbar0 = LHAPDF::xfx( PDFSet, x0, pdfQ_, -3 );

  double fg1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  0 );
  double fu1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  2 );
  double fd1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  1 );
  double fc1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  4 );
  double fs1    = LHAPDF::xfx( PDFSet, x1, pdfQ_,  3 );
  double fubar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -2 );
  double fdbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -1 );
  double fcbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -4 );
  double fsbar1 = LHAPDF::xfx( PDFSet, x1, pdfQ_, -3 );

  //gg->ttbar
  current_ME2 = getME_gg_ttbar( p0 );
  pdf_sum     = fg0*fg1;
  ME2        += pdf_sum * current_ME2;


  //qqbar->ttbar
  current_ME2 = getME_qqbar_ttbar( p0 );
  pdf_sum     = fu0*fubar1+fd0*fdbar1+fc0*fcbar1+fs0*fsbar1;
  ME2        += pdf_sum * current_ME2;
  //
  current_ME2 = getME_qqbar_ttbar( p1 );
  pdf_sum     = fu1*fubar0+fd1*fdbar0+fc1*fcbar0+fs1*fsbar0;
  ME2        += pdf_sum * current_ME2;
 
  return ME2 / (x0 * x1);

}
