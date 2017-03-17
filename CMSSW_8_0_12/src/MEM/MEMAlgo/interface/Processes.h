/* 
 * File:   Processes.h
 * Author: grasseau
 *
 * Created on 1 décembre 2014, 19:15
 */

#ifndef PROCESSES_H
#define	PROCESSES_H

//# include "TLorentzVector.h"
# include "MEMDataFormats/MEMDataFormats/interface/RunConfig.h"

# include <vector>
using namespace std;

static const int PDFSet = 1;

double getME_ttH_gg( vector<double*> &p );
double getME_ttH_qqbar( vector<double*> &p );

double getME_ttZ_gg( vector<double*> &p );
double getME_ttZ_ddbar( vector<double*> &p );
double getME_ttZ_uubar( vector<double*> &p );

double getME_ttZ_Zonly_gg( vector<double*> &p );
double getME_ttZ_Zonly_uubar( vector<double*> &p );

double getME_ttZ_Zonly_Zll_gg( vector<double*> &p );
double getME_ttZ_Zonly_Zll_uubar( vector<double*> &p );


double getME_ttWm_dubar( vector<double*> &p );

double getME_ttWmJets_gd( vector<double*> &p );
double getME_ttWmJets_gubar( vector<double*> &p );
double getME_ttWmJets_dubar( vector<double*> &p );

double getME_gg_ttg( vector<double*> &p );
double getME_gq_ttq( vector<double*> &p );
double getME_qqbar_ttg( vector<double*> &p );

double getME_gg_ttbar( vector<double*> &p );
double getME_qqbar_ttbar( vector<double*> &p );

// ME
void init_ttH_HTauTauMEProcesses( const char *configname );
void init_ttH_HTauTauMEProcesses( RunConfig cfg );
void init_ttH_HTauTauMEProcesses( const char *fileName, double pdfQ );
void init_ttH_HTauTauMEProcesses( string fileName, double pdfQ );

double get_ttHWeightedME( double x0, double x1, const vector<double*> &p );
double get_ttZWeightedME( double x0, double x1, const vector<double*> &p );
double get_ttZZonlyWeightedME( double x0, double x1, const vector<double*> &p );
double get_ttZZonlyZllWeightedME( double x0, double x1, const vector<double*> &p );
double get_ttWmWeightedME( double x0, double x1, const vector<double*> &p );
double get_ttWpWeightedME( double x0, double x1, const vector<double*> &p );
double get_ttWmJetsWeightedME( double x0, double x1, const vector<double*> &p );
double get_ttWpJetsWeightedME( double x0, double x1, const vector<double*> &p );
double get_ttjetsWeightedME( double x0, double x1, const vector<double*> &p );
double get_ttbarWeightedME( double x0, double x1, const vector<double*> &p );

#endif	/* PROCESSES_H */

