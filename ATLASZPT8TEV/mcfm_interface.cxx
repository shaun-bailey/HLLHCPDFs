//
//   mcfm_interface.cxx
//
//
//   Copyright (C) 2007 M.Sutton (sutt@cern.ch)
//
//   $Id: mcfm_interfce.cxx, v   Fri  8 Nov 2013 09:07:01 CET sutt


#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <sys/stat.h>

#include <cstdlib>
#include <sys/time.h>

#include "TFile.h"
#include "TH1D.h"
#include "TVectorT.h"
#include "TString.h"

#include "mcfm_grid.h"

bool file_exists( const std::string& filename ) {
  struct stat sb;
  if ( stat( filename.c_str(), &sb)==0 ) return true; // && S_ISREG(sb.st_mode ))
  else return false;
}

static const int mxpart = 14;    // mcfm parameter : max number of partons in event record.
                                 // Defined in Inc/constants.f

static const int _Ngrids = 6;
static       int  Ngrids = 1;
appl::mcfm_grid* mygrid[_Ngrids];

// The below line (and others where values are assigned to this variable) changed to compile without err 
static       char *gridFiles[_Ngrids];

static       double Observable[_Ngrids];
int                 nObsBins[_Ngrids];


 // --------------- eta binning --------------------------//

 // Default

 static const double eta[] =
   {
     -7.0 , -6.5 , -6.0 , -5.5 , -5.0 ,
     -4.5 , -4.0 , -3.5 , -3.0 , -2.5 ,
     -2.0 , -1.5 , -1.0 , -0.5 ,  0.0 ,
      0.5 ,  1.0 ,  1.5 ,  2.0 ,  2.5 ,
      3.0 ,  3.5 ,  4.0 ,  4.5 ,  5.0 ,
      5.5 ,  6.0 ,  6.5 ,  7.0
   };

 // ----------------  pt binning ------------------- //

 // Default

 static const double pt[] =  {
   0.0,    5.0,   10.0,   15.0,   20.0,
   25.0,  30.0,   35.0,   40.0,   45.0,
   50.0,  55.0,   60.0,   65.0,   70.0,
   75.0,  80.0,   85.0,   90.0,   95.0,
   100.0, 105.0,  110.0,  115.0,  120.0,
   125.0, 130.0,  135.0,  140.0,  145.0,
   150.0, 175.0,  200.0,  225.0,  250.0,
   275.0, 300.0,  350.0,  400.0,  450.0,
   500.0, 750.0, 1000.0, 1500.0, 2000.0,
   2500.0
 };

 // ------------ Invariant mass binning ------------ //

 // Default

 static const double mass[] =  {
    5.0,   7.5,  10.0,   12.5,   15.0,
   20.0,  25.0,  30.0,   40.0,   50.0,
   75.0, 100.0, 125.0,  150.0,  200.0,
  300.0, 500.0, 750.0, 1000.0, 1500.0
 };

 long unsigned int runs  =  0;
 bool isBooked           =  false;
 std::string glabel      =  "";

 void getObservable( const double evt[][mxpart] );

 int  cuts(int);

 std::string date() {
   time_t _t;
   time(&_t);
   return ctime(&_t);
 }

 void book_grid()  // inital grid booking
 {
   if (isBooked) return;

   std::cout<<" ***********************************************"<<std::endl;
   std::cout<<" booking the grids " << date() << std::endl;

   // binning information for the grid constructor
   double xLow    = 1.0e-9, xUp = 1.0;
   int    nXbins  = 40;
   int    xorder  = 6;
   double q2Low   = 1.0*1.0;
   double q2Up    = 13000*13000;
   int    nQ2bins = 15;
   int    qorder  = 3;
   // set transform2 value
   double apramval=5.;
   appl::grid::transformvar(apramval);

   // lowest order in alphas
   int lowest_order = 0;
   // how many loops
   int nloops = 1;

   // Initialise the array of gridFiles
   
   for(int c = 0; c < _Ngrids; c++)
   {  gridFiles[c]=(char*)malloc(80*sizeof(char)); gridFiles[c][0]=" "; }

   // number of observables and binning for observables
   const double *obsBins[_Ngrids] = { pt };

   std::string pdf_function;

   glabel = "grid-40/6-15/3";

   const char* DatasetID = std::getenv("DatasetID");
   if ( DatasetID && std::string(DatasetID)!="" ) glabel = DatasetID;
   std::cout << "Dataset ID: " << DatasetID << std::endl;
  //  std::cout << "q2low " << q2lower << "\tq2up " << q2upper << std::endl;
  //  std::cout << "Process : " << nproc_.nproc << std::endl;

  if ( glabel == "ATLASWPT31PB_WP" )
    {
      std::cout << "Atlas Wpt measurement (31 pb^-1) - W+" << std::endl;
      pdf_function = "mcfmwp.config";
      q2Low   = 6463.83, q2Up = 6463.85;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 2;
      strcpy(gridFiles[0],"_pt34.root");
      strcpy(gridFiles[1],"_tot.root");

      nObsBins[0] = 11;
      nObsBins[1] =  1;

      double _pt[12] = { 0.0, 8.0, 23.0, 38.0, 55.0, 75.0,  95.0, 120.0, 145.0, 175.0, 210.0, 300.0 };
      double _tot[2] = { 0.0, 300.0 };

      obsBins[0] = _pt;
      obsBins[1] = _tot;
    }

  else if ( glabel == "ATLASWPT31PB_WM" )
    {
      std::cout << "Atlas Wpt measurement (31 pb^-1) - W-" << std::endl;
      pdf_function = "mcfmwm.config";
      q2Low   = 6463.83, q2Up = 6463.85;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 2;
      strcpy(gridFiles[0],"_pt34.root");
      strcpy(gridFiles[1],"_tot.root");

      nObsBins[0] = 11;
      nObsBins[1] =  1;

      double _pt[12] = { 0.0, 8.0, 23.0, 38.0, 55.0, 75.0,  95.0, 120.0, 145.0, 175.0, 210.0, 300.0 };
      double _tot[2] = { 0.0, 300.0 };

      obsBins[0] = _pt;
      obsBins[1] = _tot;
    }

  else if ( glabel == "ATLASPHOTONS11" )
    {
      std::cout << "ATLAS Inclusive isolated prompt photons, 7 TeV, 4.6 fb^-1" << std::endl;
      //pdf_function = "photonLO.config:photonNLO.config";
	pdf_function = "basic";
      Ngrids  = 3;
      strcpy(gridFiles[0],"_etagamma.root");
      strcpy(gridFiles[1],"_Etctr.root");
      strcpy(gridFiles[2],"_Etfwd.root");

      nObsBins[0] = 12;
      nObsBins[1] = 13;
      nObsBins[2] = 10;

      double _etag[13] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.37, 1.52, 1.8, 2.0, 2.2, 2.37 };
      double _Etc[14]  = { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000 };
      double _Etf[11]  = { 100, 125, 150, 175, 200, 250, 300, 350, 400, 500, 600 };

      obsBins[0] = _etag;
      obsBins[1] = _Etc;
      obsBins[2] = _Etf;

    }

  else if (glabel == "ATLASPHT12" ) 
     {
      std::cout << "ATLAS Inclusive isolated prompt photons, 8 TeV, 20 fb^-1" << std::endl;
//    pdf_function = "photonLO.config:photonNLO.config";
      pdf_function = "basic";
      Ngrids = 5;
      strcpy(gridFiles[0],"_Etagamma.root");
      strcpy(gridFiles[1],"_Et_1bin.root");
      strcpy(gridFiles[2],"_Et_2bin.root");
      strcpy(gridFiles[3],"_Et_3bin.root");
      strcpy(gridFiles[4],"_Et_4bin.root");

      nObsBins[0] = 12; 
      nObsBins[1] = 18;
      nObsBins[2] = 17;
      nObsBins[3] = 14;
      nObsBins[4] = 14;
  
      double _Et1[19] = {65,75,85,105,125,150,175,200,250,300,350,400,470,550,650,750,900,1100,1500};
      double _Et2[18] = {65,75,85,105,125,150,175,200,250,300,350,400,470,550,650,750,900,1100}; 
      double _Et3[15] = {65,75,85,105,125,150,175,200,250,300,350,400,470,550,650};
      double _Et4[15] = {65,75,85,105,125,150,175,200,250,300,350,400,470,550,650};
      double _etag[13] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.37, 1.56, 1.8, 2.0, 2.2, 2.37};

      obsBins[0] = _etag;	
      obsBins[1] = _Et1;
      obsBins[2] = _Et2;
      obsBins[3] = _Et3;
      obsBins[4] = _Et4;	

    }


  else if (glabel== "ATLASPHT15")
    {

      std::cout << "ATLAS Inclusive isolated prompt photons, 13 TeV, 3.2 fb^-1" << std::endl;
//    pdf_function = "photonLO.config:photonNLO.config";
      pdf_function = "basic";
      Ngrids = 5;

      strcpy(gridFiles[0],"_Etagamma.root");
      strcpy(gridFiles[1],"_Et_1bin.root");
      strcpy(gridFiles[2],"_Et_2bin.root");
      strcpy(gridFiles[3],"_Et_3bin.root");
      strcpy(gridFiles[4],"_Et_4bin.root");

      nObsBins[0] = 12;
      nObsBins[1] = 14;
      nObsBins[2] = 14;
      nObsBins[3] = 13;
      nObsBins[4] = 12;

      double _Et1[19] = {125,150,175,200,250,300,350,400,470,550,650,750,900,1100,1500};
      double _Et2[18] = {125,150,175,200,250,300,350,400,470,550,650,750,900,1100,1500};
      double _Et3[15] = {125,150,175,200,250,300,350,400,470,550,650,750,900,1100};
      double _Et4[15] = {125,150,175,200,250,300,350,400,470,550,650,750,900};
      double _etag[13] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.37, 1.56, 1.8, 2.0, 2.2, 2.37};


      obsBins[0] = _etag;
      obsBins[1] = _Et1;
      obsBins[2] = _Et2;
      obsBins[3] = _Et3;
      obsBins[4] = _Et4;

    }

  else if ( glabel == "ATLASTTB11" )
    {
      std::cout << "Atlas ttbar differential distributions (4.6 fb^-1)" << std::endl;
      pdf_function = "mcfm-TT";
      lowest_order = 2;
      q2Low   = std::pow(172.499,2), q2Up = std::pow(172.501,2);
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 5;
      strcpy(gridFiles[0],"_yttb.root");
      strcpy(gridFiles[1],"_pTt.root");
      strcpy(gridFiles[2],"_pTttb.root");
      strcpy(gridFiles[3],"mttb.root");
      strcpy(gridFiles[4],"_tot.root");

      nObsBins[0] = 3;
      double _y[4] = { 0.0 , 0.5 , 1.0 , 2.5 };
      obsBins[0] = _y;

      nObsBins[1] = 7;
      double _pt[8] = { 0. , 50. , 100. , 150. , 200. , 250. , 350. , 800. };
      obsBins[1] = _pt;

      nObsBins[2] = 4;
      double _pttb[5] = { 0. , 40. , 170. , 340. , 1000. };
      obsBins[2] = _pttb;

      nObsBins[3] = 5;
      double _mttb[6] = { 250. , 450. , 550. , 700. , 950. , 2700. };
      obsBins[3] = _mttb;

      nObsBins[4] = 1;
      double _tot[2] = { -20. , 20. };
      obsBins[4] = _tot;
    }

  else if ( glabel == "ATLASWZRAP36PB_Z" )
    {
      std::cout << "Atlas Z rapidity measurement (36 pb^-1)" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = 8315.17, q2Up = 8315.19;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_y34.root");

      nObsBins[0] = 8;

      double _y[9] = { 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.6 };

      obsBins[0] = _y;
    }

  else if ( glabel == "ATLASWZRAP36PB_WP" )
    {
      std::cout << "Atlas W+ rapidity measurement (36 pb^-1)" << std::endl;
      pdf_function = "mcfmwp.config";
      q2Low   = 6463.83, q2Up = 6463.85;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_etalept.root");

      nObsBins[0] = 11;

      double _y[12] = { 0.00, 0.21, 0.42, 0.63, 0.84, 1.05, 1.37, 1.52, 1.74, 1.95, 2.18, 2.50 };

      obsBins[0] = _y;
    }

  else if ( glabel == "ATLASWZRAP36PB_WM" )
    {
      std::cout << "Atlas W- rapidity measurement (36 pb^-1)" << std::endl;
      pdf_function = "mcfmwm.config";
      q2Low   = 6463.83, q2Up = 6463.85;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_etalept.root");

      nObsBins[0] = 11;

      double _y[12] = { 0.00, 0.21, 0.42, 0.63, 0.84, 1.05, 1.37, 1.52, 1.74, 1.95, 2.18, 2.50 };

      obsBins[0] = _y;
    }

  else if ( glabel == "ATLASLOMASSDY11" )
    {
      std::cout << "Atlas low-mass Drell-Yan, 7 TeV, 'nominal' analysis" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(26.,2), q2Up = std::pow(66.,2);
      nQ2bins = 15;
      qorder  = 3;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_m34.root");

      nObsBins[0] = 8;

      double _m34[9] = { 26., 31., 36., 41., 46., 51., 56., 61., 66. };

      obsBins[0] = _m34;
    }

  else if ( glabel == "ATLASLOMASSDY11EXT" )
    {
      std::cout << "Atlas low-mass Drell-Yan, 7 TeV, 'extended' analysis" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(12.,2), q2Up = std::pow(66.,2);
      nQ2bins = 15;
      qorder  = 3;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_m34.root");

      nObsBins[0] = 6;

      double _m34[7] = { 12., 17., 22., 28., 36., 46., 66. };

      obsBins[0] = _m34;
    }


  else if ( glabel == "ATLASZHIGHMASS49FB" )
    {
      std::cout << "Atlas Hi-Mass Drell-Yan measurement" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(116.,2), q2Up = std::pow(1500.,2);
      nQ2bins = 15;
      qorder  = 3;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_m34.root");

      nObsBins[0] = 13;

      double _m34[14] = { 116.0, 130.0, 150.0, 170.0, 190.0, 210.0, 230.0, 250.0, 300.0, 400.0,
			  500.0, 700.0, 1000.0, 1500.0 };
      obsBins[0] = _m34;
    }

 else if ( glabel == "ATLASZPT47FB" )
    {
      std::cout << "Atlas Zpt measurement (4.7 fb^-1) - Z" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = 8315.17, q2Up = 8315.19;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 3;
      strcpy(gridFiles[0],"_pt.root");
      strcpy(gridFiles[1],"_tot_pt.root");
      strcpy(gridFiles[2],"_tot_y.root");

      nObsBins[0] = 26;
      nObsBins[1] =  1;
      nObsBins[2] =  1;

      double _pt[27] = {  0.0,  2.0,   4.0,   6.0,   8.0,  10.0,  12.0, 14.0, 16.0, 18.0,
			 22.0, 26.0,  30.0,  34.0,  38.0,  42.0,  46.0, 50.0, 54.0, 60.0,
			 70.0, 80.0, 100.0, 150.0, 200.0, 300.0, 800.0 };
      double _tot_pt[2] = { 0.0, 800.0 };
      double _tot_y[2] = { -4.9, 4.9 };

      obsBins[0] = _pt;
      obsBins[1] = _tot_pt;
      obsBins[2] = _tot_y;
    }

 else if ( glabel == "ATLASZPT47FB-NLO" )
    {
      std::cout << "Atlas Zpt measurement (4.7 fb^-1) - Z, done using Z+j trick" << std::endl;
      pdf_function = "basic";
      q2Low   = 8315.17, q2Up = 8315.19;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 3;
      strcpy(gridFiles[0],"_pt.root");
      strcpy(gridFiles[1],"_tot_pt.root");
      strcpy(gridFiles[2],"_tot_y.root");

      nObsBins[0] = 26;
      nObsBins[1] =  1;
      nObsBins[2] =  1;

      double _pt[27] = {  0.0,  2.0,   4.0,   6.0,   8.0,  10.0,  12.0, 14.0, 16.0, 18.0,
			 22.0, 26.0,  30.0,  34.0,  38.0,  42.0,  46.0, 50.0, 54.0, 60.0,
			 70.0, 80.0, 100.0, 150.0, 200.0, 300.0, 800.0 };
      double _tot_pt[2] = { 0.0, 800.0 };
      double _tot_y[2] = { -4.9, 4.9 };

      obsBins[0] = _pt;
      obsBins[1] = _tot_pt;
      obsBins[2] = _tot_y;
    }

  else if ( glabel == "CMSDY2D11-BIN1" )
    {
      std::cout << "CMS double differential Drell-Yan measurement (2011) - Bin 1" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(24.99,2), q2Up = std::pow(25.01,2);
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_rapidity.root");

      nObsBins[0] = 24;

      double _y34[25] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
			  1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
			  2.0, 2.1, 2.2, 2.3, 2.4  };
      obsBins[0] = _y34;
    }

  else if ( glabel == "CMSDY2D11-BIN2" )
    {
      std::cout << "CMS double differential Drell-Yan measurement (2011) - Bin 2" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(37.49,2), q2Up = std::pow(37.51,2);
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_rapidity.root");

      nObsBins[0] = 24;

      double _y34[25] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
			  1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
			  2.0, 2.1, 2.2, 2.3, 2.4  };
      obsBins[0] = _y34;
    }

  else if ( glabel == "CMSDY2D11-BIN3" )
    {
      std::cout << "CMS double differential Drell-Yan measurement (2011) - Bin 3" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(52.49,2), q2Up = std::pow(52.51,2);
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_rapidity.root");

      nObsBins[0] = 24;

      double _y34[25] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
			  1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
			  2.0, 2.1, 2.2, 2.3, 2.4  };
      obsBins[0] = _y34;
    }

  else if ( glabel == "CMSDY2D11-BIN4" )
    {
      std::cout << "CMS double differential Drell-Yan measurement (2011) - Bin 4" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(89.99,2), q2Up = std::pow(90.01,2);
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_rapidity.root");

      nObsBins[0] = 24;

      double _y34[25] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
			  1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
			  2.0, 2.1, 2.2, 2.3, 2.4  };
      obsBins[0] = _y34;
    }

  else if ( glabel == "CMSDY2D11-BIN5" )
    {
      std::cout << "CMS double differential Drell-Yan measurement (2011) - Bin 5" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(159.99,2), q2Up = std::pow(160.01,2);
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_rapidity.root");

      nObsBins[0] = 24;

      double _y34[25] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
			  1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
			  2.0, 2.1, 2.2, 2.3, 2.4  };
      obsBins[0] = _y34;
    }

  else if ( glabel == "CMSDY2D11-BIN6" )
    {
      std::cout << "CMS double differential Drell-Yan measurement (2011) - Bin 6" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(849.99,2), q2Up = std::pow(850.01,2);
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[1],"_rapidity.root");

      nObsBins[0] = 12;

      double _y34[13] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
			  1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4  };

      obsBins[0] = _y34;
    }

  else if ( glabel == "CMSDY2D12-BIN1" )
    {
      std::cout << "CMS double differential Drell-Yan measurement (2011) - Bin 1" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(24.99,2), q2Up = std::pow(25.01,2);
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_rapidity.root");

      nObsBins[0] = 24;

      double _y34[25] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
			  1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
			  2.0, 2.1, 2.2, 2.3, 2.4  };
      obsBins[0] = _y34;
    }

  else if ( glabel == "CMSDY2D12-BIN2" )
    {
      std::cout << "CMS double differential Drell-Yan measurement (2011) - Bin 2" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(37.49,2), q2Up = std::pow(37.51,2);
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_rapidity.root");

      nObsBins[0] = 24;

      double _y34[25] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
			  1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
			  2.0, 2.1, 2.2, 2.3, 2.4  };
      obsBins[0] = _y34;
    }

  else if ( glabel == "CMSDY2D12-BIN3" )
    {
      std::cout << "CMS double differential Drell-Yan measurement (2011) - Bin 3" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(52.49,2), q2Up = std::pow(52.51,2);
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_rapidity.root");

      nObsBins[0] = 24;

      double _y34[25] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
			  1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
			  2.0, 2.1, 2.2, 2.3, 2.4  };
      obsBins[0] = _y34;
    }

  else if ( glabel == "CMSDY2D12-BIN4" )
    {
      std::cout << "CMS double differential Drell-Yan measurement (2011) - Bin 4" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(89.99,2), q2Up = std::pow(90.01,2);
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_rapidity.root");

      nObsBins[0] = 24;

      double _y34[25] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
			  1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
			  2.0, 2.1, 2.2, 2.3, 2.4  };
      obsBins[0] = _y34;
    }

  else if ( glabel == "CMSDY2D12-BIN5" )
    {
      std::cout << "CMS double differential Drell-Yan measurement (2011) - Bin 5" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(159.99,2), q2Up = std::pow(160.01,2);
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_rapidity.root");

      nObsBins[0] = 24;

      double _y34[25] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
			  1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
			  2.0, 2.1, 2.2, 2.3, 2.4  };
      obsBins[0] = _y34;
    }

  else if ( glabel == "CMSDY2D12-BIN6" )
    {
      std::cout << "CMS double differential Drell-Yan measurement (2011) - Bin 6" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = std::pow(849.99,2), q2Up = std::pow(850.01,2);
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_rapidity.root");

      nObsBins[0] = 12;

      double _y34[13] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
			  1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4  };

      obsBins[0] = _y34;
    }

  else if ( glabel == "CMSWCHARM_WPCB" )
    {
      std::cout << "CMS W+Charm production, 7 TeV - W+ C~" << std::endl;
      pdf_function = "mcfm-wpc";
      lowest_order = 1;
      q2Low   = 6463.83, q2Up =  6463.85;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_leptrap.root");

      nObsBins[0] = 5;

      double _y4[6] = { 0.0, 0.35, 0.7, 1.1, 1.6, 2.1 };

      obsBins[0] = _y4;
    }

  else if ( glabel == "CMSWCHARM_WMC" )
    {
      std::cout << "CMS W+Charm production, 7 TeV - W- C" << std::endl;
      pdf_function = "mcfm-wmc";
      lowest_order = 1;
      q2Low   = 6463.83, q2Up =  6463.85;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_leptrap.root");

      nObsBins[0] = 5;

      double _y3[6] = { 0.0, 0.35, 0.7, 1.1, 1.6, 2.1 };

      obsBins[0] = _y3;
    }

  else if ( glabel == "CMSWEASY840PB_WP" )
    {
      std::cout << "CMS W electron asymmety, 7 TeV (840 pb^-1)- W+" << std::endl;
      pdf_function = "mcfmwp.config";
      q2Low   = 6463.83, q2Up =  6463.85;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_leptrap.root");

      nObsBins[0] = 12;

      double _y4[13] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
			1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};

      obsBins[0] = _y4;
    }

  else if ( glabel == "CMSWEASY840PB_WM" )
    {
      std::cout << "CMS W electron asymmety, 7 TeV (840 pb^-1)- W-" << std::endl;
      pdf_function = "mcfmwm.config";
      q2Low   = 6463.83, q2Up =  6463.85;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_leptrap.root");

      nObsBins[0] = 12;

      double _y3[13] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
			1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};

      obsBins[0] = _y3;
    }

  else if ( glabel == "CMSWMASY47FB_WP" )
    {
      std::cout << "CMS W electron asymmety, 7 TeV (4.7 fb^-1)- W+" << std::endl;
      pdf_function = "mcfmwp.config";
      q2Low   = 6463.83, q2Up =  6463.85;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_leptrap.root");

      nObsBins[0] = 11;

      double _y4[12] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
			1.2, 1.4, 1.6, 1.85, 2.1, 2.4};
      obsBins[0] = _y4;
    }

  else if ( glabel == "CMSWMASY47FB_WM" )
    {
      std::cout << "CMS W electron asymmetry, 7 TeV (4.7 fb^-1)- W-" << std::endl;
      pdf_function = "mcfmwm.config";
      q2Low   = 6463.83, q2Up =  6463.85;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_leptrap.root");

      nObsBins[0] = 11;

      double _y3[12] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
			1.2, 1.4, 1.6, 1.85, 2.1, 2.4};
      obsBins[0] = _y3;
    }

  else if ( glabel == "CMSTTB12" )
    {
      std::cout << "CMS ttbar differential distributions (8 TeV, 19.7 fb^-1)" << std::endl;
      pdf_function = "mcfm-TT";
      lowest_order = 2;
      q2Low   = std::pow(173.299,2), q2Up = std::pow(173.301,2);
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 6;
      strcpy(gridFiles[0],"_yt");
      strcpy(gridFiles[1],"_yttb.root");
      strcpy(gridFiles[2],"_pTt.root");
      strcpy(gridFiles[3],"_pTttb.root");
      strcpy(gridFiles[4],"mttb.root");
      strcpy(gridFiles[5],"_tot.root");

      nObsBins[0] = 10;
      double _yt[11] = { -2.5 , -1.6 , -1.2 , -0.8 , -0.4 , 0.0 , 0.4 , 0.8 , 1.2 , 1.6 , 2.5 };
      obsBins[0] = _yt;

      nObsBins[1] = 10;
      double _yttb[11] = { -2.5 , -1.3, -0.9 , -0.6, -0.3 , 0.0 , 0.3, 0.6 , 0.9 , 1.3 , 2.5 };
      obsBins[1] = _yttb;

      nObsBins[2] = 8;
      double _pt[9] = { 0. , 60. , 100. , 150. , 200. , 260. , 320. , 400. , 500. };
      obsBins[2] = _pt;

      nObsBins[3] = 6;
      double _pttb[7] = { 0. , 20. , 45. , 75. , 120. , 190. , 300. };
      obsBins[3] = _pttb;

      nObsBins[4] = 7;
      double _mttb[8] = { 345. , 400. , 470. , 550. , 650. , 800. , 1100. , 1600. };
      obsBins[4] = _mttb;

      nObsBins[5] = 1;
      double _tot[2] = { -20. , 20. };
      obsBins[5] = _tot;
    }

  else if ( glabel.substr(0,10) == "CMSZDIFF12" ) 
    {
      std::cout << "CMS double differential Z pt and y distribution, 8 TeV (19.7 fb^-1)" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = 8315.17, q2Up = 8315.19;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids    = 2;
      strcpy(gridFiles[0],"_ptZ.root");
      strcpy(gridFiles[1],"_yZ.root");

      nObsBins[0] = 10;
      nObsBins[1] =  5;
	
      double _ptZ[11] = { 0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 170.0, 200.0, 1000.0 };
      obsBins[0] = _ptZ;

      double _yZ[6] = { 0.0, 0.4, 0.8, 1.2, 1.6, 2.0 };
      obsBins[1] = _yZ;
      
    }

  else if ( glabel == "LHCBWZ36PB_WP" )
    {
      std::cout << "LHCb W+ lepton rapidity distribution, 7 TeV (36 pb^-1)" << std::endl;
      pdf_function = "mcfmwp.config";
      q2Low   = 6463.83, q2Up =  6463.85;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_leptrap.root");

      nObsBins[0] = 5;

      double _y4[6] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5};
      obsBins[0] = _y4;
    }

  else if ( glabel == "LHCBWZ36PB_WM" )
    {
      std::cout << "LHCb W- lepton rapidity distribution, 7 TeV (36 pb^-1)" << std::endl;
      pdf_function = "mcfmwm.config";
      q2Low   = 6463.83, q2Up =  6463.85;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_leptrap.root");

      nObsBins[0] = 5;
      double _y3[6] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5};

      obsBins[0] = _y3;
    }

  else if ( glabel == "LHCBWZ940PB" )
    {
      std::cout << "LHCb Z rapidity distribution, 7 TeV (940 pb^-1)" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = 8315.17, q2Up = 8315.19;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_y34.root");

      nObsBins[0] = 10;
      double _y[11] = { 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5 };

      obsBins[0] = _y;
    }

    else if ( glabel == "LHCBWZMU7TEV_Z" )
      {
        std::cout << "LHCb Z -> mu mu rapidity distribution, 1 fb-1, 7 TeV" << std::endl;
        pdf_function = "mcfm-z";
        q2Low   = 8315.17, q2Up = 8315.19;
        nQ2bins = 3;
        qorder  = 1;

        Ngrids  = 1;
        strcpy(gridFiles[0],"_yZ.root");

        nObsBins[0] = 17;
        double _y[18] = { 2.000, 2.125, 2.250, 2.375, 2.500, 2.625, 2.750, 2.875,
                          3.000, 3.125, 3.250, 3.375, 3.500, 3.625, 3.750, 3.875,
                          4.000, 4.250 };

        obsBins[0] = _y;
      }

    else if ( glabel == "LHCBWZMU7TEV_WP" )
      {
        std::cout << "LHCb W+ -> mu+ nu rapidity distribution, 1 fb-1, 7 TeV" << std::endl;
        pdf_function = "mcfmwp.config";
        q2Low   = 6463.83, q2Up = 6463.85;
        nQ2bins = 3;
        qorder  = 1;

        Ngrids  = 1;
        strcpy(gridFiles[0],"_leptrap.root");

        nObsBins[0] = 8;
        double _y[9] = { 2.000, 2.250, 2.500, 2.750,
                         3.000, 3.250, 3.500, 4.000, 4.500 };

        obsBins[0] = _y;
      }

    else if ( glabel == "LHCBWZMU7TEV_WM" )
      {
        std::cout << "LHCb W- -> mu- nu~ rapidity distribution, 1 fb-1, 7 TeV" << std::endl;
        pdf_function = "mcfmwm.config";
        q2Low   = 6463.83, q2Up = 6463.85;
        nQ2bins = 3;
        qorder  = 1;

        Ngrids  = 1;
        strcpy(gridFiles[0],"_leptrap.root");

        nObsBins[0] = 8;
        double _y[9] = { 2.000, 2.250, 2.500, 2.750,
                         3.000, 3.250, 3.500, 4.000, 4.500 };

        obsBins[0] = _y;
      }

    else if ( glabel == "LHCBWZMU8TEV_Z" )
      {
        std::cout << "LHCb Z -> mu mu rapidity distribution, 2 fb-1, 8 TeV" << std::endl;
        pdf_function = "mcfm-z";
        q2Low   = 8315.17, q2Up = 8315.19;
        nQ2bins = 3;
        qorder  = 1;

        Ngrids  = 1;
        strcpy(gridFiles[0],"_yZ.root");

        nObsBins[0] = 18;
        double _y[19] = { 2.000, 2.125, 2.250, 2.375, 2.500, 2.625, 2.750, 2.875,
                          3.000, 3.125, 3.250, 3.375, 3.500, 3.625, 3.750, 3.875,
                          4.000, 4.250, 4.500 };

        obsBins[0] = _y;
      }

    else if ( glabel == "LHCBWZMU8TEV_WP" )
      {
        std::cout << "LHCb W+ -> mu+ nu rapidity distribution, 2 fb-1, 8 TeV" << std::endl;
        pdf_function = "mcfmwp.config";
        q2Low   = 6463.83, q2Up = 6463.85;
        nQ2bins = 3;
        qorder  = 1;

        Ngrids  = 1;
        strcpy(gridFiles[0],"_leptrap.root");

        nObsBins[0] = 8;
        double _y[9] = { 2.000, 2.250, 2.500, 2.750,
                         3.000, 3.250, 3.500, 4.000, 4.500 };

        obsBins[0] = _y;
      }

    else if ( glabel == "LHCBWZMU8TEV_WM" )
      {
        std::cout << "LHCb W- -> mu- nu~ rapidity distribution, 2 fb-1, 8 TeV" << std::endl;
        pdf_function = "mcfmwm.config";
        q2Low   = 6463.83, q2Up = 6463.85;
        nQ2bins = 3;
        qorder  = 1;

        Ngrids  = 1;
        strcpy(gridFiles[0],"_leptrap.root");

        nObsBins[0] = 8;
        double _y[9] = { 2.000, 2.250, 2.500, 2.750,
                         3.000, 3.250, 3.500, 4.000, 4.500 };

        obsBins[0] = _y;
      }


  else if ( glabel == "CDFZRAP" )
    {
      std::cout << "CDF Z rapidity distribution, 1.96 TeV (2.1 fb-1)" << std::endl;
      pdf_function = "basic";
      q2Low   = 8315.17, q2Up = 8315.19;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_y34.root");

      nObsBins[0] = 29;
      double _y[30] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                        1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
			2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9};

      obsBins[0] = _y;
    }

  else if ( glabel == "D0ZRAP" )
    {
      std::cout << "D0 Z rapidity distribution, 1.96 TeV (0.4 fb-1)" << std::endl;
      pdf_function = "basic";
      q2Low   = 8315.17, q2Up = 8315.19;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids  = 1;
      strcpy(gridFiles[0],"_y34.root");

      nObsBins[0] = 28;
      double _y[29] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                        1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
			2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8};

      obsBins[0] = _y;
    }

//  Pseudo data for theoretical studies

  else if ( glabel == "ZCHARM8TEV" || glabel == "ZCHARM13TEV" )
    {
      std::cout << "Associated production Z + Charm" << std::endl;
      pdf_function ="basic";

      Ngrids  = 2;
      strcpy(gridFiles[0],"_yZ.root");
      strcpy(gridFiles[1],"_ptZ.root");

      nObsBins[0] = 18;
      double _yZ[19] = 
	{ 
	  0.00, 0.25, 0.50, 0.75, 1.00, 1.25,
	  1.50, 1.75, 2.00, 2.25, 2.50, 2.75,
	  3.00, 3.25, 3.50, 3.75, 4.00, 4.25,
	  4.50 
	};
      
      nObsBins[1] = 21;
      double _ptZ[22] =
	{
	  50.,  75., 100., 125.,  150., 175.,
	  200., 225., 250., 275.,  300., 350.,
	  400., 450., 500., 550.,  600., 650.,
	  700., 800., 900., 1000.
	};
	
      obsBins[0] = _yZ;
      obsBins[1] = _ptZ;

    }

  else if ( glabel.substr(0,10) == "GammaC8TEV" || glabel.substr(0,10) == "GammaC13TEV" )
    {
      std::cout << "Associated production Gamma + Charm "  << std::endl;
      pdf_function ="basic";

      Ngrids  = 3;
      strcpy(gridFiles[0],"_yc.root");
      strcpy(gridFiles[1],"_yA.root");
      strcpy(gridFiles[2],"_ptA.root");

      nObsBins[0] = 18;
      double _yC[19] = 
	{ 
	  0.00, 0.25, 0.50, 0.75, 1.00, 1.25,
          1.50, 1.75, 2.00, 2.25, 2.50, 2.75,
	  3.00, 3.25, 3.50, 3.75, 4.00, 4.25,
	  4.50
	};

      nObsBins[1] = 18;
      double _yA[19] = 
	{ 
	  0.00, 0.25, 0.50, 0.75, 1.00, 1.25,
	  1.50, 1.75, 2.00, 2.25, 2.50, 2.75,
	  3.00, 3.25, 3.50, 3.75, 4.00, 4.25,
	  4.50
	};

      nObsBins[2] = 21;
      double _ptA[22] = 
	{ 
	    50.,  75., 100., 125.,  150., 175.,
	   200., 225., 250., 275.,  300., 350.,
	   400., 450., 500., 550.,  600., 650.,
	   700., 800., 900., 1000.
	};

      obsBins[0] = _yC;
      obsBins[1] = _yA;
      obsBins[2] = _ptA;
    }

  else if ( glabel == "CCbar8TEV" || glabel == "CCbar13TEV" )
    {
      std::cout << "CCbar production" << std::endl;
      pdf_function ="basic";

      Ngrids  = 4;
      strcpy(gridFiles[0],"_yc.root");
      strcpy(gridFiles[1],"_ptc.root");
      strcpy(gridFiles[2],"_ycbar.root");
      strcpy(gridFiles[3],"_ptcbar.root");

      nObsBins[0] = 36;
      double _yc[37] =
        {
	  -4.50, -4.25, -4.00, -3.75, -3.50, -3.25,
	  -3.00, -2.75, -2.50, -2.25, -2.00, -1.75,
	  -1.50, -1.25, -1.00, -0.75, -0.50, -0.25,
	   0.00,  0.25,  0.50,  0.75,  1.00,  1.25,
           1.50,  1.75,  2.00,  2.25,  2.50,  2.75,
           3.00,  3.25,  3.50,  3.75,  4.00,  4.25,
           4.50
        };

      nObsBins[1] = 15;
      double _ptc[16] =
	{
	   0.,  2.,  4.,  8.,  10.,  15.,  20.,  25., 30.,
	  40., 50., 60., 80., 100., 150., 200. 
	};

      nObsBins[2] = 36;
      double _ycbar[37] =
        {
	  -4.50, -4.25, -4.00, -3.75, -3.50, -3.25,
	  -3.00, -2.75, -2.50, -2.25, -2.00, -1.75,
	  -1.50, -1.25, -1.00, -0.75, -0.50, -0.25,
	   0.00,  0.25,  0.50,  0.75,  1.00,  1.25,
           1.50,  1.75,  2.00,  2.25,  2.50,  2.75,
           3.00,  3.25,  3.50,  3.75,  4.00,  4.25,
           4.50
        };

      nObsBins[3] = 15;
      double _ptcbar[16] =
	{
	   0.,  2.,  4.,  8.,  10.,  15.,  20.,  25., 30.,
	  40., 50., 60., 80., 100., 150., 200. 
	};

      obsBins[0] = _yc;
      obsBins[1] = _ptc;
      obsBins[2] = _ycbar;
      obsBins[3] = _ptcbar;

    }
  else if ( glabel == "ATLASZPT3AB" )
    {
      std::cout << "Atlas Zpt PseudoData (3 ab^-1) - Z" << std::endl;
      pdf_function = "mcfm-z";
      q2Low   = 8315.17, q2Up = 8315.19;
      nQ2bins = 3;
      qorder  = 1;

      Ngrids=7;//Ngrids  = 3;
      
      strcpy(gridFiles[0],"_pt34_1.root");
      strcpy(gridFiles[1],"_pt34_2.root");
      strcpy(gridFiles[2],"_pt34_3.root");
      strcpy(gridFiles[3],"_pt34_4.root");
      strcpy(gridFiles[4],"_pt34_5.root");
      strcpy(gridFiles[5],"_pt34_6.root");
      strcpy(gridFiles[6],"_m34.root");

      nObsBins[0] = 18;
      nObsBins[1] = 6;

      double _pt34[19] = {10., 13., 16., 20., 25., 30., 37., 45., 55., 65., 75., 85., 105., 150., 200., 900., 2000., 5000., 13000.};
      double _m34[7]= {12., 20., 30., 46., 66., 116., 150.};
  
      obsBins[0] = _pt34;
      obsBins[1] = _m34;
}
  else
    {
      std::cerr << "Process unknown" << std::endl;
      std::exit(-1);
    }
  std::cout << std::endl;


  /// Read the ckm matrix from mcfm to store in the grid automatically
  /// NB: we store 13 x 13 ckm matrix - mcfm only stores 11 x 11 so we
  ///     must add 1 to each index to keep them aligned

  std::vector< std::vector<double> > ckm_vsq( 13, std::vector<double>( 13, 0 ) );

  for ( int ic=0 ; ic<__nf2__ ; ic++ ) {
    for ( int ic1=0 ; ic1<__nf2__ ; ic1++ ) ckm_vsq[ic+1][ic1+1] = ckm_.vsq[ic][ic1];
  }

  std::vector<std::vector<double> > __ckm( 3, std::vector<double>(3, 0) );
  __ckm[0][0] = cabib_.Vud;
  __ckm[0][1] = cabib_.Vus;
  __ckm[0][2] = cabib_.Vub;
  __ckm[1][0] = cabib_.Vcd;
  __ckm[1][1] = cabib_.Vcs;
  __ckm[1][2] = cabib_.Vcb;


  for(int igrid=0; igrid < Ngrids; igrid++)
    {

      bool create_new = false;

      // if the file does not exist, create a new grid...
      if ( !file_exists(glabel+gridFiles[igrid]) )  create_new = true;

      // or if it does exists but root file is a zombie...
      if ( !create_new ) {
	TFile testFile( (glabel+gridFiles[igrid]).c_str() );
	if ( testFile.IsZombie() ) create_new = true;
	testFile.Close();
      }

      if ( create_new )
	{
	  std::cout << "Creating NEW grid... " << std::endl;

	  std::cout << "grid interpolation: "
		    << "\tQ2 " << nQ2bins << " " <<  q2Low << " " <<  q2Up << " " <<  qorder
		    << "\tx "  <<  nXbins << " " <<   xLow << " " <<   xUp << " " <<  xorder
		    << std::endl;



	  mygrid[igrid] = new appl::mcfm_grid( nObsBins[igrid], obsBins[igrid],      // obs bins
					       nQ2bins, q2Low, q2Up, qorder,         // Q2 bins and interpolation order
					       nXbins,   xLow,  xUp, xorder,         // x bins and interpolation order
					       pdf_function, lowest_order, nloops );
	  /// try reweighting for a bit
	  mygrid[igrid]->reweight(true);
	  mygrid[igrid]->setCMSScale( energy_.sqrts );


	  /// store the ckm matrix
	  //	  mygrid[igrid]->setckm2( ckm_vsq );

	  mygrid[igrid]->setckm( __ckm );

	  //	  grid_.nSubProcess = mygrid[igrid]->subProcesses();

	  std::cout << "reference histo name = "
	       << mygrid[igrid]->getReference()->GetName() << std::endl;

	  std::cout<<*mygrid[igrid]<<std::endl;
	}
      else
	{
	  std::cout << "Using existing grid file " << (glabel+gridFiles[igrid]) << std::endl;

	  mygrid[igrid] = new appl::mcfm_grid(glabel+gridFiles[igrid]); //optimise grid x,Q2 bins
	  //       grid_.nSubProcess = mygrid[igrid]->subProcesses();
	  mygrid[igrid]->getReference()->Reset();
	  mygrid[igrid]->optimise(nQ2bins, nXbins);

	  std::cout<<*(mygrid[igrid])<<std::endl;
	}
      // CTEQ like reweighting
      //      mygrid[igrid]->reweight( false );
    }

  runs = 0;
  isBooked = true;
  std::cout<<" ***********************************************"<<std::endl;
}


void fill_grid( const double evt[][mxpart] )
{

  static unsigned evtcounter = 1;
  if ( evtcounter%1000000==0 ) std::cout << "fill_grid() filled " << evtcounter << " weights " << date();
  evtcounter++;

  if (!isBooked)
    {
      book_grid();
      return;
    }

  getObservable( evt );

  for(int igrid = 0; igrid < Ngrids; igrid++)
    if(cuts(igrid)){

      mygrid[igrid]->fillMCFM( Observable[igrid] );
    }

  runs++; // counter of number of events
}


//
// just normalise to bin width
//
void Normalise(TH1D* h)
{
  if (h->GetNbinsX() == 1) return;

  for ( int ibin=1 ; ibin<=h->GetNbinsX() ; ibin++ )
    {
      double width = h->GetBinLowEdge(ibin+1) - h->GetBinLowEdge(ibin);
      h->SetBinContent( ibin, h->GetBinContent(ibin)/width );
    }
  return;
}



void write_grid(double& xstotal)   // writes out grid after some events
{
  std::cout<<"Write out grids ..."<<std::endl;

  for(int igrid = 0; igrid < Ngrids; igrid++)
    {
      std::cout << "saving grid N=" << igrid+1 << "\tof " << Ngrids << "\t";

      std::system("sleep 1");


      mygrid[igrid]->setNormalised( false );
      mygrid[igrid]->run() = (iterat_.ncall2)*(iterat_.itmx2);

      mygrid[igrid]->untrim();
      int untrim_size = mygrid[igrid]->size();

      mygrid[igrid]->trim();
      int trim_size = mygrid[igrid]->size();

      /// scale up by number of weights
      (*mygrid[igrid]) *= mygrid[igrid]->run();

      // normalise the reference histogram by bin width
      Normalise( mygrid[igrid]->getReference() );

      /// now scale *down* the reference histogram because we've just
      /// scaled it up ...
      //      mygrid[igrid]->getReference()->Scale( 1/mygrid[igrid]->run() );

      std::string filename = glabel+gridFiles[igrid];

#if 0

      std::string newpdfname = "";

      /// automatically optimise subprocesses - doesn't quite work yet
      /// when it does it may be moved into the grid itself
      if ( mygrid[igrid]->getGenpdf().find("basic")!=std::string::npos ) {

       	std::stringstream ss;
	ss << "proc" << nproc_.nproc;
	newpdfname = ss.str();

	std::cout << "appl::grid::Write() " << newpdfname << std::endl;

	mygrid[igrid]->Write( filename, "grid", newpdfname );
      }
      else {
	mygrid[igrid]->Write( glabel+gridFiles[igrid] );
      }

#else

      mygrid[igrid]->Write( filename );

#endif


      std::cout << "size(untrimmed)=" << untrim_size
		<< "\tsize(trimmed)=" << trim_size
		<< "\tfraction="      << 100.*trim_size/untrim_size << " %" << std::endl;

      //      int nsub = mygrid[igrid]->subProcesses();

      delete mygrid[igrid];

    }

  time_t _t;
  time(&_t);

  std::cout<<" ***********************************************"<<std::endl;
  std::cout<<" saved grids " << ctime(&_t);
  std::cout<<" ***********************************************"<<std::endl;
}



//
// ----------------------------------------------
//    analysis
// ----------------------------------------------
//

void getObservable(const double evt[][mxpart])
{
  // evt[momentum][particle number-1]
  // momentum[0,1,2,3] = (x,y,z,E)
  //

  // calculate observables
  for(int igrid = 0; igrid < Ngrids; igrid++)Observable[igrid] = 0.0; // initialize

  double p3[4] = {evt[3][2],evt[0][2],evt[1][2],evt[2][2]}; // (E,x,y,z)
  double p4[4] = {evt[3][3],evt[0][3],evt[1][3],evt[2][3]};
  double p5[4] = {evt[3][4],evt[0][4],evt[1][4],evt[2][4]};
  double p6[4] = {evt[3][5],evt[0][5],evt[1][5],evt[2][5]};
  double p7[4] = {evt[3][6],evt[0][6],evt[1][6],evt[2][6]};
  double p8[4] = {evt[3][7],evt[0][7],evt[1][7],evt[2][7]};

  double rapidity3 = 0.0;
  rapidity3 = (p3[0] + p3[3])/(p3[0] - p3[3]);
  (rapidity3 < 1e-13) ? rapidity3 = 100.0 : rapidity3 = 0.5*std::log(rapidity3);

  double rapidity4 = 0.0;
  rapidity4 = (p4[0] + p4[3])/(p4[0] - p4[3]);
  (rapidity4 < 1e-13) ? rapidity4 = 100.0 : rapidity4 = 0.5*std::log(rapidity4);

  double rapidity5 = 0.0;
  rapidity5 = (p5[0] + p5[3])/(p5[0] - p5[3]);
  (rapidity5 < 1e-13) ? rapidity3 = 100.0 : rapidity5 = 0.5*std::log(rapidity5);

  double rapidity6 = 0.0;
  rapidity6 = (p6[0] + p6[3])/(p6[0] - p6[3]);
  (rapidity6 < 1e-13) ? rapidity6 = 100.0 : rapidity6 = 0.5*std::log(rapidity6);

  double rapidity7 = 0.0;
  rapidity7 = (p7[0] + p7[3])/(p7[0] - p7[3]);
  (rapidity7 < 1e-13) ? rapidity7 = 100.0 : rapidity7 = 0.5*std::log(rapidity7);

  double rapidity8 = 0.0;
  rapidity8 = (p8[0] + p8[3])/(p8[0] - p8[3]);
  (rapidity8 < 1e-13) ? rapidity8 = 100.0 : rapidity8 = 0.5*std::log(rapidity8);

  double rapidity34 = 0.0;                      // rapidity of particle (3+4) in event record
  rapidity34  = (p3[0] + p4[0]) + (p3[3] + p4[3]);
  rapidity34 /= (p3[0] + p4[0]) - (p3[3] + p4[3]);
  (rapidity34 < 1e-13) ? rapidity34 = 100.0 : rapidity34 = 0.5*std::log(rapidity34);

  double rapidity56 = 0.0;                      // rapidity of particle (5+6) in event record
  rapidity56  = (p5[0] + p6[0]) + (p5[3] + p6[3]);
  rapidity56 /= (p5[0] + p6[0]) - (p5[3] + p6[3]);
  (rapidity56 < 1e-13) ? rapidity56 = 100.0 : rapidity56 = 0.5*std::log(rapidity56);

  double rapidity345 = 0.0;                      // rapidity of particle (3+4+5) in event record
  rapidity345  = (p3[0] + p4[0] + p5[0]) + (p3[3] + p4[3] + p5[3]);
  rapidity345 /= (p3[0] + p4[0] + p5[0]) - (p3[3] + p4[3] + p5[3]);
  (rapidity345 < 1e-13) ? rapidity345 = 100.0 : rapidity345 = 0.5*std::log(rapidity345);

  double rapidity678 = 0.0;                      // rapidity of particle (7+6+8) in event record
  rapidity678  = (p6[0] + p7[0] + p8[0]) + (p6[3] + p7[3] + p8[3]);
  rapidity678 /= (p6[0] + p7[0] + p8[0]) - (p6[3] + p7[3] + p8[3]);
  (rapidity678 < 1e-13) ? rapidity678 = 100.0 : rapidity678 = 0.5*std::log(rapidity678);

  double pt3 = 0;
  pt3 = std::sqrt( p3[1]*p3[1] + p3[2]*p3[2] );

  double pt4 = 0;
  pt4 = std::sqrt( p4[1]*p4[1] + p4[2]*p4[2] );

  double pt5 = 0;
  pt5 = std::sqrt( p5[1]*p5[1] + p5[2]*p5[2] );

  double pt6 = 0;
  pt6 = std::sqrt( p6[1]*p6[1] + p6[2]*p6[2] );

  double pt7 = 0;
  pt7 = std::sqrt( p7[1]*p7[1] + p7[2]*p7[2] );

  double pt8 = 0;
  pt8 = std::sqrt( p8[1]*p8[1] + p8[2]*p8[2] );

  double pt34 = 0;
  pt34 = std::sqrt( std::pow(p3[1] + p4[1],2) + std::pow(p3[2] + p4[2],2) );

  double pt56 = 0;
  pt56 = std::sqrt( std::pow(p5[1] + p6[1],2) + std::pow(p5[2] + p6[2],2) );

  double pt345 = 0;
  pt345 = std::sqrt( std::pow(p3[1] + p4[1] + p5[1],2) + std::pow(p3[2] + p4[2] + p5[2],2) );

  double pt678 = 0;
  pt678 = std::sqrt( std::pow(p6[1] + p7[1] + p8[1],2) + std::pow(p6[2] + p7[2] + p8[2],2) );

  double mass34 = std::sqrt(   (p3[0] + p4[0]) * (p3[0] + p4[0])
			     - (p3[1] + p4[1]) * (p3[1] + p4[1])
			     - (p3[2] + p4[2]) * (p3[2] + p4[2])
			     - (p3[3] + p4[3]) * (p3[3] + p4[3]) );

  if ( glabel == "ATLASWPT31PB_WP" ||
       glabel == "ATLASWPT31PB_WM" )
    {
      Observable [ 0 ] = pt34;
      Observable [ 1 ] = pt34;
    }

  else if ( glabel == "ATLASZPT47FB"     || 
	    glabel == "ATLASZPT47FB-NLO" )
    {
      Observable [ 0 ] = pt34;
      Observable [ 1 ] = pt34;
      Observable [ 2 ] = rapidity34;
    }

  else if ( glabel.substr(0,10) == "CMSZDIFF12" )
    {
      Observable [ 0 ] = pt34;
      Observable [ 1 ] = rapidity34;
    }

  else if ( glabel == "ATLASWZRAP36PB_Z" ||
	    glabel == "LHCBWZ940PB"      ||
	    glabel == "LHCBZERAP11"      ||
	    glabel == "LHCBWZMU7TEV_Z"   ||
	    glabel == "LHCBWZMU8TEV_Z"   ||
	    glabel == "CDFZRAP"          ||
	    glabel == "D0ZRAP"             )
    {
      Observable [ 0 ] = rapidity34;
    }

  else if ( glabel == "ZCHARM8TEV"  ||
	    glabel == "ZCHARM13TEV"   )
    {
      Observable [ 0 ] = rapidity34;
      Observable [ 1 ] = pt34;
    }
  
  else if ( glabel.substr(0,10) == "GammaC8TEV"   ||
	    glabel.substr(0,10) == "GammaC13TEV"    )
    {
      Observable [ 0 ] = rapidity4;
      Observable [ 1 ] = rapidity3;
      Observable [ 2 ] = pt3;
    }

  else if ( glabel == "ATLASWZRAP36PB_WP" ||
	    glabel == "CMSWEASY840PB_WP"  ||
	    glabel == "CMSWMASY47FB_WP"   ||
	    glabel == "LHCBWZ36PB_WP"     || 
	    glabel == "LHCBWZMU7TEV_WP"   ||
	    glabel == "LHCBWZMU8TEV_WP"     )
    {
      Observable [ 0 ] = rapidity4;
    }
  else if ( glabel == "ATLASWZRAP36PB_WM" ||
	    glabel == "CMSWEASY840PB_WM"  ||
	    glabel == "CMSWMASY47FB_WM"   ||
	    glabel == "LHCBWZ36PB_WM"     ||
	    glabel == "LHCBWZMU7TEV_WM"   ||
	    glabel == "LHCBWZMU8TEV_WM"     )
    {
      Observable [ 0 ] = rapidity3;
    }
  else if ( glabel == "ATLASZHIGHMASS49FB" ||
	    glabel.substr(0,13) == "ATLASLOMASSDY" )
    {
      Observable [ 0 ] = mass34;
    }
  else if ( glabel.substr(0,9) == "CMSDY2D11" ||
	    glabel.substr(0,9) == "CMSDY2D12" )
    {
      Observable [ 0 ] = std::fabs(rapidity34);
    }
  else if ( glabel == "CMSWCHARM_WPCB")
    {
      Observable [ 0 ] = std::fabs(rapidity4);
    }
  else if ( glabel == "CMSWCHARM_WMC")
    {
      Observable [ 0 ] = std::fabs(rapidity3);
    }

  else if ( glabel == "ATLASTTB11" )
    {
      Observable [ 0 ] = std::fabs(rapidity34);
      Observable [ 1 ] = pt3 + pt4;
      Observable [ 1 ] = pt34;
      Observable [ 3 ] = mass34;
      Observable [ 4 ] = rapidity34;
    }

  else if ( glabel == "CMSTTB12" )
    {
      Observable [ 0 ] = rapidity3;
      Observable [ 1 ] = rapidity34;
      Observable [ 2 ] = pt3;
      Observable [ 3 ] = pt34;
      Observable [ 4 ] = mass34;
      Observable [ 5 ] = rapidity34;
    }

  else if ( glabel == "ATLASPHOTONS11" )
    {
      Observable [ 0 ] = std::fabs(rapidity3);
      Observable [ 1 ] = pt3;
      Observable [ 2 ] = pt3;
    }

  else if ( glabel == "ATLASPHT12" )
    {
      Observable [ 0 ] = std::fabs(rapidity3);
      Observable [ 1 ] = pt3;
      Observable [ 2 ] = pt3;
      Observable [ 3 ] = pt3;
      Observable [ 4 ]  = pt3;
    }

  else if ( glabel == "ATLASPHT15" )
    {
      Observable [ 0 ] = std::fabs(rapidity3);
      Observable [ 1 ] = pt3;
      Observable [ 2 ] = pt3;
      Observable [ 3 ] = pt3;
      Observable [ 4 ]  = pt3;
    } 

 else if ( glabel == "HiggsDiff" )
      {
        Observable [ 0 ] = rapidity34;
      }
  else if ( glabel == "ppZ_nlo_8TeV"  || glabel == "ppZ_nlo_13TeV"  ||
	    glabel == "ppWp_nlo_8TeV" || glabel == "ppWp_nlo_13TeV" ||
	    glabel == "ppWm_nlo_8TeV" || glabel == "ppWm_nlo_13TeV" ||
	    glabel == "pptt_nlo_8TeV" || glabel == "pptt_nlo_13TeV" )
    {
      Observable [ 0 ] = pt34;
    }

  else if ( glabel == "CCbar8TEV" || glabel == "CCbar13TEV" )
    {
      Observable[0] = rapidity3;
      Observable[1] = pt3;
      Observable[2] = rapidity4;
      Observable[3] = pt4;
    }
  else if ( glabel == "ATLASZPT3AB")
    {
      Observable [ 0 ] = pt34;
      Observable [ 1 ] = mass34;
      // Observable [ 1 ] = pt34;
      //Observable [ 2 ] = rapidity34;
    }

  else
    {
      std::cout << "What the heck?!?!?!?! - DAMN" << std::endl;
    }

}

int cuts(int igrid)
{
  int fill = 0;
  switch(igrid)
    {
    case(0):
      if ( glabel == "CMSZDIFF12-BIN1" ) 
	{
	  if ( std::fabs(Observable[1])<0.4 ) fill = 1;
	}      
      else if ( glabel == "CMSZDIFF12-BIN2" )
        {
          if ( std::fabs(Observable[1])<0.8 ) fill = 1;
        }
      else if ( glabel == "CMSZDIFF12-BIN3" )
        {
          if ( std::fabs(Observable[1])<1.2 ) fill = 1;
        }
      else if ( glabel == "CMSZDIFF12-BIN4" )
        {
          if ( std::fabs(Observable[1])<1.6 ) fill = 1;
        }
      else if ( glabel == "CMSZDIFF12-BIN1" )
        {
          if ( std::fabs(Observable[1])<2.0 ) fill = 1;
        }
 	else 
	fill = 1;
      break;

    case(1):
       if ( glabel == "ATLASPHOTONS11" ) {
      	if ( std::fabs(Observable[0])<1.37 ) fill = 1;
      }
      else if (glabel == "ATLASPHT12") {
		if (std::fabs(Observable[0])<0.6 ) fill = 1;
	}

      else if (glabel == "ATLASPHT15") {
                if (std::fabs(Observable[0])<0.6 ) fill = 1;
        }
	else
      fill = 1;
      break;
    case(2):
      if ( glabel == "ATLASPHOTONS11" ) {
      	if ( std::fabs(Observable[0])>=1.52 && std::fabs(Observable[0])<2.37 ) fill = 1;
       }
      else if (glabel == "ATLASPHT12") {
		if (std::fabs(Observable[0])>=0.6 && std::fabs(Observable[0])<1.37) fill = 1;
	}

      else if (glabel == "ATLASPHT15") {
                if (std::fabs(Observable[0])>=0.6 && std::fabs(Observable[0])<1.37) fill = 1;
        }
	else
      fill = 1;
      break;
    case(3):
	  if (glabel == "ATLASPHT12") {
                if (std::fabs(Observable[0])>=1.56 && std::fabs(Observable[0])<1.81) fill = 1;
        }
	else if  (glabel == "ATLASPHT15") {
                if (std::fabs(Observable[0])>=1.56 && std::fabs(Observable[0])<1.81) fill = 1;
	}	
	else
      fill = 1;
      break;
    case(4):
        if ( glabel == "ATLASPHT12") {
                if (std::fabs(Observable[0])>=1.82 && std::fabs(Observable[0])<2.37) fill = 1;
        }
        else if ( glabel == "ATLASPHT15") {
                if (std::fabs(Observable[0])>=1.82 && std::fabs(Observable[0])<2.37) fill = 1;
        }
        else
      fill = 1;
      break;
    case(5):
      fill = 1;
      break;
    case(6):
      if (glabel == "ATLASZPT3AB") { //[12., 20., 30., 46., 66., 116., 150.};
	if (Observable[1]>=12. && Observable[1]<20.) fill = 1;
      }
      else
	fill = 1;
      break;
    case(7):
      if (glabel == "ATLASZPT3AB") {
        if (Observable[1]>=20. && Observable[1]<30.) fill = 1;
      }
      else
        fill = 1;
      break;
    case(8):
      if (glabel == "ATLASZPT3AB") {
        if (Observable[1]>=30. && Observable[1]<46.) fill = 1;
      }
      else
        fill = 1;
      break;
    case(9):
      if (glabel == "ATLASZPT3AB") {
        if (Observable[1]>=46. && Observable[1]<66.) fill = 1;
      }
      else
        fill = 1;
      break;
    case(10):
      if (glabel == "ATLASZPT3AB") {
        if (Observable[1]>=66. && Observable[1]<116.) fill = 1;
      }
      else
        fill = 1;
      break;
    case(11):
      if (glabel == "ATLASZPT3AB") {
        if (Observable[1]>=116. && Observable[1]<150.) fill = 1;
      }
      else
        fill = 1;
      break;
    default:
      std::cerr<<" In gridwrap.cpp::cuts(int). No such process : "<<igrid<<std::endl;
      std::exit(-1);
    }
  return fill;
}



// namespace mcfm_bridge;

/// function pointer hooks - set to 0 when no functions defined and applgrid not linked
extern void (*book_gridptr)();
extern void (*fill_gridptr)(const double evt[][mxpart] );
extern void (*write_gridptr)(double& );


extern "C" bool setup_mcfmbridge() {
  std::cout << "setup_mcfmbridge()" << std::endl;
  book_gridptr  = book_grid;
  fill_gridptr  = fill_grid;
  write_gridptr = write_grid;
  return true;
}

extern "C" bool setup_mcfmbridge_() {
  std::cout << "setup_mcfmbridge()" << std::endl;
  return setup_mcfmbridge();
}


bool mcfm_bridge_status = setup_mcfmbridge();
