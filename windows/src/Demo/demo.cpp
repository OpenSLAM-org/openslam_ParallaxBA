//==============================================================================================================
// University of Technology, Sydney, Australia
// 
// Authors:  Liang Zhao         -- liang.zhao@imperial.ac.uk 
// 		  Shoudong Huang        -- Shoudong.Huang@uts.edu.au
//        Yanbiao Sun           -- syb51@pku.edu.cn
// 		  Gamini Dissanayake    -- Gamini.Dissanayake@uts.edu.au
// 
// 		  Centre for Autonomous Systems
// 
// 		  Faculty of Engineering and Information Technology
// 
// 		  University of Technology, Sydney
// 
// 		  NSW 2007, Australia
// 
// 		  License
// 
// 		  ParallaxBA by Liang Zhao, Shoudong Huang, Yanbiao Sun, Gamini Dissanayake is licensed under a 
// 
// 		  Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// 
// 		  Please contact Yanbiao Sun {syb51@pku.edu.cn} if you have any questions/comments about the code.
//==============================================================================================================

#include "stdafx.h"
#include "ParallaxBA.h"

#include <vector>
#include <map>
using namespace std;

#pragma  comment( lib, "../../bin/libamd.lib")
#pragma  comment( lib, "../../bin/libcamd.lib")	
#pragma  comment( lib, "../../bin/libccolamd.lib")
#pragma  comment( lib, "../../bin/libcholmod.lib")
#pragma  comment( lib, "../../bin/libcolamd.lib")
#pragma  comment( lib, "../../bin/libmetis_CHOLMOD.lib")
#pragma  comment( lib, "../../bin/libgoto2.lib")
#pragma  comment( lib, "../../bin/ParallaxBA.lib")

int _tmain(int argc, char* argv[] )
{
	IParallaxBA* pBA = newParallaxBA();

	//========================================================================================
	//Malaga dataset
 	char*  szCam = "../../../data/Malaga/Cam170.txt";
 	char*  szFea = "../../../data/Malaga/Feature170.txt";
 	char*  szCalib = "../../../data/Malaga/cal170.txt";
 	char*  szReport = "../../../data/Malaga/report.txt";
 
 	bool bLM     = true;   //true is Levenberg-Marquardt
   	pBA->pba_run( false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//========================================================================================
	//Village dataset
// 	  char*  szCam = "../../../data/Village/Cam90.txt";
// 	  char*  szFea = "../../../data/Village/Feature90.txt";
// 	  char*  szCalib = "../../../data/Village/cal90.txt";
// 	  char*  szReport = "../../../data/Village/report.txt";
// 
// 	  bool bLM     = false;   //true is Levenberg-Marquardt
//    pBA->pba_run( false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//========================================================================================
	//College dataset
// 	  char*  szCam = "../../../data/College/Cam468.txt";
// 	  char*  szFea = "../../../data/College/Feature468.txt";
// 	  char*  szCalib = "../../../data/College/cal468.txt";
// 	  char*  szReport = "../../../data/College/report.txt";
// 
// 	  bool bLM     = false;
//    pBA->pba_run( false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//========================================================================================
	//NewCollege dataset
// 	  char*  szCam = "../../../data/NewCollege/Cam3500.txt";
// 	  char*  szFea = "../../../data/NewCollege/Feature3500.txt";
// 	  char*  szCalib = "../../../data/NewCollege/cal3500.txt";
// 	  char*  szXYZ =  "../../../data/NewCollege/XYZ3500.txt";
// 	  char*  szReport = "../../../data/NewCollege/report.txt";
// 
// 	  bool bLM     = false;
//    pBA->pba_run( false, bLM, 100, szCam, szFea, szXYZ, szCalib, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//========================================================================================
	//Venice dataset
//    char*  szCam = "../../../data/Venice/Cam871.txt";
//    char*  szFea = "../../../data/Venice/Feature871.txt";
//    char*  szXYZ =  "../../../data/Venice/XYZ871.txt";
//    char*  szReport = "../../../data/Venice/report.txt";
//    
//    bool bLM     = true;
//    pBA->pba_run( false, bLM, 100, szCam, szFea, szXYZ, NULL, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//========================================================================================
	//Toronto dataset
// 	  char*  szCam = "../../../data/Toronto/Cam13.txt";
// 	  char*  szFea = "../../../data/Toronto/Feature13.txt";
// 	  char*  szCalib = "../../../data/Toronto/cal13.txt";
// 	  char*  szReport = "../../../data/Toronto/report.txt";
// 
// 	  bool bLM     = false;   //true is Levenberg-Marquardt
//    pBA->pba_run( false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//========================================================================================
	//Vaihingen dataset
//    char*  szCam = "../../../data/Vaihingen/Cam20.txt";
//    char*  szFea = "../../../data/Vaihingen/Feature20.txt";
//    char*  szCalib = "../../../data/Vaihingen/cal20.txt";
//    char*  szReport = "../../../data/Vaihingen/report.txt";
//    
//    bool bLM     = false;   //true is Levenberg-Marquardt
//    pBA->pba_run( false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//========================================================================================
	//DunHuan dataset
//    char*  szCam = "../../../data/DunHuan/Cam63.txt";
//    char*  szFea = "../../../data/DunHuan/Feature63.txt";
//    char*  szCalib = "../../../data/DunHuan/cal63.txt";
//    char*  szReport = "../../../data/DunHuan/report.txt";
//    
//    bool bLM     = false;   //true is Levenberg-Marquardt
//    pBA->pba_run( false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//========================================================================================
	//Jinan dataset
//    char*  szCam = "../../../data/Jinan/Cam76.txt";
//    char*  szFea = "../../../data/Jinan/Feature76.txt";
//    char*  szCalib = "../../../data/Jinan/cal76.txt";
//    char*  szReport = "../../../data/Jinan/report.txt";
//     
//    bool bLM     = true;   //true is Levenberg-Marquardt
//    pBA->pba_run( false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL, 1E-3 );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	//========================================================================================
	//Taian dataset
// 	char*  szCam = "../../../data/Taian/Cam737.txt";
// 	char*  szFea = "../../../data/Taian/Feature737.txt";
// 	char*  szCalib = "../../../data/Taian/cal737.txt";
// 	char*  szReport = "../../../data/Taian/report.txt";
// 	
// 	bool bLM     = false;   //true is Levenberg-Marquardt
// 	pBA->pba_run( false, bLM, 100, szCam, szFea, NULL, szCalib, szReport, NULL, NULL, 1E-3 );
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	freeParallaxBA( pBA );
	system( "pause" );

	return 0;
}

