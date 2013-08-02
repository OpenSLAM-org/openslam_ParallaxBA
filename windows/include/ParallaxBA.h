//==============================================================================================================
// University of Technology, Sydney, Australia
// 
// Authors:  Liang Zhao         -- Liang.Zhao-1@uts.edu.au 
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

#ifndef _PARALLAXBA_H
#define _PARALLAXBA_H

#ifdef PARALLAXBA_EXPORTS
	#define ParallaxBAapi  __declspec(dllexport)
#else
	#define ParallaxBAapi __declspec(dllimport)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/StdVector>
#include <Eigen/Cholesky>

#include "cholmod/cholmod.h"

#define PI  3.1415926535898 
#define MAXARCHOR 0.5

typedef enum PBA_FixType
{	
	PBA_FixDefault = -1, //fix the longest axis  
	PBA_FixX = 0 ,	    // fix X axis	
	PBA_FixY = 1 ,	    // fix Y axis
	PBA_FixZ = 2 ,	    // fix Z axis
} FixType ; 
	
class ParallaxBAapi IParallaxBA
{
public:
	virtual bool pba_run( bool bRobust, bool bLM, int nMaxIter, char* szCam, char* szFea, char* szXYZ = NULL, char* szCalib = NULL, char* szReport = NULL, 
		char* szPose = NULL, char* sz3D = NULL, double Tau = 1E-6 ) = 0;

	virtual	bool pba_run( int argc, char** argv ) = 0;

	virtual bool pba_initialize( char* szCamera, char* szFeature, char* szCalib =  NULL, char* szXYZ = NULL ) = 0;

	virtual bool pba_motstr_levmar( ) = 0;

	virtual bool pba_motstr_gn( FixType ft = PBA_FixDefault ) = 0;
};

typedef IParallaxBA IParallaxBAPtr;
ParallaxBAapi IParallaxBA*	newParallaxBA();
ParallaxBAapi void    freeParallaxBA( IParallaxBA* ptr );

#endif