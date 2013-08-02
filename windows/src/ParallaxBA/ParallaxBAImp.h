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

#pragma once
#include "ParallaxBA.h"
#include <map>
#include <vector>
#include <set>
using namespace std;
using namespace Eigen;

//use sba_crsm structure from SBA (http://www.ics.forth.gr/~lourakis/sba/) to store sparse matrix
struct sba_crsm	
{
    int nr, nc;   
    int nnz;      
    int *val;     
    int *colidx;  
    int *rowptr;  
};

class CParallaxBA : public IParallaxBA
{
public:
	CParallaxBA(void);
	~CParallaxBA(void);
	virtual bool CParallaxBA::pba_run( bool bRobust, bool bLM, int nMaxIter, char* szCam, char* szFea, char* szXYZ = NULL, char* szCalib = NULL, char* szReport = NULL, 
		char* szPose = NULL, char* sz3D = NULL, double Tau = 1E-6 );

	virtual	bool pba_run( int argc, char** argv );

	virtual bool pba_initialize( char* szCamera, char* szFeature, char* szCalib = NULL, char* szXYZ = NULL );

	virtual bool pba_motstr_levmar( );

	virtual bool pba_motstr_gn(  FixType ft = PBA_FixDefault  );

private:
	//initialize feature. You can provide xyz or system also provide them by itself;
	bool    pba_initializeMainArchor( double* imgpts, double* camera,double* K,double* feature, int nP, int FID, double* KR );
	bool    pba_initializeAssoArchor( double* imgpts, int* photo, double* camera,double* K,double* feature,int nMI, int nAI, int FID, bool bLast );
	bool	pba_initializeOtheArchors( double* imgpts, int* photo, double* camera,double* K,double* feature,int* archorSort,int nfeacout, int nOI, int FID );

	//update KR KdA KdB KdG
	void	pba_updateKR( double *KR, double *KdA, double *KdB, double *KdG, double *K, double *p );

	//compute reprojection error
	void	pba_cost(double *p, double *hx, int* archor );

	//compute Jacobian, not save Jp, Jc etc
	void	pba_jacobian(double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature );
	void	pba_jacobian_RobustKernel( double *p, int* archor,struct sba_crsm *Uidxij, double*e,double* U,double* ea, double* V, double* eb, double *W, int n, int m, int numprojs, int *nphoto, int *nfeature );

	//construct S matrix, S = U - W*V^-1*W^T
	void	pba_constructSLM( double* S, double* E, double* U, double* V, double* W, double* ea, double* eb, sba_crsm& Sidxij, double mu );
	void	pba_constructSGN( double* S, double* E, double* U, double* V, double* W, double* ea, double* eb, sba_crsm& Sidxij );

	//prepare for ordering of CHOLMOD
	void	pba_constructAuxCSSLM( int *Ap, int *Aii );
	void	pba_constructAuxCSSGN( int *Ap, int *Aii );

	//set CSS format
	void	pba_constructCSSLM( int* Si, int* Sp, double* Sx, double* S, cholmod_sparse* m_cholSparseS, sba_crsm& Sidxij, bool init);
	void	pba_constructCSSGN( int* Si, int* Sp, double* Sx, double* S, cholmod_sparse* m_cholSparseS, sba_crsm& Sidxij, bool init, int nft);

	//solve sparse system using CHOLMOD
	bool	pba_solveCholmodLM( int* Ap, int* Aii, bool init, bool ordering);
	bool	pba_solveCholmodGN( int* Ap, int* Aii, bool init, bool ordering);

	//compute feature parameters
	void	pba_solveFeatures( double *W, double *IV, double *ea, double *eb, double *dpa, double *dpb);

private:
	//compute smask of S matrix
	int		pba_ConstructSmask( sba_crsm& Sidxij, sba_crsm& Uidxij );

	//compute inverse V matrix
	void	pba_inverseVLM( double* V, double* IV, sba_crsm& Uidxij, double mu );
	void	pba_inverseVGN( double* U, double* V, sba_crsm& Uidxij );

	//compute initial mu of LM
	double	pba_computeInitialmu( double* U, double* V, sba_crsm& Uidxij, double tau, int nvars );

	//transform angle into XYZ
	int		pba_angle2xytGN( double *p );
	int		pba_angle2xytLM( double *p );
	void	pba_saveXYZ( const char* camera, const char* sz3Dpt, double *p, bool gn = true );

	bool	pba_parseArgs( int argc, char* argv[] );
	void	pba_printHelp();

private:
	void	pba_readAndInitialize( char *camsfname, char *ptsfname, int *ncams, int *n3Dpts, int *n2Dprojs,double **motstruct, 
		double **imgpts, int **archor, char **vmask, char **umask,int **nphoto, int** nfeature, int** archorSort );

	void	pba_readCameraPoseration( char *fname, double* ical );

	void	pba_readProjectionAndInitilizeFeature(	FILE *fp, double *params, double *projs, char *vmask, int ncams, 
		int *archor,char* umask,int* nphoto, int* nfeature, int* archorSort );

	void	pba_readCameraPose(FILE *fp, double *params);

	void	pba_readCablibration(FILE* fp, double *K );

private:
	int		findNcameras(FILE *fp);
	int		countNDoubles(FILE *fp);
	int		readNInts(FILE *fp, int *vals, int nvals);
	int		readNDoubles(FILE *fp, double *vals, int nvals);
	int		skipNDoubles(FILE *fp, int nvals);	
	void	readNpointsAndNprojections(FILE *fp, int *n3Dpts, int pnp, int *nprojs, int mnp);
	double	nrmL2xmy(double *const e, const double *const x, const double *const y, const int n);
	void	readNpointsAndNprojectionsFromProj(FILE *fp, int &n3Dpts, int &nprojs);
	void	readPointProjections(FILE *fp,double *imgpts, int *photo,int* imgptsSum, int n3Dpts, int n2Dprojs );
	void    readImagePts( const char* szProj, double **imgpts, int **photo,int** imgptsSum, int &n3Dpts, int &n2Dprojs );

private:
	int		m_ncams, m_n3Dpts, m_n2Dprojs, m_nS;  //number of camera, 3D points, 2D projection points, non-zero element of S matrix
	int		*m_archor;							  
	int		*m_photo, *m_feature;	  			  
	double	*m_motstruct, *m_imgpts;			  //6 camera pose and 3 feature parameters
	double	*m_XYZ;								  //initial XYZ provided 	
	double  *m_K;								  //calibration parameters

	char	*m_vmask, *m_umask, *m_smask;		  
	int		*m_imgptsSum, *m_struct, *m_pnt2main, *m_archorSort;
	double	*m_KR, *m_KdA, *m_KdB, *m_KdG;		 

	bool    m_bProvideXYZ, m_bFocal;	
	char*	m_szCameraInit;
	char*	m_szFeatures;
	char*   m_szCalibration;
	char*	m_szXYZ;
	char*   m_szCamePose;
	char*   m_sz3Dpts;
	char*   m_szReport;
	int		m_nMaxIter;
	double  m_Tau, m_e1, m_e2, m_e3, m_e4;
	bool	m_bRobustKernel;
	bool	m_bsolverLM;
	bool	m_bsolverGN;
	int		m_nRobustType;
	double  m_delt;

private:
	//Solve Sparse Matrix using CHOLMOD (http://www.cise.ufl.edu/research/sparse/SuiteSparse/) 
	cholmod_sparse *m_cholSparseS;				
	cholmod_factor *m_cholFactorS; 
	cholmod_common m_cS; 
	cholmod_dense  *m_cholSparseR, *m_cholSparseE;
};

