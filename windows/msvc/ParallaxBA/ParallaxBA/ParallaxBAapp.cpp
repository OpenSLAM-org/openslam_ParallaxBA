// ParallaxBA.cpp : Defines the exported functions for the DLL application.

#include "stdafx.h"
#include "ParallaxBA.h"
#include "../../../src/ParallaxBA/ParallaxBAImp.h"

ParallaxBAapi IParallaxBA*	newParallaxBA()
{
	return (IParallaxBA*)(new CParallaxBA );
}

ParallaxBAapi void    freeParallaxBA( IParallaxBA* ptr )
{
	free( ptr );
}