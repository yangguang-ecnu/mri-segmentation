/*=========================================================================
  Program:   MATITK: Extending MATLAB with ITK
  Version:   2.4.02
  Language:  C++

  Copyright (c) Vincent Chu and Ghassan Hamarneh, Medical Image Analysis 
  Lab, Simon Fraser University. All rights reserved. 
  See http://www.cs.sfu.ca/~hamarneh/software/matitk/copyright.txt for
  details.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notices for more information.
 =========================================================================*/

#ifndef MATITK_H
#define MATITK_H
#include "itkImage.h"
#include "matrix.h"
#include "itkImportImageFilter.h"
#include <math.h>
#include "mex.h"
#define DIMENSION 3
#define LRESULTINDEX 0
#define ROPINDEX 0
#define RSEEDINDEX 4
#define RINPUTAINDEX 2
#define RINPUTBINDEX 3
#define RPARAMINDEX 1
#define IMPORTFILTERA 0
#define IMPORTFILTERB 1
#define RSPACINGINDEX 5

typedef double MATPARAMTYPE;
typedef double MATSEEDTYPE;

class FunctionCall{
	typedef void (*pt2Function)();
public:
	char* OpCode;
	char* OpName;
	pt2Function ptrFcn;
	FunctionCall(char* pstrzCode, char* pstrzName, pt2Function fcn){
		OpCode=pstrzCode;
		OpName=pstrzName;
		ptrFcn=fcn;
	}
};

template <class ITKPIXELTYPE>
class ITKPIXELTYPEArray{
	#include "typedefs.inl"
public:
	ITKPIXELTYPE* theArray;
	unsigned long numElem;
	bool needTranspose;
	ITKPIXELTYPEArray(ITKPIXELTYPE* ptrDoubleArray, unsigned long lnumberOfElements, bool bneedTranspose=true){
		theArray=ptrDoubleArray;
		numElem=lnumberOfElements;
		needTranspose = bneedTranspose;
	}
};
#endif