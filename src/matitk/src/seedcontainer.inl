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

template <class ITKPIXELTYPE>
class SeedContainer{
	#include "typedefs.inl"
private:
	int nSeeds;
	typename InternalImageType::IndexType** ppSeedList;
public:
	SeedContainer(){
		nSeeds=0;
	}
	void set
		(ITKPIXELTYPE* pfseeds, int arraySize, int* volumeDimension)
	{
		nSeeds=(arraySize/3);
		ppSeedList=new InternalImageType::IndexType*[nSeeds];
		for (int i=0; i<arraySize;i++){
			if (pfseeds[i]<=0) mexErrMsgTxt("Please note that array in matlab starts from 1");
		}
		for (int i=0; i<nSeeds; i++){
			ppSeedList[i]=new InternalImageType::IndexType();
			//this is done to correct order, and the fact that matlab starts from 1, and c starts from 0
			(*(ppSeedList[i]))[0]=pfseeds[i*3+1]-1;
			(*(ppSeedList[i]))[1]=pfseeds[i*3]-1;
			(*(ppSeedList[i]))[2]=pfseeds[i*3+2]-1;
			if (pfseeds[i*3+1]>=volumeDimension[0] || 
				pfseeds[i*3]>=volumeDimension[1] || pfseeds[i*3+2]>=volumeDimension[2])
				mexErrMsgTxt("Location of seed outside volume A image space");
		}
	}
	int getNumberOfSeeds(){
		return nSeeds;
	}
	typename InternalImageType::IndexType& getIndex(int i){
		if (i>=nSeeds) mexErrMsgTxt("You did not specifiy enough seed points");
		return *(ppSeedList[i]);
	}

	~SeedContainer(){
		for (int i=0; i<nSeeds; i++) delete ppSeedList[i];
		delete [] ppSeedList;
	}
};