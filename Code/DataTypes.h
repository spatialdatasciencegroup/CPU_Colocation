#pragma once
#pragma once
#include "includes.h"
typedef unsigned int Integer;
typedef int SInteger;
typedef float Real;
struct sInstance {
	Integer featureType;
	Integer instanceID;
	//Integer cellID;
};
struct stuple {
	Integer size;
	vector<struct sInstance> ituple; //bounded by size tuple of instance
									 //vector<Integer> cellSortedID;
};
struct sColocations {
	Integer size;
	vector<string> featureList; //bounded by size
	vector<Integer> featureIndex;  //Index based on feature type
	vector<struct stuple> colocInstances; //instances of colocation
	vector<Integer> participation;
};

struct sNeighbor {
	vector<struct sInstance> neighbors;
	//vector<Integer> neigbborsIDList;
};

struct sFeatureStats {
	Integer start;
	Integer end;
	Integer count;
};

struct sLocations {
	Real x;
	Real y;
};


struct sColocationStats {
	Integer Degree;
	Integer start;
	Integer end;
};

struct sCellStats {
	Integer start;
	Integer end;
};

struct MyComparator
{
	const vector<Integer> & value_vector;

	MyComparator(const vector<Integer> & val_vec) :
		value_vector(val_vec) {}

	bool operator()(Integer i1, Integer i2)
	{
		return value_vector[i1] < value_vector[i2];
	}
};

struct gridDomain {
	Integer width;
	Integer height;
	Integer numCells;
	Integer MaxCellID;

};

struct param {
	Real thresholdDistance;
	Real PIthreshold;
	Integer FilterON_OFF; 
};

struct scalar2 {
	Real x;
	Real y;
};

struct Integer2 {
	Integer x;
	Integer y;

};

struct cBox {
	scalar2 minBound;
	scalar2 maxBound;

};
struct gridBox {
	cBox computeZone;
	Integer2 gridSize;
	scalar2 offsets;
	Integer totalCells;
	Real cellWidth;
};

