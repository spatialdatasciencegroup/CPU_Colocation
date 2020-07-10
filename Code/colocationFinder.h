#pragma once
#include "DataTypes.h"
#include <limits>
#include <string>
#include <chrono>
#include <math.h>
#include <iomanip>
#include <ctime>
using  sec = chrono::seconds;
using get_time = chrono::steady_clock;
class colocationFinder {

public:
	colocationFinder(void);

public:
	~colocationFinder(void);

private:
	//Host Variables
	size_t maxInstancesNum;
	size_t maxFeaturesNum;
	size_t maxCellNum;
	std::string location;
	struct gridBox gridStructure;
	struct gridDomain grid;
	vector <string> featureTypes;   //done
	Integer * instanceList;         //done
	//vector<struct sInstance> instanceList;  //sort by cell id after degree 2 and think on for each sorted groups sort them by type
	vector<struct sInstance> cellSortedinstanceList;
	//vector <struct sFeatureStats> featureStats;
	Integer *featureInstanceStart;  //done
	Integer * featureInstanceEnd;	//done
	Integer * featureInstanceCount; //done
	Integer* instanceCellIDs;  
	Real * instanceLocationX;    //done
	Real * instanceLocationY;    //done
	vector<Integer> candiColocations;
	//Integer *prevalantColocations;
	//Integer **prevalentInstanceTables;
	//Integer *prevalentInstanceTableSize;
	//Integer **candidateInstanceTables;

	vector<vector<Integer>> h_prevalentInstanceTables;
	vector<Integer> h_prevalantColocations;
	Integer h_prevelentColocationCount;
	vector<Integer> h_prevalentInstanceTableSize;

	//Integer** h_prevalentInstanceTable2;
	vector<vector<Integer>> h_prevalentInstanceTables2;
	vector<Integer> h_prevalantColocations2;
	Integer h_prevelentColocationCount2;
	vector<Integer> h_prevalentInstanceTableSize2;

	Integer candiColocCounter;
	struct param parameter;
	//for subset
	vector<Integer> subset;
	vector<vector<Integer>> subsetList;

	//for MRF
	vector<Integer> coarseInstanceStart;  //done
	vector<Integer> coarseInstanceEnd;	//done
	vector<Integer> coarseInstanceList;


	//vector<Integer> coarsePrevalance;
	//Integer ** coarsePrevalentInstanceTables;
	//Integer *coarsePrevalentInstanceTableSize;
	vector<vector<Integer>> h_coarsePrevalentInstanceTables;
	vector<Integer> h_coarsePrevalentInstanceTableSize;

	vector<vector<Integer>> h_coarsePrevalentInstanceTables2;
	vector<Integer> h_coarsePrevalentInstanceTableSize2;

	Integer* cellEventInstanceCount;
	vector<Integer> coarseInstanceCount;

	//vector<vector<Integer>> h_coarsePrevalentInstanceTables2;
	//Integer** h_prevalentCoarseInstanceTable2;
	//vector<Integer> h_prevalentCoarseInstanceTableSize2;



	//added in version 1_2
	SInteger* cellEventInstanceStart;
	SInteger* cellEventInstanceEnd;
	Integer* cellBasedSortedIndex;

	bool needInstanceTable;


	//LOG VARIABLES
	std::string logFileName;
	std::ofstream timeLogStream;
	std::string compLogFileName;
	std::ofstream compLogStream;

	std::string outputLocation;
	std::string outputFile;
	Integer totalCandidatePatterns;
	Integer totalFilteredPatterns;
	Real totalFilterTime;
	Real totalRefineTime;
	Real degree2FilterTime;
	Real degree2RefineTime;
	Real memoryTracker;
public:
	void Begin(int argc, char**argv);
	void candidateColocationGeneral(Integer degree);
	void initializeHash();
	void populateData(std::string datasetFilename);
	void setFeatureStats(Integer start, Integer end);
	struct sColocations helper2degCandidateColoc(string f1, string f2, Integer i, Integer j);
	struct sInstance getInstance(string value, string cellid);
	struct sNeighbor populateInstanceNeighborList(string nlist, string cellids);
	Integer instanceIDLookup(struct sInstance ins);
	void degree2CandidateColocationGenerator();
	void colocationInstanceGeneration();
	bool generateDegreeNCandidateColoc();
	void setColocStats(Integer degree, Integer start, Integer end);
	bool generateandCheckSubset(std::vector<Integer> &members);
	Integer getColocationIndex(std::vector<Integer> &features, Integer degreeIndex);
	std::vector<Integer> getNeighborCells(Integer cellID);
	Real DistanceInMeters(Integer idx1, Integer idx2);
	bool generateandCheckSubsets(vector<Integer> &inter, Integer degree);
	void subsetGen(vector<Integer> &inter, int k, int n, int idx);
	bool checkSubset(vector<Integer> subsetElem);
	void generateInstanceTableGeneral(Integer degree);
	void degree2Processing();
	void extractPrevalentColocations(Integer degree);
	void clearSubsetVectors();
	Integer getIndex(vector <Integer> inner,Integer degree);
	void resetPrevalentData();
	void resetCandidateData();
	void resetBitmapData();
	void savetoFile(Integer degree);
	void LoadParameter(std::string parameterFileName);
		void setGrid(Real maxX, Real maxY, Real minX, Real minY);
	Integer getCellID(Real lon, Real lat);
	Integer getTypeID(Integer i);
	bool checkifNeighbors(Integer cell1, Integer cell2);
	void extractCoarsePrevalentColocations(Integer degree);
	SInteger getNeighborCellID(Integer x, Integer y, Integer TID);
	//log related functions
	void timeLog(Integer degree, std::string function, std::string eventType, Integer duration);
	void loggingFileOpenFunction();
	void loggingFileCloseFunction();
	void compLog(Integer i, Integer degree, Integer candicolocNumber, std::string Remark, std::string PIType, Real PI);
	void compLog(Integer degree, std::string function);
	void compLog(Integer degree, Integer i, std::string totalmem, Integer candicolocNumber);
	Real getUpperbound(vector<Integer> &coarseBitmap, Integer degree, vector<Integer> &currentColoc);
	Real getParticipationIndex(vector<Integer> &bitmap, Integer degree, vector<Integer> &currentColoc);
	void copyPrevalentColocations(Integer degree);
	void clearMemory();
	void savetoFileGen(Integer degree);

	void generateInstanceTableGeneral2(Integer degree);
	void tableGenRequired(Integer degree);
	bool generateandCheckSubsets2(vector<Integer> &inter, Integer degree);
	bool checkSubset2(vector<Integer> subsetElem);

};
