#pragma warning(disable : 4996)
#include "includes.h"
#include "DataTypes.h"
#include "colocationFinder.h"
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/compare.hpp>
#include<iterator>
#include<algorithm>
#include <numeric>

using namespace boost;
using namespace std;

colocationFinder::colocationFinder(void) {
	totalCandidatePatterns = 0;
	totalFilteredPatterns = 0;
	totalFilterTime = 0.0;
	totalRefineTime = 0.0;
	degree2FilterTime = 0.0;
	degree2RefineTime = 0.0;
}

void colocationFinder::LoadParameter(std::string parameterFileName) {
	std::string FileName= /*location+*/ parameterFileName;
	std::ifstream ifs(FileName);
	if (!ifs)
	{
		std::cout << "Parameter file missing..." << std::endl;
		exit(0);
		//return CCT_FILEERR;
	}
	ifs >> parameter.thresholdDistance;
	ifs >> parameter.PIthreshold;
	ifs >> parameter.FilterON_OFF;
	std::cout<<"Distacne Threshold  = "<< parameter.thresholdDistance<< endl<< "Participation Index Threshold = "<<parameter.PIthreshold <<endl <<"Filter Switch = "<< parameter.FilterON_OFF <<endl;
	ifs.close();
}
void colocationFinder::degree2CandidateColocationGenerator() {
	Integer counter = 0;
	for (int i = 0; i < featureTypes.size() - 1; i++) {
		for (int j = i + 1; j < featureTypes.size(); j++) {
			//candiColocations[counter++] = i;
			//candiColocations[counter++] = j;
			candiColocations.push_back(i);
			candiColocations.push_back(j);
		}
	}
}
int mode(int a[], int n) {
	int maxValue = 0, maxCount = 0, i, j;

	for (i = 0; i < n; ++i) {
		int count = 0;

		for (j = 0; j < n; ++j) {
			if (a[j] == a[i])
				++count;
		}

		if (count > maxCount) {
			maxCount = count;
			maxValue = a[i];
		}
	}

	return maxValue;
}
Integer colocationFinder::getTypeID(Integer i) {
	Integer Typeid = -1;
	for (size_t j = 0; j < maxFeaturesNum; j++) {
		if (i >= featureInstanceStart[j] && i <= featureInstanceEnd[j]) {
			Typeid = j;
			break;
		}
	}
	return Typeid;
}


void colocationFinder::setGrid(Real maxX, Real maxY, Real minX, Real minY) {
	gridStructure.computeZone.minBound.x = minX;
	gridStructure.computeZone.minBound.y = minY;
	gridStructure.computeZone.maxBound.x = maxX;
	gridStructure.computeZone.maxBound.y = maxY + parameter.thresholdDistance;
	gridStructure.cellWidth = parameter.thresholdDistance;
	scalar2 gs;
	gs.x = (gridStructure.computeZone.maxBound.x - gridStructure.computeZone.minBound.x) / gridStructure.cellWidth;
	gs.y = (gridStructure.computeZone.maxBound.y - gridStructure.computeZone.minBound.y) / gridStructure.cellWidth;
	gridStructure.gridSize.x = (Integer)ceil(gs.x);
	gridStructure.gridSize.y = (Integer)ceil(gs.y);
	gridStructure.totalCells = gridStructure.gridSize.x * gridStructure.gridSize.y;
	cellEventInstanceCount = new Integer[gridStructure.totalCells*maxFeaturesNum];
	cellBasedSortedIndex = new Integer[maxInstancesNum];
	instanceCellIDs=  new Integer[maxInstancesNum];
	for (size_t i = 0; i < gridStructure.totalCells*maxFeaturesNum; i++) {
		cellEventInstanceCount[i] = 0;
	}
	vector< pair<Integer, Integer> >v;
	for (size_t i = 0; i < maxInstancesNum; i++) {
		Integer typeID = instanceList[i];
		Integer cellID = getCellID(instanceLocationX[i], instanceLocationY[i]);
		instanceCellIDs[i] = cellID;
		Integer value = cellEventInstanceCount[cellID*maxFeaturesNum + typeID];
		value++;
		cellEventInstanceCount[cellID*maxFeaturesNum + typeID] = value;
		v.push_back(make_pair(cellID, i));
	}
	std::sort(std::begin(v), std::end(v));
	Integer cnt = 0;
	for (Integer t = 0; t < v.size();t++) {
		cellBasedSortedIndex[cnt] = v[t].second;
		cnt++;
	}
	v.clear();
	v.shrink_to_fit();
		
	size_t posCounter = 0;
	size_t maxEnd = 0;
	cellEventInstanceStart = new SInteger[gridStructure.totalCells*maxFeaturesNum];
	cellEventInstanceEnd = new SInteger[gridStructure.totalCells*maxFeaturesNum];
	for (size_t CID = 0; CID < gridStructure.totalCells; CID++) {
		for (size_t FID = 0; FID < maxFeaturesNum; FID++) {
			Integer countValue = cellEventInstanceCount[CID*maxFeaturesNum + FID];
			if (countValue > 0) {
				cellEventInstanceStart[CID*maxFeaturesNum + FID] = posCounter;
				cellEventInstanceEnd[CID*maxFeaturesNum + FID] = posCounter + countValue - 1;
				maxEnd = posCounter + countValue - 1;
				posCounter = posCounter + countValue;
			}
			else {
				cellEventInstanceStart[CID*maxFeaturesNum + FID] = -1;
				cellEventInstanceEnd[CID*maxFeaturesNum + FID] = -1;
			}
		}
	}
	size_t insCounter = 0;
	size_t insPerTypeCounter = 0;
	for (size_t t = 0; t < maxFeaturesNum; t++) {
		coarseInstanceStart.push_back(insCounter);
		size_t insPerTypeCounter = 0;
		for (size_t k = 0; k < gridStructure.totalCells; k++) {
			if (cellEventInstanceCount[k*maxFeaturesNum + t] > 0) {
				coarseInstanceList.push_back(k);
				insCounter++;
				insPerTypeCounter++;
			}
		}
		coarseInstanceEnd.push_back(insCounter -1);
		coarseInstanceCount.push_back(insPerTypeCounter);
	}

}
Integer colocationFinder::getCellID(Real x, Real y) {
	Real RelativeX = x - gridStructure.computeZone.minBound.x;
	Real RelativeY = y - gridStructure.computeZone.minBound.y;

	RelativeX /= gridStructure.cellWidth;
	RelativeY /= gridStructure.cellWidth;
	Integer i = (Integer)RelativeX;
	Integer j = (Integer)RelativeY;
	const Integer CellID = (gridStructure.gridSize.x * j) + i;
	return CellID;
}
void colocationFinder::populateData(std::string datasetFilename) { //pass data.csv as argument

	//variable decleration and defination
	vector<string> vec;
	string line;
	Integer entryCounter = 0;
	Integer eventTypeCounter = 0;
	Integer instanceCounter = 0;
	Integer lastEventID = -1;
	vector<struct sFeatureStats> featureStats;
	Real maxX, maxY, minX, minY;
	string data(/*location+*/datasetFilename);

	//File Reading
	ifstream in(data.c_str());
	typedef tokenizer< escaped_list_separator<char> > Tokenizer;
	//getline(in, line); // for header
	getline(in, line); // for instance count
	cout<< "ok";
	maxInstancesNum = stoi(line);
	cout<<"here";
	instanceLocationX = new Real[maxInstancesNum];
	instanceLocationY = new Real[maxInstancesNum];
	instanceList = new Integer[maxInstancesNum];
	getline(in, line);
	while (!line.empty())
	{
		Tokenizer tok(line);
		vec.assign(tok.begin(), tok.end());
		if (lastEventID == -1) {
			lastEventID = stoi(vec[0]);
			featureTypes.push_back(vec[0]);
			eventTypeCounter++;
		}
		else if (stoi(vec[0]) != lastEventID) {
			struct sFeatureStats tempStat;
			tempStat.start = entryCounter - instanceCounter;
			tempStat.end = entryCounter - 1;
			tempStat.count = instanceCounter;
			featureStats.push_back(tempStat);
			lastEventID = stoi(vec[0]);
			instanceCounter = 0;
			featureTypes.push_back(vec[0]);
			eventTypeCounter++;
		}
		instanceList[entryCounter] = eventTypeCounter - 1;
		instanceLocationX[entryCounter] = stof(vec[1]);
		instanceLocationY[entryCounter] = stof(vec[2]);
		entryCounter++;
		instanceCounter++;
		getline(in, line);
	}
	getline(in, line);
	Tokenizer tok(line);
	vec.assign(tok.begin(), tok.end());
	minX = stof(vec[0]);
	minY = stof(vec[1]);
	maxX = stof(vec[2]);
	maxY = stof(vec[3]);
	//for last feature
	struct sFeatureStats tempStat;
	tempStat.start = entryCounter - instanceCounter;
	tempStat.end = entryCounter - 1;
	tempStat.count = instanceCounter;
	featureStats.push_back(tempStat);

	//converting vector to array
	Integer featureSize = featureStats.size();
	featureInstanceStart = new Integer[featureSize];
	featureInstanceEnd = new Integer[featureSize];
	featureInstanceCount = new Integer[featureSize];
	for (size_t i = 0; i < featureStats.size(); i++) {
		featureInstanceStart[i] = featureStats[i].start;
		featureInstanceEnd[i] = featureStats[i].end;
		featureInstanceCount[i] = featureStats[i].count;
	}
	maxFeaturesNum = featureSize;
	in.close();
	setGrid(maxX, maxY, minX, minY);
}
Real colocationFinder::DistanceInMeters(Integer idx1, Integer idx2) {
	Real x2 = instanceLocationX[idx2];
	Real x1 = instanceLocationX[idx1];
	Real y2 = instanceLocationY[idx2];
	Real y1 = instanceLocationY[idx1];
	Real distance = sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2));
	return distance;
}
bool colocationFinder::checkifNeighbors(Integer cell1, Integer cell2) {  //for MRF
	//if same cells
	if (cell1 == cell2) {
		return true;
	}
	SInteger cell1Y = cell1 / gridStructure.gridSize.x;
	SInteger cell1X = cell1 - (gridStructure.gridSize.x*cell1Y);
	SInteger Range = 1;
	for (SInteger posj = -Range; posj <= Range; posj++) {
		for (SInteger posi = -Range; posi <= Range; posi++) {
			SInteger NCID = getNeighborCellID(cell1X + posi, cell1Y + posj, cell1);
			if (NCID < 0 || NCID >= gridStructure.totalCells)
			{
				continue;
			}
			if (NCID == cell2) {
				return true;
			}
		}
	}
	return false;
}

SInteger colocationFinder::getNeighborCellID(Integer x, Integer y, Integer TID) {
	SInteger NCID;
	NCID = y*gridStructure.gridSize.x + x;
	if (TID%gridStructure.gridSize.x == 0) {
		if ((NCID + 1) % gridStructure.gridSize.x == 0) {
			return -1;
		}
		return NCID;
	}
	else if ((TID + 1) % gridStructure.gridSize.x == 0) {
		if (NCID % gridStructure.gridSize.x == 0) {
			return -1;
		}
		return NCID;
	}
	else {
		return NCID;
	}
}
Real colocationFinder::getUpperbound(vector<Integer> &coarseBitmap, Integer degree, vector<Integer> &currentColoc) {
	Real PI = 1.0;
	Integer from = 0;
	Integer to = 0;
	for (Integer j = 0; j < degree; j++) {
		Integer index = currentColoc[j];
		if (j == 0) {
			from = 0;
			to = gridStructure.totalCells;
		}
		else {
			from = to - 1;
			to = from + gridStructure.totalCells + 1;
		}
		//size_t bitSum = thrust::reduce(d_bitmapPtr + from, d_bitmapPtr + to - 1);
		size_t bitSum = std::accumulate(coarseBitmap.begin() + from, coarseBitmap.begin() + to - 1,0);
		Real pr = bitSum / (Real)featureInstanceCount[index];
		if (pr < PI) {
			PI = pr;
		}
	}
	//thrust::device_free(d_bitmapPtr);
	return PI;
}
Real colocationFinder::getParticipationIndex(vector<Integer> &bitmap, Integer degree, vector<Integer> &currentColoc) {
	Real PI = 1.0;
	Integer from;
	Integer to;
	for (Integer j = 0; j < degree; j++) {
		Integer index = currentColoc[j];
		Integer totalInstance = featureInstanceCount[index];
		if (j == 0) {
			from = 0;
			to = totalInstance;
		}
		else {
			from = to - 1;
			to = from + totalInstance + 1;
		}
		size_t bitSum = std::accumulate(bitmap.begin() + from, bitmap.begin() + to - 1, 0);
		//size_t bitSum = thrust::reduce(d_bitmapPtr + from, d_bitmapPtr + to - 1);
		Real pr = bitSum / (Real)totalInstance;
		if (pr < PI) {
			PI = pr;
		}
	}
	//thrust::device_free(d_bitmapPtr);
	return PI;
}

void colocationFinder::degree2Processing() {
	memoryTracker = 0.0;
	std::cout << "\n\n\nFor degree 2......" << std::endl;
	auto start = get_time::now();
	//timeLog(2, "candidate colocgen", "start",0);
	Integer degree = 2;
	h_prevelentColocationCount = 0;
	candiColocCounter = featureTypes.size() * (featureTypes.size() - 1) / 2;
	totalCandidatePatterns += candiColocCounter;
	degree2CandidateColocationGenerator();
	auto end = get_time::now();
	auto diff = end - start;
	//timeLog(2, "candidate colocgen " + std::to_string(candiColocCounter), "end", chrono::duration_cast<sec>(diff).count());
	std::cout << "degree2CandidateColocationGenerator() took  " << chrono::duration_cast<sec>(diff).count() << " seconds " << endl;
	std::cout << "MultiResolution Pruning Starts... " << "Degree: " << degree << "candidate Colocations: " << candiColocCounter << endl;
	//start = get_time::now();
	//timeLog(2, "multiResolution pruning ", "start", 0);
	//compLog(2, "MultiResolution Filter Phase");
	vector<Integer> coarseIntermediate;
	Integer filterpruneCounter = 0;
	Integer totalCoarseInstances = 0;
	Integer totalInstances = 0;
	vector<Real> upperboundList;
	vector<Real> PIList;
	for (Integer i = 0; i < candiColocCounter; i++) {
		std::cout << "MRF Degree:" << degree << " Working on: " << i + 1 << "/" << candiColocCounter << std::endl;
		size_t firstFeatureIndex = candiColocations[i*degree];
		size_t secondFeatureIndex = candiColocations[i*degree + 1];
		vector<Integer> CurrentColoc;
		CurrentColoc.push_back(firstFeatureIndex);
		CurrentColoc.push_back(secondFeatureIndex);
		//MultiResolution Filter start
		if (parameter.FilterON_OFF == 1) {
			auto fstart = get_time::now();
			coarseIntermediate.clear();
			coarseIntermediate.shrink_to_fit();
			vector<Integer>coarseBitmap(degree*gridStructure.totalCells, 0);
			size_t start = coarseInstanceStart[firstFeatureIndex];
			size_t end = coarseInstanceEnd[firstFeatureIndex];
			size_t secondstart = coarseInstanceStart[secondFeatureIndex];
			size_t secondend = coarseInstanceEnd[secondFeatureIndex];
			Integer counter = 0;
			//Integer Ps = secondstart;
			for (size_t idx1 = start; idx1 <= end; idx1++) {
				Integer coarseInstance1 = coarseInstanceList[idx1];
				SInteger cell1Y = coarseInstance1 / gridStructure.gridSize.x;
				SInteger cell1X = coarseInstance1 - (gridStructure.gridSize.x*cell1Y);
				SInteger Range = 1;
				for (SInteger posj = -Range; posj <= Range; posj++) {
					for (SInteger posi = -Range; posi <= Range; posi++) {
						SInteger NCID = getNeighborCellID(cell1X + posi, cell1Y + posj, coarseInstance1);
						if (NCID < 0 || NCID >= gridStructure.totalCells)
						{
							continue;
						}
						size_t indx = NCID*maxFeaturesNum + secondFeatureIndex;
						size_t count = cellEventInstanceCount[indx];
						if (count > 0) {
							coarseIntermediate.push_back(coarseInstance1);
							coarseIntermediate.push_back(NCID);
							size_t bitIndex1 = coarseInstance1;
							size_t bitIndex2 = NCID + gridStructure.totalCells;
							Integer index1 = coarseInstance1*maxFeaturesNum + firstFeatureIndex;
							Integer index2 = NCID*maxFeaturesNum + secondFeatureIndex;
							Integer type1Count = cellEventInstanceCount[index1];
							Integer type2Count = cellEventInstanceCount[index2];
							coarseBitmap[bitIndex1] = type1Count;
							coarseBitmap[bitIndex2] = type2Count;
						}
					}
				}
			}
			Real upperbound = getUpperbound(coarseBitmap, degree, CurrentColoc);
			upperboundList.push_back(upperbound);
			totalCoarseInstances = coarseIntermediate.size() / degree;
			coarseBitmap.clear();
			coarseBitmap.shrink_to_fit();
			if (upperbound < parameter.PIthreshold) {
				filterpruneCounter++;
				coarseIntermediate.clear();
				coarseIntermediate.shrink_to_fit();
				auto fend = get_time::now();
				auto fdiff = fend - fstart;
				Real msecs = chrono::duration_cast<chrono::milliseconds>(fdiff).count();
				totalFilterTime += (msecs / 1000.0);
				continue;
			}
			auto fend = get_time::now();
			auto fdiff = fend - fstart;
			Real msecs = chrono::duration_cast<chrono::milliseconds>(fdiff).count();
			totalFilterTime += (msecs / 1000.0);
		}
	//MultiResolution Filter end


	//***********C++ Implementation: instance table generation ***************
		memoryTracker = 0.0;

		auto rStart = get_time::now();
		vector<Integer> intermediate;
		std::cout << "Table Gen Degree:" << degree << " Working on: " << i + 1 << "/" << candiColocCounter << std::endl;
		size_t start = featureInstanceStart[firstFeatureIndex];
		size_t end = featureInstanceEnd[firstFeatureIndex];
		size_t secondstart = featureInstanceStart[secondFeatureIndex];
		size_t secondend = featureInstanceEnd[secondFeatureIndex];
		size_t counter = 0;
		size_t bitmapSize = featureInstanceCount[firstFeatureIndex] + featureInstanceCount[secondFeatureIndex];
		vector<Integer> bitmap(bitmapSize,0);
		for (size_t idx1 = start; idx1 <= end; idx1++) {
			if (idx1 == 425) {
				std::cout<<std::endl;
			}
			Integer currentCell = instanceCellIDs[idx1];
			vector<Integer> nList;
			Integer CIDY = currentCell / gridStructure.gridSize.x;
			Integer CIDX = currentCell - (gridStructure.gridSize.x*CIDY);
			SInteger Range = 1;
			for (SInteger j = -1; j <= Range; ++j)
			{
				for (SInteger i = -1; i <= Range; ++i)
				{
					SInteger NCID = getNeighborCellID(CIDX + i, CIDY + j, currentCell);
					if (NCID < 0 || NCID >= gridStructure.totalCells)
					{
						continue;
					}
					size_t indx = NCID*maxFeaturesNum + secondFeatureIndex;
					size_t count = cellEventInstanceCount[indx];
					if (count > 0) {
						size_t startIndex = cellEventInstanceStart[indx];
						size_t endIndex = cellEventInstanceEnd[indx];
						for (size_t pos = startIndex; pos <= endIndex; pos++) {
							Integer instanceID = cellBasedSortedIndex[pos];
							Real dist = DistanceInMeters(idx1, instanceID);
							if (dist <= parameter.thresholdDistance) {
								nList.push_back(instanceID);
								size_t firstcount = end - start;
								size_t bitIndex1 = idx1-start;
								size_t bitIndex2 = instanceID - secondstart + firstcount;
								bitmap[bitIndex1] = 1;
								bitmap[bitIndex2] = 1;
							}
						}
					}
				}
			}
			std::sort(nList.begin(),nList.end());
			for (int indexList = 0; indexList < nList.size(); indexList++) {
				intermediate.push_back(idx1);
				intermediate.push_back(nList[indexList]);
			}
			nList.clear();
			nList.shrink_to_fit();
		}
		Real participationIndex = getParticipationIndex(bitmap, degree, CurrentColoc);
		PIList.push_back(participationIndex);
		if (participationIndex < parameter.PIthreshold) {
			if (parameter.FilterON_OFF == 1) {
				coarseIntermediate.clear();
				coarseIntermediate.shrink_to_fit();
			}
			auto rend = get_time::now();
			auto rdiff = rend - rStart;
			Real mrsecs = chrono::duration_cast<chrono::milliseconds>(rdiff).count();
			totalRefineTime += (mrsecs / 1000.0);
			continue;
		}
		intermediate.shrink_to_fit();
		totalInstances = intermediate.size() / degree;
		h_prevalentInstanceTables.push_back(intermediate);
		h_prevalentInstanceTableSize.push_back(intermediate.size() / degree);
		h_prevalantColocations.push_back(firstFeatureIndex);
		h_prevalantColocations.push_back(secondFeatureIndex);
		h_prevelentColocationCount++;
		if (parameter.FilterON_OFF == 1) {
			h_coarsePrevalentInstanceTableSize.push_back(coarseIntermediate.size()/degree);
			coarseIntermediate.shrink_to_fit();
			h_coarsePrevalentInstanceTables.push_back(coarseIntermediate);
			coarseIntermediate.clear();
			coarseIntermediate.shrink_to_fit();
		}
		memoryTracker += (Real)(intermediate.size()*4) / (1024 * 1024 * 1024);
		std::string msg = std::to_string(i + 1) + "/" + std::to_string(candiColocCounter) + "Instance Table Size = " + std::to_string(intermediate.size());
		//compLog(2, msg);
		intermediate.clear();
		intermediate.shrink_to_fit();
		auto rend = get_time::now();
		auto rdiff = rend - rStart;
		Real mrsecs = chrono::duration_cast<chrono::milliseconds>(rdiff).count();
		totalRefineTime += (mrsecs / 1000.0);
		std::cout << "Degree " << degree << ": Candidate pattern " << i + 1 << "/" << candiColocCounter << " end" << std::endl;
	}
	std::cout << "Degree " << degree << " Instance table Memory Requirement " << memoryTracker << " GB" << std::endl;
	delete[] cellEventInstanceEnd;
	delete[] cellEventInstanceStart;
	delete[] cellBasedSortedIndex;
	degree2FilterTime = totalFilterTime;
	degree2RefineTime = totalRefineTime;
	totalFilteredPatterns += filterpruneCounter;
}

void colocationFinder::clearSubsetVectors() {
	subset.clear();
	subset.shrink_to_fit();
	subsetList.clear();
	subsetList.shrink_to_fit();
}

void colocationFinder::tableGenRequired(Integer degree) {
	Integer lastIndex = degree - 2;
	bool flag;
	needInstanceTable = false;
	for (size_t i = 0; i < candiColocCounter; i++) {
		flag = true;
		for (size_t j = i + 1; j < candiColocCounter; j++) {
			for (size_t k = 0; k < lastIndex; k++) {
				if (candiColocations[i*(degree - 1) + k] != candiColocations[j*(degree - 1) + k]) {
					flag = false;
					break;
				}
			}
			if (flag && candiColocations[i*(degree - 1) + lastIndex]< candiColocations[j*(degree - 1) + lastIndex]) {
				vector<Integer> inter;
				for (Integer t = 0; t < degree - 1; t++) {
					inter.push_back(candiColocations[i*(degree - 1) + t]);
				}
				inter.push_back(candiColocations[j*(degree - 1) + lastIndex]);
				//module to generate k-1 degree subset
				flag = generateandCheckSubsets2(inter, degree);
				clearSubsetVectors();
				if (!flag) {
					break;
				}
				else {
					needInstanceTable = true;
					return;
				}
			}
			else {
				break;
			}
		}
	}
	totalCandidatePatterns += candiColocCounter;
}

void colocationFinder::candidateColocationGeneral(Integer degree) {
	Integer lastIndex = degree - 2;
	bool flag;
	candiColocations.clear();
	candiColocations.shrink_to_fit();
	candiColocCounter = 0;
	for (size_t i = 0; i < h_prevelentColocationCount; i++) {
		flag = true;
		for (size_t j = i + 1; j < h_prevelentColocationCount; j++) {
			for (size_t k= 0; k < lastIndex; k++) {
				if (h_prevalantColocations[i*(degree - 1) + k] != h_prevalantColocations[j*(degree - 1) + k]) {
					flag = false;
					break;
				}
			}
			if (flag && h_prevalantColocations[i*(degree - 1) + lastIndex]< h_prevalantColocations[j*(degree - 1) + lastIndex]) {
				vector<Integer> inter;
				for (Integer t = 0; t < degree-1; t++) {
					inter.push_back(h_prevalantColocations[i*(degree - 1) + t]);
				}
				inter.push_back(h_prevalantColocations[j*(degree - 1) + lastIndex]);
				//module to generate k-1 degree subset
				flag = generateandCheckSubsets(inter, degree);  
				clearSubsetVectors();
				if (!flag) {
					break;
				}
				else {
					for (Integer l = 0; l < degree; l++) {
						candiColocations.push_back(inter[l]);
					}
					candiColocCounter++;
				}
			}
			else {
				break;
			}
		}
	}
	totalCandidatePatterns += candiColocCounter;
}

void colocationFinder::subsetGen(vector<Integer> &inter, int k, int n, int idx) {
	if (idx == n)
		return;

	if (k == 1) {
		for (int i = idx; i<n; i++)
		{
			subset.push_back(inter[i]);
			subsetList.push_back(subset);
			subset.pop_back();
		}
	}

	for (int j = idx; j<n; j++) {
		subset.push_back(inter[j]);
		subsetGen(inter, k - 1, n, j + 1);
		subset.pop_back();
	}
}

bool colocationFinder::checkSubset(vector<Integer> subsetElem) {
	Integer degree = subsetElem.size();
	Integer flag = true;
	for (size_t i = 0; i < h_prevelentColocationCount; i++) {
		for (size_t j = 0; j < degree; j++) {
			if (h_prevalantColocations[i*degree + j] != subsetElem[j]) {
				flag = false;
				break;
			}
			flag = true;
		}
		if (flag) {
			return flag;
		}
	}
	return flag;
}

bool colocationFinder::checkSubset2(vector<Integer> subsetElem) {
	Integer degree = subsetElem.size();
	Integer flag = true;
	for (size_t i = 0; i < candiColocCounter; i++) {
		for (size_t j = 0; j < degree; j++) {
			if (candiColocations[i*degree + j] != subsetElem[j]) {
				flag = false;
				break;
			}
			flag = true;
		}
		if (flag) {
			return flag;
		}
	}
	return flag;
}

bool colocationFinder::generateandCheckSubsets(vector<Integer> &inter, Integer degree) {
	subsetGen(inter, degree - 1, degree, 0);
	bool flag = true;
	for (Integer i = 2; i < degree; i++) {
		flag = checkSubset(subsetList[i]);
		if (!flag) {
			return flag;
		}
	}
	return flag;
}

bool colocationFinder::generateandCheckSubsets2(vector<Integer> &inter, Integer degree) {
	subsetGen(inter, degree - 1, degree, 0);
	bool flag = true;
	for (size_t i = 2; i < degree; i++) {
		flag = checkSubset2(subsetList[i]);
		if (!flag) {
			return flag;
		}
	}
	return flag;
}


Integer colocationFinder::getIndex(vector<Integer> inner,Integer degree) {
	bool flag = false;
	for (Integer i = 0; i < h_prevelentColocationCount; i++) {
		for (Integer j = 0; j < degree; j++) {
			if (inner[j] != h_prevalantColocations[i*degree + j]) {
				flag = false;
				break;
			}
			flag = true;
		}
		if (flag) {
			return i;
		}
	}
	return -1;
}

void colocationFinder::generateInstanceTableGeneral(Integer degree) {
	Integer lastIndex = degree - 2;
	memoryTracker = 0.0;
	//MultiResolution Filter starts
	//timeLog(2, "multiResolution pruning ", "start", 0);
	vector<Integer> coarseIntermediate;
	Integer filterpruneCounter = 0;
	vector<Real> upperboundList;
	vector<Real> PIList;
	for (size_t i = 0; i < candiColocCounter; i++) {
		std::cout << "MRF Degree:" << degree << " Working on: " << i + 1 << "/" << candiColocCounter << std::endl;
		vector <Integer> coloc;
		vector<Integer> combiColoc1;
		vector<Integer> combiColoc2;
		for (Integer t = 0; t < degree; t++) {
			coloc.push_back(candiColocations[i*degree + t]);
			if (t < degree - 2) {
				combiColoc1.push_back(candiColocations[i*degree + t]);
				combiColoc2.push_back(candiColocations[i*degree + t]);
			}
		}
		combiColoc1.push_back(candiColocations[i*degree + (degree - 2)]);
		combiColoc2.push_back(candiColocations[i*degree + (degree - 1)]);

		size_t table1Index = getIndex(combiColoc1, degree - 1);
		size_t table2Index = getIndex(combiColoc2, degree - 1);
		if (parameter.FilterON_OFF == 1) {
			auto fstart = get_time::now();
			vector<Integer> coarseBitmap(degree*gridStructure.totalCells, 0);
			size_t instaceCountTable1 = h_coarsePrevalentInstanceTableSize[table1Index];
			size_t instaceCountTable2 = h_coarsePrevalentInstanceTableSize[table2Index];
			
			vector<Integer> CoarseTable1_ChunkStartArray;
			vector<Integer> CoarseTable1_ChunkEndArray;
			vector<Integer> CoarseTable2_ChunkStartArray;
			vector<Integer> CoarseTable2_ChunkEndArray;

			for (Integer CAIndex1 = 0; CAIndex1 < instaceCountTable1; CAIndex1++) {
				if (CAIndex1 == 0) {
					CoarseTable1_ChunkStartArray.push_back(0);
					continue;
				}
				for (Integer degIdx = 0; degIdx < degree - 2; degIdx++) {
					if (h_coarsePrevalentInstanceTables[table1Index][CAIndex1*(degree - 1) + degIdx] != h_coarsePrevalentInstanceTables[table1Index][(CAIndex1 - 1)*(degree - 1) + degIdx]) {
						CoarseTable1_ChunkEndArray.push_back(CAIndex1 - 1);
						CoarseTable1_ChunkStartArray.push_back(CAIndex1);
						break;
					}
				}
			}
			CoarseTable1_ChunkEndArray.push_back(instaceCountTable1 - 1);
			for (Integer CAIndex2 = 0; CAIndex2 < instaceCountTable2; CAIndex2++) {
				if (CAIndex2 == 0) {
					CoarseTable2_ChunkStartArray.push_back(0);
					continue;
				}
				for (Integer degIdx = 0; degIdx < degree - 2; degIdx++) {
					if (h_coarsePrevalentInstanceTables[table2Index][CAIndex2*(degree - 1) + degIdx] != h_coarsePrevalentInstanceTables[table2Index][(CAIndex2 - 1)*(degree - 1) + degIdx]) {
						CoarseTable2_ChunkEndArray.push_back(CAIndex2 - 1);
						CoarseTable2_ChunkStartArray.push_back(CAIndex2);
						break;
					}
				}
			}
			CoarseTable2_ChunkEndArray.push_back(instaceCountTable2 - 1);
			Integer p = 0, q = 0;
			Integer flag = 1;
			//vector<Integer> intermediate;
			while (p < CoarseTable1_ChunkStartArray.size() && q< CoarseTable2_ChunkStartArray.size()) {
				flag = 1;
				Integer pStartIndex = CoarseTable1_ChunkStartArray[p];
				Integer qStartIndex = CoarseTable2_ChunkStartArray[q];
				for (Integer degIdx = 0; degIdx < degree - 2; degIdx++) {
					if (h_coarsePrevalentInstanceTables[table1Index][pStartIndex*(degree - 1) + degIdx] < h_coarsePrevalentInstanceTables[table2Index][qStartIndex*(degree - 1) + degIdx]) {
						p++;
						flag = 0;
						break;
					}
					if (h_coarsePrevalentInstanceTables[table1Index][pStartIndex*(degree - 1) + degIdx] > h_coarsePrevalentInstanceTables[table2Index][qStartIndex*(degree - 1) + degIdx]) {
						q++;
						flag = 0;
						break;
					}
				}
				if (flag == 0) {
					continue;
				}
				for (size_t idx1 = CoarseTable1_ChunkStartArray[p]; idx1 <= CoarseTable1_ChunkEndArray[p]; idx1++) {
					for (size_t idx2 = CoarseTable2_ChunkStartArray[q]; idx2 <= CoarseTable2_ChunkEndArray[q]; idx2++) {
						bool ifNeighbour = checkifNeighbors(h_coarsePrevalentInstanceTables[table1Index][idx1*(degree - 1) + lastIndex], h_coarsePrevalentInstanceTables[table2Index][idx2*(degree - 1) + lastIndex]);
						if (ifNeighbour) {
							size_t offset = 0;
							size_t bitIndex = 0;
							Integer fType;
							Integer cellID;
							for (Integer idx4 = 0; idx4 < degree - 1; idx4++) {
								coarseIntermediate.push_back(h_coarsePrevalentInstanceTables[table1Index][idx1*(degree - 1) + idx4]);
								Integer value = h_coarsePrevalentInstanceTables[table1Index][idx1*(degree - 1) + idx4];
								fType = coloc[idx4];
								size_t fstart = featureInstanceStart[fType];
								size_t fcount = featureInstanceCount[fType];
								Integer bitIndex = value + offset;
								offset += gridStructure.totalCells;
								Integer index1 = value*maxFeaturesNum + fType;
								Integer typeCount = cellEventInstanceCount[index1];
								coarseBitmap[bitIndex] = typeCount;
							}
							coarseIntermediate.push_back(h_coarsePrevalentInstanceTables[table2Index][idx2*(degree - 1) + lastIndex]);
							cellID = h_coarsePrevalentInstanceTables[table2Index][idx2*(degree - 1) + lastIndex];
							fType = coloc[degree-1];
							bitIndex = cellID + offset;
							Integer index2 = cellID*maxFeaturesNum + fType;
							Integer typeCount2 = cellEventInstanceCount[index2];
							coarseBitmap[bitIndex] = typeCount2;
						}
					}
				}
				p++;
			}
			Real upperbound = getUpperbound(coarseBitmap, degree, coloc);
			upperboundList.push_back(upperbound);
			if (upperbound < parameter.PIthreshold) {
				coarseIntermediate.clear();
				coarseIntermediate.shrink_to_fit();
				auto fend = get_time::now();
				auto fdiff = fend - fstart;
				Real msecs = chrono::duration_cast<chrono::milliseconds>(fdiff).count();
				totalFilterTime += (msecs / 1000.0);
				filterpruneCounter++;
				continue;
			}
			memoryTracker += (Real)(coarseIntermediate.size()*4) / (1024 * 1024 * 1024);
			std::string msg1 = std::to_string(i + 1) + "/" + std::to_string(candiColocCounter) + "coarse Table Size = " + std::to_string(coarseIntermediate.size());
			//compLog(degree, msg1);

			CoarseTable1_ChunkStartArray.clear();
			CoarseTable1_ChunkStartArray.shrink_to_fit();
			CoarseTable1_ChunkEndArray.clear();
			CoarseTable1_ChunkEndArray.shrink_to_fit();

			CoarseTable2_ChunkStartArray.clear();
			CoarseTable2_ChunkStartArray.shrink_to_fit();
			CoarseTable2_ChunkEndArray.clear();
			CoarseTable2_ChunkEndArray.shrink_to_fit();
			auto fend = get_time::now();
			auto fdiff = fend - fstart;
			//std::cout << "Degree " << degree << " MultiResolution Pruning Memory Requirement " << memoryTracker << " GB" << std::endl;
			Real msecs = chrono::duration_cast<chrono::milliseconds>(fdiff).count();
			totalFilterTime += (msecs / 1000.0);
		}
		//MultiResoultion Filter ends
		memoryTracker = 0.0;
		auto rStart = get_time::now();
		auto start = get_time::now();
		//timelog(degree, "generateInstanceTableGeneral ", "start",0);
		std::cout << "Table Gen Degree:" << degree << " Working on: " << i + 1 << "/" << candiColocCounter << std::endl;
		Integer bitmapSize = 0;
		for (size_t Idx = 0; Idx < degree; Idx++) {
			Integer fType = coloc[Idx];
			bitmapSize += featureInstanceCount[fType];
		}
		vector<Integer> bitmap(bitmapSize, 0);
		size_t instaceCountTable1 = h_prevalentInstanceTableSize[table1Index]/*/(degree-1)*/;
		size_t instaceCountTable2 = h_prevalentInstanceTableSize[table2Index] /*/ (degree - 1)*/;
		
		size_t sizeofFirstFeatuere = featureInstanceCount[coloc[0]];
		size_t firstEventStartsFrom = featureInstanceStart[coloc[0]];

		vector<Integer> Table1_ChunkStartArray;
		vector<Integer> Table1_ChunkEndArray;
		vector<Integer> Table2_ChunkStartArray;
		vector<Integer> Table2_ChunkEndArray;

		for (Integer CAIndex1 = 0; CAIndex1 < instaceCountTable1; CAIndex1++) {
			if (CAIndex1 == 0) {
				Table1_ChunkStartArray.push_back(0);
				continue;
			}
			for (Integer degIdx = 0; degIdx < degree - 2; degIdx++) {
				if (h_prevalentInstanceTables[table1Index][CAIndex1*(degree - 1) + degIdx] != h_prevalentInstanceTables[table1Index][(CAIndex1 - 1)*(degree - 1) + degIdx]) {
					Table1_ChunkEndArray.push_back(CAIndex1 - 1);
					Table1_ChunkStartArray.push_back(CAIndex1);
					break;
				}
			}
		
		}
		Table1_ChunkEndArray.push_back(instaceCountTable1-1);
		for (Integer CAIndex2 = 0; CAIndex2 < instaceCountTable2; CAIndex2++) {
			if (CAIndex2 == 0) {
				Table2_ChunkStartArray.push_back(0);
				continue;
			}
			for (Integer degIdx = 0; degIdx < degree - 2; degIdx++) {
				if (h_prevalentInstanceTables[table2Index][CAIndex2*(degree - 1) + degIdx] != h_prevalentInstanceTables[table2Index][(CAIndex2 - 1)*(degree - 1) + degIdx]) {
					Table2_ChunkEndArray.push_back(CAIndex2 - 1);
					Table2_ChunkStartArray.push_back(CAIndex2);
					break;
				}
			}
		}
		Table2_ChunkEndArray.push_back(instaceCountTable2 - 1);
		Integer p = 0, q = 0;
		Integer flag = 1;
		vector<Integer> intermediate;
		while (p < Table1_ChunkStartArray.size() &&	q< Table2_ChunkStartArray.size()) {
			flag = 1;
			Integer pStartIndex = Table1_ChunkStartArray[p];
			Integer qStartIndex = Table2_ChunkStartArray[q];
			for (Integer degIdx = 0; degIdx < degree - 2; degIdx++) {
				if (h_prevalentInstanceTables[table1Index][pStartIndex*(degree - 1) + degIdx] < h_prevalentInstanceTables[table2Index][qStartIndex*(degree - 1) + degIdx]) {
					p++;
					flag = 0;
					break;
				}
				if (h_prevalentInstanceTables[table1Index][pStartIndex*(degree - 1) + degIdx] > h_prevalentInstanceTables[table2Index][qStartIndex*(degree - 1) + degIdx]) {
					q++;
					flag = 0;
					break;
				}
			}
			if (flag == 0) {
				continue;
			}
			for (size_t idx1 = Table1_ChunkStartArray[p]; idx1 <= Table1_ChunkEndArray[p]; idx1++) {
				for (size_t idx2 = Table2_ChunkStartArray[q]; idx2 <= Table2_ChunkEndArray[q]; idx2++) {
					Real dist = DistanceInMeters(h_prevalentInstanceTables[table1Index][idx1*(degree - 1) + lastIndex], h_prevalentInstanceTables[table2Index][idx2*(degree - 1) + lastIndex]);
					if (dist <= parameter.thresholdDistance) {
						size_t offset = 0;
						size_t bitIndex = 0;
						Integer fType;
						Integer value = 0;
						for (Integer idx4 = 0; idx4 < degree - 1; idx4++) {
							value = h_prevalentInstanceTables[table1Index][idx1*(degree - 1) + idx4];
							intermediate.push_back(value);
							fType = coloc[idx4];
							size_t fstart = featureInstanceStart[fType];
							size_t fcount = featureInstanceCount[fType];
							bitIndex = value - fstart + offset;
							offset += fcount;
							bitmap[bitIndex] = 1;
						}
						value = h_prevalentInstanceTables[table2Index][idx2*(degree - 1) + lastIndex];
						fType = coloc[degree-1];
						size_t lastStart = featureInstanceStart[fType];
						bitIndex = value - lastStart + offset;
						bitmap[bitIndex] = 1;
						intermediate.push_back(value);
					}
				}
			}
			p++;
		}

		Real participationIndex = getParticipationIndex(bitmap, degree, coloc);
		PIList.push_back(participationIndex);
		if (participationIndex < parameter.PIthreshold) {
			if (parameter.FilterON_OFF == 1)
			{
				coarseIntermediate.clear();
				coarseIntermediate.shrink_to_fit();
			}
			auto rend = get_time::now();
			auto rdiff = rend - rStart;
			Real mrsecs = chrono::duration_cast<chrono::milliseconds>(rdiff).count();
			totalRefineTime += (mrsecs / 1000.0);
			continue;
		}
		intermediate.shrink_to_fit();
		if (parameter.FilterON_OFF == 1) {
			h_coarsePrevalentInstanceTableSize2.push_back(coarseIntermediate.size()/degree);
			coarseIntermediate.shrink_to_fit();
			h_coarsePrevalentInstanceTables2.push_back(coarseIntermediate);
			coarseIntermediate.clear();
			coarseIntermediate.shrink_to_fit();
		}
		for (Integer Idx = 0; Idx < degree; Idx++) {
			h_prevalantColocations2.push_back(coloc[Idx]);
		}
		h_prevelentColocationCount2++;
		h_prevalentInstanceTableSize2.push_back(intermediate.size()/degree);
		intermediate.shrink_to_fit();
		h_prevalentInstanceTables2.push_back(intermediate);
		memoryTracker += (Real)(intermediate.size() * 4) / (1024 * 1024 * 1024);
		std::string msg = std::to_string(i + 1) + "/" + std::to_string(candiColocCounter) + "Instance Table Size = " + std::to_string(intermediate.size());
		//compLog(degree, msg);
		intermediate.clear();
		intermediate.shrink_to_fit();
		coloc.clear();
		coloc.shrink_to_fit();
		combiColoc1.clear();
		combiColoc1.shrink_to_fit();
		combiColoc2.clear();
		combiColoc2.shrink_to_fit();
		Table1_ChunkStartArray.clear();
		Table1_ChunkStartArray.shrink_to_fit();
		Table1_ChunkEndArray.clear();
		Table1_ChunkEndArray.shrink_to_fit();
		
		Table2_ChunkStartArray.clear();
		Table2_ChunkStartArray.shrink_to_fit();
		Table2_ChunkEndArray.clear();
		Table2_ChunkEndArray.shrink_to_fit();
		auto rend = get_time::now();
		auto rdiff = rend - rStart;
		Real mrsecs = chrono::duration_cast<chrono::milliseconds>(rdiff).count();
		totalRefineTime += (mrsecs / 1000.0);
	}
	//std::cout << "Degree " << degree << " Instance table Memory Requirement " << memoryTracker << " GB" << std::endl;
	totalFilteredPatterns += filterpruneCounter;
}
void colocationFinder::generateInstanceTableGeneral2(Integer degree) {
	Integer lastIndex = degree - 2;
	memoryTracker = 0.0;
	//MultiResolution Filter starts
	//if (parameter.FilterON_OFF == 1) {
	//auto start = get_time::now();
	//timeLog(2, "multiResolution pruning ", "start", 0);
	vector<Integer> coarseIntermediate;
	Integer filterpruneCounter = 0;
	for (size_t i = 0; i < candiColocCounter; i++) {
		std::cout << "MRF Degree:" << degree << " Working on: " << i + 1 << "/" << candiColocCounter << std::endl;
		vector <Integer> coloc;
		vector<Integer> combiColoc1;
		vector<Integer> combiColoc2;
		for (Integer t = 0; t < degree; t++) {
			coloc.push_back(candiColocations[i*degree + t]);
			if (t < degree - 2) {
				combiColoc1.push_back(candiColocations[i*degree + t]);
				combiColoc2.push_back(candiColocations[i*degree + t]);
			}
		}
		combiColoc1.push_back(candiColocations[i*degree + (degree - 2)]);
		combiColoc2.push_back(candiColocations[i*degree + (degree - 1)]);

		size_t table1Index = getIndex(combiColoc1, degree - 1);
		size_t table2Index = getIndex(combiColoc2, degree - 1);
		if (parameter.FilterON_OFF == 1) {
			auto fstart = get_time::now();
			vector<Integer> coarseBitmap(degree*gridStructure.totalCells, 0);
			//size_t instaceCountTable1 = h_prevalentCoarseInstanceTableSize[table1Index] / (degree - 1);
			//size_t instaceCountTable2 = h_prevalentCoarseInstanceTableSize[table2Index] / (degree - 1);
			size_t instaceCountTable1 = h_coarsePrevalentInstanceTableSize[table1Index];
			size_t instaceCountTable2 = h_coarsePrevalentInstanceTableSize[table2Index];

			vector<Integer> CoarseTable1_ChunkStartArray;
			vector<Integer> CoarseTable1_ChunkEndArray;
			vector<Integer> CoarseTable2_ChunkStartArray;
			vector<Integer> CoarseTable2_ChunkEndArray;

			for (Integer CAIndex1 = 0; CAIndex1 < instaceCountTable1; CAIndex1++) {
				if (CAIndex1 == 0) {
					CoarseTable1_ChunkStartArray.push_back(0);
					continue;
				}
				for (Integer degIdx = 0; degIdx < degree - 2; degIdx++) {
					if (h_coarsePrevalentInstanceTables[table1Index][CAIndex1*(degree - 1) + degIdx] != h_coarsePrevalentInstanceTables[table1Index][(CAIndex1 - 1)*(degree - 1) + degIdx]) {
						CoarseTable1_ChunkEndArray.push_back(CAIndex1 - 1);
						CoarseTable1_ChunkStartArray.push_back(CAIndex1);
						break;
					}
				}
			}
			CoarseTable1_ChunkEndArray.push_back(instaceCountTable1 - 1);
			for (Integer CAIndex2 = 0; CAIndex2 < instaceCountTable2; CAIndex2++) {
				if (CAIndex2 == 0) {
					CoarseTable2_ChunkStartArray.push_back(0);
					continue;
				}
				for (Integer degIdx = 0; degIdx < degree - 2; degIdx++) {
					if (h_coarsePrevalentInstanceTables[table2Index][CAIndex2*(degree - 1) + degIdx] != h_coarsePrevalentInstanceTables[table2Index][(CAIndex2 - 1)*(degree - 1) + degIdx]) {
						CoarseTable2_ChunkEndArray.push_back(CAIndex2 - 1);
						CoarseTable2_ChunkStartArray.push_back(CAIndex2);
						break;
					}
				}
			}
			CoarseTable2_ChunkEndArray.push_back(instaceCountTable2 - 1);
			Integer p = 0, q = 0;
			Integer flag = 1;
			//vector<Integer> intermediate;
			while (p < CoarseTable1_ChunkStartArray.size() && q< CoarseTable2_ChunkStartArray.size()) {
				flag = 1;
				Integer pStartIndex = CoarseTable1_ChunkStartArray[p];
				Integer qStartIndex = CoarseTable2_ChunkStartArray[q];
				for (Integer degIdx = 0; degIdx < degree - 2; degIdx++) {
					if (h_coarsePrevalentInstanceTables[table1Index][pStartIndex*(degree - 1) + degIdx] < h_coarsePrevalentInstanceTables[table2Index][qStartIndex*(degree - 1) + degIdx]) {
						p++;
						flag = 0;
						break;
					}
					if (h_coarsePrevalentInstanceTables[table1Index][pStartIndex*(degree - 1) + degIdx] > h_coarsePrevalentInstanceTables[table2Index][qStartIndex*(degree - 1) + degIdx]) {
						q++;
						flag = 0;
						break;
					}
				}
				if (flag == 0) {
					continue;
				}
				for (size_t idx1 = CoarseTable1_ChunkStartArray[p]; idx1 <= CoarseTable1_ChunkEndArray[p]; idx1++) {
					for (size_t idx2 = CoarseTable2_ChunkStartArray[q]; idx2 <= CoarseTable2_ChunkEndArray[q]; idx2++) {
						bool ifNeighbour = checkifNeighbors(h_coarsePrevalentInstanceTables[table1Index][idx1*(degree - 1) + lastIndex], h_coarsePrevalentInstanceTables[table2Index][idx2*(degree - 1) + lastIndex]);
						if (ifNeighbour) {
							size_t offset = 0;
							size_t bitIndex = 0;
							Integer fType;
							Integer cellID;
							for (Integer idx4 = 0; idx4 < degree - 1; idx4++) {
								Integer value = h_coarsePrevalentInstanceTables[table1Index][idx1*(degree - 1) + idx4];
								fType = coloc[idx4];
								size_t fstart = featureInstanceStart[fType];
								size_t fcount = featureInstanceCount[fType];
								Integer bitIndex = value + offset;
								offset += gridStructure.totalCells;
								Integer index1 = value*maxFeaturesNum + fType;
								Integer typeCount = cellEventInstanceCount[index1];
								coarseBitmap[bitIndex] = typeCount;
							}
							cellID = h_coarsePrevalentInstanceTables[table2Index][idx2*(degree - 1) + lastIndex];
							fType = coloc[degree - 1];
							bitIndex = cellID + offset;
							Integer index2 = cellID*maxFeaturesNum + fType;
							Integer typeCount2 = cellEventInstanceCount[index2];
							coarseBitmap[bitIndex] = typeCount2;
						}
					}
				}
				p++;
			}
			Real upperbound = getUpperbound(coarseBitmap, degree, coloc);
			if (upperbound < parameter.PIthreshold) {
				auto fend = get_time::now();
				auto fdiff = fend - fstart;
				Real msecs = chrono::duration_cast<chrono::milliseconds>(fdiff).count();
				totalFilterTime += (msecs / 1000.0);
				filterpruneCounter++;
				continue;
			}

			CoarseTable1_ChunkStartArray.clear();
			CoarseTable1_ChunkStartArray.shrink_to_fit();
			CoarseTable1_ChunkEndArray.clear();
			CoarseTable1_ChunkEndArray.shrink_to_fit();

			CoarseTable2_ChunkStartArray.clear();
			CoarseTable2_ChunkStartArray.shrink_to_fit();
			CoarseTable2_ChunkEndArray.clear();
			CoarseTable2_ChunkEndArray.shrink_to_fit();
			auto fend = get_time::now();
			auto fdiff = fend - fstart;
		//	std::cout << "Degree " << degree << " MultiResolution Pruning Memory Requirement " << memoryTracker << " GB" << std::endl;
			Real msecs = chrono::duration_cast<chrono::milliseconds>(fdiff).count();
			totalFilterTime += (msecs / 1000.0);
		}
		//MultiResoultion Filter ends
		memoryTracker = 0.0;
		auto rStart = get_time::now();
		auto start = get_time::now();
		//timelog(degree, "generateInstanceTableGeneral ", "start", 0);
		std::cout << "Table Gen Degree:" << degree << " Working on: " << i + 1 << "/" << candiColocCounter << std::endl;
		Integer bitmapSize = 0;
		for (size_t Idx = 0; Idx < degree; Idx++) {
			Integer fType = coloc[Idx];
			bitmapSize += featureInstanceCount[fType];
		}
		vector<Integer> bitmap(bitmapSize, 0);
		size_t instaceCountTable1 = h_prevalentInstanceTableSize[table1Index]/*/(degree-1)*/;
		size_t instaceCountTable2 = h_prevalentInstanceTableSize[table2Index] /*/ (degree - 1)*/;

		size_t sizeofFirstFeatuere = featureInstanceCount[coloc[0]];
		size_t firstEventStartsFrom = featureInstanceStart[coloc[0]];

		vector<Integer> Table1_ChunkStartArray;
		vector<Integer> Table1_ChunkEndArray;
		vector<Integer> Table2_ChunkStartArray;
		vector<Integer> Table2_ChunkEndArray;

		for (Integer CAIndex1 = 0; CAIndex1 < instaceCountTable1; CAIndex1++) {
			if (CAIndex1 == 0) {
				Table1_ChunkStartArray.push_back(0);
				continue;
			}
			for (Integer degIdx = 0; degIdx < degree - 2; degIdx++) {
				if (h_prevalentInstanceTables[table1Index][CAIndex1*(degree - 1) + degIdx] != h_prevalentInstanceTables[table1Index][(CAIndex1 - 1)*(degree - 1) + degIdx]) {
					Table1_ChunkEndArray.push_back(CAIndex1 - 1);
					Table1_ChunkStartArray.push_back(CAIndex1);
					break;
				}
			}

		}
		Table1_ChunkEndArray.push_back(instaceCountTable1 - 1);
		for (Integer CAIndex2 = 0; CAIndex2 < instaceCountTable2; CAIndex2++) {
			if (CAIndex2 == 0) {
				Table2_ChunkStartArray.push_back(0);
				continue;
			}
			for (Integer degIdx = 0; degIdx < degree - 2; degIdx++) {
				if (h_prevalentInstanceTables[table2Index][CAIndex2*(degree - 1) + degIdx] != h_prevalentInstanceTables[table2Index][(CAIndex2 - 1)*(degree - 1) + degIdx]) {
					Table2_ChunkEndArray.push_back(CAIndex2 - 1);
					Table2_ChunkStartArray.push_back(CAIndex2);
					break;
				}
			}
		}
		Table2_ChunkEndArray.push_back(instaceCountTable2 - 1);
		Integer p = 0, q = 0;
		Integer flag = 1;
		vector<Integer> intermediate;
		while (p < Table1_ChunkStartArray.size() && q< Table2_ChunkStartArray.size()) {
			flag = 1;
			Integer pStartIndex = Table1_ChunkStartArray[p];
			Integer qStartIndex = Table2_ChunkStartArray[q];
			for (Integer degIdx = 0; degIdx < degree - 2; degIdx++) {
				if (h_prevalentInstanceTables[table1Index][pStartIndex*(degree - 1) + degIdx] < h_prevalentInstanceTables[table2Index][qStartIndex*(degree - 1) + degIdx]) {
					p++;
					flag = 0;
					break;
				}
				if (h_prevalentInstanceTables[table1Index][pStartIndex*(degree - 1) + degIdx] > h_prevalentInstanceTables[table2Index][qStartIndex*(degree - 1) + degIdx]) {
					q++;
					flag = 0;
					break;
				}
			}
			if (flag == 0) {
				continue;
			}
			for (size_t idx1 = Table1_ChunkStartArray[p]; idx1 <= Table1_ChunkEndArray[p]; idx1++) {
				for (size_t idx2 = Table2_ChunkStartArray[q]; idx2 <= Table2_ChunkEndArray[q]; idx2++) {
					Real dist = DistanceInMeters(h_prevalentInstanceTables[table1Index][idx1*(degree - 1) + lastIndex], h_prevalentInstanceTables[table2Index][idx2*(degree - 1) + lastIndex]);
					if (dist <= parameter.thresholdDistance) {
						size_t offset = 0;
						size_t bitIndex = 0;
						Integer fType;
						Integer value = 0;
						for (Integer idx4 = 0; idx4 < degree - 1; idx4++) {
							value = h_prevalentInstanceTables[table1Index][idx1*(degree - 1) + idx4];
							fType = coloc[idx4];
							size_t fstart = featureInstanceStart[fType];
							size_t fcount = featureInstanceCount[fType];
							bitIndex = value - fstart + offset;
							offset += fcount;
							bitmap[bitIndex] = 1;
						}
						value = h_prevalentInstanceTables[table2Index][idx2*(degree - 1) + lastIndex];
						fType = coloc[degree - 1];
						size_t lastStart = featureInstanceStart[fType];
						bitIndex = value - lastStart + offset;
						bitmap[bitIndex] = 1;
					}
				}
			}
			p++;
		}

		Real participationIndex = getParticipationIndex(bitmap, degree, coloc);
		if (participationIndex < parameter.PIthreshold) {
			if (parameter.FilterON_OFF == 1)
			{
				coarseIntermediate.clear();
				coarseIntermediate.shrink_to_fit();
			}
			auto rend = get_time::now();
			auto rdiff = rend - rStart;
			Real mrsecs = chrono::duration_cast<chrono::milliseconds>(rdiff).count();
			totalRefineTime += (mrsecs / 1000.0);
			continue;
		}
		intermediate.shrink_to_fit();
		if (parameter.FilterON_OFF == 1) {
			h_coarsePrevalentInstanceTableSize2.push_back(coarseIntermediate.size() / degree);
			coarseIntermediate.shrink_to_fit();
			h_coarsePrevalentInstanceTables2.push_back(coarseIntermediate);
			coarseIntermediate.clear();
			coarseIntermediate.shrink_to_fit();
		}
		for (Integer Idx = 0; Idx < degree; Idx++) {
			h_prevalantColocations2.push_back(coloc[Idx]);
		}
		h_prevelentColocationCount2++;
		coloc.clear();
		coloc.shrink_to_fit();
		combiColoc1.clear();
		combiColoc1.shrink_to_fit();
		combiColoc2.clear();
		combiColoc2.shrink_to_fit();
		Table1_ChunkStartArray.clear();
		Table1_ChunkStartArray.shrink_to_fit();
		Table1_ChunkEndArray.clear();
		Table1_ChunkEndArray.shrink_to_fit();

		Table2_ChunkStartArray.clear();
		Table2_ChunkStartArray.shrink_to_fit();
		Table2_ChunkEndArray.clear();
		Table2_ChunkEndArray.shrink_to_fit();
		auto rend = get_time::now();
		auto rdiff = rend - rStart;
		Real mrsecs = chrono::duration_cast<chrono::milliseconds>(rdiff).count();
		totalRefineTime += (mrsecs / 1000.0);
	}
	//std::cout << "Degree " << degree << " Instance table Memory Requirement " << memoryTracker << " GB" << std::endl;
	totalFilteredPatterns += filterpruneCounter;
}

void colocationFinder::copyPrevalentColocations(Integer degree) {
	h_prevalantColocations.clear();
	h_prevalantColocations.shrink_to_fit();
	h_prevalantColocations.swap(h_prevalantColocations2);
	h_prevalantColocations2.clear();
	h_prevalantColocations2.shrink_to_fit();
	h_prevalentInstanceTables.clear();
	h_prevalentInstanceTables.shrink_to_fit();
	h_prevalentInstanceTables.swap(h_prevalentInstanceTables2);
	h_prevalentInstanceTables2.clear();
	h_prevalentInstanceTables2.shrink_to_fit();
	h_prevelentColocationCount = h_prevelentColocationCount2;
	h_prevelentColocationCount2 = 0;
	h_prevalentInstanceTableSize.clear();
	h_prevalentInstanceTableSize.shrink_to_fit();
	h_prevalentInstanceTableSize.swap(h_prevalentInstanceTableSize2);
	h_prevalentInstanceTableSize2.clear();
	h_prevalentInstanceTableSize2.shrink_to_fit();
	if (parameter.FilterON_OFF == 1) {
		h_coarsePrevalentInstanceTables.clear();
		h_coarsePrevalentInstanceTables.shrink_to_fit();
		h_coarsePrevalentInstanceTables.swap(h_coarsePrevalentInstanceTables2);
		h_coarsePrevalentInstanceTables2.clear();
		h_coarsePrevalentInstanceTables2.shrink_to_fit();
		h_coarsePrevalentInstanceTableSize.clear();
		h_coarsePrevalentInstanceTableSize.shrink_to_fit();
		h_coarsePrevalentInstanceTableSize.swap(h_coarsePrevalentInstanceTableSize2);
		h_coarsePrevalentInstanceTableSize2.clear();
		h_coarsePrevalentInstanceTableSize2.shrink_to_fit();
	}
	
}
void colocationFinder::resetCandidateData() {
	//free memory
	candiColocations.clear();
	candiColocations.shrink_to_fit();
}

void colocationFinder::savetoFile(Integer degree) {
	std::string degStr = std::to_string(degree);
	std::string FileName = location + "ColocationDegree_" + degStr + ".txt";
	std::ofstream fout(FileName.c_str());
	for (Integer i = 0; i < h_prevelentColocationCount; i++) {
		fout << "( ";
		for (Integer j = 0; j < degree; j++) {
			Integer index = h_prevalantColocations[i*degree + j];
			fout << featureTypes[index] << " ";
			if (j != degree - 1) {
				fout << "| ";
			}
			else {
				fout << ")" << "\n";
			}
		}
	}
	fout.close();
}
void colocationFinder::savetoFileGen(Integer degree) {
	std::string degStr = std::to_string(degree);
	std::string FileName = location + "ColocationDegree_" + degStr + ".txt";
	std::ofstream fout(FileName.c_str());
	for (size_t i = 0; i < h_prevelentColocationCount2; i++) {
		fout << "( ";
		for (size_t j = 0; j < degree; j++) {
			size_t index = h_prevalantColocations2[i*degree + j];
			fout << featureTypes[index] << " ";
			if (j != degree - 1) {
				fout << "| ";
			}
			else {
				fout << ")" << "\n";
			}
		}
	}
	fout.close();
}

void colocationFinder::timeLog(Integer degree, std::string function, std::string eventType, Integer duration) {
	time_t currentTime;
	struct tm *localTime;
	time(&currentTime);                   // Get the current time
	localTime = localtime(&currentTime);  // Convert the current time to the local time
	timeLogStream.open(logFileName.c_str(), ofstream::app);
	timeLogStream << degree << " " << function << " " << eventType << " " << localTime->tm_year + 1900 << "\\" << localTime->tm_mon + 1 << "\\" << localTime->tm_mday << " " << localTime->tm_hour << ":" << localTime->tm_min << ":" << localTime->tm_sec << "\n";
	if (eventType == "end") {
		timeLogStream << degree << " " << function << " " << "Duration: " << " " << duration << "seconds \n";
	}
	timeLogStream.close();
}
void colocationFinder::compLog(Integer degree, std::string function) {
	compLogStream.open(compLogFileName.c_str(), ofstream::app);
	compLogStream << "Degree: " << degree << ", " << function << "\n";
	compLogStream.close();
}
void colocationFinder::compLog(Integer i, Integer degree, Integer candicolocNumber, std::string Remark, std::string PIType, Real PI) {
	compLogStream.open(compLogFileName.c_str(), ofstream::app);
	if (Remark != "") {
		compLogStream << "Degree: " << degree << ", Candidate " << i << "/" << candicolocNumber << PIType << PI << " (" << Remark << ") \n";
	}
	else {
		compLogStream << "Degree: " << degree << ", Candidate " << i << "/" << candicolocNumber << PIType << PI << "\n";
	}
	compLogStream.close();
}
void colocationFinder::compLog(Integer degree, Integer i, std::string totalmem, Integer candiColocNum) {
	compLogStream.open(compLogFileName.c_str(), ofstream::app);
	compLogStream << "Degree: " << degree << " pattern No.: " << i << "/" << candiColocNum << " Total memRequired =" << totalmem << "GB\n";
	compLogStream.close();
}
void colocationFinder::loggingFileOpenFunction() {
	logFileName = location + "CPU_TimeLog.txt";
	compLogFileName = location + "CPU_ComputationLog.txt";
}
void colocationFinder::loggingFileCloseFunction() {
	timeLogStream.open(logFileName.c_str(), ofstream::app);
	timeLogStream << "distance Threshold: " << parameter.thresholdDistance << "\n";
	timeLogStream << "prevalance Threshold: " << parameter.PIthreshold << "\n";
	timeLogStream.close();
	//compLogStream.close();
}

void colocationFinder::clearMemory() {
	h_prevalentInstanceTables.clear();
	h_prevalentInstanceTables.shrink_to_fit();

	h_prevalentInstanceTables2.clear();
	h_prevalentInstanceTables2.shrink_to_fit();

	h_prevalantColocations.clear();
	h_prevalantColocations.shrink_to_fit();

	h_prevalantColocations2.clear();
	h_prevalantColocations2.shrink_to_fit();

	h_prevalentInstanceTableSize.clear();
	h_prevalentInstanceTableSize.shrink_to_fit();

	h_prevalentInstanceTableSize2.clear();
	h_prevalentInstanceTableSize2.shrink_to_fit();
}

void colocationFinder::Begin(int argc, char**argv) {
	auto begin = get_time::now();
	std::string datasetFilename;
	if (argc > 1) {
		std::string configurationFile = argv[1];
		ifstream in(configurationFile.c_str());
		std::string line;
		std::string parameterFilename;
		getline(in, line);
		location = line;
		getline(in, line);
		parameterFilename = line;
		getline(in, line);
		datasetFilename = line;
		getline(in, line);
		outputLocation = line;
		getline(in, line);
		outputFile = line;
		in.close();
		loggingFileOpenFunction();
		LoadParameter(parameterFilename);
		auto start = get_time::now();
		//timeLog(1, "populate Data", "start",0);
		populateData(datasetFilename);
		auto end = get_time::now();
		auto diff = end - start;
		//timeLog(1, "populate Data", "end", chrono::duration_cast<sec>(diff).count());
		}
	else {
		std::cout << "Configuration file missing...";
		exit(0);
	}
	degree2Processing();
	//std::cout << "Saving degree 2 colcation data..."<<std::endl;
	savetoFile(2);
	//std::cout << "\n\nDegree 2 processing ends...\n" << std::endl;
	Integer degree = 3;
	Integer featureCount = featureTypes.size();
	std::cout << "\n\nDegree " << degree << " processing starts..." << std::endl;
	auto start = get_time::now();
	//timelog(degree, "candidateColocationGeneral ", "start", 0);
	candidateColocationGeneral(degree);
	auto end = get_time::now();
	auto diff = end - start;
	//timelog(degree, "candidateColocationGeneral " + std::to_string(candiColocCounter), "end", chrono::duration_cast<sec>(diff).count());
	//std::cout << "candidateColocationGeneral() took  " << chrono::duration_cast<sec>(diff).count() << " seconds " << endl;
	tableGenRequired(degree+1);
	while (candiColocCounter > 0 && degree <= featureCount) {
		std::cout << "\n\nDegree " << degree << " processing ...\n" << std::endl;
		if (needInstanceTable) {
			generateInstanceTableGeneral(degree);
		}
		else {
			generateInstanceTableGeneral2(degree);
		}
			
		savetoFileGen(degree);
		//std::cout << "\n\nDegree " << degree << " processing ends...\n" << std::endl;
		//std::cout << std::endl;

		copyPrevalentColocations(degree);

		//candidate colocation pattern generator 
		//std::cout << "\n\nDegree " << degree << " processing starts..." << std::endl;
		start = get_time::now();
		//timelog(degree, "candidateColocationGeneral ", "start", 0);
		candidateColocationGeneral(degree+1);
		tableGenRequired(degree+2);
		end = get_time::now();
		diff = end - start;
		//timeLog(degree , "candidateColocationGeneral " + std::to_string(candiColocCounter), "end", chrono::duration_cast<sec>(diff).count());
	//	std::cout << "candidateColocationGeneral() took  " << chrono::duration_cast<sec>(diff).count() << " seconds " << endl;
		if (candiColocCounter == 0) {
			clearMemory();//clean memory
			break;
		}
		degree++;
	}
	auto last = get_time::now();
	auto difference = last - begin;
	//timelog(degree, "End of Program.. ", "end", chrono::duration_cast<sec>(difference).count());
	std::cout << "End of program :  " << chrono::duration_cast<sec>(difference).count() << " seconds " << endl;
	loggingFileCloseFunction();
	std::ofstream expResult;
	std::string expFile = outputLocation + outputFile;
	expResult.open(expFile.c_str(), ofstream::app);
	std::string exe = argv[0];
	int pos = 0;
	for (int i = 0; i < exe.length(); i++) {
		if (exe[i] == '\\') {
			pos = i;
		}
	}
	exe.erase(0, pos + 1);
	int pos2 = 0;
	for (int i = 0; i < datasetFilename.length(); i++) {
		if (datasetFilename[i] == '\\') {
			pos2 = i;
		}
	}
	datasetFilename.erase(0, pos2 + 1);
	expResult << datasetFilename << "," << exe << "," << maxInstancesNum << "," << chrono::duration_cast<sec>(difference).count() << "," << parameter.PIthreshold << "," << parameter.thresholdDistance << "," << totalFilteredPatterns << "," << totalCandidatePatterns << "," << (Real)totalFilteredPatterns / totalCandidatePatterns << "," << degree2FilterTime << "," << totalFilterTime << "," << totalRefineTime << "," << totalFilterTime + totalRefineTime << "\n";
	expResult.close();
}
