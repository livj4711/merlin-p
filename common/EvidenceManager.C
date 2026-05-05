#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include <algorithm>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"
#include "Evidence.H"
#include "EvidenceManager.H"

EvidenceManager::EvidenceManager()
{
	foldCnt=1;
	preRandomizeSplit=false;
	randseed=0;
}

Error::ErrorCode
EvidenceManager::loadEvidenceFromFile(const char* inFName)
{
	ifstream inFile(inFName);
	char* buffer=NULL;
	string buffstr;
	int bufflen=0;

	// skip the first line (gene headers)
	if(inFile.good())
	{
		getline(inFile,buffstr);
	}

	while(inFile.good())
	{
		getline(inFile,buffstr);

		if(buffstr.length()<=0)
		{
			continue;
		}
		if(bufflen<=buffstr.length())
		{
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			bufflen=buffstr.length()+1;
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());

		//All the evidences for each variable are stored in a map, indexed by the varId
		EMAP* evidMap=new EMAP;
		char* tok=strtok(buffer,"\t");

		// The first token in each row is the cell/sample name //L
		if(tok==NULL) //L
		{ //L 
			continue; //L 
		} //L
		cellNames.push_back(tok); //L
		tok=strtok(NULL,"\t"); //L

		// cout << "Reading evidence for sample " << cellNames.back() << endl; //L
		int vId = 0;
		while(tok!=NULL)
		{
			Evidence* evid = new Evidence;
			evid->assocVariable(vId);
			//double varVal=log(atof(tok));
			double varVal=atof(tok);
			if(isinf(varVal) || isnan(varVal))
			{
				//cout <<"Found nan! " << tok << endl;
				cerr << "Please remove NaNs from the expression data or check the data format. Not a valid number: " << tok << endl;
				exit(-1);
			}
			evid->setEvidVal(varVal);
			(*evidMap)[vId]=evid;
			tok=strtok(NULL,"\t");
			vId++;
		}
		evidenceSet.push_back(evidMap);
	}

	inFile.close();

	cout <<"Number of samples read: " << evidenceSet.size() << endl;

	return Error::SUCCESS;
}



// Reads a two-column tab-delimited file (cell_name \t pseudotime_value) into
// cellPseudotime. Also clears any previously built pseudotimeOrder so it will
// be lazily rebuilt by buildPseudotimeOrder() the next time splitData() runs.
// Returns Error::SUCCESS on clean parse, Error::UNKNOWN on any I/O or format error.
Error::ErrorCode
EvidenceManager::readPseudotime(const char* inFName) //L
{
	ifstream inFile(inFName);
	if(!inFile.is_open())
	{
		cerr << "Error: pseudotime file path incorrect or file cannot be opened: " << inFName << endl;
		return Error::UNKNOWN;
	}

	cellPseudotime.clear();
	pseudotimeOrder.clear();

	string buffstr;
	while(getline(inFile,buffstr))
	{
		if(buffstr.length()<=0) // skip blank lines
		{
			continue;
		}
		char buffer[1024];
		strncpy(buffer,buffstr.c_str(),1023);
		buffer[1023]='\0';

		char* tok=strtok(buffer,"\t"); // extract the first tab-delimited token (cell name)
		if(tok==NULL) // if line contained no tab at all, treat as empty/corrupt
		{
			continue;
		}
		string cellName(tok); //L string cellName = first token in pseudotime file row
		tok=strtok(NULL,"\t");
		if(tok==NULL) // if second token missing, error
		{
			cerr << "Error: malformed pseudotime row for cell " << cellName << endl;
			return Error::UNKNOWN;
		}
		cellPseudotime[cellName]=atof(tok); //L cellPseudotime[cellName] = double(rank)
	}

	inFile.close();
	return Error::SUCCESS;
}

//We create a matrix of randomized evidence, where each evidence has some value
//for the random variables. We populate the matrix one random variable at a time.
//We first generate a vector of permuted indices, in the range of 0 to the total
//number of evidences. Then we populate the part of the matrix associated with
//this variable by querying values from the original matrix in the order specified
//by the permuted indices

int
EvidenceManager::randomizeEvidence(gsl_rng* r, VariableManager* vMgr)
{
	//First create all the evidence sets
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=new EMAP;
		randEvidenceSet.push_back(evidMap);
	}
	//Populate variable wise
	VSET& variableSet=vMgr->getVariableSet();
	int* randInds=new int[trainIndex.size()];
	for(VSET_ITER vIter=variableSet.begin();vIter!=variableSet.end();vIter++)
	{
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		populateRandIntegers(r,randInds,trainIndex,trainIndex.size());
		int j=0;
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=NULL;
			if(trainIndex.find(i)!=trainIndex.end())
			{
				evidMap=evidenceSet[randInds[j]];
				j++;
			}
			else
			{
				evidMap=evidenceSet[i];
			}
			EMAP* randEvidMap=randEvidenceSet[i];
			Evidence* evid=(*evidMap)[vIter->first];
			(*randEvidMap)[vIter->first]=evid;
		}
	}
	return 0;
}

EMAP*
EvidenceManager::getEvidenceAt(int evId)
{
	if((evId>=evidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return evidenceSet[evId];
}

EMAP*
EvidenceManager::getRandomEvidenceAt(int evId)
{
	if((evId>=randEvidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return randEvidenceSet[evId];
}

// Builds pseudotimeOrder: a list of expression row indices sorted by
// ascending pseudotime value. After this runs, pseudotimeOrder[k] gives the
// original row index of the k-th earliest cell, letting splitData() assign
// cross-validation folds in temporal order.
// Returns 0 on success, -1 if any validation check fails.
int
EvidenceManager::buildPseudotimeOrder() //L
{
	// if(cellPseudotime.empty()) //L it should never be, it only gets called if cellPseudotime is nonempty
	// {
	// 	pseudotimeOrder.clear();
	// 	return 0;
	// }

	if(cellPseudotime.size() != cellNames.size())
	{
		cerr << "Error: pseudotime entries do not match the number of expression rows" << endl;
		return -1;
	}

	// Have rankedRows[i] = (pseudotime_value, original_row_index) for cell i
	vector<pair<double,int>> rankedRows;
	rankedRows.reserve(cellNames.size()); // pre-allocate to avoid repeated heap reallocations in the loop
	for(int i=0;i<cellNames.size();i++)
	{
		map<string,double>::iterator pIter = cellPseudotime.find(cellNames[i]);
		if(pIter==cellPseudotime.end()) // if cell exists in the expression matrix but not in the pseudotime map
		{
			cerr << "Error: missing pseudotime for cell " << cellNames[i] << endl;
			return -1;
		}
		rankedRows.push_back(make_pair(pIter->second,i)); // record (pseudotime, row_index) for this cell
	}

	// Sort rankedRows ascending by pseudotime so index 0 = earliest cell, last = latest
	// Ties are broken by original row index
	sort(rankedRows.begin(), rankedRows.end(),
		[](const pair<double,int>& a, const pair<double,int>& b)
		{
			if(a.first == b.first)
			{
				return a.second < b.second; // row index tie-breaker: earlier row wins
			}
			return a.first < b.first;
		});

	// pseudotimeOrder[k] = original row index of k-th earliest cell
	pseudotimeOrder.clear();
	for(int i=0;i<rankedRows.size();i++)
	{
		pseudotimeOrder.push_back(rankedRows[i].second); // append the row index at sorted position i
	}

	return 0;
}

int
EvidenceManager::setFoldCnt(int f)
{
	foldCnt=f;
	return 0;
}

const vector<string>&
EvidenceManager::getCellNames() const //L
{
	return cellNames;
}

int
EvidenceManager::setPreRandomizeSplit()
{
	preRandomizeSplit=true;
	return 0;
}

int
EvidenceManager::splitData(int s)
{
	//L -->
	// If pseudotime values exist but ordering hasn’t been computed, build it
	if(!cellPseudotime.empty() && pseudotimeOrder.empty()) 
	{
		if(buildPseudotimeOrder()!=0) //L check if pseudotime order can be built successfully
		{
			return -1;
		}
	}
	if(!pseudotimeOrder.empty() && preRandomizeSplit)
	{
		cerr << "Error: preRandomizeSplit cannot be used with pseudotime-ordered folds" << endl;
		return -1;
	}

	const vector<int>* splitOrder = nullptr;

	// Use pseudotime order if available, otherwise use original order
	if(!pseudotimeOrder.empty())
	{
		splitOrder = &pseudotimeOrder;
	}
	int sampleCount = splitOrder ? splitOrder->size() : evidenceSet.size(); 
	//L <--

	// Define fold boundaries
	int testSetSize=sampleCount/foldCnt;
	int testStartIndex=s*testSetSize;
	int testEndIndex=(s+1)*testSetSize;
	if(s==foldCnt-1)
	{
		testEndIndex=sampleCount;
	}
	if(foldCnt==1)
	{
		testStartIndex=-1;
		testEndIndex=-1;
	}
	trainIndex.clear();
	testIndex.clear();

	// Optional randomization
	int* randInds=NULL;
	if(preRandomizeSplit)
	{
		randInds=new int[sampleCount];
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
		randseed=getpid();
		gsl_rng_set(r,randseed);
		populateRandIntegers(r,randInds,sampleCount);
		gsl_rng_free(r);
		cout <<"Random seed " << randseed << endl;
	}

	// Fold assignment loop
	for(int i=0;i<sampleCount;i++)
	{
		int eInd=splitOrder ? (*splitOrder)[i] : i; //L get index from pseudotime order if available, otherwise use original index
		if(randInds!=NULL)
		{
			eInd=randInds[i];
		}
		if((i>=testStartIndex) && (i<testEndIndex))
		{
			testIndex[eInd]=0;
		}
		else
		{
			trainIndex[eInd]=0;
		}
	}
	if(preRandomizeSplit)
	{
		delete[] randInds;
	}
	return 0;
}

INTINTMAP&
EvidenceManager::getTrainingSet()
{
	return trainIndex;
}

INTINTMAP&
EvidenceManager::getTestSet()
{
	return testIndex;
}

int
EvidenceManager::populateRandIntegers(gsl_rng* r, int* randInds,int size)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	for(int i=0;i<size;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds[i]=rind;
	}
	usedInit.clear();
	return 0;
}

int
EvidenceManager::populateRandIntegers(gsl_rng* r, int* randInds, INTINTMAP& populateFrom, int size)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	map<int,int> temp;
	INTINTMAP_ITER tIter=populateFrom.begin();
	for(int i=0;i<size;i++)
	{
		int tid=tIter->first;
		temp[i]=tid;
		tIter++;
	}
	for(int i=0;i<size;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds[i]=temp[rind];
	}
	usedInit.clear();
	return 0;
}

int
EvidenceManager::populateRandIntegers(gsl_rng* r, vector<int>& randInds,int size, int subsetsize)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	for(int i=0;i<subsetsize;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds.push_back(rind);
	}
	return 0;
}

