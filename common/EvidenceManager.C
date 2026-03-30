#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
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

EvidenceManager::~EvidenceManager()
{
}

int
EvidenceManager::setVariableManager(VariableManager* aPtr)
{
	vMgr=aPtr;
	return 0;
}


Error::ErrorCode
EvidenceManager::loadEvidenceFromFile(const char* inFName)
{
	ifstream inFile(inFName);
	char* buffer=NULL;
	string buffstr;
	int bufflen=0;
	int lineNo=0;

	// skip the first line (gene headers)
	if(inFile.good())
	{
		getline(inFile,buffstr);
	}

	while(inFile.good()) //L for each line (cell) in the data file
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
		char* tok=strtok(buffer,"\t"); //L all the gene expression values for this cell, as strings
		//The toks take the form of varid and value

		int vId = 0;
		while(tok!=NULL) //L for each gene expression value in the line (cell)
		{
			Evidence* evid = new Evidence;
			evid->assocVariable(vId);
			//double varVal=log(atof(tok));
			double varVal=atof(tok);
			if(isinf(varVal) || isnan(varVal))
			{
				//cout <<"Found nan! " << tok << endl;
				cerr << "Please remove NaNs from the expression data or check the data format, this is not a valid number (" << tok << ") " << endl; //L say NaN instead of zero
				exit(-1);	
			}
			evid->setEvidVal(varVal);
			(*evidMap)[vId]=evid;
			tok=strtok(NULL,"\t");
			vId++;
		}
		evidenceSet.push_back(evidMap); //L add the evidence map for this cell to the evidenceSet vector, which holds the evidence maps for all the cells
		lineNo++;
	}

	inFile.close();

	cout <<"Read " << evidenceSet.size() << " different cells " << endl; //L say cells instead of datapoints

	return Error::SUCCESS;
}

//We create a matrix of randomized evidence, where each evidence has some value
//for the random variables. We populate the matrix one random variable at a time.
//We first generate a vector of permuted indices, in the range of 0 to the total
//number of evidences. Then we populate the part of the matrix associated with
//this variable by querying values from the original matrix in the order specified
//by the permuted indices

int
EvidenceManager::randomizeEvidence(gsl_rng* r)
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
		string& geneName=(string&)vIter->second->getName();
		if((strcmp(geneName.c_str(),"FBgn0002631")==0) 
		|| (strcmp(geneName.c_str(),"FBgn0000411")==0) 
		|| (strcmp(geneName.c_str(),"FBgn0004915")==0) 
		|| (strcmp(geneName.c_str(), "FBgn0002573")==0)
		|| (strcmp(geneName.c_str(),"FBgn0005596")==0)  
		|| (strcmp(geneName.c_str(),"FBgn0035769")==0) 
		|| (strcmp(geneName.c_str(),"FBgn0011655")==0)
		|| (strcmp(geneName.c_str(),"FBgn0000576")==0))
		{
			cout <<geneName<<"IDs";
			for(int i=0;i<trainIndex.size();i++)
			{
				cout <<"\t" <<randInds[i];
			}
			cout << endl;
			cout <<geneName;
			for(INTINTMAP_ITER tIter=trainIndex.begin();tIter!=trainIndex.end();tIter++)
			{
				EMAP* emap=randEvidenceSet[tIter->first];
				cout <<"\t" << (*emap)[vIter->first]->getEvidVal();
			}
			cout << endl;
			cout <<geneName;
			for(INTINTMAP_ITER tIter=trainIndex.begin();tIter!=trainIndex.end();tIter++)
			{
				EMAP* emap=evidenceSet[tIter->first];
				cout <<"\t" << (*emap)[vIter->first]->getEvidVal();
			}
			cout << endl;
		}
	}
	return 0;
}

int 
EvidenceManager::getNumberOfEvidences()
{
	return evidenceSet.size();
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


int
EvidenceManager::addEvidence(EMAP* evidSet)
{
	evidenceSet.push_back(evidSet);
	return 0;
}


int
EvidenceManager::addToEvidence(int eSetID, int vId, INTDBLMAP& evidData)
{
	EMAP* emap=evidenceSet[eSetID];
	Evidence* evid=new Evidence;
	evid->setData(evidData);
	(*emap)[vId]=evid;
	return 0;
}

int
EvidenceManager::dumpEvidenceSet(ostream& oFile)
{
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=evidenceSet[i];
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			Evidence* evid=eIter->second;
			if(eIter!=evidMap->begin())
			{
				oFile<<"\t";
			}
			evid->dumpEvidence(oFile);
		}
		oFile << endl;
	}
	return 0;
}

int
EvidenceManager::getMLSettings(ostream& oFile)
{
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=evidenceSet[i];
		if(i==0)
		{
			for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
			{
				if(eIter!=evidMap->begin())
				{
					oFile<<"\t";
				}
				oFile<< eIter->first;
			}
			oFile << endl;
		}

		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			if(eIter!=evidMap->begin())
			{
				oFile<<"\t";
			}
			Evidence* evid=eIter->second;
			oFile << evid->getMLVal();
		}
		oFile << endl;
	}
	return 0;
}


int 
EvidenceManager::setFoldCnt(int f)
{
	foldCnt=f;
	return 0;
}

int
EvidenceManager::generateValidationSet(const char* vFName, int vSetSize,gsl_rng* r)
{
	ifstream inFile(vFName);
	if(inFile.good())
	{
		char buffer[256];
		while(inFile.good())
		{
			inFile.getline(buffer,255);
			if(strlen(buffer)<=0)
			{
				continue;
			}
			int dId=atoi(buffer);
			validationIndex[dId]=0;
		}
		inFile.close();
	}
	else
	{
		populateRandIntegers(r,validationIndex,evidenceSet.size(),vSetSize);
		ofstream oFile(vFName);
		for(INTINTMAP_ITER vIter=validationIndex.begin();vIter!=validationIndex.end();vIter++)
		{
			oFile << vIter->first << endl;
		}
		oFile.close();
	}
	return 0;
}


int 
EvidenceManager::setPreRandomizeSplit()
{
	preRandomizeSplit=true;
	return 0;
}

int 
EvidenceManager::splitData(int s) //L split data for k-fold cross-validation
{
	int testSetSize=(evidenceSet.size()-validationIndex.size())/foldCnt; //L evidenceSet.size()=number of cells. validationIndex NEVER GETS CALLED, so its size is 0
	int testStartIndex=s*testSetSize;
	int testEndIndex=(s+1)*testSetSize;
	if(s==foldCnt-1) //L if this is the last fold
	{
		testEndIndex=evidenceSet.size()-validationIndex.size(); //L set the last index to the leftover cells
	}
	if(foldCnt==1) //L if we aren't doing cross-fold validation on the data
	{
		testStartIndex=-1;
		testEndIndex=-1;
	}
	trainIndex.clear();
	testIndex.clear();
	int m=0;
	int* randInds=NULL;
	if(preRandomizeSplit) //L preRandomizeSplit is always false
	{ //L generates random permutation of the dataset indices so that later data partitioning happens in random order instead of sequential order.
		randInds=new int[evidenceSet.size()];
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
		randseed=getpid();
		gsl_rng_set(r,randseed);
		populateRandIntegers(r,randInds,evidenceSet.size()); //L fills the array with a random permutation of integers
		gsl_rng_free(r);
		cout <<"Random seed " << randseed << endl;
	}
	for(int i=0;i<evidenceSet.size();i++)
	{
		int eInd=i;
		if(randInds!=NULL)
		{
			eInd=randInds[i];
		}
		if(validationIndex.find(eInd)!=validationIndex.end()) //L if this cell exists in validationIndex, skip it. but validationIndex is always empty
		{
			continue;
		}
		if((m>=testStartIndex) && (m<testEndIndex)) //L if it's a 'test' cell, add the index to the testIndex
		{
			testIndex[eInd]=0;
		}
		else //L if it's not a 'test' cell, add the index to the trainIndex
		{
			trainIndex[eInd]=0;
		}
		m++; //L m increments with the cell
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

INTINTMAP&
EvidenceManager::getValidationSet()
{	
	return validationIndex;
}


int 
EvidenceManager::standardizeData()
{
	VSET& varSet=vMgr->getVariableSet();
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		double mean=0;
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=evidenceSet[i];
			mean=mean+(*evidMap)[vIter->first]->getEvidVal();
		}
		mean=mean/evidenceSet.size();
		double std=0;
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=evidenceSet[i];
			double diff=mean-(*evidMap)[vIter->first]->getEvidVal();
			std=std+(diff*diff);
		}
		std=sqrt(std/(evidenceSet.size()-1));
		//Standardize
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=evidenceSet[i];
			Evidence* evid=(*evidMap)[vIter->first];
			double tval=evid->getEvidVal();
			double sval=(tval-mean)/std;
			evid->setEvidVal(sval);
		}
	}
	return 0;
}

//J remove partitionData()

int 
EvidenceManager::populateEvidence(Evidence** evid,const char* evidStr)
{
	//first check for validity of evidStr
	if(strchr(evidStr,'=')==NULL)
	{
		return -1;
	}
	*evid=new Evidence;
	
	INTDBLMAP evidData;
	int currInd=0;
	int ttInd=0;
	int tokId=0;
	char tempTok[256];
	while(evidStr[currInd]!='\0')
	{
		if((evidStr[currInd]=='=') || 
		   (evidStr[currInd]==']') ||
		   (evidStr[currInd]==',')
		  )
		{
			tempTok[ttInd]='\0';
			ttInd=0;
			if(tokId==0)
			{
				//This is the variable
				int vId=atoi(tempTok);
				Variable* var=vMgr->getVariableAt(vId);
				var->initEvidence(evidData);
				(*evid)->assocVariable(vId);
			}
			else
			{
				char* pos=strchr(tempTok,'|');
				//Hard evidence
				if(pos==NULL)
				{
					int varVal=atoi(tempTok);
					if(evidData.find(varVal)==evidData.end())
					{
						cout <<"No value "<< varVal << " in the domain of  a variable" << endl;
						return -1;
					}
					evidData[varVal]=1.0;
					(*evid)->setType(Evidence::HARD);
				}
				else
				{
					*pos='\0';
					int varVal=atoi(tempTok);
					double varValProb=atof(pos+1);
					if(evidData.find(varVal)==evidData.end())
					{
						cout <<"No value "<< varVal << " in the domain of  a variable" << endl;
						return -1;
					}
					evidData[varVal]=varValProb;
					//Will be setting it multiple times but its ok for now.
					(*evid)->setType(Evidence::SOFT);
				}
			}
			tokId++;
		}
		else if(evidStr[currInd]!='[')
		{
			tempTok[ttInd]=evidStr[currInd];
			ttInd++;
		}
		currInd++;
	}
	(*evid)->setData(evidData);
	return 0;
}

int
EvidenceManager::dumpSummaryStat(ostream& oFile)
{
	//Need an evidence like object but to store the frequency over all evidences
	//rather than a single evidence
	map<int,INTDBLMAP*> summary;
	map<int,double> normFactors;
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=evidenceSet[i];
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			INTDBLMAP* evCnt=NULL;
			if(summary.find(eIter->first)==summary.end())
			{
				evCnt=new INTDBLMAP;
				summary[eIter->first]=evCnt;
			}
			else
			{
				evCnt=summary[eIter->first];
			}
			//Get data and add to evCnt
			INTDBLMAP& data=eIter->second->getData();
			for(INTDBLMAP_ITER idIter=data.begin();idIter!=data.end();idIter++)
			{
				if(evCnt->find(idIter->first)==evCnt->end())
				{
					(*evCnt)[idIter->first]=idIter->second;
				}
				else
				{
					(*evCnt)[idIter->first]=(*evCnt)[idIter->first]+idIter->second;
				}
				//Add the normalization factor for all the freq or exp. freq cnts
				if(normFactors.find(eIter->first)==normFactors.end())
				{
					normFactors[eIter->first]=idIter->second;
				}
				else
				{
					normFactors[eIter->first]=normFactors[eIter->first]+idIter->second;
				}
			}
		}
	}

	//Now iterate over the evidence summary, normalize and display values
	for(map<int,INTDBLMAP*>::iterator aIter=summary.begin();aIter!=summary.end();aIter++)
	{
		double normConst=normFactors[aIter->first];
		INTDBLMAP* evCnt=aIter->second;
		oFile <<"Distribution of "<< aIter->first;
		for(INTDBLMAP_ITER idIter=evCnt->begin();idIter!=evCnt->end();idIter++)
		{
			oFile << " " << idIter->first<<"=" << idIter->second/normConst;
		}
		oFile << endl;
	}
	return 0;
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
EvidenceManager::populateRandIntegers(gsl_rng* r, INTINTMAP& randInds,int size, int subsetsize)
{
	double step=1.0/(double)size;
	for(int i=0;i<subsetsize;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(randInds.find(rind)!=randInds.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		randInds[rind]=0;
	}
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

