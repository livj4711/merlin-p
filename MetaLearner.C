#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <time.h>

#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"

#include "Evidence.H"
#include "EvidenceManager.H"

#include "Potential.H"
#include "SlimFactor.H"
#include "PotentialManager.H"

#include "FactorGraph.H"
#include "MetaMove.H"
#include "HierarchicalClusterNode.H"
#include "HierarchicalCluster.H"
#include "HyperGeomPval.H"
#include "MetaLearner.H"

MetaLearner::MetaLearner()
{
	restrictedFName[0]='\0';
	trueGraphFName[0]='\0';
	preRandomizeSplit=false;
	random=false;
	lambda=0;
	clusterThreshold=0.5;
	specificFold=-1;
	convThreshold=1e-3;
	factorGraph=nullptr;
	currPLL=nullptr;
	correlationDistances=nullptr;
}

MetaLearner::~MetaLearner()
{
}

int
MetaLearner::setMaxFactorSize_Approx(int aVal)
{
	maxFactorSizeApprox=aVal;
	return 0;
}

// int 
// MetaLearner::setPenalty(double aVal)
// {
// 	penalty=aVal;
// 	return 0;
// }

int
MetaLearner::setBeta1(double aval)
{
	beta1=aval;
	return 0;
}

int
MetaLearner::initEdgePriorMeta_All() //L reverse of setPriorGraph_All...
{
	for(map<string,map<string,map<string,double>*>*>::iterator gIter=priorgraphmap.begin();gIter!=priorgraphmap.end();gIter++) //L for each prior graph in the priorgraphmap
	{
		map<string,map<string,double>*>* priorgraph = gIter->second; //L get the actual prior graph (edge information) for this prior graph
		map<int,INTDBLMAP*>* edgeprior = new map<int,INTDBLMAP*>();
		edgepriormap[gIter->first] = edgeprior; //L initialize a new empty graph into edgepriormap for this prior graph
		initEdgePriorMeta(gIter->first,*priorgraph,*edgeprior); //L i edited to send  prior name as well. This initializes edgeprior as the "reverse" of priorgraph, which is keyed by the gene whose values are its regulators in the prior (which have edge weights)
	}
	return 0;
}

int
MetaLearner::setPriorGraph_All(const char* aFName)
{
	ifstream inFile(aFName);
    if (!inFile.is_open())
    {
        std::cerr << "Error: Prior config file path incorrect or file cannot be opened: " << aFName << std::endl;
    }
	
	char buffer[1024];
	while(inFile.good()) //L for each prior network provided in the prior config file
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string gname;
		string fname;
		double gbeta;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0) //L first tok = prior graph name
			{
				gname.append(tok);
			}
			else if(tokCnt==1) //L second tok = prior graph file location/name
			{
				fname.append(tok);
			}
			else if(tokCnt==2) //L third tok = beta/confidence for this prior graph
			{
				gbeta=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		betamap[gname] = gbeta;
		map<string,map<string,double>*>* priorGraph = new map<string,map<string,double>*>();
		setPriorGraph(fname.c_str(),*priorGraph); //L read all edge information of this prior network into priorGraph
		priorgraphmap[gname] = priorGraph;
	}
	inFile.close();
	return 0;
}

int
MetaLearner::setBeta_Motif(double aval)
{
	beta_motif=aval;
	return 0;
}

int
MetaLearner::setLambda(double l)
{
	lambda=l;
	return 0;
}

int 
MetaLearner::setConvergenceThreshold(double aVal)
{
	convThreshold=aVal;
	return 0;
}

int
MetaLearner::setRestrictedList(const char* aFName)
{
	strcpy(restrictedFName,aFName);
	ifstream inFile(restrictedFName);
	string buffer;

	int count = 0; // counter for number of restricted regulators

	while(inFile.good()) //L for each line (regulator) in the restricted regulator list file
	{
		getline(inFile,buffer);
		if(buffer.length()<=0)
		{
			continue;
		}
		restrictedVarList[buffer]=0;
		count++;
	}
	inFile.close();
	std::cout << "Number of regulators read: " << count << std::endl;
	return 0;
}


int 
MetaLearner::setPreRandomizeSplit()
{
	preRandomizeSplit=true;
	return 0;
}

int
MetaLearner::setGlobalEvidenceManager(EvidenceManager* anEvMgr)
{
	evidenceManager=anEvMgr;
	return 0;
}

int 
MetaLearner::setVariableManager(VariableManager* aPtr)
{
	varManager=aPtr;
	return 0;
}

int
MetaLearner::setOutputDirName(const char* dirPath)
{
	strcpy(outputDirName,dirPath);
	return 0;
}

int 
MetaLearner::setClusteringThreshold(double aVal)
{
	clusterThreshold=aVal;
	return 0;
}
 

int
MetaLearner::setSpecificFold(int fid)
{
	specificFold=fid;
	return 0;
}

int 
MetaLearner::setPriorGraph(const char* aFName, map<string,map<string,double>*>& priorGraph)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good()) //L for each line (edge) in the prior graph file
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string tfName;
		string tgtName;
		double edgeStrength;
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0) //L first tok = TF name
			{
				tfName.append(tok);
			}
			else if(tokCnt==1) //L second tok = target gene name
			{
				tgtName.append(tok);
			}
			else if(tokCnt==2) //L third tok = edge strength
			{
				edgeStrength=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		map<string,double>* tgtSet=NULL;
		if(priorGraph.find(tfName)==priorGraph.end()) //L if this TF is not already in the prior graph, add it with an empty target set
		{
			tgtSet=new map<string,double>;
			priorGraph[tfName]=tgtSet;
		}
		else //L else if this TF is already in the prior graph, get its target set
		{
			tgtSet=priorGraph[tfName];
		}
		(*tgtSet)[tgtName]=edgeStrength;
	}
	inFile.close();
	return 0;
}

int
MetaLearner::setRandom(bool flag)
{
	random=flag;
	return 0;
}

int
MetaLearner::readModuleMembership(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good()) //L for each line [gene, moduleID] in the module membership file
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		string geneName;
		int moduleID;
		int tokCnt=0;
		char* tok=strtok(buffer,"\t");
		while(tok!=NULL)
		{
			if(tokCnt==0) //L first tok = gene name
			{
				geneName.append(tok);
			}
			else if(tokCnt==1) //L second tok = module ID
			{
				moduleID=atoi(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		map<string,int>* geneSet=NULL;
		if(moduleGeneSet.find(moduleID)==moduleGeneSet.end())
		{
			geneSet=new map<string,int>;
			moduleGeneSet[moduleID]=geneSet;
		}
		else
		{
			geneSet=moduleGeneSet[moduleID];
		}
		(*geneSet)[geneName]=0;
		geneModuleID[geneName]=moduleID;
	}
	inFile.close();
	return 0;
}

int
MetaLearner::setDefaultModuleMembership()
{
	VSET& varSet=varManager->getVariableSet();
	int vCnt=varSet.size();
	int moduleCnt=(int) sqrt(vCnt/2);
	if(moduleCnt>30)
	{
		moduleCnt=30;
	}
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	//Randomly partition the variables into clusterassignments
	vector<int> randIndex;
	double step=1.0/(double)vCnt;
	map<int,int> usedInit;
	for(int i=0;i<vCnt;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randIndex.push_back(rind);
	}
	//For each partition estimate the mean and covariance
	int clusterSize=vCnt/moduleCnt;
	for(int e=0;e<moduleCnt;e++)
	{
		int startInd=e*clusterSize;
		int endInd=(e+1)*clusterSize;
		if(e==moduleCnt-1)
		{
			endInd=clusterSize;
		}
		map<string,int>* geneSet=NULL;
		geneSet=new map<string,int>;
		moduleGeneSet[e]=geneSet;
		for(int i=startInd;i<endInd;i++)
		{
			int dataId=randIndex[i];
			Variable* v=varSet[dataId];
			(*geneSet)[v->getName()]=0;
			geneModuleID[v->getName()]=e;
		}
	}
	randIndex.clear();
	usedInit.clear();
	return 0;
}

//J completely remove initPartitions()

int
MetaLearner::initEdgePriorMeta(const string& priorName, map<string,map<string,double>*>& graph, map<int,INTDBLMAP*>& edgePriors)
{
	VSET& varSet=varManager->getVariableSet();
	cout << "Initializing prior: \"" << priorName << "\" " << endl;
	for(map<string,int>::iterator rIter=restrictedVarList.begin();rIter!=restrictedVarList.end();rIter++) //L for each restricted regulator
	{
		int regId=varManager->getVarID(rIter->first.c_str());
		if(regId==-1) //Lif regulator wasn't in expression matrix, skip
		{
			continue;
		}
		if(graph.find(rIter->first)==graph.end()) //L if this regulator is not in the prior graph, skip
		{
			continue;
		}
		int tfhit=0;
		map<string,double>* tgtSet=graph[rIter->first]; //L tgtSet is the set of targets for this regulator in the prior graph
		for(map<string,double>::iterator vIter=tgtSet->begin();vIter!=tgtSet->end();vIter++) //L for each target of this regulator
		{
			INTDBLMAP* edgePriorGene=NULL;
			int tgtId=varManager->getVarID(vIter->first.c_str());
			if(tgtId==-1)
			{
				continue;
			}
			if(edgePriors.find(tgtId)==edgePriors.end()) //L if this target gene is not already in edgePriors, add it with an empty set of regulators
			{
				edgePriorGene=new INTDBLMAP;
				edgePriors[tgtId]=edgePriorGene;
			}
			else
			{
				edgePriorGene=edgePriors[tgtId];
			}
			double ewt=fabs(vIter->second); //L take the absolute value of the edge strength as the edge weight
			if(edgePriorGene->find(regId)==edgePriorGene->end()) //L if this regulator is not already a regulator of this target in edgePriors, add it with the edge weight as its value
			{
				//(*edgePriorGene)[regId]=vIter->second;
				(*edgePriorGene)[regId]=ewt; //l * just points to this gene key in edgePriors
			}
			else //L else if this regulator is already a regulator of this target in edgePriors, add the edge weight to the existing value (in case there are multiple edges from this regulator to this target in the prior graph?)
			{
				//(*edgePriorGene)[regId]=(*edgePriorGene)[regId]+vIter->second;
				(*edgePriorGene)[regId]=(*edgePriorGene)[regId]+ewt;
			}
			tfhit++;
		}
		//L cout << " Regulator "<< rIter->first << " has " << tfhit << " targets" << endl; //L comment this to reduce verbosity
	}

	return 0;
}

int
MetaLearner::doCrossValidation(int foldCnt)
{
	gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
	rnd=gsl_rng_alloc(gsl_rng_default);

	evidenceManager->setFoldCnt(foldCnt);
	evidenceManager->splitData(0); //L I think this splitData is redundant...

	potManager = new PotentialManager;

	//The first key is for the fold number
	//For each fold we have a trained model. For each trained model we have the likelihood on 
	//all the test sets, including the self test.
	int foldBegin=0;
	int foldEnd=foldCnt;
	if(specificFold>-1)
	{
		foldBegin=specificFold;
		foldEnd=specificFold+1;
	}
	for(int f=foldBegin;f<foldEnd;f++) //L loop over cross-validation folds
	{
		//J remove loop through evMgrSet (set of all conditions)	
		evidenceManager->splitData(f);
		if(random)
		{	
			evidenceManager->randomizeEvidence(r, varManager);
		}

		vector<int> regIDs;
		for (map<string,int>::iterator iter = restrictedVarList.begin(); iter != restrictedVarList.end(); iter++)
		{
			int regID = varManager->getVarID(iter->first.c_str());
			regIDs.push_back(regID);
		}

		potManager->init(evidenceManager, random, regIDs);

		factorGraph = new FactorGraph(varManager);

		//L Prepare output directory for this fold
		char outputDir[1024];
		sprintf(outputDir,"%s/fold%d",outputDirName,f);
		char foldOutputDirCmd[1024];
		sprintf(foldOutputDirCmd,"mkdir -p %s",outputDir);
		system(foldOutputDirCmd);

		//L Begin identifying regulators/inferring modules for this fold
		start(f);

		getPredictionError_CrossValid(f);
		clearFoldSpecData(); //L clear fold-specific structures before the next fold
	}
	gsl_rng_free(r);

	gsl_rng_free(rnd);
	return 0;
}

int
MetaLearner::start(int f)
{
	currFold=f;
	sprintf(foldoutDirName,"%s/fold%d",outputDirName,f); //L set the output directory for this fold ("examples/out_dir/fold0")
	//L int maxMBSizeApprox=maxFactorSizeApprox-1; //L remove to make this more intuitive
	int maxNumRegs = maxFactorSizeApprox-1; // max num of regulators a gene can have
	rnd=gsl_rng_alloc(gsl_rng_default);
	int rseed=getpid();
	gsl_rng_set(rnd,rseed);
	cout << "Random seed: " << rseed << endl;
	initEdgePriorMeta_All(); //L initialize "edgeprior" for each prior network, which is the "reverse" of priorgraph (key=genes, values=regulators/edgewt dict)
	initEdgeSet();
	initPhysicalDegree();

	if(strlen(trueGraphFName)!=0)
	{
		return 0;
	}

	VSET& varSet=varManager->getVariableSet();

	for (VSET_ITER vIter=varSet.begin(); vIter != varSet.end(); vIter++)
	{
		Variable *var = vIter->second;
		variableStatus[var->getName()] = 0;
	}

	double currGlobalScore=getInitPLLScore();
	double initScore=getInitPrior();
	int showid=0;
	int moduleiter=0;
	bool notConvergedTop=true;
	while(moduleiter<1 && notConvergedTop)
	{
		int iter=0;
		bool notConverged=true;
		while(notConverged && iter<50)
		{
			int attemptedMoves=0;
			int subiter=0;
			double scorePremodule=currGlobalScore;
			while(subiter<varSet.size())
			{
				int vID=subiter;
				Variable* v=varSet[vID];

				// If 5 iterations have passed without finding a score improving parent, then skip.
				int lastiter = variableStatus[v->getName()];
				if((iter - lastiter) >= 5)
				{
					cout <<"Skipping " << v->getName() << endl;
					subiter++;
					continue;
				}

				MetaMove* nextMove = getNextMove(maxNumRegs, vID);
				if (nextMove == nullptr)
				{
					subiter++;
					continue;
				}

				makeMove(nextMove, iter);
				delete nextMove;

				currGlobalScore=getPLLScore();

				subiter++;
				showid++;
				attemptedMoves++;
			}
			if((currGlobalScore-scorePremodule)<=convThreshold)
			{
				notConverged=false;
			}
			else
			{
				redefineModules();
			}
			iter++;
			scorePremodule=currGlobalScore;
			dumpAllGraphs(maxNumRegs,f,iter);
		}
		moduleiter++;
	}
	cout <<"Final Score " << currGlobalScore << endl;
	finalScores[f]=currGlobalScore;
	return 0;
}

double
MetaLearner::getInitPLLScore()
{
	double initScore=0;
	VSET& varSet=varManager->getVariableSet();
	//Initially we just sum up the marginal likelihoods
	currPLL=new INTDBLMAP;
	for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		if(varNeighborhoodPrior.find(vIter->first)==varNeighborhoodPrior.end())
		{
			continue;
		}
		Variable* var=varSet[vIter->first];
		double newPLL_s=getInitPLLScore(vIter->first);
		double priorScore=varNeighborhoodPrior[vIter->first];
		(*currPLL)[vIter->first]=newPLL_s+priorScore;
		initScore=initScore+(*currPLL)[vIter->first];
	}
	return initScore;
}

double
MetaLearner::getPLLScore()
{
	double gScore=0;
	for(INTDBLMAP_ITER dIter=currPLL->begin();dIter!=currPLL->end();dIter++)
	{
		if(isnan(gScore) || isinf(gScore))
		{
			cout << "Found nan/inf for variable " << dIter->first << endl;
		}
		gScore=gScore+dIter->second;
	}
	return gScore;
}

double 
MetaLearner::getInitPrior()
{
	double graphPrior=0;
	double edgePresence=1/(1+exp(-1*beta1));
	for(map<string,double>::iterator aIter=edgePresenceProb.begin();aIter!=edgePresenceProb.end();aIter++)
	{
		//graphPrior=graphPrior+log(1-edgePresence);
		graphPrior=graphPrior+log(1-aIter->second);
		if(isinf(graphPrior)|| isnan(graphPrior))
		{
			cout <<"Graph prior is "<< graphPrior << " after " << aIter->first << " for " << aIter->second << endl;
		}
	}
	return graphPrior;
}

int
MetaLearner::clearFoldSpecData()
{
	if (factorGraph != nullptr)
	{
		delete factorGraph;
		factorGraph = nullptr;
	}
	edgeMap.clear();
	if (currPLL != nullptr)
	{
		delete currPLL;
		currPLL = nullptr;
	}
	return 0;
}

int
MetaLearner::initEdgeSet()
{
	VSET& varSet=varManager->getVariableSet();
	for(VSET_ITER uIter=varSet.begin();uIter!=varSet.end();uIter++) //L for each reg u
	{
		Variable* u=varSet[uIter->first];
		if((restrictedVarList.size()>0) && (restrictedVarList.find(u->getName())==restrictedVarList.end())) //L if there is a reg list and regulator u is not in it, skip
		{
			continue;
		}

		for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++) //L for each gene v
		{
			if(uIter->first==vIter->first) //L skip self loops u -> v
			{
				continue;
			}
			Variable* v=varSet[vIter->first];
			if(geneModuleID.find(v->getName())==geneModuleID.end()) //L if gene isnt associated with a module (it alwasy should be), skip.
			{	
				continue;
			}
			string edgeKey;
			//This is going to be a directed graph. edgeKey looks like "reg_name\tgene_name"
			edgeKey.append(u->getName().c_str());
			edgeKey.append("\t");
			edgeKey.append(v->getName().c_str());
			//L edgeMap is an empty map default-initialized when MetaLearner is constructed. we will be a binary indicator map
			edgeMap[edgeKey]=0; //initialize this edge to absent for the future LEARNED graph

			double initPrior=getEdgePrior(uIter->first,vIter->first); //L beta_1 + sum(beta_g * edgewt) across all prior networks. beta_1 is sparsity param (default -5), beta_g is confidence in a prior.
			initPrior = 1/(1+exp(-1*initPrior)); //L convert this log-odds edge prior to probability between 0 and 1
			if(initPrior<1e-6)
			{
				initPrior=1e-6;
			}
			if(initPrior==1)
			{
				initPrior=1-1e-6;
			}
			edgePresenceProb[edgeKey]=initPrior;

			//L varNeighborhoodPrior is a map of gene ID -> sum of log(1-initPrior) for all edges pointing to this gene.
			if(varNeighborhoodPrior.find(vIter->first)==varNeighborhoodPrior.end()) //L varNeighborhoodPrior default-initialized to an empty map when MetaLearner is constructed
			{ 
				varNeighborhoodPrior[vIter->first]=log(1-initPrior);
			}
			else
			{
				varNeighborhoodPrior[vIter->first]=varNeighborhoodPrior[vIter->first]+log(1-initPrior);
			}
		}
	}
	//L cout <<"Restricted varlist size: " << restrictedVarList.size() << endl; //L we print the resticted regulators list size previously
	int n=varSet.size();
	int r=restrictedVarList.size();
	//L int expEdgeCnt=((r*(r-1))/2) + (r*(n-r)) ; //L this is the expected edge count for undirected edges
	int expEdgeCnt=r*(n-1); //L this is the expected edge count for directed edges 
	//L cout <<"Initialized " << edgeMap.size() << " edges. Expected " << expEdgeCnt << endl; //L not necessary to print

	// Init the potentials
	for(int f=0;f<factorGraph->getFactorCnt();f++)
	{
		SlimFactor* sFactor=factorGraph->getFactorAt(f);
		sFactor->potFunc=potManager->createPotential(sFactor->fId);
	}

	return 0;
}

int
MetaLearner::getPredictionError_CrossValid(int foldid)
{
	VSET& varSet=varManager->getVariableSet();
	char foldoutDirName[1024];
	sprintf(foldoutDirName,"%s/fold%d",outputDirName,foldid);
	INTINTMAP& testSet=evidenceManager->getTestSet();
	map<int,double> varPLL;
	for(INTINTMAP_ITER dIter=testSet.begin();dIter!=testSet.end();dIter++)
	{
		//for each gc, get the expected value of this datapoint
		EMAP* evidMap=evidenceManager->getEvidenceAt(dIter->first);

		for(map<string,int>::iterator vIter=geneModuleID.begin();vIter!=geneModuleID.end();vIter++)
		{
			int vId=varManager->getVarID(vIter->first.c_str());
			if(vId==-1)
			{
				continue;
			}
			Variable* v=varSet[vId];
			SlimFactor* sFactor=factorGraph->getFactorAt(vId);
			Potential* sPot=sFactor->potFunc;
			if(sPot==NULL)
			{
				cout <<"Found null for factor="<< sFactor->fId
					<< " variable=" <<varSet[sFactor->fId]->getName() << endl;
			}
			double pval=sPot->evaluateProbabilityDensity(evidMap);
			if(pval<1e-50)
			{
				pval=1e-50;
			}
			if(isinf(pval) || isnan(pval))
			{
				cout <<"Stop here. Found nan/inf for " << vIter->first << " dtpt "<< dIter->first << endl;
			}
			double cll=log(pval);
			if(varPLL.find(vId)==varPLL.end())
			{
				varPLL[vId]=cll;
			}
			else
			{
				varPLL[vId]=varPLL[vId]+cll;
			}
		}
	}
	/*
	for(map<int,double>::iterator pIter=varPLL.begin();pIter!=varPLL.end();pIter++)
	{
		oFile << varSet[pIter->first]->getName() << "\t" << pIter->second << endl;
	}
	pFile << "\tRMSE\tNormRMSE\tCoeff_Det_aka_R^2\tCC"<< endl;
	*/
	vector<double> truevect;
	vector<double> predvect;
	for(map<string,int>::iterator vIter=geneModuleID.begin();vIter!=geneModuleID.end();vIter++)
	{
		int vId=varManager->getVarID(vIter->first.c_str());
		if(vId==-1)
		{
			continue;
		}
		//pFile <<vIter->first;
		double norm=0;
		double maxval=-100000;
		double minval=1000000;
		double totalvar=0;
		double truemean=0;
		truevect.clear();
		predvect.clear();

		for(INTINTMAP_ITER dIter=testSet.begin();dIter!=testSet.end();dIter++)
		{
			EMAP* evidMap=evidenceManager->getEvidenceAt(dIter->first);
			Evidence* evid=(*evidMap)[vId];
			double trueval=evid->getEvidVal();
			truemean=truemean+trueval;
			truevect.push_back(trueval);
		}

		truemean=truemean/((double)testSet.size());

		//First the predicted time course
		SlimFactor* sFactor=factorGraph->getFactorAt(vId);
		Potential* sPot=sFactor->potFunc;
		for(INTINTMAP_ITER dIter=testSet.begin();dIter!=testSet.end();dIter++)
		{
			EMAP* evidMap=evidenceManager->getEvidenceAt(dIter->first);
			double predval=sPot->getExpectation(evidMap);
			Evidence* evid=(*evidMap)[vId];
			double trueval=evid->getEvidVal();
			totalvar=totalvar+((trueval-truemean)*(trueval-truemean));
			//also called residuals
			predvect.push_back(predval);
			//norm=norm+(trueval*trueval);
			norm=norm+1;
			if(trueval>maxval)
			{
				maxval=trueval;
			}
			if(trueval<minval)
			{
				minval=trueval;
			}
		}
	}
	//oFile.close();
	//pFile.close();
	varPLL.clear();
	return 0;
}

MetaMove*
MetaLearner::getNextMove(int maxNumRegs, int vID)
{
	VSET& varSet=varManager->getVariableSet();
	Variable* v = varSet[vID];

	if(geneModuleID.find(v->getName()) == geneModuleID.end())
	{
		return nullptr;
	}

	// If v already has the max number of parents, dont test adding another.
	SlimFactor* dFactor = factorGraph->getFactorAt(vID);
	if(dFactor->mergedMB.size() >= maxNumRegs)
	{
		return nullptr;
	}

	// Collect the new set of parents for v
	vector<int> parentIDs;
	for (INTINTMAP_ITER iter = dFactor->mergedMB.begin(); iter != dFactor->mergedMB.end(); iter++)
	{
		parentIDs.push_back(iter->first);
	}

	int moduleID=geneModuleID[v->getName()];
	map<string,int>* moduleMembers=moduleGeneSet[moduleID];

	double bestTargetScore=0;
	double bestScoreImprovement=0;
	Variable* bestu=NULL;
	Potential* bestPot=NULL;

	for(map<string,int>::iterator uIter=restrictedVarList.begin();uIter!=restrictedVarList.end();uIter++)
	{
		int regID=varManager->getVarID(uIter->first.c_str());

		// Ensure we can find the regulator, and that it isnt the same node as the target.
		if(regID==-1 || vID==regID)
		{
			continue;
		}

		Variable* u=varSet[regID];

		string edgeKey;
		edgeKey.append(u->getName().c_str());
		edgeKey.append("\t");
		edgeKey.append(v->getName().c_str());

		// If the edge already exists, no need to test adding it.
		int edgeValue=edgeMap[edgeKey];
		if(edgeValue==1)
		{
			continue;
		}

		double improvement = 0;
		double score = 0;
		Potential* aPot = NULL;

		parentIDs.push_back(u->getID());

		getNewPLLScore(u, v, parentIDs, edgeKey, score, improvement, &aPot);

		parentIDs.pop_back();

		bool betterMoveExists = (bestu != NULL) && (bestScoreImprovement >= improvement);

		// If there is no score improvement, cleanup aPot and continue.
		if (improvement < 0 || betterMoveExists)
		{
			if (aPot != NULL)
			{
				delete aPot;
			}
			continue;
		}

		if(bestPot != NULL)
		{
			delete bestPot;
		}

		bestu = u;
		bestTargetScore = score;
		bestScoreImprovement = improvement;
		bestPot = aPot;
	}

	// We could not find a parent to add to v that would improve the score.
	if((bestu == NULL) || (bestScoreImprovement <= 0))
	{
		return nullptr;
	}

	MetaMove* nextMove = new MetaMove;
	nextMove->setSrcVertex(bestu->getID());
	nextMove->setTargetVertex(v->getID());
	nextMove->setTargetMBScore(bestTargetScore);
	nextMove->setScoreImprovement(bestScoreImprovement);
	nextMove->setDestPot(bestPot);
	return nextMove;
}

void
MetaLearner::getNewPLLScore(Variable* u, Variable* v, vector<int>& parentIDs, string& edgeKey, double& mbScore, double& scoreImprovement, Potential** newdPot)
{
	int factorID = v->getID();
	VSET& varSet = varManager->getVariableSet();

	double plus = 0;
	double minus = 0;
	for (vector<int>::iterator iter = parentIDs.begin(); iter != parentIDs.end(); iter++)
	{
<<<<<<< Updated upstream
		Variable* parentVar = varSet[*iter];
		double eprior = getEdgePrior(*iter, factorID);
		double moduleContrib = getModuleContribLogistic((string&)v->getName(), (string&)parentVar->getName());
		double edgeProb = 1 / (1 + exp(-1 * (eprior + moduleContrib)));
		double edgeProbOld = 1 / (1 + exp(-1 * eprior));
		minus += log(1 - edgeProbOld);
		plus += log(edgeProb);
=======
		int nID=varManager->getVarID(mIter->first.c_str());
		if(nID==vID)
		{
			continue;
		}

		if(nID==uID)
		{
			continue;
		}
		SlimFactor* sFactor=factorGraph->getFactorAt(nID);
		if(sFactor->mergedMB.find(uID)!=sFactor->mergedMB.end())
		{
			continue;
		}
		double scoreImprovement;
		double mbScore;
		Potential* dPot=NULL;
		string edgeKey(u->getName());
		edgeKey.append("\t");
		edgeKey.append(mIter->first.c_str());
		getNewPLLScore(u,varSet[nID],edgeKey,mbScore,scoreImprovement,&dPot);
		if(dPot!=NULL)
		{
			delete dPot;
		}
		moduleScore=moduleScore+scoreImprovement;
		neighborCnt++;
	}
	double score=moduleScore/neighborCnt;
	return score;
}


int
MetaLearner::getNewPLLScore(Variable* u, Variable* v, string& edgeKey, double& mbScore, double& scoreImprovement, Potential** newdPot)
{ //L simulate adding edge (u, v), rebuild local conditional model, and compute the prior-side score delta
	VSET& varSet=varManager->getVariableSet();
	bool dPotDel=true;
	scoreImprovement=0;
	double currPrior=varNeighborhoodPrior[v->getID()];
	double plus=0;
	double minus=0;

	SlimFactor* dFactor=factorGraph->getFactorAt(v->getID());

	// If u is already a parent of v, then we dont need to remove it at the end of this function.
	if(dFactor->mergedMB.find(u->getID())!=dFactor->mergedMB.end())
	{
		dPotDel=false;
	}
	dFactor->mergedMB[u->getID()]=0;

	Potential *dPot=new Potential;
	dPot->setAssocVariable(varSet[dFactor->fId],Potential::FACTOR);
	for(INTINTMAP_ITER mIter=dFactor->mergedMB.begin();mIter!=dFactor->mergedMB.end();mIter++) //L for each current parent of v (including new reg u), add to the potential
	{
		Variable* aVar=varSet[mIter->first];
		dPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
		double eprior=getEdgePrior(mIter->first,v->getID());
		double moduleContrib=getModuleContribLogistic((string&)v->getName(),(string&)aVar->getName()); //L adding the module prior to the edge priors. I change from u to aVar since we want to consider the contribution of all current parents of v, not just the new reg u
		double edgeProb=1/(1+exp(-1*(eprior+moduleContrib)));
		double edgeProbOld=1/(1+exp(-1*(eprior)));
		minus=minus+log(1-edgeProbOld);
		plus=plus+log(edgeProb);
	}
	dPot->potZeroInit();
	dPot->setCondBias(dFactor->potFunc->getCondBias());
	dPot->setCondVariance(dFactor->potFunc->getCondVariance());
	dPot->setCondWeight(dFactor->potFunc->getCondWeight());

	currPrior=currPrior+plus-minus;
	double newPLL_d=0;
	*newdPot=dPot;
	potManager->populatePotential(*newdPot);
	(*newdPot)->initMBCovMean();

	if((dPot->getCondVariance()<0) || (isnan(dPot->getCondVariance())) || (isinf(dPot->getCondVariance())))
	{
		scoreImprovement=-1;
>>>>>>> Stashed changes
	}

	INTINTMAP* tSet = &evidenceManager->getTrainingSet();
	int datasize = tSet->size();

	double currPrior = varNeighborhoodPrior[factorID] + plus - minus;
	double condLL = potManager->computeLL(factorID, parentIDs, datasize, newdPot);

	double varCnt = (double)parentIDs.size() + 1;
	double paramCnt = 2 * varCnt + varCnt * (varCnt - 1) / 2;
	double complexityPrior = -lambda * paramCnt * log(datasize);

	mbScore = condLL + complexityPrior + currPrior;
	scoreImprovement = mbScore - (*currPLL)[factorID];
}

double
MetaLearner::getInitPLLScore(int vId)
{
	SlimFactor* sFactor=factorGraph->getFactorAt(vId);
	Potential* sPot=sFactor->potFunc;

	double pll=0; 

	INTINTMAP* tSet=&evidenceManager->getTrainingSet();
	for(INTINTMAP_ITER eIter=tSet->begin();eIter!=tSet->end();eIter++)
	{
		EMAP* evidMap=evidenceManager->getEvidenceAt(eIter->first);
		double pval=sPot->evaluateProbabilityDensity(evidMap);
		if(isnan(pval))
		{
			cout <<"Pval is nan for datapoint " << eIter->first << endl;
		}
		if(pval<1e-50)
		{
			pval=1e-50;
		}
		pll += log(pval);
	}

	// The initial graph has no edges, meaning is variable is univariate
	// gaussian, with just 2 params (mean, variance).
	double complexityPrior = lambda * 2 * log(tSet->size());
	pll -= complexityPrior;
	return pll;
}

double 
MetaLearner::getEdgePrior(int tfID, int targetID)
{
	INTDBLMAP* regPriors=NULL;
	double prior; //L removed initializing this to beta1
	double fwt = 0;
	for (map<string,map<int,INTDBLMAP*>*>::iterator pItr=edgepriormap.begin(); pItr!=edgepriormap.end(); pItr++) //L for each prior network
	{
		double eweight=0;
		double gbeta = 0;
		map<int,INTDBLMAP*>* edgeprior = pItr->second;
		if(edgeprior->find(targetID)!=edgeprior->end()) //L if the target gene is in the prior
		{
			regPriors=(*edgeprior)[targetID]; //L get all regualtors of the gene in the prior
			if(regPriors->find(tfID)!=regPriors->end()) //L if the reguator is a regulator of the gene in the prior
			{
				eweight=(*regPriors)[tfID]; 
				gbeta = betamap[pItr->first]; //L get the weight of this prior network
				fwt = fwt + gbeta*eweight; 
			}
		}
	}
	prior=beta1+fwt; //L beta1 is the sparsity param, default = -5
	//if(prior<1e-6)
	//{
	//	prior=1e-6;
	//}
	//if(prior==1)
	//{
	//	prior=1-1e-6;
	//}
	return prior;
}

void
MetaLearner::makeMove(MetaMove* nextMove, int currIteration)
{
	VSET& varSet = varManager->getVariableSet();
	Variable* u = varSet[nextMove->getSrcVertex()];
	Variable* v = varSet[nextMove->getTargetVertex()];

	string edgeKey;
	edgeKey.append(u->getName().c_str());
	edgeKey.append("\t");
	edgeKey.append(v->getName().c_str());

	SlimFactor* dFactor = factorGraph->getFactorAt(nextMove->getTargetVertex());

	// Clean up the old potential
	delete dFactor->potFunc;

	// Add the new parent and update the potential
	dFactor->mergedMB[nextMove->getSrcVertex()] = 0;
	dFactor->potFunc = nextMove->getDestPot();

	// Update the current score for this factor
	(*currPLL)[dFactor->fId] = nextMove->getTargetMBScore();

	int mID = geneModuleID[v->getName()];

	// Get or create an indegree map for this module
	map<string,int>* currIndegree=NULL;
	if(moduleIndegree.find(mID) == moduleIndegree.end())
	{
		currIndegree = new map<string,int>;
		moduleIndegree[mID] = currIndegree;
	}
	else
	{
		currIndegree = moduleIndegree[mID];
	}

	// Increment the count of edges from u to the module of v
	if(currIndegree->find(u->getName()) == currIndegree->end())
	{
		(*currIndegree)[u->getName()] = 1;
	}
	else
	{	
		(*currIndegree)[u->getName()] += 1;
	}

	// Increment the count of edges from u
	if(regulatorModuleOutdegree.find(u->getName()) == regulatorModuleOutdegree.end())
	{
		regulatorModuleOutdegree[u->getName()] = 1;
	}
	else
	{
		regulatorModuleOutdegree[u->getName()] += 1;
	}

	edgeMap[edgeKey] = 1;

	variableStatus[v->getName()] = currIteration;
}

int
MetaLearner::dumpAllGraphs(int maxNumRegs,int foldid,int iter)
{
	VSET& varSet=varManager->getVariableSet();
	char aFName[1024];
	sprintf(aFName,"%s/prediction_k%d.txt",foldoutDirName,maxNumRegs+1);
	ofstream oFile(aFName);
	factorGraph->dumpVarMB(oFile,varSet);
	oFile.close();
	return 0;
}

int
MetaLearner::initPhysicalDegree()
{
	for(map<int,map<string,int>*>::iterator mIter=moduleGeneSet.begin();mIter!=moduleGeneSet.end();mIter++)
	{
		map<string,int>* indegree=NULL;
		map<string,map<string,int>*> innet;
		map<string,int>* geneSet=mIter->second;
		for(map<string,map<string,map<string,double>*>*>::iterator gpIter=priorgraphmap.begin();gpIter!=priorgraphmap.end();gpIter++)
		{
			map<string,int>enrichedTFs;
			map<string,map<string,double>*>* priorgraph = gpIter->second;
			getEnrichedTFs(enrichedTFs,geneSet,*priorgraph);
			for(map<string,int>::iterator tfIter=enrichedTFs.begin();tfIter!=enrichedTFs.end();tfIter++)
			{
				map<string,double>* motiftgts=(*priorgraph)[tfIter->first];
				map<string,int>* ttgts;
				if (innet.find(tfIter->first) == innet.end())
				{
					ttgts = new map<string,int>();
					innet[tfIter->first] = ttgts;
				}
				else
				{
					ttgts = innet[tfIter->first];
				}
				for(map<string,double>::iterator gIter=motiftgts->begin();gIter!=motiftgts->end();gIter++)
				{
					if(geneSet->find(gIter->first)==geneSet->end())
					{
						continue;
					}
					(*ttgts)[gIter->first] = 0;
				}
			}
		}
		for (map<string,map<string,int>*>::iterator tItr=innet.begin();tItr!=innet.end();tItr++)
		{
			string tf = tItr->first;
			map<string,int>* ttgts = tItr->second;
			if (ttgts->size() == 0)
			{
				continue;
			}
			if (indegree == NULL)
			{
				indegree=new map<string,int>;
			}
			(*indegree)[tf] = ttgts->size();
		}
		
		if(indegree!=NULL)
		{
			moduleIndegree[mIter->first]=indegree;
			cout <<"Physical_indegree for module " << mIter->first << endl;
			for(map<string,int>::iterator dIter=indegree->begin();dIter!=indegree->end();dIter++)
			{
				cout << dIter->first <<"\t" << dIter->second << endl;
				if(regulatorModuleOutdegree.find(dIter->first)==regulatorModuleOutdegree.end())
				{
					regulatorModuleOutdegree[dIter->first]=dIter->second;
				}
				else
				{
					regulatorModuleOutdegree[dIter->first]=regulatorModuleOutdegree[dIter->first]+dIter->second;
				}
			}
		}
	}
	return 0;
}

int
MetaLearner::getEnrichedTFs(map<string,int>& tfSet,map<string,int>* genes,map<string,map<string,double>*>& edgeSet)
{
	VSET& varSet=varManager->getVariableSet();
	int total=varSet.size();
	int k=genes->size();
	HyperGeomPval hgp;
	for(map<string,map<string,double>*>::iterator fIter=edgeSet.begin();fIter!=edgeSet.end();fIter++)
	{
		int vID=varManager->getVarID(fIter->first.c_str());
		if(vID<0)
		{
			continue;
		}
		map<string,double>* tgtSet=fIter->second;
		//int n=tgtSet->size();
		int n=0;
		int hit=0;
		for(map<string,double>::iterator gIter=tgtSet->begin();gIter!=tgtSet->end();gIter++)
		//for(map<string,int>::iterator gIter=genes->begin();gIter!=genes->end();gIter++)
		{
			
			int vID=varManager->getVarID(gIter->first.c_str());
                        if(vID<0)
                        {
                                continue;
                        }
                        n++;
			//if(tgtSet->find(gIter->first)==tgtSet->end())
                	if(genes->find(gIter->first)==genes->end())
			{
				continue;
			}
			hit++;
		}
		double enpval=hgp.getOverRepPval(k,hit,n,total-n);
		if(enpval<0.05 && hit>4)
		//if(hit>0)
		{
			tfSet[fIter->first]=hit;
		}
	}
	return 0;
}

double
MetaLearner::getModuleContribLogistic(string& tgtName, string& tfName)
{
	if(geneModuleID.find(tgtName)==geneModuleID.end())
	{
		return 0;
	}

	int moduleID=geneModuleID[tgtName];
	if(moduleIndegree.find(moduleID)==moduleIndegree.end())
	{
		return 0;
	}

	map<string,int>* moddegree=moduleIndegree[moduleID];
	if(moddegree->find(tfName)==moddegree->end())
	{
		return 0;
	}

	int degree=(*moddegree)[tfName];

	int regDegree=0;
	if(regulatorModuleOutdegree.find(tfName)!=regulatorModuleOutdegree.end())
	{
		regDegree=regulatorModuleOutdegree[tfName];
	}

	double contrib=((double) degree)/((double) regDegree);
	return beta_motif*contrib;
}

//To redefine the modules we will start with the original set of modules 
//For each original module, find for every gene its pairwise similarity to every other
//gene. merge two nodes that have the greatest pairwise similarity. replace by the merged
//regulatory program. recompute similarity of all nodes to this merged node. repeat with
//finding the next most similar pair of nodes.

int
MetaLearner::redefineModules()
{
	INTINTMAP& tSet=evidenceManager->getTrainingSet();

	if (correlationDistances == nullptr)
	{
		initCorrelationDistances();
	}

	map<string,int> genesWithNoNeighbors;

	// Create a node for each member of each module
	for(map<int,map<string,int>*>::iterator gIter=moduleGeneSet.begin();gIter!=moduleGeneSet.end();gIter++)
	{
		map<string,int>* moduleMembers=gIter->second;
		for(map<string,int>::iterator mIter=moduleMembers->begin();mIter!=moduleMembers->end();mIter++)
		{
			int mID=varManager->getVarID(mIter->first.c_str());
			if(mID<0)
			{
				continue;
			}
			SlimFactor* mFactor=factorGraph->getFactorAt(mID);

			// If a gene has no neighbors, we dont include it in the clustering algorithm.
			INTINTMAP& mbvars1=mFactor->mergedMB;
			if(mbvars1.size()==0)
			{
				genesWithNoNeighbors[mIter->first]=0;
				continue;
			}

			// Create a node for this gene
			HierarchicalClusterNode* node = hc.getNode(mIter->first);
			if (node == nullptr)
			{
				node = new HierarchicalClusterNode;
				node->nodeName.append(mIter->first);
				node->varID = mID;
				hc.addNode(node);

				// Add expression data on the new node
				for(INTINTMAP_ITER eIter=tSet.begin();eIter!=tSet.end();eIter++)
				{
					EMAP* evidMap=evidenceManager->getEvidenceAt(eIter->first);
					Evidence* evid=(*evidMap)[mID];
					double v=evid->getEvidVal();
					node->expr.push_back(v);
				}
			}

			// Add weights for incoming edges onto the node
			INTDBLMAP& regWts=mFactor->potFunc->getWeights();
			for(INTDBLMAP_ITER bIter=regWts.begin();bIter!=regWts.end();bIter++)
			{
				node->attrib[bIter->first]=bIter->second;
			}
		}
	}

	// Perform the new clustering
	map<int,map<string,int>*> newModules;
	hc.cluster(newModules, clusterThreshold, correlationDistances);

	// Clear out any data representing the old module assignments
	moduleGeneSet.clear();
	geneModuleID.clear();
	regulatorModuleOutdegree.clear();
	for(map<int,map<string,int>*>::iterator mIter=moduleIndegree.begin();mIter!=moduleIndegree.end();mIter++)
	{
		mIter->second->clear();
		delete mIter->second;
	}
	moduleIndegree.clear();

	char moduleFName[1024];
	sprintf(moduleFName,"%s/fold%d/modules.txt",outputDirName,currFold);
	ofstream modFile(moduleFName);

	// Read in the new module assignments
	int largestModuleID=0;
	VSET& varSet=varManager->getVariableSet();
	for(map<int,map<string,int>*>::iterator mIter=newModules.begin();mIter!=newModules.end();mIter++)
	{
		moduleGeneSet[mIter->first]=mIter->second;
		map<string,int>* geneSet=mIter->second;
		map<string,int>* indegree=new map<string,int>;
		for(map<string,int>::iterator gIter=geneSet->begin();gIter!=geneSet->end();gIter++)
		{
			modFile << gIter->first <<"\t" << mIter->first << endl;
			geneModuleID[gIter->first]=mIter->first;
			int mID=varManager->getVarID(gIter->first.c_str());
			SlimFactor* mFactor=factorGraph->getFactorAt(mID);
			INTINTMAP& mbvars1=mFactor->mergedMB;

			for(INTINTMAP_ITER nIter=mbvars1.begin();nIter!=mbvars1.end();nIter++)
			{
				// Count incoming edges to this module per regulator
				Variable* var=varSet[nIter->first];
				if(indegree->find(var->getName())==indegree->end())
				{
					(*indegree)[var->getName()]=1;
				}
				else
				{
					(*indegree)[var->getName()]=(*indegree)[var->getName()]+1;
				}
				// Count outgoing edges from regulator to any module
				if(regulatorModuleOutdegree.find(var->getName())==regulatorModuleOutdegree.end())
				{
					regulatorModuleOutdegree[var->getName()]=1;
				}	
				else
				{
					regulatorModuleOutdegree[var->getName()]=regulatorModuleOutdegree[var->getName()]+1;
				}
			}
		}
		moduleIndegree[mIter->first]=indegree;
		largestModuleID=mIter->first;
	}
	modFile.close();

	// For any genes with no neighbors, create single gene modules
	for(map<string,int>::iterator gIter=genesWithNoNeighbors.begin();gIter!=genesWithNoNeighbors.end();gIter++)
	{
		largestModuleID++;
		map<string,int>* newmodule=new map<string,int>;
		(*newmodule)[gIter->first]=0;
		moduleGeneSet[largestModuleID]=newmodule;
		geneModuleID[gIter->first]=largestModuleID;
	}
	genesWithNoNeighbors.clear();

	return 0;
}

void
MetaLearner::initCorrelationDistances()
{
	INTINTMAP& samples = evidenceManager->getTrainingSet();
	VSET& varSet = varManager->getVariableSet();

	int varCount = varSet.size();
	int sampleCount = samples.size();

	vector<double> means(varCount, 0);

	for (INTINTMAP_ITER iter = samples.begin(); iter != samples.end(); iter++)
	{
		EMAP* evidMap = evidenceManager->getEvidenceAt(iter->first);
		for (int i = 0; i < varCount; i++)
		{
			Evidence* evid=(*evidMap)[i];
			means[i] += evid->getEvidVal();
		}
	}

	for (int i = 0; i < means.size(); i++) {
		means[i] /= sampleCount;
	}

	vector<double> ssd(varCount, 0);
	vector<vector<double>> deviations(varCount, vector<double>(sampleCount, 0));

	int sampleIndex = 0;
	for (INTINTMAP_ITER iter = samples.begin(); iter != samples.end(); iter++)
	{
		EMAP* evidMap = evidenceManager->getEvidenceAt(iter->first);
		for (int i = 0; i < varSet.size(); i++)
		{
			double deviation = (*evidMap)[i]->getEvidVal() - means[i];
			deviations[i][sampleIndex] = deviation;
			ssd[i] += deviation * deviation;
		}
		sampleIndex++;
	}

	correlationDistances = new Matrix(varCount, varCount);

	double threshold = sampleCount / 2.0;

	for (int i = 0; i < varCount; i++)
	{
		double xx = ssd[i];
		double* dev_i = deviations[i].data();

		for (int j = i; j < varCount; j++)
		{
			double* dev_j = deviations[j].data();
			double xy = 0;
			int oppRel = 0;

			for(int k = 0; k < sampleCount; k++)
			{
				double diff1 = dev_i[k];
				double diff2 = dev_j[k];
				double val = diff1 * diff2;
				xy += val;
				oppRel += (val < 0);
			}

			double yy = ssd[j];
			double cc = abs(xy) / sqrt(xx * yy);

			if(oppRel > threshold)
			{
				cc *= -1;
			}

			cc = 0.5 * (1 - cc);

			correlationDistances->setValue(cc, i, j);
			correlationDistances->setValue(cc, j, i);
		}
	}
}
