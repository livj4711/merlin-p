#include <iostream>
#include <cstring>
#include <math.h>

#include "CommonTypes.H"
#include "Error.H"
#include "Variable.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "SlimFactor.H"
#include "PotentialManager.H"

PotentialManager::PotentialManager()
{
	data=NULL;
	meanMat=NULL;
	covMat=NULL;
	ludecomp=NULL;
	perm=NULL;
	lambda=0;
	randomData=false;
}

PotentialManager::~PotentialManager()
{
	if(ludecomp!=NULL)
	{
		gsl_matrix_free(ludecomp);
	}
	if(perm!=NULL)
	{
		gsl_permutation_free(perm);
	}
	if(data!=NULL)
	{
		delete data;
	}
	if(meanMat!=NULL)
	{
		delete meanMat;
	}
	if(covMat!=NULL)
	{
		delete covMat;
	}
}

int 
PotentialManager::setEvidenceManager(EvidenceManager* aPtr)
{
	evMgr=aPtr;
	return 0;
}

int 
PotentialManager::setRestrictedNeighborSet(map<int,Variable*>& rSet)
{
	for(map<int,Variable*>::iterator vIter=rSet.begin();vIter!=rSet.end();vIter++)
	{
		restrictedNeighborSet[vIter->first]=vIter->second;
	}
	return 0;
}

int
PotentialManager::setLambda(double lVal)
{
	lambda=lVal;
	return 0;
}

int
PotentialManager::setRandom(bool flag)
{
	randomData=flag;
	return 0;
}

int
PotentialManager::init()
{
	initData(); //L initialize the expression matrix and the mean "vector" (mean of gene expression) but set cov matrix to all -1!
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}

int
PotentialManager::reset()
{
	if(ludecomp!=NULL)
	{
		gsl_matrix_free(ludecomp);
		ludecomp=NULL;
	}
	if(perm!=NULL)
	{
		gsl_permutation_free(perm);
		perm=NULL;
	}
	if(data!=NULL)
	{
		delete data;
		data=NULL;
	}
	if(meanMat!=NULL)
	{
		delete meanMat;
		meanMat=NULL;
	}
	if(covMat!=NULL)
	{
		delete covMat;
		covMat=NULL;
	}
	return 0;
}

int
PotentialManager::initData()
{
	INTINTMAP& trainEvidSet=evMgr->getTrainingSet();
	EMAP* evidMap=evMgr->getEvidenceAt(trainEvidSet.begin()->first); //L get the map of cells for the first training sample, just to get the number of genes
	int varCnt=evidMap->size(); //L number of genes in the data

	// data is the data matrix which will have the variable (gene) by sample (cell) information
	if(data==NULL) //L if the data matrix has not been initialized yet, initialize it with the number of genes and number of training samples
	{
		data=new Matrix(varCnt,trainEvidSet.size()); //L num_genes x num_training_cells matrix
		meanMat=new Matrix(varCnt,1); //L num_genes x 1 matrix of mean gene expression
		meanMat->setAllValues(0); //L initialize the mean matrix to all 0ss
		covMat=new Matrix(varCnt,varCnt); //L num_genes x num_genes covariance matrix
		covMat->setAllValues(-1); //L initialize the covariance matrix to all -1s
	}

	// Copy all the samples into the data matrix
	int sampleIndex = 0;
	for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		EMAP* evidMap=NULL;
		if(randomData)
		{
			evidMap=evMgr->getRandomEvidenceAt(eIter->first);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(eIter->first); //L get this cell's gene expression map from this PotentialManager's evMgr
		}
		for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++) //L for each gene in this cell
		{
			int vId=vIter->first;
			Evidence* evid=vIter->second;
			double val=evid->getEvidVal();
			data->setValue(val,vId,sampleIndex);
		}
		sampleIndex++;
	}

	// Done copying. Now we can go over the rows of data and get the means
	for(int i=0;i<varCnt;i++) //L for each gene (row) in the data matrix
	{
		double sampleSum=0;
		for(int j=0;j<data->getColCnt();j++) //L for each cell (column) in the data matrix
		{
			sampleSum += data->getValue(i,j); //L sum up the expression values for this gene across all the training cells
		}
		double sampleSize=(double) data->getColCnt(); //L number of training cells
		meanMat->setValue(sampleSum/sampleSize,i,0); //L set the mean expression for this gene in the meanMat
	}

	return 0;
}

int
PotentialManager::estimateCovariance(int uId, int vId)
{
	double ssd=0; //L sum of squared differences, this is what we will use to calculate the covariance for this pair of genes
	for(int i=0;i<data->getColCnt();i++)
	{
		double vval=data->getValue(vId,i);
		double vmean=meanMat->getValue(vId,0);
		double uval=data->getValue(uId,i);
		double umean=meanMat->getValue(uId,0);
		ssd += (vval-vmean)*(uval-umean);
	}
	if(uId==vId)
	{
		ssd += 0.001; //L add a small value to avoid division by zero
	}
	double var = ssd/((double)(data->getColCnt()-1));
	covMat->setValue(var,uId,vId);
	covMat->setValue(var,vId,uId);
	return 0;
}

Error::ErrorCode
PotentialManager::populatePotentialsSlimFactors(map<int,SlimFactor*>& factorSet,VSET& varSet)
{
	//The set of flags to keep status of the potentials that have been calculated
	map<int,bool> doneFlag;
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++) //L for each factor in the factorSet. factorSet contains all the SlimFactors (one per gene)
	{
		doneFlag[fIter->first]=false; //L initialize doneFlag to false
	}
	int popFId=0;
	for(map<int,SlimFactor*>::reverse_iterator rIter=factorSet.rbegin();rIter!=factorSet.rend();rIter++) //L Iterates over all slim factors in reverse key order
	{
		//If we have computed the potential for this flag move one
		if(doneFlag[rIter->first])
		{
			popFId++;
			continue;
		}
		SlimFactor* sFactor=rIter->second;

		//Otherwise create the potential
		Potential* aPotFunc=new Potential; //L a potential FOR EACH SLIMFACTOR (GENE)
		for(int j=0;j<sFactor->vCnt;j++) //L for each gene in this SlimFactor (should just be one)
		{
			Variable* aVar=varSet[sFactor->vIds[j]]; //L set aVar to be the gene ID of this gene (variable)

			aPotFunc->setAssocVariable(aVar,Potential::FACTOR); //L there is only one gene in this factor, so we set the association type to FACTOR (rather than MAKOV_BNKT)
			// if(j==sFactor->vCnt-1) 
			// {
			// 	aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
			// }
			// else
			// {
			// 	aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
			// }
		}
		aPotFunc->potZeroInit(); //L initialize mean vector and covariance matrix for this one gene (and its markov blanket but that's empty right now) to 0 
		populatePotential(aPotFunc); //L init this potential with the actual mean and covariance value (covar = variance, since theres only one gene in this potential) based on the entire training data, and also the normalzing factor for the gaussian pdf
		aPotFunc->calculateJointEntropy(); //L get the joint entropy for this potential based on the covariance
		sFactor->jointEntropy=aPotFunc->getJointEntropy(); //L set this slimFactor's joint entropy
		if(sFactor->jointEntropy<0)
		{
		//	sFactor->jointEntropy=0;
		//	cout <<"Negative entropy for " << sFactor->fId << endl;
		}
		doneFlag[rIter->first]=true;
		delete aPotFunc;
		popFId++;
		if(popFId%100000==0) //L print progress for every 100k factors (genes)
		{
			cout <<"Done with " << factorSet.size()-popFId << " factors " << endl;
		}
	}
	return Error::SUCCESS;
}

int
PotentialManager::populatePotential(Potential* aPot)
{
	VSET& potVars=aPot->getAssocVariables(); //L the gene and all its regulators (of the markov blanket)
	for(VSET_ITER vIter=potVars.begin();vIter!=potVars.end(); vIter++) //L for each gene or reg v in this potential
	{
		double mean=meanMat->getValue(vIter->first,0); 
		aPot->updateMean(vIter->first,mean); //L set the mean value for this gene to its mean of the entire training data

		for(VSET_ITER uIter=vIter;uIter!=potVars.end();uIter++) //L for each other gene or reg u in this porential
		{
			double cval=covMat->getValue(vIter->first,uIter->first); //L get the covar for this pair from the entire training data
			if(cval==-1) //L not initialized yet
			{
				estimateCovariance(uIter->first,vIter->first); //L estimate the covariance for this pair based on the entire training data and set it in the covMat
				cval=covMat->getValue(vIter->first,uIter->first); //L extract that covariance
			}
			aPot->updateCovariance(vIter->first,uIter->first,cval); //L update that pairs overall covariance into this potential
			aPot->updateCovariance(uIter->first,vIter->first,cval); //L the symmetric value too
		}
	}
	aPot->makeValidJPD(ludecomp,perm);
	return 0;
}

double
PotentialManager::getLikelihood(SlimFactor* sFactor,VSET& varSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

double
PotentialManager::getLikelihood(SlimFactor* sFactor,VSET& varSet,map<int,int>& visitedVertices )
{
	double dll=0;
	cout <<"Not implemented" << endl;
	return dll;
}

int
PotentialManager::estimateConditionalPotential(SlimFactor* sFactor,VSET& varSet,Potential** pot, STRDBLMAP& counts)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int 
PotentialManager::estimateCanonicalPotential(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int 
PotentialManager::estimateCanonicalPotential_Abbeel(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int 
PotentialManager::estimateCanonicalPotential_Approximate(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int
PotentialManager::resetPotFuncs()
{
	for(map<int,Potential*>::iterator pIter=potFuncs.begin();pIter!=potFuncs.end();pIter++)
	{
		delete pIter->second;
	}
	potFuncs.clear();
	return 0;
}

int 
PotentialManager::estimateCanonicalPotential_Joint(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

Potential*
PotentialManager::getPotential(int fId)
{
	if(potFuncs.find(fId)==potFuncs.end())
	{
		return NULL;
	}
	return potFuncs[fId];
}

double 
PotentialManager::getSampleLikelihood(map<int,SlimFactor*>& factorSet, VSET& varSet, INTINTMAP* sample)
{
	double sampleLL=0;
	cout <<"Not implemented " <<endl;
	return sampleLL;
}

int 
PotentialManager::getVariableSample(INTINTMAP& jointConf,VSET& varSet,int vId,SlimFactor* sFactor, gsl_rng* r)
{
	Potential* pot=NULL;
	if(potFuncs.find(vId)==potFuncs.end())
	{
		pot=new Potential;
		STRDBLMAP counts;
		estimateConditionalPotential(sFactor,varSet,&pot,counts);
		potFuncs[vId]=pot;
	}
	else
	{
		pot=potFuncs[vId];
	}
	//int sample=pot->generateSample(jointConf,vId,r);
	int sample=-1;
	return sample;
}

int
PotentialManager::clearJointEntropies()
{
	jointEntropies.clear();
	return 0;
}
