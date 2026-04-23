#include <fstream>
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"


VariableManager::VariableManager()
{
}

VariableManager::~VariableManager()
{
}

//Reads the schema of the variables

Error::ErrorCode
VariableManager::readVariables(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[400000];
	//L int nodeCnt=0; //unused

	if(inFile.good()) 
	{
		inFile.getline(buffer,400000); //L read the first line, the header, containing the gene names

		if(strlen(buffer)<=0)
		{
			cout <<"Error: gene expression header is empty" << endl; //L change to be more informative
			return Error::VARSCHEMA_ERR;
		}

		char* tok=strtok(buffer,"\t");  //L split the header into tokens by tabs
		int tokCnt=0;  //L will become variable id, so gene1 has ID 0, gene2 has ID

		while(tok!=NULL) //L loop through all gene names
		{
			Variable* var=new Variable;
			var->setID(tokCnt);
			var->setName(tok);
			variableSet[tokCnt]=var;
			
			string varKey(tok);
			varNameIDMap[varKey]=tokCnt;  //L create a name --> ID lookup
			tokCnt++;
			tok=strtok(NULL,"\t");
		}
	}

	inFile.close();

	cout <<"Number of genes read: " << variableSet.size() << endl; //L say genes instead of variables

	return Error::SUCCESS;
}

int
VariableManager::getVarID(const char* varName)
{
	string varKey(varName);
	if(varNameIDMap.find(varKey)==varNameIDMap.end())
	{
		return -1;
	}
	int vId=varNameIDMap[varKey];
	return vId;
}

bool 
VariableManager::isValid(int varID,int varVal)
{
	Variable* rVar=variableSet[varID];
	return rVar->isValidValue(varVal);
}

map<int,Variable*>& 
VariableManager::getVariableSet()
{
	return variableSet;
}


Variable* 
VariableManager::getVariableAt(int vId)
{
	if(variableSet.find(vId)==variableSet.end())
	{
		cout << "Illegal variable id " << vId << endl;
		return NULL;
	}
	return variableSet[vId];
}
