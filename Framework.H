/*The class which brings it all together*/
#ifndef _FRAMEWORK
#define _FRAMEWORK

//All the managers need forward declarations
class MetaLearner;
class VariableManager;
class EvidenceManager;
class Framework
{
	public:
		Framework();
		~Framework();
		Error::ErrorCode init(int, char**);
		int start();

	private:
		MetaLearner metaLearner;
		VariableManager varManager;
		EvidenceManager evManager;
		EvidenceManager holdoutEvManager;
		//L double epsThreshold; // no longer used in current pipeline
		//L char actionStr[256]; //no longer used in current pipeline
		int cvCnt;
		//L map<int,int> idVidMap; //no longer used in current pipeline
};
#endif
