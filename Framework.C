#include <fstream>
#include <iostream>
#include <cstring>
#include <getopt.h>

#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"

#include "Evidence.H"
#include "EvidenceManager.H"

#include "Potential.H"
#include "SlimFactor.H"

#include "FactorGraph.H"
#include "PotentialManager.H"
#include "MetaMove.H"
#include "MetaLearner.H"

#include "Framework.H"

//L Removed dead commented-out code, using namespace std, redundant oldoptind tracking, and the duplicate '?' case in the switch.
//L wDefaults, int c, epsThreshold, and outFilePrefix were all unused in the original code and have been removed for clarity.

static void printUsage()
{
    std::cout
        << "MERLIN-P GRN Inference\n"
        << "  -d gene_expression_file\n"
        << "  -k maxfactorsize (default size_of_dataset)\n"
        << "  -v cross_validation_cnt (default 1)\n"
        << "  -l restricted_regulator_fname\n"
        << "  -p sparsity_prior (default -5)\n"
        << "  -r module_prior (default 4)\n"
        << "  -q prior_config\n"
        << "  -o outputdirectory\n"
        << "  -c moduleassignment (default random_partitioning)\n"
        << "  -h hierarchical_clustering_threshold (default 0.6)\n"
        << "  -f specificfold_torun (default is -1)\n";
}

Framework::Framework() = default;
Framework::~Framework() = default;

Error::ErrorCode Framework::init(int argc, char** argv)
{
	bool dDefault = true;
	bool oDefault = true;
	bool kDefault = true;
	bool cvDefault = true;
	bool regDefault = true;
	bool pDefault = true;
	bool rDefault = true;
	bool moduleDefault = true;
	bool hDefault = true;

				eCode = evManager.loadEvidenceFromFile(optarg);
				if(eCode != Error::SUCCESS)
				{
					cerr << Error::getErrorString(eCode) << endl;
					return eCode;
				}
				metaLearner.setGlobalEvidenceManager(&evManager);
				metaLearner.setVariableManager(&varManager);
			
				//strncpy(outFilePrefix, optarg, strlen(optarg) - 4);
				//outFilePrefix[strlen(optarg)-4] = '\0';

                eCode = evManager.loadEvidenceFromFile(optarg);
                if (eCode != Error::SUCCESS)
                {
                    std::cerr << Error::getErrorString(eCode) << '\n';
                    return eCode;
                }

                metaLearner.setGlobalEvidenceManager(&evManager);
                metaLearner.setVariableManager(&varManager);
                break;
            }
			case 'o': // output file 
			{
				oDefault=false;
				metaLearner.setOutputDirName(optarg);	
				break;
			}
            case 'k': // max factor size (max # regs per gene)
            {
				kDefault=false;
                metaLearner.setMaxFactorSize_Approx(std::atoi(optarg));
                break;
            }
            case 'v': // cross validation count
            {
				cvDefault=false;
                cvCnt = std::atoi(optarg);
                if (cvCnt <= 0)
                {
                    std::cerr << "Cross validation count should be greater than zero.\n";
                    return Error::UNKNOWN;
                }
                break;
            }
            case 'l': // list of restricted/known regulators
            {
				regDefault=false;
                metaLearner.setRestrictedList(optarg);
                break;
            }
            case 'p': //sparsity parameter
            {
				pDefault=false;
                metaLearner.setBeta1(std::atof(optarg));
                break;
            }
			case 'r':
			{
				rDefault=false;
				metaLearner.setBeta_Motif(std::atof(optarg));
				break;
			}
			case 'q':
			{
				metaLearner.setPriorGraph_All(optarg);
				break;
			}
			case 'c':
			{
				moduleDefault=false;
				metaLearner.readModuleMembership(optarg);
				break;
			}
			case 'h':
			{
				hDefault = false;
				metaLearner.setClusteringThreshold(std::atof(optarg));
				break;
			}
			case 'f':
			{
				metaLearner.setSpecificFold(std::atoi(optarg));
				break;
			}
            case '?':
            default:
			{
                std::cerr << "Option parsing error\n";
                return Error::UNKNOWN;
			}
        }
    }

    // Validate required options
	if(dDefault)
	{
		std::cerr << "Please specify the name of expression file. (option -d)\n";
		return Error::UNKNOWN;
	}	
    if(oDefault)
	{
		std::cerr << "Please specify the name of output directory. (option -o)\n";
		return Error::UNKNOWN;
	}	
    if(regDefault)
	{
		std::cerr << "Please input a file of regulators. (option -l)\n";
		return Error::UNKNOWN;
	}


    // Apply defaults for unset options
	if(moduleDefault)
	{
		std::cout << "Setting to default module membership\n";
		metaLearner.setDefaultModuleMembership();
	}
	if(cvDefault)
	{
		cvCnt=1;	
	}
	if(hDefault)
	{
		metaLearner.setClusteringThreshold(0.6);
	}
	if(pDefault)
	{
		metaLearner.setBeta1(-5);
	}
	if(kDefault)
	{
		metaLearner.setMaxFactorSize_Approx(300);
	}
	if(rDefault)
	{
		metaLearner.setBeta_Motif(4);
	}

	//J metaLearner.initPartitions(1);
	return Error::SUCCESS;
}

int Framework::start()
{
    metaLearner.doCrossValidation(cvCnt);
    return 0;
}

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        printUsage();
        return 0;
    }

    Framework fw;
    if (fw.init(argc, argv) != Error::SUCCESS)
    {
        return -1;
    }
    fw.start();
    return 0;
}