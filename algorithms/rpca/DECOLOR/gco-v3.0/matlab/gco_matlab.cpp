#include <mex.h>
#include <stdio.h>
#include <map>
#include <string>
#include "../GCoptimization.h"

#if !defined(MX_API_VER) || MX_API_VER < 0x07030000
//commented for MATLAB R2016a
//typedef int mwSize;
//typedef int mwIndex;
#endif

// MATLAB 2014a: mxCreateReference is NO LONGER AVAILABLE :(
//extern "C" mxArray *mxCreateReference(const mxArray*); // undocumented mex function

#define GCO_EXPORT(func) \
	extern "C" void func(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]); \
	FuncRegistry::Entry regentry_##func(#func,func); \
	extern "C" void func(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

#define MATLAB_ASSERT(expr,msg) if (!(expr)) { throw MatlabError(msg); }
#define MATLAB_ASSERT_ARGCOUNT(nout, nin) \
	MATLAB_ASSERT(nlhs >= nout, "Not enough output arguments, expected " #nout); \
	MATLAB_ASSERT(nlhs <= nout, "Too many output arguments, expected " #nout); \
	MATLAB_ASSERT(nrhs >= nin,  "Not enough input arguments, expected " #nin); \
	MATLAB_ASSERT(nrhs <= nin,  "Too many input arguments, expected " #nin);
#define MATLAB_ASSERT_INTYPE(arg, type) \
	MATLAB_ASSERT(mxGetClassID(prhs[arg-1]) == c##type##ClassID, "Expected argument " #arg " to be of type " #type);
#define MATLAB_ASSERT_INCLASS(arg, classid) \
	MATLAB_ASSERT(mxGetClassID(prhs[arg-1]) == classid, "Expected argument " #arg " to be of type " #classid);
#define MATLAB_ASSERT_HANDLE(arg) \
	MATLAB_ASSERT(mxGetClassID(prhs[arg-1]) == mxINT32_CLASS, "Expected argument " #arg " to be a valid GCO instance handle");

struct MatlabError {
	MatlabError(const char* msg): msg(msg) { }
	const char* msg;
};

struct FuncRegistry {
	typedef void (*Func)(int, mxArray*[], int, const mxArray*[]);
	typedef std::map<std::string,Func> LookupTable;
	static LookupTable sLookup;
	struct Entry {
		Entry(const char* name, Func ptr) { sLookup[name] = ptr; }
	};
};
FuncRegistry::LookupTable FuncRegistry::sLookup;


struct GCInstanceInfo {
	GCInstanceInfo(): gco(0), grid(false), dc(0), sc(0) { }
	~GCInstanceInfo() {
		if (sc) mxDestroyArray(sc);
		if (dc) mxDestroyArray(dc);
		if (gco) delete gco;
	}
	GCoptimization* gco;
	bool grid;
	mxArray* dc;
	mxArray* sc;
private:
};

template <typename ctype> struct ctype2mx { };
template <> struct ctype2mx<bool>            { static const mxClassID classid = mxLOGICAL_CLASS;   static const char* classname() { return "logical"; } };
template <> struct ctype2mx<float>           { static const mxClassID classid = mxSINGLE_CLASS;    static const char* classname() { return "single"; }  };
template <> struct ctype2mx<double>          { static const mxClassID classid = mxDOUBLE_CLASS;    static const char* classname() { return "double"; }  };
template <> struct ctype2mx<char>            { static const mxClassID classid = mxINT8_CLASS;      static const char* classname() { return "int8"; }    };
template <> struct ctype2mx<unsigned char>   { static const mxClassID classid = mxUINT8_CLASS;     static const char* classname() { return "uint8"; }   };
template <> struct ctype2mx<short>           { static const mxClassID classid = mxINT16_CLASS;     static const char* classname() { return "int16"; }   };
template <> struct ctype2mx<unsigned short>  { static const mxClassID classid = mxUINT16_CLASS;    static const char* classname() { return "uint16"; }  };
template <> struct ctype2mx<int>             { static const mxClassID classid = mxINT32_CLASS;     static const char* classname() { return "int32"; }   };
template <> struct ctype2mx<unsigned int>    { static const mxClassID classid = mxUINT32_CLASS;    static const char* classname() { return "uint32"; }  };
template <> struct ctype2mx<long long>       { static const mxClassID classid = mxINT64_CLASS;     static const char* classname() { return "int64"; }   };
template <> struct ctype2mx<unsigned long long> { static const mxClassID classid = mxUINT64_CLASS; static const char* classname() { return "uint64"; }  };


typedef GCoptimization::LabelID LabelID;
typedef GCoptimization::SiteID SiteID;
typedef GCoptimization::EnergyType EnergyType;
typedef GCoptimization::EnergyTermType EnergyTermType;
mxClassID cLabelClassID      = ctype2mx<LabelID>::classid;
mxClassID cSiteClassID       = ctype2mx<SiteID>::classid;
mxClassID cEnergyTermClassID = ctype2mx<EnergyTermType>::classid;
mxClassID cEnergyClassID     = ctype2mx<EnergyType>::classid;

typedef std::map<int,GCInstanceInfo> GCInstanceMap;

static int gNextInstanceID = 10001; // some start id for the first GC object
static GCInstanceMap gInstanceMap;

GCInstanceMap::mapped_type& sGetGCInstance(int id) {
	GCInstanceMap::iterator it = gInstanceMap.find(id);
	MATLAB_ASSERT(it != gInstanceMap.end(), "Invalid handle; no such GCoptimization object");
	return it->second;
}

/*
// Not currently used because it was too slow to shuffle per-edge smoothcosts between MATLAB and C++
void sBatchSmoothCostFn(SiteID size, const SiteID sites[][2], const LabelID labels[][2], EnergyTermType* energies, void *smoothFn)
{
	mxArray* sitesArr = mxCreateNumericMatrix(2, size, cSiteClassID, mxREAL);
	mxArray* labelsArr = mxCreateNumericMatrix(2, size, cLabelClassID, mxREAL);
	LabelID* sitesArrData = (LabelID*)mxGetData(sitesArr);
	LabelID* labelsArrData = (LabelID*)mxGetData(labelsArr);

	// Convert indices from 0..N-1 to 1..N for Matlab callback
	for (SiteID i = 0; i < size*2; ++i)
		sitesArrData[i] = ((SiteID*)sites)[i]+1;
	for (SiteID i = 0; i < size*2; ++i)
		labelsArrData[i] = ((LabelID*)labels)[i]+1;
	

	// Call the user-specified smoothcost callback, which is either a string or a function_handle
	mxArray* plhs[1] = { 0 };
	mxArray* prhs[3] = { (mxArray*)smoothFn, sitesArr, labelsArr };
	mexCallMATLAB(1, plhs, 3, prhs, "feval");

	// Copy the result into the output 'energies' array
	mxArray* result = plhs[0];
	MATLAB_ASSERT(result, "Batch smoothcost function failed to return a result");
	MATLAB_ASSERT(mxGetClassID(result) == cEnergyTermClassID, "Batch smoothcost function must return numeric array of type GCoptimization::EnergyTermType");
	MATLAB_ASSERT(mxGetNumberOfElements(result) >= size, "Batch smoothcost function returned too few energy terms");
	MATLAB_ASSERT(mxGetNumberOfElements(result) <= size, "Batch smoothcost function returned too many energy terms");
	memcpy(energies, mxGetData(result), size*sizeof(EnergyTermType));

	// Clean up
	mxDestroyArray(sitesArr);
	mxDestroyArray(labelsArr);
	mxDestroyArray(result);
}
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs == 0 || !mxIsChar(prhs[0]))
		mexErrMsgTxt("Do not use gco_matlab() directly, instead use the GCO functions such as GCO_Create"); 
	mwSize nameLen = mxGetN(prhs[0])*sizeof(mxChar)+1;
	char funcName[512];
	mxGetString(prhs[0], funcName, nameLen); 
	FuncRegistry::LookupTable::const_iterator it = FuncRegistry::sLookup.find(funcName);
	if (it == FuncRegistry::sLookup.end())
		mexErrMsgTxt("Specified function does not exist within gco_matlab module"); 
	try {
		it->second(nlhs, plhs, nrhs-1, prhs+1);
	} catch (GCException err) {
		mexErrMsgTxt(err.message);
	} catch (MatlabError err) {
		mexErrMsgTxt(err.msg);
	}
}

GCO_EXPORT(gco_get_energyterm_class)
{
	// This function is provided so that the GCO_SetDataCost/SetSmoothCost/SetLabelCost
	// matlab functions can determine if they should automatically convert float input 
	// matrices to int32 before passing them
	MATLAB_ASSERT_ARGCOUNT(1,0);
	char* str = (char*)mxCalloc(32, sizeof(char));
	strcpy(str, ctype2mx<EnergyTermType>::classname());
	plhs[0] = mxCreateString(str);
}

GCO_EXPORT(gco_create_general)
{
	int instanceID = 0;
	try {
		MATLAB_ASSERT_ARGCOUNT(1,2);
		MATLAB_ASSERT_INTYPE(1,Site);
		MATLAB_ASSERT_INTYPE(2,Label);
		SiteID  numSites  = *(SiteID* )mxGetData(prhs[0]); MATLAB_ASSERT(numSites  >= 1, "Number of sites must be positive");
		LabelID numLabels = *(LabelID*)mxGetData(prhs[1]); MATLAB_ASSERT(numLabels >= 2, "Number of labels must be positive");
		instanceID = gNextInstanceID++;
		GCInstanceInfo& gcinstance = gInstanceMap[instanceID];
		gcinstance.gco = new GCoptimizationGeneralGraph(numSites, numLabels);
		gcinstance.grid = false;
		mwSize outSize = 1;
		plhs[0] = mxCreateNumericArray(1, &outSize, mxINT32_CLASS, mxREAL);
		*(int*)mxGetData(plhs[0]) = instanceID;
	} catch (MatlabError) {
		if (instanceID) 
			gInstanceMap.erase(instanceID);
		throw;
	}
}

GCO_EXPORT(gco_delete)
{
	MATLAB_ASSERT_HANDLE(1);
	const int* instanceIDs = (int*)mxGetData(prhs[0]);
	MATLAB_ASSERT(mxGetN(prhs[0]) == 1 || mxGetM(prhs[0]) == 1, "Input must be a scalar or a vector");
	mwIndex count = mxGetNumberOfElements(prhs[0]);
	for (mwIndex i = 0; i < count; ++i) {
		MATLAB_ASSERT(gInstanceMap.find(instanceIDs[i]) != gInstanceMap.end(), "Invalid handle (no such GCoptimization object)");
		gInstanceMap.erase(instanceIDs[i]);
	}
}

GCO_EXPORT(gco_listhandles)
{
	MATLAB_ASSERT_ARGCOUNT(1,0);
	mwSize outSize = (mwSize)gInstanceMap.size();
	plhs[0] = mxCreateNumericArray(1, &outSize, mxINT32_CLASS, mxREAL);
	int* instanceIDs = (int*)mxGetData(plhs[0]);
	for (GCInstanceMap::const_iterator i = gInstanceMap.begin(); i != gInstanceMap.end(); ++i)
		*(instanceIDs++) = i->first;
}


GCO_EXPORT(gco_setdatacost)
{
	MATLAB_ASSERT_ARGCOUNT(0,3);
	MATLAB_ASSERT_HANDLE(1);
	MATLAB_ASSERT_INTYPE(2,EnergyTerm);
	MATLAB_ASSERT_INTYPE(3,Label);
	GCInstanceInfo& gcinstance = sGetGCInstance(*(int*)mxGetData(prhs[0]));
	const mxArray* dc = prhs[1];
	LabelID label = *(LabelID*)mxGetData(prhs[2]);
	if (label == 0) {
		// Dense data costs
		MATLAB_ASSERT(mxGetN(dc) == gcinstance.gco->numSites() && mxGetM(dc) == gcinstance.gco->numLabels(),
					  "Numeric data cost must be NumLabels x NumSites in size");
        // Make a copy since MATLAB no longer supports incrementing reference 
        // count from mex extensions... I miss you mxCreateReference :(
        mxArray* dc_copy = mxDuplicateArray(dc);
        mexMakeArrayPersistent(dc_copy);
            
		gcinstance.gco->setDataCost((EnergyTermType*)mxGetData(dc_copy));
		if (gcinstance.dc)
			mxDestroyArray(gcinstance.dc);
		gcinstance.dc = dc_copy;
	} else {
		// Sparse data costs
		MATLAB_ASSERT(mxGetM(dc) == 2, "Sparse data cost must have 2 rows");
		MATLAB_ASSERT(sizeof(SiteID) == sizeof(EnergyTermType), 
			"Sparse data costs cannot be used because GCoptimization was compiled "
			"such that sizeof(SiteID) != sizeof(EnergyTermType)");
		MATLAB_ASSERT(cSiteClassID == cEnergyTermClassID, "Sparse data costs only supported when compiled for int32 energy terms");
		SiteID count = (SiteID)mxGetN(dc);
		GCoptimization::SparseDataCost* dcmem = (GCoptimization::SparseDataCost*)mxGetData(dc);
		try {
			for (LabelID i = 0; i < count; ++i)
				dcmem[i].site--; // from 1-based to 0-based index
			gcinstance.gco->setDataCost(label-1,dcmem,count);
			if (gcinstance.dc) {
				mxDestroyArray(gcinstance.dc); // the gco object forgets about the pointer, so we let MATLAB know
				gcinstance.dc = 0;
			}
			for (LabelID i = 0; i < count; ++i)
				dcmem[i].site++; // from 0-based to 1-based index
		} catch (...) {
			for (LabelID i = 0; i < count; ++i)
				dcmem[i].site++; // from 0-based to 1-based index
			throw;
		}
	}
}

GCO_EXPORT(gco_setsmoothcost)
{
	MATLAB_ASSERT_ARGCOUNT(0,2);
	MATLAB_ASSERT_HANDLE(1);
	MATLAB_ASSERT_INTYPE(2,EnergyTerm);
	GCInstanceInfo& gcinstance = sGetGCInstance(*(int*)mxGetData(prhs[0]));
	const mxArray* sc = prhs[1];
	// Smooth costs provided as numeric array, and is applied to all neighbouring variables.
	MATLAB_ASSERT(mxGetN(sc) == gcinstance.gco->numLabels() && mxGetM(sc) == gcinstance.gco->numLabels(),
	              "Numeric smooth cost must be NumLabels x NumLabels in size");
    // Make a copy since MATLAB no longer supports incrementing reference 
    // count from mex extensions... I miss you mxCreateReference :(
    mxArray* sc_copy = mxDuplicateArray(sc);
    mexMakeArrayPersistent(sc_copy);
	gcinstance.gco->setSmoothCost((EnergyTermType*)mxGetData(sc_copy));
	if (gcinstance.sc) mxDestroyArray(gcinstance.sc);
	gcinstance.sc = sc_copy;
}

GCO_EXPORT(gco_setlabelcost)
{
	LabelID* subset = 0;
	LabelID  subsetSize = 0;
	try {
		MATLAB_ASSERT(nlhs == 0, "Too many output arguments, expected 0");
		MATLAB_ASSERT(nrhs >= 2, "Not enough input arguments, expected at least 2");
		MATLAB_ASSERT(nrhs <= 3, "Too many input arguments, expected at most 3");
		MATLAB_ASSERT_HANDLE(1);
		MATLAB_ASSERT_INTYPE(2,EnergyTerm);
		GCInstanceInfo& gcinstance = sGetGCInstance(*(int*)mxGetData(prhs[0]));
		const mxArray* lc = prhs[1];
		if (mxGetN(lc) == 1 && mxGetM(lc) == 1) {
			if (nrhs == 2) {
				// Single cost used independently for all labels
				gcinstance.gco->setLabelCost(*(EnergyTermType*)mxGetData(lc));
			} else {
				// Cost for using an element from a specific subset of labels
				MATLAB_ASSERT_INTYPE(3,Label);
				subset = (LabelID*)mxGetData(prhs[2]);
				subsetSize = (LabelID)mxGetNumberOfElements(prhs[2]);
				for (LabelID i = 0; i < subsetSize; ++i)
					subset[i]--; // from 1-based to 0-based index
				gcinstance.gco->setLabelSubsetCost(subset,subsetSize,*(EnergyTermType*)mxGetData(lc));
				for (LabelID i = 0; i < subsetSize; ++i)
					subset[i]++; // from 0-based to 1-based index
			}
		} else {
			// Label costs provided as complete numeric array.
			MATLAB_ASSERT(mxGetNumberOfElements(lc) == gcinstance.gco->numLabels(),
			              "Numeric label cost must have either one element, or NumLabels elements");
			gcinstance.gco->setLabelCost((EnergyTermType*)mxGetData(lc));
		}
	} catch (...) {
		for (LabelID i = 0; i < subsetSize; ++i)
			subset[i]++; // from 0-based to 1-based index
		throw;
	}
}

GCO_EXPORT(gco_setneighbors)
{
	MATLAB_ASSERT_ARGCOUNT(0,2);
	MATLAB_ASSERT_HANDLE(1);
	MATLAB_ASSERT_INCLASS(2,mxDOUBLE_CLASS);
	GCInstanceInfo& gcinstance = sGetGCInstance(*(int*)mxGetData(prhs[0]));
	MATLAB_ASSERT(gcinstance.grid == false, "SetNeighbors can only be called on general graphs");
	GCoptimizationGeneralGraph* gco = static_cast<GCoptimizationGeneralGraph*>(gcinstance.gco);
	const mxArray* nb = prhs[1];
	MATLAB_ASSERT(mxIsSparse(nb), "Expected sparse array for neighbours");
	MATLAB_ASSERT(mxGetN(nb) == gcinstance.gco->numSites() && mxGetM(nb) == gcinstance.gco->numSites(),
	              "Neighbours array must be NumSites x NumSites in size");
	bool warned_rounding = false;
	bool warned_nonupper = false;
	mwIndex n = (mwIndex)mxGetN(nb);
	const mwIndex* ir = mxGetIr(nb);
	const mwIndex* jc = mxGetJc(nb);
	double*        pr = mxGetPr(nb);
	mwIndex count = 0;
	for (mwIndex c = 0; c < n; ++c) {
		mwIndex rowStart = jc[c]; 
		mwIndex rowEnd   = jc[c+1]; 
		for (mwIndex ri = rowStart; ri < rowEnd; ++ri)  {
			mwIndex r = ir[ri];
			MATLAB_ASSERT(r != c, "A site cannot neighbor itself; make sure diagonal is all zero");

			double dw = pr[count++];
			if (ctype2mx<EnergyTermType>::classid == mxINT32_CLASS) {
				int w = (int)dw;
				if ((double)w != dw && !warned_rounding) {
					mexWarnMsgTxt("Non-integer weight detected; rounding to int32");
					warned_rounding = true;
				}
			}
			if (r < c)
				gco->setNeighbors((SiteID)r, (SiteID)c, (EnergyTermType)dw);
			else {
				if (!warned_nonupper) {
					mexWarnMsgTxt("Neighbours array should be upper-triangular; entries below the diagnonal will be ignored");
					warned_nonupper = true;
				}
			}

		}
	}
}

GCO_EXPORT(gco_setlabelorder)
{
	MATLAB_ASSERT_ARGCOUNT(0,2);
	MATLAB_ASSERT_HANDLE(1);
	MATLAB_ASSERT_INTYPE(2,Label);
	GCInstanceInfo& gcinstance = sGetGCInstance(*(int*)mxGetData(prhs[0]));
	LabelID* order = (LabelID*)mxGetData(prhs[1]);
	LabelID size = (LabelID)mxGetNumberOfElements(prhs[1]);
	try {
		for (LabelID i = 0; i < size; ++i)
			order[i]--; // from 1-based to 0-based index
		gcinstance.gco->setLabelOrder(order,size);
		for (LabelID i = 0; i < size; ++i)
			order[i]++; // from 0-based to 1-based index
	} catch (...) {
		for (LabelID i = 0; i < size; ++i)
			order[i]++; // from 0-based to 1-based index
		throw;
	}
}

GCO_EXPORT(gco_setverbosity)
{
	MATLAB_ASSERT_ARGCOUNT(0,2);
	MATLAB_ASSERT_HANDLE(1);
	MATLAB_ASSERT_INCLASS(2,mxINT32_CLASS);
	GCInstanceInfo& gcinstance = sGetGCInstance(*(int*)mxGetData(prhs[0]));
	int level = *(int*)mxGetData(prhs[1]);
	MATLAB_ASSERT(level >= 0 && level <= 2,"Level must be in range 0..2");
	gcinstance.gco->setVerbosity(level);
}

GCO_EXPORT(gco_expansion)
{
	MATLAB_ASSERT_ARGCOUNT(1,2);
	MATLAB_ASSERT_HANDLE(1);
	MATLAB_ASSERT_INCLASS(2,mxINT32_CLASS);
	GCInstanceInfo& gcinstance = sGetGCInstance(*(int*)mxGetData(prhs[0]));
	int maxIter = *(int*)mxGetData(prhs[1]);
	EnergyType energy = gcinstance.gco->expansion(maxIter);
	mwSize outdim = 1;
	plhs[0] = mxCreateNumericArray(1, &outdim, cEnergyClassID, mxREAL);
	*(EnergyType*)mxGetData(plhs[0]) = energy;
}

GCO_EXPORT(gco_swap)
{
	MATLAB_ASSERT_ARGCOUNT(1,2);
	MATLAB_ASSERT_HANDLE(1);
	MATLAB_ASSERT_INCLASS(2,mxINT32_CLASS);
	GCInstanceInfo& gcinstance = sGetGCInstance(*(int*)mxGetData(prhs[0]));
	int maxIter = *(int*)mxGetData(prhs[1]);
	EnergyType energy = gcinstance.gco->swap();
	mwSize outdim = 1;
	plhs[0] = mxCreateNumericArray(1, &outdim, cEnergyClassID, mxREAL);
	*(EnergyType*)mxGetData(plhs[0]) = energy;
}

GCO_EXPORT(gco_alphaexpansion)
{
	MATLAB_ASSERT_ARGCOUNT(0,2);
	MATLAB_ASSERT_HANDLE(1);
	MATLAB_ASSERT_INTYPE(2,Label);
	GCInstanceInfo& gcinstance = sGetGCInstance(*(int*)mxGetData(prhs[0]));
	LabelID alpha = *(LabelID*)mxGetData(prhs[1])-1;
	gcinstance.gco->alpha_expansion(alpha);
}

GCO_EXPORT(gco_computeenergy)
{
	MATLAB_ASSERT_ARGCOUNT(4,1);
	MATLAB_ASSERT_HANDLE(1);
	GCInstanceInfo& gcinstance = sGetGCInstance(*(int*)mxGetData(prhs[0]));
	EnergyType energy = gcinstance.gco->compute_energy();
	mwSize outSize = 1;
	plhs[0] = mxCreateNumericArray(1, &outSize, cEnergyClassID, mxREAL);
	plhs[1] = mxCreateNumericArray(1, &outSize, cEnergyClassID, mxREAL);
	plhs[2] = mxCreateNumericArray(1, &outSize, cEnergyClassID, mxREAL);
	plhs[3] = mxCreateNumericArray(1, &outSize, cEnergyClassID, mxREAL);
	*(EnergyType*)mxGetData(plhs[3]) = gcinstance.gco->giveLabelEnergy();
	*(EnergyType*)mxGetData(plhs[2]) = gcinstance.gco->giveSmoothEnergy();
	*(EnergyType*)mxGetData(plhs[1]) = gcinstance.gco->giveDataEnergy();
	*(EnergyType*)mxGetData(plhs[0]) = *(EnergyType*)mxGetData(plhs[1])+*(EnergyType*)mxGetData(plhs[2])+*(EnergyType*)mxGetData(plhs[3]);
}

GCO_EXPORT(gco_setlabeling)
{
	MATLAB_ASSERT_ARGCOUNT(0,2);
	MATLAB_ASSERT_HANDLE(1);
	MATLAB_ASSERT_INTYPE(2,Label);
	GCInstanceInfo& gcinstance = sGetGCInstance(*(int*)mxGetData(prhs[0]));
	const mxArray* labeling = prhs[1];
	MATLAB_ASSERT(mxGetN(labeling) == gcinstance.gco->numSites() || mxGetM(labeling) == gcinstance.gco->numSites(),
	              "Labeling must be of length NumSites");
	LabelID* labeldata = (LabelID*)mxGetData(labeling);
	for (mwIndex i = 0; i < gcinstance.gco->numSites(); ++i) 
		MATLAB_ASSERT(labeldata[i] >= 1 && labeldata[i] <= gcinstance.gco->numLabels(),
		             "Labeling must be in range 1..NumLabels");
	for (mwIndex i = 0; i < gcinstance.gco->numSites(); ++i)
		gcinstance.gco->setLabel((SiteID)i, labeldata[i]-1);
}

GCO_EXPORT(gco_getlabeling)
{
	MATLAB_ASSERT_ARGCOUNT(1,3);
	MATLAB_ASSERT_HANDLE(1);
	MATLAB_ASSERT_INTYPE(2,Site);
	MATLAB_ASSERT_INTYPE(3,Site);
	GCInstanceInfo& gcinstance = sGetGCInstance(*(int*)mxGetData(prhs[0]));
	MATLAB_ASSERT(*(SiteID*)mxGetData(prhs[1]) > 0, "Start index must be in range 1..NumSites");
	SiteID start = *(SiteID*)mxGetData(prhs[1])-1;
	SiteID count = *(SiteID*)mxGetData(prhs[2]);
	MATLAB_ASSERT(start+count <= gcinstance.gco->numSites(), "End index must be in range 1..NumSites");
	mwSize mlcount = (mwSize)count;
	plhs[0] = mxCreateNumericArray(1, &mlcount, cLabelClassID, mxREAL);
	LabelID* labeling = (LabelID*)mxGetData(plhs[0]);
	gcinstance.gco->whatLabel(start, count, labeling);
	for ( SiteID i = 0; i < count; ++i )
		labeling[i]++; // convert C index to Matlab index
}

GCO_EXPORT(gco_getnumsites)
{
	MATLAB_ASSERT_ARGCOUNT(1,1);
	MATLAB_ASSERT_HANDLE(1);
	GCInstanceInfo& gcinstance = sGetGCInstance(*(int*)mxGetData(prhs[0]));
	mwSize outdim = 1;
	plhs[0] = mxCreateNumericArray(1, &outdim, cSiteClassID, mxREAL);
	*(SiteID*)mxGetData(plhs[0]) = gcinstance.gco->numSites();
}

GCO_EXPORT(gco_getnumlabels)
{
	MATLAB_ASSERT_ARGCOUNT(1,1);
	MATLAB_ASSERT_HANDLE(1);
	GCInstanceInfo& gcinstance = sGetGCInstance(*(int*)mxGetData(prhs[0]));
	mwSize outSize = 1;
	plhs[0] = mxCreateNumericArray(1, &outSize, cLabelClassID, mxREAL);
	*(LabelID*)mxGetData(plhs[0]) = gcinstance.gco->numLabels();
}
