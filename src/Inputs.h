/* 
 * File:   Inputs.h
 * Author: shailesh
 *
 * Created on 12 July, 2014, 9:54 PM
 */

#ifndef INPUTS_H
#define	INPUTS_H

#include <string>
#include <vector>

#include "configcentre.h"
#include "Centre.h"

class Inputs {
private:
	bool isgood_;
	std::string inpfile;

	class Control {
	private:
		hbin_t nsteps_;
		bool run_real2int_;
		bool calc_entropy_;
		bool gen_convgdata_;
		u_int dscrinfofreq_;
		u_int entcinfofreq_;
		u_int cachedscrtdimspercpu_;
		u_int cacheentcdimspercpu_;
		std::string outfilepath_;
		std::string infilepath_;
		std::string infofile_;
	public:
		const hbin_t getNsteps() const {
			return nsteps_;
		}

		const bool isDiscetize() const {
			return run_real2int_;
		}

		const bool isCalcEntropy() const {
			return calc_entropy_;
		}

		const bool isGenConvgdata() const {
			return gen_convgdata_;
		}

		const std::string getOutfilepath() const {
			return outfilepath_;
		}

		const std::string getInfilepath() const {
			return infilepath_;
		}

		void setNsteps(hbin_t v) {
			nsteps_ = v;
		}

		void setDiscretize(bool v) {
			run_real2int_ = v;
		}

		void setCalcEntropy(bool v) {
			calc_entropy_ = v;
		}

		void setGenConvgdata(bool v) {
			gen_convgdata_ = v;
		}

		void setOutfilepath(std::string v) {
			outfilepath_ = v;
		}

		void setInfilepath(std::string v) {
			infilepath_ = v;
		}

		u_int getDscrinfofreq() const {
			return dscrinfofreq_;
		}

		void setDscrinfofreq(u_int dscrinfofreq) {
			dscrinfofreq_ = dscrinfofreq;
		}

		u_int getEntcinfofreq() const {
			return entcinfofreq_;
		}

		void setEntcinfofreq(u_int entcinfofreq) {
			entcinfofreq_ = entcinfofreq;
		}

		const std::string& getInfofile() const {
			return infofile_;
		}

		void setInfofile(const std::string &infofile) {
			infofile_ = infofile;
		}

		u_int getCachedscrtdimspercpu() const {
			return cachedscrtdimspercpu_;
		}

		void setCachedscrtdimspercpu(u_int cachedscrtdimspercpu) {
			cachedscrtdimspercpu_ = cachedscrtdimspercpu;
		}

		u_int getCacheentcdimspercpu() const {
			return cacheentcdimspercpu_;
		}

		void setCacheentcdimspercpu(u_int cacheentcdimspercpu) {
			cacheentcdimspercpu_ = cacheentcdimspercpu;
		}
	};

	class Entropy {
	private:
		ull_int nframes_eff;
		ull_int startframe;
		ull_int strideframe;
		ull_int numframe;
		hbin_t nfiles;
		BATSet workset;
		bool shuffle_blocks;
		bool use_subset;
		bool use_neighbor;
		bool jacobian;
		u_int shuffle_block_times;
		std::string subsetfile;
		std::string neighborfile;
		std::string convrgfileext;
	public:
		BATSet getWorkSet() const {
			return workset;
		}

		void setWorkSet(BATSet ws) {
			this->workset = ws;
		}

		bool isJacobian() const {
			return jacobian;
		}

		void setJacobian(bool jacobian) {
			this->jacobian = jacobian;
		}

		const std::string& getNeighborfile() const {
			return neighborfile;
		}

		void setConvergenceFileExt(const std::string &cfext) {
			this->convrgfileext = cfext;
		}

		const std::string& getConvergenceFileExt() const {
			return convrgfileext;
		}
		void setNeighborfile(const std::string &neighborfile) {
			this->neighborfile = neighborfile;
		}

		hbin_t getNfiles() const {
			return nfiles;
		}

		void setNfiles(hbin_t nfiles) {
			this->nfiles = nfiles;
		}

		ull_int getNframesEff() const {
			return nframes_eff;
		}

		void setNframesEff(ull_int nframesEff) {
			nframes_eff = nframesEff;
		}

		ull_int getNumframe() const {
			return numframe;
		}

		void setNumframe(ull_int numframe) {
			this->numframe = numframe;
		}

		u_int getShuffleBlockTimes() const {
			return shuffle_block_times;
		}

		void setShuffleBlockTimes(u_int shuffleBlockTimes) {
			shuffle_block_times = shuffleBlockTimes;
		}

		bool isShuffleBlocks() const {
			return shuffle_blocks;
		}

		void setShuffleBlocks(bool shuffleBlocks) {
			shuffle_blocks = shuffleBlocks;
		}

		ull_int getStartframe() const {
			return startframe;
		}

		void setStartframe(ull_int startframe) {
			this->startframe = startframe;
		}

		ull_int getStrideframe() const {
			return strideframe;
		}

		void setStrideframe(ull_int strideframe) {
			this->strideframe = strideframe;
		}

		const std::string& getSubsetfile() const {
			return subsetfile;
		}

		void setSubsetfile(const std::string &subsetfile) {
			this->subsetfile = subsetfile;
		}

		bool isUseNeighbor() const {
			return use_neighbor;
		}

		void setUseNeighbor(bool useNeighbor) {
			use_neighbor = useNeighbor;
		}

		bool isUseSubset() const {
			return use_subset;
		}

		void setUseSubset(bool useSubset) {
			use_subset = useSubset;
		}
	};

	class Discretize {
	private:
		u_int nbond;
		u_int nangle;
		u_int ndihed;
		bool shuffle_frames;
		bool optimizedih;
		bool shuffle_dofs;
		ull_int randseed;
		PDFMethod pdfmethod;
		std::vector<std::string> fnames;
		std::vector<ull_int> nframes;
	public:
		const std::vector<std::string>& getFnames() const {
			return fnames;
		}

		void addFname(const std::string fn) {
			this->fnames.push_back(fn);
		}
		void setFnames(const std::vector<std::string> &fnames) {
			this->fnames = fnames;
		}

		u_int getNangle() const {
			return nangle;
		}

		void setNangle(u_int nangle) {
			this->nangle = nangle;
		}

		u_int getNbond() const {
			return nbond;
		}

		void setNbond(u_int nbond) {
			this->nbond = nbond;
		}

		u_int getNdihed() const {
			return ndihed;
		}

		void setNdihed(u_int ndihed) {
			this->ndihed = ndihed;
		}

		const std::vector<ull_int>& getNframes() const {
			return nframes;
		}

		void setNframes(const std::vector<ull_int> &nframes) {
			this->nframes = nframes;
		}

		void setNframe(const ull_int nframe, u_int index) {
			this->nframes[index] = nframe;
		}

		void addNframe(const ull_int nf) {
			this->nframes.push_back(nf);
		}

		bool isOptimizedih() const {
			return optimizedih;
		}

		void setOptimizedih(bool optimizedih) {
			this->optimizedih = optimizedih;
		}

		ull_int getRandseed() const {
			return randseed;
		}

		void setRandseed(ull_int randseed) {
			this->randseed = randseed;
		}

		bool isShuffleFrames() const {
			return shuffle_frames;
		}

		void setShuffleFrames(bool shuffle) {
			this->shuffle_frames = shuffle;
		}

		bool isShuffleDofs() const {
			return shuffle_dofs;
		}

		void setShuffleDofs(bool shuffleDofs) {
			shuffle_dofs = shuffleDofs;
		}

		PDFMethod getPdfmethod() const {
			return pdfmethod;
		}

		void setPdfmethod(PDFMethod pdfmethod) {
			this->pdfmethod = pdfmethod;
		}
	};

	class vMisesKDE {
	private:
		bool writefreq;
		BATSet writeset;
		hbin_t writestepstart;
		hbin_t writestepstride;
		hbin_t nmaxconf;
		double kappa_value;
		u_int sdo_nstep;
		u_int sdo_niter;
		double sdo_conv_limit;
	public:
		BATSet getWriteSet() const {
			return writeset;
		}

		void setWriteSet(BATSet ws) {
			this->writeset = ws;
		}

		hbin_t getNmaxconf() const {
			return nmaxconf;
		}

		void setNmaxconf(hbin_t nmaxconf) {
			this->nmaxconf = nmaxconf;
		}

		bool isWritefreq() const {
			return writefreq;
		}

		void setWritefreq(bool writefreq) {
			this->writefreq = writefreq;
		}

		hbin_t getWritestepstart() const {
			return writestepstart;
		}

		void setWritestepstart(hbin_t writestepstart) {
			this->writestepstart = writestepstart;
		}

		hbin_t getWritestepstride() const {
			return writestepstride;
		}

		void setWritestepstride(hbin_t writestepstride) {
			this->writestepstride = writestepstride;
		}

		double getSdoConvLimit() const {
			return sdo_conv_limit;
		}

		void setSdoConvLimit(double sdoConvLimit) {
			sdo_conv_limit = sdoConvLimit;
		}

		u_int getSdoNiter() const {
			return sdo_niter;
		}

		void setSdoNiter(u_int sdoNiter) {
			sdo_niter = sdoNiter;
		}

		u_int getSdoNstep() const {
			return sdo_nstep;
		}

		void setSdoNstep(u_int sdoNstep) {
			sdo_nstep = sdoNstep;
		}

		double getKappaValue() const {
			return kappa_value;
		}

		void setKappaValue(double kappaValue) {
			kappa_value = kappaValue;
		}
	};

	class Histogram {
	private:
		bool writefreq;
		BATSet writeset;
		hbin_t writestepstart;
		hbin_t writestepstride;
		hbin_t nbin_ref;
		std::vector<hbin_t> bin_schemes;
	public:
		BATSet getWriteSet() const {
			return writeset;
		}

		void setWriteSet(BATSet ws) {
			this->writeset = ws;
		}
		const std::vector<hbin_t>& getBinSchemes() const {
			return bin_schemes;
		}

		void setBinSchemes(const std::vector<hbin_t> &binSchemes) {
			bin_schemes = binSchemes;
		}

		void addBinScheme(hbin_t bs) {
			bin_schemes.push_back(bs);
		}

		hbin_t getNbinRef() const {
			return nbin_ref;
		}

		void setNbinRef(hbin_t nbinRef) {
			nbin_ref = nbinRef;
		}

		bool isWritefreq() const {
			return writefreq;
		}

		void setWritefreq(bool writefreq) {
			this->writefreq = writefreq;
		}

		hbin_t getWritestepstart() const {
			return writestepstart;
		}

		void setWritestepstart(hbin_t writestepstart) {
			this->writestepstart = writestepstart;
		}

		hbin_t getWritestepstride() const {
			return writestepstride;
		}

		void setWritestepstride(hbin_t writestepstride) {
			this->writestepstride = writestepstride;
		}
	};

	Control control;
	EntropyMethods method;
	std::vector<EntropyEstimators> estimators;
	Entropy4Data entropy4data;
	SubsetData subset;
	Neighbors neighbors;
	Discretize bats;
	Histogram hist;
	vMisesKDE vmkde;
	Entropy entropy;
	std::vector<std::string> requiredfiles;
	std::vector<std::string> unwantedfiles;

	int status = 0;
public:

	Inputs(const std::string &inpFileName);

	bool Init(const std::string &inpFileName);

	friend std::ostream& operator<<(std::ostream &os, const Inputs &inp);

	bool isGood() const {
		return isgood_;
	}

	const Discretize& getBats() const {
		return bats;
	}

	void setBats(const Discretize &bats) {
		this->bats = bats;
	}

	bool isCalcEntropy() const {
		return control.isCalcEntropy();
	}

	const Control& getControl() const {
		return control;
	}

	void setControl(const Control &control) {
		this->control = control;
	}

	const Entropy& getEntropy() const {
		return entropy;
	}

	void setEntropy(const Entropy &entropy) {
		this->entropy = entropy;
	}

	Entropy4Data getEntropy4data() const {
		return entropy4data;
	}

	void setEntropy4data(Entropy4Data entropy4data) {
		this->entropy4data = entropy4data;
	}

	const std::vector<EntropyEstimators> getEstimators() const {
		return estimators;
	}

	const std::vector<std::string> getRequiredfiles() const {
		return requiredfiles;
	}

	const std::vector<std::string> getUnwantedfiles() const {
		return unwantedfiles;
	}

	void addEstimator(EntropyEstimators estimator) {
		this->estimators.push_back(estimator);
	}

	bool isGenConvgdata() const {
		return control.isGenConvgdata();
	}

	void setGenConvgdata(bool genConvgdata) {
		control.setGenConvgdata(genConvgdata);
	}

	const Histogram& getHist() const {
		return hist;
	}

	void setHist(const Histogram &hist) {
		this->hist = hist;
	}

	const std::string& getInpfile() const {
		return inpfile;
	}

	void setInpfile(const std::string &inpfile) {
		this->inpfile = inpfile;
	}

	EntropyMethods getMethod() const {
		return method;
	}

	void setMethod(EntropyMethods method) {
		this->method = method;
	}

	const Neighbors& getNeighbors() const {
		return neighbors;
	}

//      void setNeighbors(Neighbors*& neighbors) {
//         this->neighbors = neighbors;
//      }

	bool isRunReal2int() const {
		return control.isDiscetize();
	}

	void setRunReal2int(bool runReal2int) {
		control.setDiscretize(runReal2int);
	}

	int getStatus() const {
		return status;
	}

	const SubsetData& getSubset() const {
		return subset;
	}

	vMisesKDE getVmkde() const {
		return vmkde;
	}

	void setVmkde(vMisesKDE vmkde) {
		this->vmkde = vmkde;
	}
};

#endif	/* INPUTS_H */

