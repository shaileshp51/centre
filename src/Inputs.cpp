#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <sstream>
#include <cctype>

#include "Utils.h"
#include "Inputs.h"
#include "Netcdf_BAT.h"

Inputs::Inputs(const std::string &inpFileName)
{
	inpfile = std::string(inpFileName.c_str());
	entropy4data = Entropy4Data::ALL;
	method = EntropyMethods::MIE;
	//subset = NULL;
	//neighbors = NULL;

	control.setDiscretize(true);
	control.setCalcEntropy(true);
	control.setGenConvgdata(true);
	control.setNsteps(20);
	control.setCachedscrtdimspercpu(5);
	control.setCacheentcdimspercpu(10);
	control.setDscrinfofreq(25);
	control.setEntcinfofreq(500);
	control.setInfofile("progress.info");

	bats.setShuffleDofs(false);
	bats.setShuffleFrames(false);
	bats.setOptimizedih(true);
	bats.setNbond(0);
	bats.setNangle(0);
	bats.setNdihed(0);
	bats.setRandseed(0);
	bats.setPdfmethod(PDFMethod::HISTOGRAM);

	hist.setWritefreq(false);
	hist.setWriteSet(BATSet::NOTHING);
	hist.setWritestepstart(1);
	hist.setWritestepstride(1);
	hist.setWritestepstart(255); // a Magic number: just to singifily unintialized value

	vmkde.setWritefreq(false);
	vmkde.setWriteSet(BATSet::NOTHING);
	vmkde.setWritestepstride(1);
	vmkde.setWritestepstart(255); // a Magic number: just to singifily unintialized value
	vmkde.setNmaxconf(5);
	vmkde.setKappaValue(1.0);
	vmkde.setSdoNstep(5);
	vmkde.setSdoNiter(1000);
	vmkde.setSdoConvLimit(0.0001);

	entropy.setNfiles(0);
	entropy.setStartframe(0);
	entropy.setStrideframe(1);
	entropy.setNumframe(0);
	entropy.setNframesEff(0);
	entropy.setShuffleBlockTimes(1);
	entropy.setWorkSet(BATSet::XY2D);
	entropy.setJacobian(true);
	entropy.setShuffleBlocks(false);
	entropy.setUseSubset(false);
	entropy.setUseNeighbor(false);
	entropy.setConvergenceFileExt("csv");
	entropy.setSubsetfile("");
	entropy.setNeighborfile("");
	method = EntropyMethods::MIST;

	isgood_ = false;
}

bool Inputs::Init(const std::string &inpFileName)
{
	inpfile = std::string(inpFileName.c_str());
	entropy4data = Entropy4Data::ALL;
	method = EntropyMethods::MIE;
	//subset = NULL;
	//neighbors = NULL;
	hist.setWritestepstride(1);
	hist.setWritestepstart(255); // a Magic number: just to represent unintialized value

	control.setDscrinfofreq(25);
	control.setEntcinfofreq(500);
	bool isOK = true;
	if (Utils::fileExists(inpFileName))
	{
		isOK = true;
		auto symbolTable = Utils::parseConfiguration(inpFileName);

		// std::cout << "printing parsed options" << std::endl;
		for (auto x : symbolTable)
		{
			// printParsedConfiguration(x);

			std::string tagO = x.first;
			if (tagO == "control.nsteps")
			{
				(*this).control.setNsteps((hbin_t)atoi(x.second[0].c_str()));
			}
			else if (tagO == "control.discretize")
			{
				(*this).control.setDiscretize(Utils::string2bool(x.second[0]));
			}
			else if (tagO == "control.calcentropy")
			{
				(*this).control.setCalcEntropy(Utils::string2bool(x.second[0]));
			}
			else if (tagO == "control.genconvgdata")
			{
				(*this).control.setGenConvgdata(
					Utils::string2bool(x.second[0]));
			}
			else if (tagO == "control.dscrinfofreq")
			{
				(*this).control.setDscrinfofreq(
					(u_int)atoi(x.second[0].c_str()));
			}
			else if (tagO == "control.entcinfofreq")
			{
				(*this).control.setEntcinfofreq(
					(u_int)atoi(x.second[0].c_str()));
			}
			else if (tagO == "control.cachedscrtdimspercpu")
			{
				(*this).control.setCachedscrtdimspercpu(
					(u_int)atoi(x.second[0].c_str()));
			}
			else if (tagO == "control.cacheentcdimspercpu")
			{
				(*this).control.setCacheentcdimspercpu(
					(u_int)atoi(x.second[0].c_str()));
			}
			else if (tagO == "control.outfilepath")
			{
				(*this).control.setOutfilepath(x.second[0]);
			}
			else if (tagO == "control.infilepath")
			{
				(*this).control.setInfilepath(x.second[0]);
			}
			else if (tagO == "control.infofile")
			{
				(*this).control.setInfofile(x.second[0]);
			}
			else if (tagO == "discretization.nbond")
			{
				bats.setNbond((u_int)atoi(x.second[0].c_str()));
			}
			else if (tagO == "discretization.nangle")
			{
				bats.setNangle((u_int)atoi(x.second[0].c_str()));
			}
			else if (tagO == "discretization.ndihed")
			{
				bats.setNdihed((u_int)atoi(x.second[0].c_str()));
			}
			else if (tagO == "discretization.fname")
			{
				bats.addFname(x.second[0].c_str());
			}
			else if (tagO == "discretization.nframe")
			{
				bats.addNframe((ull_int)atoi(x.second[0].c_str()));
			}
			else if (tagO == "discretization.randseed")
			{
				bats.setRandseed((ull_int)atoi(x.second[0].c_str()));
			}
			else if (tagO == "discretization.shuffleframes")
			{
				bats.setShuffleFrames(Utils::string2bool(x.second[0]));
			}
			else if (tagO == "discretization.shuffledofs")
			{
				bats.setShuffleDofs(Utils::string2bool(x.second[0]));
			}
			else if (tagO == "discretization.optimizedih")
			{ // optimizedih
				bats.setOptimizedih(Utils::string2bool(x.second[0]));
			}
			else if (tagO == "discretization.pdfmethod")
			{
				std::transform(x.second[0].begin(), x.second[0].end(),
							   x.second[0].begin(), ::tolower);
				if (x.second[0] == "histogram")
				{
					bats.setPdfmethod(PDFMethod::HISTOGRAM);
				}
				else if (x.second[0] == "vonmiseskde")
				{
					bats.setPdfmethod(PDFMethod::vonMisesKDE);
				}
			}
			else if (tagO == "vonmiseskde.nmaxconf")
			{
				vmkde.setNmaxconf((hbin_t)atoi(x.second[0].c_str()));
			}
			else if (tagO == "vonmiseskde.writefreq")
			{
				vmkde.setWritefreq(Utils::string2bool(x.second[0]));
			}
			else if (tagO == "vonmiseskde.writeset")
			{
				int prev_i = 0;
				for (auto tgv : x.second)
				{
					BATSet curr;
					if (tgv == "NONE")
					{
						vmkde.setWriteSet(BATSet::NOTHING);
						break;
					}
					else if (tgv == "B1D")
					{
						curr = BATSet::B1D;
					}
					else if (tgv == "A1D")
					{
						curr = BATSet::A1D;
					}
					else if (tgv == "D1D")
					{
						curr = BATSet::D1D;
					}
					else if (tgv == "1D")
					{
						curr = BATSet::X1D;
					}
					else if (tgv == "BB2D")
					{
						curr = BATSet::BB2D;
					}
					else if (tgv == "AA2D")
					{
						curr = BATSet::AA2D;
					}
					else if (tgv == "DD2D")
					{
						curr = BATSet::DD2D;
					}
					else if (tgv == "XX2D")
					{
						curr = BATSet::XX2D;
					}
					else if (tgv == "BA2D")
					{
						curr = BATSet::BA2D;
					}
					else if (tgv == "BD2D")
					{
						curr = BATSet::BD2D;
					}
					else if (tgv == "AD2D")
					{
						curr = BATSet::AD2D;
					}
					else if (tgv == "2D")
					{
						curr = BATSet::XY2D;
					}
					if (prev_i == 0)
					{
						vmkde.setWriteSet(curr);
					}
					else
					{
						BATSet prv = vmkde.getWriteSet();
						curr = static_cast<BATSet>(prv | curr);
						vmkde.setWriteSet(curr);
					}
					++prev_i;
				}
			}
			else if (tagO == "vonmiseskde.writestepstart")
			{
				vmkde.setWritestepstart((hbin_t)atoi(x.second[0].c_str()));
			}
			else if (tagO == "vonmiseskde.writestepstride")
			{
				vmkde.setWritestepstride((hbin_t)atoi(x.second[0].c_str()));
			}
			else if (tagO == "vonmiseskde.kappa")
			{
				vmkde.setKappaValue((double)atof(x.second[0].c_str()));
			}
			else if (tagO == "vonmiseskde.sdosteps")
			{
				vmkde.setSdoNstep((u_int)atoi(x.second[0].c_str()));
			}
			else if (tagO == "vonmiseskde.sdoiterations")
			{
				vmkde.setSdoNiter((u_int)atoi(x.second[0].c_str()));
			}
			else if (tagO == "vonmiseskde.sdoconvlimit")
			{
				vmkde.setSdoConvLimit((double)atof(x.second[0].c_str()));
			}
			else if (tagO == "histogram.referencenbins")
			{
				hist.setNbinRef((hbin_t)atoi(x.second[0].c_str()));
			}
			else if (tagO == "histogram.nbins")
			{
				for (auto tgv : x.second)
				{
					hist.addBinScheme((hbin_t)atoi(tgv.c_str()));
				}
			}
			else if (tagO == "histogram.writefreq")
			{
				hist.setWritefreq(Utils::string2bool(x.second[0]));
			}
			else if (tagO == "histogram.writeset")
			{
				int prev_i = 0;
				for (auto tgv : x.second)
				{
					BATSet curr;
					if (tgv == "NONE")
					{
						hist.setWriteSet(BATSet::NOTHING);
						break;
					}
					else if (tgv == "B1D")
					{
						curr = BATSet::B1D;
					}
					else if (tgv == "A1D")
					{
						curr = BATSet::A1D;
					}
					else if (tgv == "D1D")
					{
						curr = BATSet::D1D;
					}
					else if (tgv == "1D")
					{
						curr = BATSet::X1D;
					}
					else if (tgv == "BB2D")
					{
						curr = BATSet::BB2D;
					}
					else if (tgv == "AA2D")
					{
						curr = BATSet::AA2D;
					}
					else if (tgv == "DD2D")
					{
						curr = BATSet::DD2D;
					}
					else if (tgv == "XX2D")
					{
						curr = BATSet::XX2D;
					}
					else if (tgv == "BA2D")
					{
						curr = BATSet::BA2D;
					}
					else if (tgv == "BD2D")
					{
						curr = BATSet::BD2D;
					}
					else if (tgv == "AD2D")
					{
						curr = BATSet::AD2D;
					}
					else if (tgv == "2D")
					{
						curr = BATSet::XY2D;
					}
					if (prev_i == 0)
					{
						hist.setWriteSet(curr);
					}
					else
					{
						BATSet prv = hist.getWriteSet();
						curr = static_cast<BATSet>(prv | curr);
						hist.setWriteSet(curr);
					}
					++prev_i;
				}
			}
			else if (tagO == "histogram.writestepstart")
			{
				hist.setWritestepstart((hbin_t)atoi(x.second[0].c_str()));
			}
			else if (tagO == "histogram.writestepstride")
			{
				hist.setWritestepstride((hbin_t)atoi(x.second[0].c_str()));
			}
			else if (tagO == "entropy.workset")
			{
				int prev_i = 0;
				for (auto tgv : x.second)
				{
					BATSet curr;
					if (tgv == "NONE")
					{
						entropy.setWorkSet(BATSet::NOTHING);
						break;
					}
					else if (tgv == "B1D")
					{
						curr = BATSet::B1D;
					}
					else if (tgv == "A1D")
					{
						curr = BATSet::A1D;
					}
					else if (tgv == "D1D")
					{
						curr = BATSet::D1D;
					}
					else if (tgv == "1D")
					{
						curr = BATSet::X1D;
					}
					else if (tgv == "BB2D")
					{
						curr = BATSet::BB2D;
					}
					else if (tgv == "AA2D")
					{
						curr = BATSet::AA2D;
					}
					else if (tgv == "DD2D")
					{
						curr = BATSet::DD2D;
					}
					else if (tgv == "XX2D")
					{
						curr = BATSet::XX2D;
					}
					else if (tgv == "BA2D")
					{
						curr = BATSet::BA2D;
					}
					else if (tgv == "BD2D")
					{
						curr = BATSet::BD2D;
					}
					else if (tgv == "AD2D")
					{
						curr = BATSet::AD2D;
					}
					else if (tgv == "2D")
					{
						curr = BATSet::XY2D;
					}
					if (prev_i == 0)
					{
						entropy.setWorkSet(curr);
					}
					else
					{
						BATSet prv = entropy.getWorkSet();
						curr = static_cast<BATSet>(prv | curr);
						entropy.setWorkSet(curr);
					}
					++prev_i;
				}
			}
			else if (tagO == "entropy.shuffleblocks")
			{
				entropy.setShuffleBlocks(Utils::string2bool(x.second[0]));
			}
			else if (tagO == "entropy.shuffleblocktimes")
			{
				entropy.setShuffleBlockTimes(atoi(x.second[0].c_str()));
			}
			else if (tagO == "entropy.usesubset")
			{
				entropy.setUseSubset(Utils::string2bool(x.second[0]));
			}
			else if (tagO == "entropy.useneighbor")
			{
				entropy.setUseNeighbor(Utils::string2bool(x.second[0]));
			}
			else if (tagO == "entropy.jacobian")
			{
				entropy.setJacobian(Utils::string2bool(x.second[0]));
			}
			else if (tagO == "entropy.subsetfile")
			{
				entropy.setSubsetfile(x.second[0]);
			}
			else if (tagO == "entropy.neighborfile")
			{
				entropy.setNeighborfile(x.second[0]);
			}
			else if (tagO == "entropy.startframe")
			{
				entropy.setStartframe((ull_int)atoi(x.second[0].c_str()));
			}
			else if (tagO == "entropy.strideframe")
			{
				entropy.setStrideframe((ull_int)atoi(x.second[0].c_str()));
			}
			else if (tagO == "entropy.nframe")
			{
				entropy.setNumframe((ull_int)atoi(x.second[0].c_str()));
			}
			else if (tagO == "entropy.scoringmethod")
			{
				if (x.second[0] == "MIE")
				{
					(*this).method = EntropyMethods::MIE;
				}
				else if (x.second[0] == "AMIE")
				{
					(*this).method = EntropyMethods::AMIE;
				}
				else if (x.second[0] == "MIST")
				{
					(*this).method = EntropyMethods::MIST;
				}
				else if (x.second[0] == "AMIST")
				{
					(*this).method = EntropyMethods::AMIST;
				}
			}
			else if (tagO == "entropy.estimator")
			{
				for (auto tgv : x.second)
				{
					if (tgv == "ML")
					{
						(*this).addEstimator(EntropyEstimators::ML);
					}
					else if (tgv == "MM")
					{
						(*this).addEstimator(EntropyEstimators::MM);
					}
					else if (tgv == "CS")
					{
						(*this).addEstimator(EntropyEstimators::CS);
					}
					else if (tgv == "JS")
					{
						(*this).addEstimator(EntropyEstimators::JS);
					}
				}
			}
			else if (tagO == "entropy.reportfileext")
			{
				entropy.setConvergenceFileExt(x.second[0]);
			}
		}

		ull_int numfrm = 0;
		if (bats.getFnames().size() != bats.getNframes().size())
		{
			for (u_int i = bats.getNframes().size();
				 i < bats.getFnames().size(); ++i)
			{
				bats.addNframe((ull_int)0);
			}
			for (u_int i = 0; i < bats.getFnames().size(); ++i)
			{
				string fname = control.getInfilepath() + bats.getFnames()[i];
				Netcdf_BAT trajIn;
				NetcdfFile::NCTYPE fltype = trajIn.GetNetcdfConventions(
					fname.c_str());

				if (fltype == NetcdfFile::NCTYPE::NC_CENTRETRAJ)
				{
					trajIn.NC_openRead(fname);
					trajIn.setupFrameDim();
					trajIn.setupTime();
					CoordinateInfo cInfo;
					trajIn.setupRead(cInfo);
					if (bats.getNbond() == 0)
					{
						bats.setNbond(cInfo.nBonds());
					}
					else if (bats.getNbond() != (u_int)cInfo.nBonds())
					{
						mprinterr(
							"earlier nbond=%d != current nbond=%d, file=%s",
							bats.getNbond(), cInfo.nBonds(), fname.c_str());
						exit(0);
					}
					if (bats.getNangle() == 0)
					{
						bats.setNangle(cInfo.nAngles());
					}
					else if (bats.getNangle() != (u_int)cInfo.nAngles())
					{
						mprinterr(
							"earlier nangle=%d != current nangle=%d, file=%s",
							bats.getNangle(), cInfo.nAngles(),
							fname.c_str());
						exit(0);
					}
					if (bats.getNdihed() == 0)
					{
						bats.setNdihed(cInfo.nDihedrals());
					}
					else if (bats.getNdihed() != (u_int)cInfo.nDihedrals())
					{
						mprinterr(
							"earlier ndihed=%d != current ndihed=%d, file=%s",
							bats.getNdihed(), cInfo.nDihedrals(),
							fname.c_str());
						exit(0);
					}
					ull_int nframes_trjs = trajIn.Ncframe();
					if (nframes_trjs <= 0)
					{
						mprinterr(
							"file '%s' input for discretization must have  more than 0 frames",
							fname.c_str());
						exit(0);
					}
					else
					{
						bats.setNframe(nframes_trjs, i);
					}
				}
			}
		}
		for (u_int i = 0; i < bats.getNframes().size(); ++i)
		{
			numfrm += bats.getNframes()[i];
		}
		if (entropy.getNumframe() > numfrm || entropy.getNumframe() <= 0)
		{
			entropy.setNumframe(numfrm);
		}
		ull_int nframes_eff = 0;
		if (entropy.getNumframe() > 0 && (entropy.getStartframe() > 0 || entropy.getStrideframe() > 1))
		{
			nframes_eff = (entropy.getNumframe() - entropy.getStartframe() + entropy.getStrideframe() - 1) / entropy.getStrideframe();
		}
		else
		{
			nframes_eff = entropy.getNumframe();
		}
		entropy.setNframesEff(nframes_eff);
		if (isOK && !Utils::isDirectory(control.getInfilepath()))
		{
			isOK = false;
			std::cout << "InputReadError:: infilepath ("
					  << control.getInfilepath() << ") is not a directory"
					  << std::endl;
		}
		if (isOK && !Utils::isDirectory(control.getOutfilepath()))
		{
			isOK &= false;
			std::cout << "InputReadError:: outfilepath ("
					  << control.getOutfilepath() << ") is not a directory"
					  << std::endl;
		}
		if (isOK && entropy.isUseSubset())
		{
			std::string fname = control.getInfilepath() + "/" + entropy.getSubsetfile();
			if (isOK && Utils::fileExists(fname))
			{
				subset = SubsetData(fname, entropy.getWorkSet());
				entropy4data = Entropy4Data::SUBSET;
				if (subset.isGood() && entropy.isUseNeighbor())
				{
					std::string fn = control.getInfilepath() + "/" + entropy.getNeighborfile();
					if (Utils::fileExists(fn))
					{
						neighbors = Neighbors(fn, entropy.getWorkSet());
						isOK &= neighbors.isGood();
					}
					else
					{
						isOK = false;
						std::cout << "InputReadError:: neighborfile ("
								  << entropy.getNeighborfile()
								  << ") does not exists" << std::endl;
					}
				}
				else
				{
					isOK &= subset.isGood();
					neighbors = Neighbors(subset, entropy.getWorkSet());
				}
			}
			else
			{
				isOK &= false;
				std::cout << "InputReadError:: subsetfile ("
						  << entropy.getSubsetfile() << ") does not exists"
						  << std::endl;
			}
		}
		else
		{
			subset = SubsetData(bats.getNbond(), bats.getNangle(),
								bats.getNdihed(), entropy.getWorkSet());
			neighbors = Neighbors(subset, entropy.getWorkSet());
		}
		//		if (control.getNsteps() == 0) {
		//			control.setNsteps(10);
		//		}
		//		if (control.getCachedscrtdimspercpu() == 0) {
		//			control.setCachedscrtdimspercpu(5);
		//		}
		//		if (control.getCacheentcdimspercpu() == 0) {
		//			control.setCacheentcdimspercpu(10);
		//		}
		if (hist.isWritefreq() && hist.getWritestepstart() == 255)
		{
			hist.setWritestepstart(control.getNsteps() - 1);
		}
		if (vmkde.isWritefreq() && vmkde.getWritestepstart() == 255)
		{
			vmkde.setWritestepstart(control.getNsteps() - 1);
		}
		if (bats.getPdfmethod() == PDFMethod::vonMisesKDE)
		{
			hist.setNbinRef(vmkde.getNmaxconf());
			std::vector<hbin_t> tmp1;
			tmp1.push_back(vmkde.getNmaxconf());
			hist.setBinSchemes(tmp1);
			hist.setWritefreq(vmkde.isWritefreq());
			hist.setWriteSet(vmkde.getWriteSet());
		}
		if (estimators.size() == 0)
		{
			addEstimator(EntropyEstimators::MM);
		}

		if (control.isDiscetize())
		{
			for (auto fln : bats.getFnames())
			{
				requiredfiles.push_back(control.getInfilepath() + "/" + fln);
			}
			BATSet workset = entropy.getWorkSet();
			if (bats.getPdfmethod() == PDFMethod::vonMisesKDE)
			{
				if (workset & BATSet::D1D)
				{
					entropy.setWorkSet(BATSet::D1D);
				}
				if (workset & BATSet::DD2D)
				{
					entropy.setWorkSet(static_cast<BATSet>(BATSet::D1D | BATSet::DD2D));
				}
			}
			workset = entropy.getWorkSet();

			if ((workset & BATSet::B1D) || (workset & BATSet::BB2D) || (workset & BATSet::BA2D) || (workset & BATSet::BD2D))
			{
				unwantedfiles.push_back(
					control.getOutfilepath() + "/bin_bonds.nc");
			}
			if ((workset & BATSet::A1D) || (workset & BATSet::AA2D) || (workset & BATSet::BA2D) || (workset & BATSet::AD2D))
			{
				unwantedfiles.push_back(
					control.getOutfilepath() + "/bin_angles.nc");
			}
			if ((workset & BATSet::D1D) || (workset & BATSet::DD2D) || (workset & BATSet::AD2D) || (workset & BATSet::BD2D))
			{
				unwantedfiles.push_back(
					control.getOutfilepath() + "/bin_torsions.nc");
			}

			if (control.isCalcEntropy())
			{
				bool isWritefreq = hist.isWritefreq();
				BATSet histwriteset = hist.getWriteSet();

				if (bats.getPdfmethod() == PDFMethod::vonMisesKDE)
				{
					isWritefreq = vmkde.isWritefreq();
					histwriteset = vmkde.getWriteSet();
				}
				if (workset & BATSet::B1D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_bnd-1d.nc");
				}
				if (workset & BATSet::A1D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_ang-1d.nc");
				}
				if (workset & BATSet::D1D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_tor-1d.nc");
				}
				if (workset & BATSet::BB2D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_bnd-2d.nc");
				}
				if (workset & BATSet::AA2D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_ang-2d.nc");
				}
				if (workset & BATSet::DD2D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_tor-2d.nc");
				}
				if (workset & BATSet::BA2D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_ba-2d.nc");
				}
				if (workset & BATSet::BD2D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_bd-2d.nc");
				}
				if (workset & BATSet::AD2D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_ad-2d.nc");
				}

				if (isWritefreq)
				{
					if (histwriteset & workset & BATSet::B1D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_bnd-1d.nc");
					}
					if (histwriteset & workset & BATSet::A1D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_ang-1d.nc");
					}
					if (histwriteset & workset & BATSet::D1D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_tor-1d.nc");
					}
					if (histwriteset & BATSet::BB2D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_bnd-2d.nc");
					}
					if (histwriteset & workset & BATSet::AA2D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_ang-2d.nc");
					}
					if (histwriteset & workset & BATSet::DD2D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_tor-2d.nc");
					}
					if (histwriteset & workset & BATSet::BA2D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_ba-2d.nc");
					}
					if (histwriteset & workset & BATSet::BD2D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_bd-2d.nc");
					}
					if (histwriteset & workset & BATSet::AD2D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_ad-2d.nc");
					}
				}
			}
		}
		else
		{
			BATSet workset = entropy.getWorkSet();
			if (control.isCalcEntropy())
			{

				if ((workset & BATSet::B1D) || (workset & BATSet::BB2D) || (workset & BATSet::BA2D) || (workset & BATSet::BD2D))
				{
					requiredfiles.push_back(
						control.getOutfilepath() + "/bin_bonds.nc");
				}
				if ((workset & BATSet::A1D) || (workset & BATSet::AA2D) || (workset & BATSet::BA2D) || (workset & BATSet::AD2D))
				{
					requiredfiles.push_back(
						control.getOutfilepath() + "/bin_angles.nc");
				}
				if ((workset & BATSet::D1D) || (workset & BATSet::DD2D) || (workset & BATSet::AD2D) || (workset & BATSet::BD2D))
				{
					requiredfiles.push_back(
						control.getOutfilepath() + "/bin_torsions.nc");
				}

				if (workset & BATSet::B1D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_bnd-1d.nc");
				}
				if (workset & BATSet::A1D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_ang-1d.nc");
				}
				if (workset & BATSet::D1D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_tor-1d.nc");
				}
				if (workset & BATSet::BB2D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_bnd-2d.nc");
				}
				if (workset & BATSet::AA2D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_ang-2d.nc");
				}
				if (workset & BATSet::DD2D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_tor-2d.nc");
				}
				if (workset & BATSet::BA2D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_ba-2d.nc");
				}
				if (workset & BATSet::BD2D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_bd-2d.nc");
				}
				if (workset & BATSet::AD2D)
				{
					unwantedfiles.push_back(
						control.getOutfilepath() + "/entcontri_ad-2d.nc");
				}

				bool isWritefreq = hist.isWritefreq();
				BATSet histwriteset = hist.getWriteSet();

				if (bats.getPdfmethod() == PDFMethod::vonMisesKDE)
				{
					isWritefreq = vmkde.isWritefreq();
					histwriteset = vmkde.getWriteSet();
				}
				if (isWritefreq)
				{
					if (histwriteset & workset & BATSet::B1D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_bnd-1d.nc");
					}
					if (histwriteset & workset & BATSet::A1D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_ang-1d.nc");
					}
					if (histwriteset & workset & BATSet::D1D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_tor-1d.nc");
					}
					if (histwriteset & workset & BATSet::BB2D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_bnd-2d.nc");
					}
					if (histwriteset & workset & BATSet::AA2D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_ang-2d.nc");
					}
					if (histwriteset & BATSet::DD2D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_tor-2d.nc");
					}
					if (histwriteset & workset & BATSet::BA2D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_ba-2d.nc");
					}
					if (histwriteset & workset & BATSet::BD2D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_bd-2d.nc");
					}
					if (histwriteset & workset & BATSet::AD2D)
					{
						unwantedfiles.push_back(
							control.getOutfilepath() + "/hist_ad-2d.nc");
					}
				}
			}
			else if (control.isGenConvgdata())
			{

				if (workset & BATSet::B1D)
				{
					requiredfiles.push_back(
						control.getOutfilepath() + "/entcontri_bnd-1d.nc");
				}
				if (workset & BATSet::A1D)
				{
					requiredfiles.push_back(
						control.getOutfilepath() + "/entcontri_ang-1d.nc");
				}
				if (workset & BATSet::D1D)
				{
					requiredfiles.push_back(
						control.getOutfilepath() + "/entcontri_tor-1d.nc");
				}
				if (workset & BATSet::BB2D)
				{
					requiredfiles.push_back(
						control.getOutfilepath() + "/entcontri_bnd-2d.nc");
				}
				if (workset & BATSet::AA2D)
				{
					requiredfiles.push_back(
						control.getOutfilepath() + "/entcontri_ang-2d.nc");
				}
				if (workset & BATSet::DD2D)
				{
					requiredfiles.push_back(
						control.getOutfilepath() + "/entcontri_tor-2d.nc");
				}
				if (workset & BATSet::BA2D)
				{
					requiredfiles.push_back(
						control.getOutfilepath() + "/entcontri_ba-2d.nc");
				}
				if (workset & BATSet::BD2D)
				{
					requiredfiles.push_back(
						control.getOutfilepath() + "/entcontri_bd-2d.nc");
				}
				if (workset & BATSet::AD2D)
				{
					requiredfiles.push_back(
						control.getOutfilepath() + "/entcontri_ad-2d.nc");
				}
			}
		}
		this->isgood_ = isOK;
	}
	else
	{
		this->isgood_ = false;
		std::cout << "InputReadError:: inputfile (" << inpFileName
				  << ") does not exists." << std::endl;
	}
	return isOK;
}

std::ostream &operator<<(std::ostream &os, const Inputs &inp)
{
	int level = 0;
	os
		<< "@@@@@@@@@@@@@@@@@@@@@@@@@  SUPPLIED INPUT PARAMETERS  @@@@@@@@@@@@@@@@@@@@@@@@@@"
		<< std::endl;
	os << std::string(level * 3, ' ') << " [control]" << std::endl;
	level++;
	std::string st = std::string(level * 3, ' ');
	os << st << std::setw(25) << std::left << "infilepath"
	   << " = "
	   << inp.control.getInfilepath() << std::endl;
	os << st << std::setw(25) << std::left << "outfilepath"
	   << " = "
	   << inp.control.getOutfilepath() << std::endl;
	os << st << std::setw(25) << std::left << "infofile"
	   << " = "
	   << inp.control.getInfofile() << std::endl;
	os << st << std::setw(25) << std::left << "nsteps"
	   << " = "
	   << (int)inp.control.getNsteps() << std::endl;
	os << st << std::setw(25) << std::left << "discretize"
	   << " = "
	   << std::boolalpha << inp.control.isDiscetize() << std::endl;
	os << st << std::setw(25) << std::left << "calcentropy"
	   << " = "
	   << std::boolalpha << inp.control.isCalcEntropy() << std::endl;
	os << st << std::setw(25) << std::left << "genconvgdata"
	   << " = "
	   << std::boolalpha << inp.control.isGenConvgdata() << std::endl;
	os << st << std::setw(25) << std::left << "dscrinfofreq"
	   << " = "
	   << (int)inp.control.getDscrinfofreq() << std::endl;
	os << st << std::setw(25) << std::left << "entcinfofreq"
	   << " = "
	   << (int)inp.control.getEntcinfofreq() << std::endl;
	os << st << std::setw(25) << std::left << "cachedscrtdimspercpu"
	   << " = "
	   << inp.control.getCachedscrtdimspercpu() << std::endl;
	os << st << std::setw(25) << std::left << "cacheentcdimspercpu"
	   << " = "
	   << (int)inp.control.getCacheentcdimspercpu() << std::endl;
	level--;
	os << std::endl;

	os << std::string(level * 3, ' ') << " [discretization]" << std::endl;
	level++;
	st = std::string(level * 3, ' ');
	os << st << std::setw(25) << std::left << "nbond"
	   << " = "
	   << inp.getBats().getNbond() << std::endl;
	os << st << std::setw(25) << std::left << "nangle"
	   << " = "
	   << inp.getBats().getNangle() << std::endl;
	os << st << std::setw(25) << std::left << "ndihed"
	   << " = "
	   << inp.getBats().getNdihed() << std::endl;
	int f = 0;
	os << st << std::setw(25) << std::left << "nframe"
	   << " = ";
	for (std::vector<ull_int>::const_iterator it =
			 inp.getBats().getNframes().begin();
		 it != inp.getBats().getNframes().end(); ++it)
	{
		if (f != 0)
			os << ", ";
		os << *it;
		++f;
	}
	os << std::endl;
	os << st << std::setw(25) << std::left << "fname"
	   << " = ";
	f = 0;
	for (std::vector<std::string>::const_iterator it =
			 (inp).getBats().getFnames().begin();
		 it != (inp).getBats().getFnames().end(); ++it)
	{
		if (f != 0)
			os << ", ";
		os << *it;
		++f;
	}
	os << std::endl;
	os << st << std::setw(25) << std::left << "randseed"
	   << " = "
	   << inp.getBats().getRandseed() << std::endl;
	os << st << std::setw(25) << std::left << std::left << "shuffleframes"
	   << " = " << std::boolalpha << (inp).getBats().isShuffleFrames()
	   << std::endl;
	os << st << std::setw(25) << std::left << std::left << "shuffledofs"
	   << " = " << std::boolalpha << (inp).getBats().isShuffleDofs()
	   << std::endl;

	std::string mname;
	if ((inp).getBats().getPdfmethod() == PDFMethod::vonMisesKDE)
	{
		mname = "vonMisesKDE";
		os << st << std::setw(25) << std::left << std::left << "pdfmethod"
		   << " = " << mname << std::endl;
		level--;
		os << std::endl;
		os << std::string(level * 3, ' ') << " [vonmiseskde]" << std::endl;
		level++;
		st = std::string(level * 3, ' ');
		os << st << std::setw(25) << std::left << "writefreq"
		   << " = "
		   << std::boolalpha << std::boolalpha
		   << inp.getVmkde().isWritefreq() << std::endl;
		BATSet writeset = inp.getVmkde().getWriteSet();
		os << st << std::setw(25) << std::left << "writeset"
		   << " = ";
		bool hasmlt = false;
		if (writeset == BATSet::NOTHING)
		{
			os << "NONE" << std::endl;
		}
		else
		{
			if ((writeset & BATSet::XY2D) == BATSet::XY2D)
			{
				os << "2D";
				hasmlt = true;
			}
			else if ((writeset & BATSet::XX2D) == BATSet::XX2D)
			{

				if (writeset & BATSet::BA2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BA2D";
					hasmlt = true;
				}
				if (writeset & BATSet::BD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BD2D";
					hasmlt = true;
				}
				if (writeset & BATSet::AD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "AD2D";
					hasmlt = true;
				}
				if (hasmlt)
					os << ", ";
				os << "XX2D";
				hasmlt = true;
			}
			else if ((writeset & BATSet::X1D) == BATSet::X1D)
			{
				if (writeset & BATSet::BA2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BA2D";
					hasmlt = true;
				}
				if (writeset & BATSet::BD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BD2D";
					hasmlt = true;
				}
				if (writeset & BATSet::AD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "AD2D";
					hasmlt = true;
				}
				if (writeset & BATSet::BB2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BB2D";
					hasmlt = true;
				}
				if (writeset & BATSet::AA2D)
				{
					if (hasmlt)
						os << ", ";
					os << "AA2D";
					hasmlt = true;
				}
				if (writeset & BATSet::DD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "DD2D";
					hasmlt = true;
				}
				if (hasmlt)
					os << ", ";
				os << "1D";
				hasmlt = true;
			}
			else
			{
				if (writeset & BATSet::BA2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BA2D";
					hasmlt = true;
				}
				if (writeset & BATSet::BD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BD2D";
					hasmlt = true;
				}
				if (writeset & BATSet::AD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "AD2D";
					hasmlt = true;
				}
				if (writeset & BATSet::BB2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BB2D";
					hasmlt = true;
				}
				if (writeset & BATSet::AA2D)
				{
					if (hasmlt)
						os << ", ";
					os << "AA2D";
					hasmlt = true;
				}
				if (writeset & BATSet::DD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "DD2D";
					hasmlt = true;
				}
				if (writeset & BATSet::B1D)
				{
					if (hasmlt)
						os << ", ";
					os << "B1D";
					hasmlt = true;
				}
				if (writeset & BATSet::A1D)
				{
					if (hasmlt)
						os << ", ";
					os << "A1D";
					hasmlt = true;
				}
				if (writeset & BATSet::D1D)
				{
					if (hasmlt)
						os << ", ";
					os << "D1D";
					hasmlt = true;
				}
			}
			os << std::endl;
		}
		os << st << std::setw(25) << std::left << "writestepstart"
		   << " = "
		   << (int)inp.getVmkde().getWritestepstart() << std::endl;
		os << st << std::setw(25) << std::left << "writestepstride"
		   << " = "
		   << (int)inp.getVmkde().getWritestepstride() << std::endl;
		os << st << std::setw(25) << std::left << "nmaxconf"
		   << " = "
		   << (int)inp.getVmkde().getNmaxconf() << std::endl;
		os << st << std::setw(25) << std::left << "kappa"
		   << " = "
		   << (double)inp.getVmkde().getKappaValue() << std::endl;
		os << st << std::setw(25) << std::left << "sdosteps"
		   << " = "
		   << (int)inp.getVmkde().getSdoNstep() << std::endl;
		os << st << std::setw(25) << std::left << "sdoiterations"
		   << " = "
		   << inp.getVmkde().getSdoNiter() << std::endl;
		os << st << std::setw(25) << std::left << "sdoconvlimit"
		   << " = "
		   << (double)inp.getVmkde().getSdoConvLimit() << std::endl;
		level--;
		os << std::endl;
	}
	else
	{
		mname = "Histogram";
		os << st << std::setw(25) << std::left << std::left << "pdfmethod"
		   << " = " << mname << std::endl;
		os << std::endl;
		level--;
		os << std::string(level * 3, ' ') << " [histogram]" << std::endl;
		level++;
		st = std::string(level * 3, ' ');
		os << st << std::setw(25) << std::left << "writefreq"
		   << " = "
		   << std::boolalpha << inp.getHist().isWritefreq() << std::endl;
		BATSet writeseth = inp.getHist().getWriteSet();
		os << st << std::setw(25) << std::left << "writeset"
		   << " = ";
		bool hasmlt = false;
		if (writeseth == BATSet::NOTHING)
		{
			os << "NONE" << std::endl;
		}
		else
		{
			if (writeseth == BATSet::XY2D)
			{
				os << "2D";
				hasmlt = true;
			}
			else if ((writeseth & BATSet::XX2D) == BATSet::XX2D)
			{

				if (writeseth & BATSet::BA2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BA2D";
					hasmlt = true;
				}
				if (writeseth & BATSet::BD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BD2D";
					hasmlt = true;
				}
				if (writeseth & BATSet::AD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "AD2D";
					hasmlt = true;
				}
				if (hasmlt)
					os << ", ";
				os << "XX2D";
				hasmlt = true;
			}
			else if ((writeseth & BATSet::X1D) == BATSet::X1D)
			{
				if (writeseth & BATSet::BA2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BA2D";
					hasmlt = true;
				}
				if (writeseth & BATSet::BD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BD2D";
					hasmlt = true;
				}
				if (writeseth & BATSet::AD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "AD2D";
					hasmlt = true;
				}
				if (writeseth & BATSet::BB2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BB2D";
					hasmlt = true;
				}
				if (writeseth & BATSet::AA2D)
				{
					if (hasmlt)
						os << ", ";
					os << "AA2D";
					hasmlt = true;
				}
				if (writeseth & BATSet::DD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "DD2D";
					hasmlt = true;
				}
				if (hasmlt)
					os << ", ";
				os << "1D";
				hasmlt = true;
			}
			else
			{
				if (writeseth & BATSet::BA2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BA2D";
					hasmlt = true;
				}
				if (writeseth & BATSet::BD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BD2D";
					hasmlt = true;
				}
				if (writeseth & BATSet::AD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "AD2D";
					hasmlt = true;
				}
				if (writeseth & BATSet::BB2D)
				{
					if (hasmlt)
						os << ", ";
					os << "BB2D";
					hasmlt = true;
				}
				if (writeseth & BATSet::AA2D)
				{
					if (hasmlt)
						os << ", ";
					os << "AA2D";
					hasmlt = true;
				}
				if (writeseth & BATSet::DD2D)
				{
					if (hasmlt)
						os << ", ";
					os << "DD2D";
					hasmlt = true;
				}
				if (writeseth & BATSet::B1D)
				{
					if (hasmlt)
						os << ", ";
					os << "B1D";
					hasmlt = true;
				}
				if (writeseth & BATSet::A1D)
				{
					if (hasmlt)
						os << ", ";
					os << "A1D";
					hasmlt = true;
				}
				if (writeseth & BATSet::D1D)
				{
					if (hasmlt)
						os << ", ";
					os << "D1D";
					hasmlt = true;
				}
			}
			os << std::endl;
		}
		os << st << std::setw(25) << std::left << "writestepstart"
		   << " = "
		   << (int)inp.getHist().getWritestepstart() << std::endl;
		os << st << std::setw(25) << std::left << "writestepstride"
		   << " = "
		   << (int)inp.getHist().getWritestepstride() << std::endl;
		os << st << std::setw(25) << std::left << "referencenbins"
		   << " = "
		   << (int)inp.getHist().getNbinRef() << std::endl;
		os << st << std::setw(25) << std::left << "nbins"
		   << " = ";
		f = 0;
		for (std::vector<hbin_t>::const_iterator it =
				 (inp).getHist().getBinSchemes().begin();
			 it != (inp).getHist().getBinSchemes().end(); ++it)
		{
			if (f != 0)
				os << ", ";
			os << (int)(*it);
			++f;
		}
		os << std::endl;
		os << std::endl;
		level--;
	}

	os << std::string(level * 3, ' ') << " [entropy]" << std::endl;
	level++;
	st = std::string(level * 3, ' ');
	BATSet workset = inp.getEntropy().getWorkSet();
	os << st << std::setw(25) << std::left << "workset"
	   << " = ";
	bool hasmlt = false;
	if (workset == BATSet::NOTHING)
	{
		os << "NONE" << std::endl;
	}
	else
	{
		if (workset == BATSet::XY2D)
		{
			os << "2D";
			hasmlt = true;
		}
		else if ((workset & BATSet::XX2D) == BATSet::XX2D)
		{

			if (workset & BATSet::BA2D)
			{
				if (hasmlt)
					os << ", ";
				os << "BA2D";
				hasmlt = true;
			}
			if (workset & BATSet::BD2D)
			{
				if (hasmlt)
					os << ", ";
				os << "BD2D";
				hasmlt = true;
			}
			if (workset & BATSet::AD2D)
			{
				if (hasmlt)
					os << ", ";
				os << "AD2D";
				hasmlt = true;
			}
			if (hasmlt)
				os << ", ";
			os << "XX2D";
			hasmlt = true;
		}
		else if ((workset & BATSet::X1D) == BATSet::X1D)
		{
			if (workset & BATSet::BA2D)
			{
				if (hasmlt)
					os << ", ";
				os << "BA2D";
				hasmlt = true;
			}
			if (workset & BATSet::BD2D)
			{
				if (hasmlt)
					os << ", ";
				os << "BD2D";
				hasmlt = true;
			}
			if (workset & BATSet::AD2D)
			{
				if (hasmlt)
					os << ", ";
				os << "AD2D";
				hasmlt = true;
			}
			if (workset & BATSet::BB2D)
			{
				if (hasmlt)
					os << ", ";
				os << "BB2D";
				hasmlt = true;
			}
			if (workset & BATSet::AA2D)
			{
				if (hasmlt)
					os << ", ";
				os << "AA2D";
				hasmlt = true;
			}
			if (workset & BATSet::DD2D)
			{
				if (hasmlt)
					os << ", ";
				os << "DD2D";
				hasmlt = true;
			}
			if (hasmlt)
				os << ", ";
			os << "1D";
			hasmlt = true;
		}
		else
		{
			if (workset & BATSet::BA2D)
			{
				if (hasmlt)
					os << ", ";
				os << "BA2D";
				hasmlt = true;
			}
			if (workset & BATSet::BD2D)
			{
				if (hasmlt)
					os << ", ";
				os << "BD2D";
				hasmlt = true;
			}
			if (workset & BATSet::AD2D)
			{
				if (hasmlt)
					os << ", ";
				os << "AD2D";
				hasmlt = true;
			}
			if (workset & BATSet::BB2D)
			{
				if (hasmlt)
					os << ", ";
				os << "BB2D";
				hasmlt = true;
			}
			if (workset & BATSet::AA2D)
			{
				if (hasmlt)
					os << ", ";
				os << "AA2D";
				hasmlt = true;
			}
			if (workset & BATSet::DD2D)
			{
				if (hasmlt)
					os << ", ";
				os << "DD2D";
				hasmlt = true;
			}
			if (workset & BATSet::B1D)
			{
				if (hasmlt)
					os << ", ";
				os << "B1D";
				hasmlt = true;
			}
			if (workset & BATSet::A1D)
			{
				if (hasmlt)
					os << ", ";
				os << "A1D";
				hasmlt = true;
			}
			if (workset & BATSet::D1D)
			{
				if (hasmlt)
					os << ", ";
				os << "D1D";
				hasmlt = true;
			}
		}
		os << std::endl;
	}

	//	os << st << std::setw(25) << std::left << "shuffleblocks" << " = "
	//			<< std::boolalpha << inp.getEntropy().isShuffleBlocks()
	//			<< std::endl;
	//	os << st << std::setw(25) << std::left << "shuffleblocktimes" << " = "
	//			<< std::boolalpha << inp.getEntropy().getShuffleBlockTimes()
	//			<< std::endl;
	os << st << std::setw(25) << std::left << "usesubset"
	   << " = "
	   << std::boolalpha << (inp).getEntropy().isUseSubset() << std::endl;
	os << st << std::setw(25) << std::left << "useneighbor"
	   << " = "
	   << std::boolalpha << (inp).getEntropy().isUseNeighbor()
	   << std::endl;
	os << st << std::setw(25) << std::left << "jacobian"
	   << " = "
	   << std::boolalpha << inp.getEntropy().isJacobian() << std::endl;

	os << st << std::setw(25) << std::left << "subsetfile"
	   << " = "
	   << inp.getEntropy().getSubsetfile() << std::endl;
	os << st << std::setw(25) << std::left << "neighborfile"
	   << " = "
	   << inp.getEntropy().getNeighborfile() << std::endl;
	std::string tmethod = "MIE";
	if (inp.getMethod() == EntropyMethods::MIST)
	{
		tmethod = "MIST";
	}
	os << st << std::setw(25) << std::left << "scoringmethod"
	   << " = "
	   << tmethod << std::endl;
	os << st << std::setw(25) << std::left << "estimator"
	   << " = ";
	f = 0;
	string ests[4] = {"ML", "MM", "CS", "JS"};
	for (auto est = 0U; est < inp.getEstimators().size(); ++est)
	{
		if (f != 0)
			os << ", ";
		os << ests[inp.getEstimators()[est]];
		++f;
	}
	os << std::endl;
	ull_int nframe_bats = 0;
	for (auto x : inp.getBats().getNframes())
	{
		nframe_bats += x;
	}
	if (inp.getEntropy().getNumframe() != nframe_bats || (inp.getEntropy().getStrideframe() != 1) || (inp.getEntropy().getStartframe() != 0))
	{
		os << st << std::setw(25) << std::left << "startframe"
		   << " = "
		   << inp.getEntropy().getStartframe() << std::endl;
		os << st << std::setw(25) << std::left << "strideframe"
		   << " = "
		   << inp.getEntropy().getStrideframe() << std::endl;
		os << st << std::setw(25) << std::left << "numframe"
		   << " = "
		   << inp.getEntropy().getNumframe() << std::endl;
	}
	level--;
	os << std::endl;
	os << std::string(80, '@') << std::endl;
	return os;
}
