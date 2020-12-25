/* 
 * File:   main.cpp
 * Author: shailesh
 *
 * Created on 5 August, 2016, 10:36 PM
 */

#include "main_common.h"

#ifndef USE_OMPMPI

#include "Discretizer.h"
#include "EntropyCalculator.h"
#include "EntropyScorer.h"
#include <sys/types.h>
#include <unistd.h>
#include <algorithm>
#include <omp.h>

using namespace std;

int main(int argc, char **argv) {

	Timer total_time;

	bool overwritefiles = false;

	char *inpfilename;
	Inputs inp = Inputs(string(""));
	bool validate_inputs = true;
	int num_threads = omp_get_max_threads();
	// Input validation
	if (validate_inputs) {
		if (cmdOptionExists(argv, argv + argc, "-h")
				|| cmdOptionExists(argv, argv + argc, "--help")) {

				printUsage();
				exit(0);

		}
		if (cmdOptionExists(argv, argv + argc, "-s")
				|| cmdOptionExists(argv, argv + argc, "--sample")) {

				printSample();
				exit(0);

		}
		if (cmdOptionExists(argv, argv + argc, "-c")
				|| cmdOptionExists(argv, argv + argc, "--cite")) {

				printFooter();
				exit(0);

		}
		if (cmdOptionExists(argv, argv + argc, "-O")) {
			overwritefiles = true;
		}
		inpfilename = getCmdOption(argv, argv + argc, "-i");
		// TODO: All validations should be done here
		if (inpfilename) {
			inp.Init(string(inpfilename));
		} else {
			cerr << "-i inputfile is a required option" << endl;
			exit(0);
		}
		if (!inp.isGood()) {
			// Input read error, abort program
			exit(0);
		}

		total_time.start();
		printHeader();
		cout << inp << endl;

		cout << "Program started using " << num_threads << " threads.";
		cout << endl << "Master pid: " << getpid() << " ppid: " << getppid()
				<< endl << endl;
		ofstream fopid(inp.getControl().getOutfilepath() + "/.master-pid.txt",
				ios::out);
		fopid << getpid();
		fopid.flush();
		fopid.close();

		bool fileexistanceError = false;
		for (auto rqrdfile : inp.getRequiredfiles()) {
			if (!Utils::fileExists(rqrdfile)) {
				cerr << "Error:: required-file (" << rqrdfile
						<< ") doen't exist." << endl;
				fileexistanceError = true;
			}
		}
		for (auto unwntfile : inp.getUnwantedfiles()) {
			if (Utils::fileExists(unwntfile)) {
				if (!overwritefiles) {
					cerr << "Error:: file-to-write (" << unwntfile
							<< ") already exist." << endl;
					fileexistanceError = true;
				} else {
					if (remove(unwntfile.c_str())) {
						cerr << "Error removing file (" << unwntfile << ") "
								<< endl;
						exit(0);
					} else {
						cout << "deleted: " << unwntfile.c_str() << endl;
					}
				}
			}
		}

		if (fileexistanceError) {
			cerr
					<< "To force overwriting existing files use -O option with program"
					<< endl;
			exit(0);
		}
	}

	const vector<hbin_t> bin_schemes = inp.getHist().getBinSchemes();
	//hbin_t n_schemes = (hbin_t) bin_schemes.size();
	u_int bin_schemes_sum = 0;
	ull_int nbins_sum_sq = 0;

	for (size_t i = 0; i < bin_schemes.size(); ++i) {
		bin_schemes_sum += bin_schemes[i];
		nbins_sum_sq += bin_schemes[i] * bin_schemes[i];
	}

	const hbin_t nsteps = inp.getControl().getNsteps();

	ull_int tmp_frm = 0;
	for (auto v : inp.getBats().getNframes()) {
		tmp_frm += v;
	}
	ull_int nfrm4rl2int = tmp_frm; // inp.getEntropy().getNumframe();

	const u_int n_bonds = inp.getBats().getNbond();
	const u_int n_angles = inp.getBats().getNangle();
	const u_int n_diheds = inp.getBats().getNdihed();
	const u_int n_bnd_eff = inp.getSubset().bondsSize();
	const u_int n_ang_eff = inp.getSubset().anglesSize();
	const u_int n_dih_eff = inp.getSubset().torsionsSize();

	u_int n_bnd2d = 0, n_ang2d = 0, n_dih2d = 0, n_ba2d = 0, n_bd2d = 0,
			n_ad2d = 0;
	if (inp.getEntropy().getWorkSet() & BATSet::BB2D) {
		vector<u_int> bnd_keys_v;
		inp.getNeighbors().bondKeys(bnd_keys_v);
		const u_int n_bnd_keys = bnd_keys_v.size();
		for (size_t d1 = 0; d1 < n_bnd_keys; ++d1) {
			n_bnd2d += inp.getNeighbors().bondNeighSize(bnd_keys_v[d1]);
		}
	}
	if (inp.getEntropy().getWorkSet() & BATSet::AA2D) {
		vector<u_int> angKeys;
		inp.getNeighbors().angleKeys(angKeys);
		const u_int num_angKeys = angKeys.size();
		for (size_t d1 = 0; d1 < num_angKeys; ++d1) {
			n_ang2d += inp.getNeighbors().angleNeighSize(angKeys[d1]);
		}
	}
	if (inp.getEntropy().getWorkSet() & BATSet::DD2D) {
		vector<u_int> dihKeys;
		inp.getNeighbors().torsionKeys(dihKeys);
		const u_int num_dihKeys = dihKeys.size();
		for (size_t d1 = 0; d1 < num_dihKeys; ++d1) {
			n_dih2d += inp.getNeighbors().torsionNeighSize(dihKeys[d1]);
		}
	}
	if (inp.getEntropy().getWorkSet() & BATSet::BA2D) {
		vector<u_int> baKeys;
		inp.getNeighbors().bacrossKeys(baKeys);
		const u_int num_baKeys = baKeys.size();
		for (size_t d1 = 0; d1 < num_baKeys; ++d1) {
			n_ba2d += inp.getNeighbors().baNeighSize(baKeys[d1]);
		}
	}
	if (inp.getEntropy().getWorkSet() & BATSet::BD2D) {
		vector<u_int> bdKeys;
		inp.getNeighbors().bdcrossKeys(bdKeys);
		const u_int num_bdKeys = bdKeys.size();
		for (size_t d1 = 0; d1 < num_bdKeys; ++d1) {
			n_bd2d += inp.getNeighbors().bdNeighSize(bdKeys[d1]);
		}
	}
	if (inp.getEntropy().getWorkSet() & BATSet::AD2D) {
		vector<u_int> adKeys;
		inp.getNeighbors().adcrossKeys(adKeys);
		const u_int num_adKeys = adKeys.size();
		for (size_t d1 = 0; d1 < num_adKeys; ++d1) {
			n_ad2d += inp.getNeighbors().adNeighSize(adKeys[d1]);
		}
	}

	vector<ull_int> shuffle_index;
	vector<hbin_t> shuffle_blocks;
	random_device rd;
	auto seed = rd();
	if (inp.getBats().getRandseed() != 0) {
		seed = inp.getBats().getRandseed() + 1000;
	}
	mt19937 rnd_gen(seed);

	ProgState prg_state;

	if (inp.getControl().isDiscetize()) {
		prg_state.dscrt.active = true;
		if ((inp.getEntropy().getWorkSet() & BATSet::B1D)
				|| (inp.getEntropy().getWorkSet() & BATSet::BB2D)) {
			prg_state.dscrt_b.active = true;
			prg_state.dscrt_b.cstt = ExecState::WAITING;
			if (!prg_state.dscrt.active) {
				prg_state.dscrt.cstt = ExecState::WAITING;
			}
			prg_state.dscrt_b.tot_tasks = n_bonds;
			prg_state.dscrt.tot_tasks += n_bonds;
		}
		if ((inp.getEntropy().getWorkSet() & BATSet::A1D)
				|| (inp.getEntropy().getWorkSet() & BATSet::AA2D)) {
			prg_state.dscrt_a.active = true;
			prg_state.dscrt_a.cstt = ExecState::WAITING;
			if (!prg_state.dscrt.active) {
				prg_state.dscrt.cstt = ExecState::WAITING;
			}
			prg_state.dscrt_a.tot_tasks = n_angles;
			prg_state.dscrt.tot_tasks += n_angles;
		}
		if ((inp.getEntropy().getWorkSet() & BATSet::D1D)
				|| (inp.getEntropy().getWorkSet() & BATSet::DD2D)) {
			prg_state.dscrt_d.active = true;
			prg_state.dscrt_d.cstt = ExecState::WAITING;
			if (!prg_state.dscrt.active) {
				prg_state.dscrt.cstt = ExecState::WAITING;
			}
			prg_state.dscrt_d.tot_tasks = n_diheds;
			prg_state.dscrt.tot_tasks += n_diheds;
		}
	}

	if (inp.getControl().isCalcEntropy()) {
		bool entc_activate = false;
		if ((inp.getEntropy().getWorkSet() & BATSet::B1D)
				|| (inp.getEntropy().getWorkSet() & BATSet::BB2D)) {
			prg_state.entc_b1d.active = true;
			prg_state.entc_b1d.cstt = ExecState::WAITING;
			prg_state.entc_b1d.tot_tasks = n_bnd_eff;
			prg_state.entc.tot_tasks += n_bnd_eff;
			entc_activate = true;
		}
		if ((inp.getEntropy().getWorkSet() & BATSet::A1D)
				|| (inp.getEntropy().getWorkSet() & BATSet::AA2D)) {
			prg_state.entc_a1d.active = true;
			prg_state.entc_a1d.cstt = ExecState::WAITING;
			prg_state.entc_a1d.tot_tasks = n_ang_eff;
			prg_state.entc.tot_tasks += n_ang_eff;
			entc_activate = true;
		}
		if ((inp.getEntropy().getWorkSet() & BATSet::D1D)
				|| (inp.getEntropy().getWorkSet() & BATSet::DD2D)) {
			prg_state.entc_d1d.active = true;
			prg_state.entc_d1d.cstt = ExecState::WAITING;
			prg_state.entc_d1d.tot_tasks = n_dih_eff;
			prg_state.entc.tot_tasks += n_dih_eff;
			entc_activate = true;
		}

		if (inp.getEntropy().getWorkSet() & BATSet::BB2D) {
			prg_state.entc_b2d.active = true;
			prg_state.entc_b2d.cstt = ExecState::WAITING;
			prg_state.entc_b2d.tot_tasks = n_bnd2d;
			prg_state.entc.tot_tasks += n_bnd2d;
			entc_activate = true;
		}
		if (inp.getEntropy().getWorkSet() & BATSet::AA2D) {
			prg_state.entc_a2d.active = true;
			prg_state.entc_a2d.cstt = ExecState::WAITING;
			prg_state.entc_a2d.tot_tasks = n_ang2d;
			prg_state.entc.tot_tasks += n_ang2d;
			entc_activate = true;
		}
		if (inp.getEntropy().getWorkSet() & BATSet::DD2D) {
			prg_state.entc_d2d.active = true;
			prg_state.entc_d2d.cstt = ExecState::WAITING;
			prg_state.entc_d2d.tot_tasks = n_dih2d;
			prg_state.entc.tot_tasks += n_dih2d;
			entc_activate = true;
		}
		if (inp.getEntropy().getWorkSet() & BATSet::BA2D) {
			prg_state.entc_ba2d.active = true;
			prg_state.entc_ba2d.cstt = ExecState::WAITING;
			prg_state.entc_ba2d.tot_tasks = n_ba2d;
			prg_state.entc.tot_tasks += n_ba2d;
			entc_activate = true;
		}
		if (inp.getEntropy().getWorkSet() & BATSet::BD2D) {
			prg_state.entc_bd2d.active = true;
			prg_state.entc_bd2d.cstt = ExecState::WAITING;
			prg_state.entc_bd2d.tot_tasks = n_bd2d;
			prg_state.entc.tot_tasks += n_bd2d;
			entc_activate = true;
		}
		if (inp.getEntropy().getWorkSet() & BATSet::AD2D) {
			prg_state.entc_ad2d.active = true;
			prg_state.entc_ad2d.cstt = ExecState::WAITING;
			prg_state.entc_ad2d.tot_tasks = n_ad2d;
			prg_state.entc.tot_tasks += n_ad2d;
			entc_activate = true;
		}
		prg_state.entc.active = entc_activate;
		prg_state.entc.cstt = ExecState::WAITING;
	}

	if (inp.getBats().isShuffleFrames()
			|| (inp.getEntropy().isShuffleBlocks()
					&& inp.getEntropy().getShuffleBlockTimes() > 1)) {
		if (inp.getBats().isShuffleFrames()) {
			shuffle_index.resize(nfrm4rl2int);
			for (ull_int i = 0; i < nfrm4rl2int; ++i) {
				shuffle_index[i] = i;
			}
			shuffle(shuffle_index.begin(), shuffle_index.end(), rnd_gen);
			ofstream ofp(
					inp.getControl().getOutfilepath() + "shuffle_indices.bin",
					ios::out | ios::binary);
			ofp.write(reinterpret_cast<const char*>(shuffle_index.data()),
					shuffle_index.size() * sizeof(shuffle_index[0]));
			ofp.close();
		}
		if (inp.getEntropy().isShuffleBlocks()
				&& inp.getEntropy().getShuffleBlockTimes() > 1) {
			shuffle_blocks.resize(
					nsteps * inp.getEntropy().getShuffleBlockTimes());
			for (size_t j = 0; j < inp.getEntropy().getShuffleBlockTimes();
					++j) {
				vector<hbin_t> shuffle_block(nsteps);
				for (size_t i = 0; i < nsteps; ++i) {
					shuffle_block[i] = i;
				}
				shuffle(shuffle_block.begin(), shuffle_block.end(), rnd_gen);
				for (size_t i = 0; i < nsteps; ++i) {
					shuffle_blocks[j * nsteps + i] = shuffle_block[i];
				}
			}

			ofstream ofp(
					inp.getControl().getOutfilepath()
							+ "/shuffle_stepblock_indices.bin",
					ios::out | ios::binary);
			ofp.write(reinterpret_cast<const char*>(shuffle_blocks.data()),
					shuffle_blocks.size() * sizeof(shuffle_blocks[0]));
			ofp.close();
		}
	} else if (inp.getBats().isShuffleDofs()) {
		shuffle_index.resize(nfrm4rl2int);
		for (ull_int ir = 0; ir < nfrm4rl2int; ++ir) {
			shuffle_index[ir] = ir;
		}
	}

	Discretizer discretizer(inp, prg_state);
	discretizer.run(shuffle_index, num_threads);

	EntropyCalculator entcalculator(inp, prg_state);
	entcalculator.run(num_threads);

	EntropyScorer entropyscores(inp);
	entropyscores.run();

	printFooter();
	total_time.stop();
	mprintf("TIME: Total execution: %.4f seconds.\n", total_time.total());
	cout << "Execution completed successfully..." << endl;

}

#endif
