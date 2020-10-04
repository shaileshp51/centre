/* 
 * File:   main.cpp
 * Author: shailesh
 *
 * Created on 5 August, 2016, 10:36 PM
 */
#include "configcentre.h"

#ifdef USE_OMPMPI
#include "main_common.h"

#include "Discretizer.h"
#include "EntropyCalculator.h"
#include "EntropyScorer.h"
#include <sys/types.h>
#include <unistd.h>
#include <algorithm>
#include <omp.h>
#include <mpi.h>
#include "MpiCommTags.h"

using namespace std;

int main(int argc, char **argv) {

	Timer total_time;

	bool overwritefiles = false;

	char *inpfilename;
	Inputs inp = Inputs(string(""));
	bool validate_inputs = true;

	int rank, namelen;
	char hostName[MPI_MAX_PROCESSOR_NAME];
	int rc = 0; // iam = 0, np = 1
	int numprocs;
	int n_thread_perproc = 1;
	//int tmpInt;
	int thread_level_provided, thread_level_claimed;

	MPI_Init_thread(&argc, &argv , MPI_THREAD_SINGLE, &thread_level_provided);
	// MPI_Init(&argc, &argv); // , MPI_THREAD_SERIALIZED, &thread_level_provided
	MPI_Query_thread(&thread_level_claimed);
	if (thread_level_claimed != thread_level_provided) {
		mprinterr(
				"Claimed thread level=%d, but got thread level=%d, aborting\n",
				thread_level_claimed, thread_level_provided);
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Get_processor_name(hostName, &namelen);

	// Input validation
	if (validate_inputs) {
		if (cmdOptionExists(argv, argv + argc, "-h")
				|| cmdOptionExists(argv, argv + argc, "--help")) {
			if (rank == MASTER_PROC) {
				printUsage();
				MPI_Abort(MPI_COMM_WORLD, rc);
			}
		}
		if (cmdOptionExists(argv, argv + argc, "-s")
				|| cmdOptionExists(argv, argv + argc, "--sample")) {
			if (rank == MASTER_PROC) {
				printSample();
				MPI_Abort(MPI_COMM_WORLD, rc);
			}
		}
		if (cmdOptionExists(argv, argv + argc, "-O")) {
			overwritefiles = true;
		}
		inpfilename = getCmdOption(argv, argv + argc, "-i");
		// TODO: All validations should be done here
	}

	if (inpfilename) {
		inp.Init(string(inpfilename));
	} else if (rank == 0) {
		cerr << "-i inputfile is a required option" << endl;
		exit(0);
	}
	#pragma omp parallel
        #pragma omp master
        {
           n_thread_perproc = omp_get_num_threads();
        }
	MPI_Barrier(MPI_COMM_WORLD);
	cout << endl;
	cout << "MPI PROCESS: rank[" << rank << "], NodeName: " << hostName
			<< ", UsedThreadsOfHost: " << n_thread_perproc << ", TotalThreadsOnHost: "
			<< omp_get_num_procs() << endl;
	MPI_Barrier(MPI_COMM_WORLD);

#pragma omp parallel
	{
		int ID = omp_get_thread_num();
		printf(
				"Initiated process = %i of %i MPI processes on machine=%s, thread number %i\n",
				rank, numprocs, hostName, ID);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	if (rank == MASTER_PROC) {
		if (!inp.isGood()) {
			// Input read error, abort program
			MPI_Abort(MPI_COMM_WORLD, 0);
		} else {
			total_time.start();
			
			printHeader();
			cout << inp << endl;
			cout << endl << "Master pid: " << getpid() << " ppid: " << getppid()
					<< endl << endl;
			ofstream fopid(
					inp.getControl().getOutfilepath() + "/.master-pid.txt",
					ios::out);
			fopid << getpid();
			fopid.flush();
			fopid.close();

			// TODO: Check if output files already exist and overwrite flag is missing if so warn and abort
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
							MPI_Abort(MPI_COMM_WORLD, 0);
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
				MPI_Abort(MPI_COMM_WORLD, 0);
			}
		}
	}

	const vector<hbin_t> bin_schemes = inp.getHist().getBinSchemes();
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
	ull_int nfrm4rl2int = tmp_frm; 
	const u_int n_bonds = inp.getBats().getNbond();
	const u_int n_angles = inp.getBats().getNangle();
	const u_int n_diheds = inp.getBats().getNdihed();
	const u_int n_bnd_eff = inp.getSubset().bondsSize();
	const u_int n_ang_eff = inp.getSubset().anglesSize();
	const u_int n_dih_eff = inp.getSubset().torsionsSize();
	
	u_int n_bnd2d = 0, n_ang2d = 0, n_dih2d = 0, n_ba2d = 0, n_bd2d = 0,
			n_ad2d = 0;
	if (inp.getEntropy().getWorkSet() & BATSet::BB2D) {
		vector<u_int> bndKeys;
		inp.getNeighbors().bondKeys(bndKeys);
		const u_int num_bndKeys = bndKeys.size();
		for (u_int d1 = 0U; d1 < num_bndKeys; ++d1) {
			n_bnd2d += inp.getNeighbors().bondNeighSize(bndKeys[d1]);
		}
	}
	if (inp.getEntropy().getWorkSet() & BATSet::AA2D) {
		vector<u_int> angKeys;
		inp.getNeighbors().angleKeys(angKeys);
		const u_int num_angKeys = angKeys.size();
		for (u_int d1 = 0U; d1 < num_angKeys; ++d1) {
			n_ang2d += inp.getNeighbors().angleNeighSize(angKeys[d1]);
		}
	}
	if (inp.getEntropy().getWorkSet() & BATSet::DD2D) {
		vector<u_int> dihKeys;
		inp.getNeighbors().torsionKeys(dihKeys);
		const u_int num_dihKeys = dihKeys.size();
		for (u_int d1 = 0U; d1 < num_dihKeys; ++d1) {
			n_dih2d += inp.getNeighbors().torsionNeighSize(dihKeys[d1]);
		}
	}
	if (inp.getEntropy().getWorkSet() & BATSet::BA2D) {
		vector<u_int> baKeys;
		inp.getNeighbors().bacrossKeys(baKeys);
		const u_int num_baKeys = baKeys.size();
		for (u_int d1 = 0U; d1 < num_baKeys; ++d1) {
			n_ba2d += inp.getNeighbors().baNeighSize(baKeys[d1]);
		}
	}
	if (inp.getEntropy().getWorkSet() & BATSet::BD2D) {
		vector<u_int> bdKeys;
		inp.getNeighbors().bdcrossKeys(bdKeys);
		const u_int num_bdKeys = bdKeys.size();
		for (u_int d1 = 0U; d1 < num_bdKeys; ++d1) {
			n_bd2d += inp.getNeighbors().bdNeighSize(bdKeys[d1]);
		}
	}
	if (inp.getEntropy().getWorkSet() & BATSet::AD2D) {
		vector<u_int> adKeys;
		inp.getNeighbors().adcrossKeys(adKeys);
		const u_int num_adKeys = adKeys.size();
		for (u_int d1 = 0U; d1 < num_adKeys; ++d1) {
			n_ad2d += inp.getNeighbors().adNeighSize(adKeys[d1]);
		}
	}

	vector<ull_int> shuffle_index;
	vector<hbin_t> shuffle_blocks;

	ProgState prg_state;

	if (rank == MASTER_PROC) {
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
	}
	if (inp.getBats().isShuffleFrames()
			|| (inp.getEntropy().isShuffleBlocks()
					&& inp.getEntropy().getShuffleBlockTimes() > 1)) {
		if (rank == MASTER_PROC) {
			random_device rd;
			auto seed = rd();
			if (inp.getBats().getRandseed() != 0) {
				seed = inp.getBats().getRandseed() + 1000;
			}
			mt19937 rnd_gen(seed);
			if (inp.getBats().isShuffleFrames()) {
				shuffle_index.resize(nfrm4rl2int);
				for (ull_int i = 0; i < nfrm4rl2int; ++i) {
					shuffle_index[i] = i;
				}
				shuffle(shuffle_index.begin(), shuffle_index.end(), rnd_gen);
				ofstream ofp(
						inp.getControl().getOutfilepath()
								+ "shuffle_indices.bin",
						ios::out | ios::binary);
				ofp.write(reinterpret_cast<const char*>(shuffle_index.data()),
						shuffle_index.size() * sizeof(shuffle_index[0]));
				ofp.close();

				MPI_Bcast(shuffle_index.data(), nfrm4rl2int,
				MPI_UNSIGNED_LONG_LONG, MASTER_PROC, MPI_COMM_WORLD);
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
					shuffle(shuffle_block.begin(), shuffle_block.end(),
							rnd_gen);
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
				MPI_Bcast(shuffle_blocks.data(),
						nsteps * inp.getEntropy().getShuffleBlockTimes(),
						MPI_UNSIGNED_CHAR, MASTER_PROC, MPI_COMM_WORLD);
			}
		} else if (rank != 0) {
			if (inp.getBats().isShuffleFrames()) {
				shuffle_index.resize(nfrm4rl2int);
				MPI_Bcast(shuffle_index.data(), nfrm4rl2int,
				MPI_UNSIGNED_LONG_LONG, MASTER_PROC, MPI_COMM_WORLD);
			}
			if (inp.getEntropy().isShuffleBlocks()
					&& inp.getEntropy().getShuffleBlockTimes() > 1) {
				shuffle_blocks.resize(
						nsteps * inp.getEntropy().getShuffleBlockTimes());
				MPI_Bcast(shuffle_blocks.data(),
						nsteps * inp.getEntropy().getShuffleBlockTimes(),
						MPI_UNSIGNED_CHAR, MASTER_PROC, MPI_COMM_WORLD);
			}
		}

	} else if (inp.getBats().isShuffleDofs()) {
		shuffle_index.resize(nfrm4rl2int);
		for (ull_int ir = 0; ir < nfrm4rl2int; ++ir) {
			shuffle_index[ir] = ir;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	Discretizer discretizer(inp, prg_state);
	MPI_Barrier(MPI_COMM_WORLD);
	discretizer.run_mpi(shuffle_index, rank, numprocs, n_thread_perproc);

	MPI_Barrier(MPI_COMM_WORLD);
	EntropyCalculator entcalculator(inp, prg_state);
	MPI_Barrier(MPI_COMM_WORLD);

	entcalculator.run_mpi(rank, numprocs, n_thread_perproc);

	MPI_Barrier(MPI_COMM_WORLD);
	EntropyScorer entropyscores(inp);
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == MASTER_PROC)
		entropyscores.run();

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == MASTER_PROC) {
		printFooter();
		total_time.stop();
		mprintf("TIME: Total execution: %.4f seconds.\n", total_time.total());
		cout << "Execution completed successfully..." << endl;
	}
	MPI_Finalize();
}

#endif
