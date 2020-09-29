#include "main_common.h"

using namespace std;

char* getCmdOption(char **begin, char **end, const string &option) {
	char **itr = find(begin, end, option);
	if (itr != end && ++itr != end) {
		return *itr;
	}
	return 0;
}

bool cmdOptionExists(char **begin, char **end, const string &option) {
	return find(begin, end, option) != end;
}

string to_durration_string(int64_t dur, bool show_ms) {
	auto qr = div((long)dur, 1000L);
	auto ms = qr.rem;
	qr = div(qr.quot, 60L);
	auto s = qr.rem;
	qr = div(qr.quot, 60L);
	auto m = qr.rem;
	auto h = qr.quot;
	ostringstream result;
	if (h > 0) {
		result << to_string(h) << 'h' << setw(2) << setfill('0') << m << 'm'
				<< setw(2) << setfill('0') << s;
		if (show_ms) {
			result << '.' << setw(3) << setfill('0') << ms;
		}
		result << 's';
	} else if (m > 0) {
		result << m << 'm' << setw(2) << setfill('0') << s;
		if (show_ms) {
			result << '.' << setw(3) << setfill('0') << ms;
		}
		result << 's';
	} else {
		result << s;
		if (show_ms) {
			result << '.' << setw(3) << setfill('0') << ms;
		}
		result << 's';
	}
	return result.str();
}

string ProgTaskSet::toString() {
	string chr_state[5] = { "Disabled", "Waiting", "Running", "", "Completed" };
	stringstream os;
	os << right << setw(8) << this->name;
	if (this->cstt > ExecState::DISABLED) {
		os << setw(10) << chr_state[this->cstt];
		u_int nfills = 15 * 2 + 19;
#ifdef DETAILED_TIMING
	      nfills = 15 * 6 + 19;
	#endif

		long time_left;
		if (this->cstt > ExecState::WAITING) {
			char tmp_time_str[24];
			time_t strt_c = Clock::to_time_t(
					this->start);
			strftime(tmp_time_str, sizeof(tmp_time_str), "%F %T",
					localtime(&strt_c));
			os << " " << setw(19) << tmp_time_str;
#ifdef DETAILED_TIMING
	         auto rd_time = chrono::duration_cast<chrono::milliseconds>(
	        		 this->curr_read - this->strt_read).count();
	         auto ex_time = chrono::duration_cast<chrono::milliseconds>(
	        		 this->curr_comp - this->strt_comp).count();
	         auto cm_time = chrono::duration_cast<chrono::milliseconds>(
	        		 this->curr_comm - this->strt_comm).count();
	         auto wr_time = chrono::duration_cast<chrono::milliseconds>(
	        		 this->curr_write - this->strt_write).count();

	         os << setw(15) << to_durration_string(rd_time) << setw(15)
	               << to_durration_string(ex_time) << setw(15) << to_durration_string(cm_time)
	               << setw(15) << to_durration_string(wr_time);
	#endif
			auto tot_time = chrono::duration_cast<chrono::milliseconds>(
					this->current - this->start).count();
			if (this->done_tasks > 0) {
				auto tsk_time = tot_time / this->done_tasks;
				os << setw(15) << to_durration_string(tsk_time);
				auto left_time = (this->tot_tasks - this->done_tasks)
						* tsk_time;
				time_left = left_time;
			} else {
				os << " ";
			}
			os << setw(15) << to_durration_string(tot_time);
		}

		if (this->cstt == ExecState::WAITING) {
			os << setw(nfills + 1) << setfill(' ') << " ";
		}
		os << setw(12) << this->done_tasks << setw(12) << this->tot_tasks;
		float perc = this->done_tasks * 100.0 / this->tot_tasks;
		os << right << fixed << setw(8) << setprecision(2) << perc;

		if (this->cstt > ExecState::WAITING) {
			os << setw(13) << to_durration_string(time_left, false);
		}
		os << endl;
	}

	return os.str();
}

ostream& operator<<(ostream &os, const ProgTaskSet &elm) {
	string chr_state[4] = { "Disabled", "Waiting", "Running", "Completed" };

	os << right << setw(8) << elm.name;
	if (elm.cstt > ExecState::DISABLED) {
		os << setw(10) << chr_state[elm.cstt];
		u_int nfills = 15 * 2 + 19;
#ifdef DETAILED_TIMING
      nfills = 15 * 6 + 19;
#endif

		long time_left;
		if (elm.cstt > ExecState::WAITING) {
			char tmp_time_str[24];
			time_t strt_c = Clock::to_time_t(elm.start);
			strftime(tmp_time_str, sizeof(tmp_time_str), "%F %T",
					localtime(&strt_c));
			os << " " << setw(19) << tmp_time_str;
#ifdef DETAILED_TIMING
         auto rd_time = chrono::duration_cast<chrono::milliseconds>(
               elm.curr_read - elm.strt_read).count();
         auto ex_time = chrono::duration_cast<chrono::milliseconds>(
               elm.curr_comp - elm.strt_comp).count();
         auto cm_time = chrono::duration_cast<chrono::milliseconds>(
               elm.curr_comm - elm.strt_comm).count();
         auto wr_time = chrono::duration_cast<chrono::milliseconds>(
               elm.curr_write - elm.strt_write).count();

         os << setw(15) << to_durration_string(rd_time) << setw(15)
               << to_durration_string(ex_time) << setw(15) << to_durration_string(cm_time)
               << setw(15) << to_durration_string(wr_time);
#endif
			auto tot_time = chrono::duration_cast<chrono::milliseconds>(
					elm.current - elm.start).count();
			if (elm.done_tasks > 0) {
				auto tsk_time = tot_time / elm.done_tasks;
				os << setw(15) << to_durration_string(tsk_time);
				auto left_time = (elm.tot_tasks - elm.done_tasks) * tsk_time;
				time_left = left_time;
			} else {
				os << " ";
			}
			os << setw(15) << to_durration_string(tot_time);
		}

		if (elm.cstt == ExecState::WAITING) {
			os << setw(nfills + 1) << setfill(' ') << " ";
		}
		os << setw(12) << elm.done_tasks << setw(12) << elm.tot_tasks;
		float perc = elm.done_tasks * 100.0 / elm.tot_tasks;
		os << right << fixed << setw(8) << setprecision(2) << perc;

		if (elm.cstt > ExecState::WAITING) {
			os << setw(13) << to_durration_string(time_left, false);
		}
		os << endl;
	}

	return os;
}

string ProgState::toString() {
	stringstream os;
	if (this->dscrt.active || this->entc.active) {
		os << right << setw(8) << "TaskSet" << setw(10) << "State";
		os << setw(20) << "Started On";
#ifdef DETAILED_TIMING
	   os << setw(15) << "ReadTime" << setw(15) << "ExecuteTime";
	   os << setw(15) << "I.P.C.Time" << setw(15) << "WriteTime";
#endif
		os << setw(15) << "AvgTimeTask" << setw(15) << "ElapsedTime" << setw(12)
				<< "DoneTasks" << setw(12) << "TotalTasks" << setw(8)
				<< "Done[%]" << setw(13) << "ExpcTimeLeft" << endl;
	}
	if (this->dscrt.active) {
		os << this->dscrt.toString();
		if (this->dscrt_b.active) {
			os << this->dscrt_b.toString();
		}
		if (this->dscrt_a.active) {
			os << this->dscrt_a.toString();
		}
		if (this->dscrt_d.active) {
			os << this->dscrt_d.toString();
		}
	}
	if (this->entc.active) {
		os << this->entc.toString();
		if (this->entc_b1d.active) {
			os << this->entc_b1d.toString();
		}
		if (this->entc_a1d.active) {
			os << this->entc_a1d.toString();
		}
		if (this->entc_d1d.active) {
			os << this->entc_d1d.toString();
		}
		if (this->entc_b2d.active) {
			os << this->entc_b2d.toString();
		}
		if (this->entc_a2d.active) {
			os << this->entc_a2d.toString();
		}
		if (this->entc_d2d.active) {
			os << this->entc_d2d.toString();
		}
		if (this->entc_ba2d.active) {
			os << this->entc_ba2d.toString();
		}
		if (this->entc_bd2d.active) {
			os << this->entc_bd2d.toString();
		}
		if (this->entc_ad2d.active) {
			os << this->entc_ad2d.toString();
		}
	}
	os << endl;
	return os.str();
}

ostream& operator<<(ostream &os, const ProgState &stt) {

	if (stt.dscrt.active || stt.entc.active) {
		os << right << setw(8) << "TaskSet" << setw(10) << "State";
		os << setw(20) << "Started On";
#ifdef DETAILED_TIMING
	   os << setw(15) << "ReadTime" << setw(15) << "ExecuteTime";
	   os << setw(15) << "I.P.C.Time" << setw(15) << "WriteTime";
#endif
		os << setw(15) << "AvgTimeTask" << setw(15) << "ElapsedTime" << setw(12)
				<< "DoneTasks" << setw(12) << "TotalTasks" << setw(8)
				<< "Done[%]" << setw(13) << "ExpcTimeLeft" << endl;
	}
	if (stt.dscrt.active) {
		os << stt.dscrt;
		if (stt.dscrt_b.active) {
			os << stt.dscrt_b;
		}
		if (stt.dscrt_a.active) {
			os << stt.dscrt_a;
		}
		if (stt.dscrt_d.active) {
			os << stt.dscrt_d;
		}
	}
	if (stt.entc.active) {
		os << stt.entc;
		if (stt.entc_b1d.active) {
			os << stt.entc_b1d;
		}
		if (stt.entc_a1d.active) {
			os << stt.entc_a1d;
		}
		if (stt.entc_d1d.active) {
			os << stt.entc_d1d;
		}
		if (stt.entc_b2d.active) {
			os << stt.entc_b2d;
		}
		if (stt.entc_a2d.active) {
			os << stt.entc_a2d;
		}
		if (stt.entc_d2d.active) {
			os << stt.entc_d2d;
		}
		if (stt.entc_ba2d.active) {
			os << stt.entc_ba2d;
		}
		if (stt.entc_bd2d.active) {
			os << stt.entc_bd2d;
		}
		if (stt.entc_ad2d.active) {
			os << stt.entc_ad2d;
		}
	}
	os << endl;
	return os;
}

bool cmpNodeMIST(nodeMIST &lh, nodeMIST &rh) {
	return lh.MI > rh.MI;
}

bool cmpEdgeMIST(const Edge &lh, const Edge &rh) {
	return lh.weight > rh.weight;
}

void printHeader() {
	chrono::time_point<chrono::system_clock> start;
	start = chrono::system_clock::now();
	time_t start_time = chrono::system_clock::to_time_t(start);
	cout << endl << endl
			<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
			<< endl
			<< "|                            CENTRE " CENTRE_VERSION_STRING "                                     |"
			<< endl
			<< "|       Configurational Entropy Estimation from Molecular Dynamics Data       |"
			<< endl
			<< "|                                       by                                    |"
			<< endl
			<< "|                          Information Theoretic Methods                      |"
			<< endl
			<< "|                              MIE/MIST/AMIE/AMIST                            |"
			<< endl
			<< "|                                      and                                    |"
			<< endl
			<< "|                          ML/MM/CS/JS Entropy Estimator                      |"
			<< endl
			<< "|                                                                             |"
			<< endl << "|             Date Time: " << ctime(&start_time)
			<< "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"
			<< endl << endl;
}

void printUsage() {

#ifdef USE_OMPMPI
	cout << "USAGE: " << endl
			<< "mpirun -n $NUM_PROC centre_OMPMPI [OPTION] -i <input-file>"
			<< endl;
#else
	cout << "USAGE: " << endl << "centre_OMP [OPTION] -i <input-file>" << endl;
#endif
	cout << "Optional arguments for the CENTRE " << CENTRE_VERSION_STRING
			<< endl;
	cout << "    -h  --help  		print this help message and exit" << endl;
	cout
			<< "    -s  --sample 		print a sample input file with description of parameters and exit"
			<< endl;
	cout
			<< "    -O 		        overwite any already existing files in output directory"
			<< endl;
}

void printSample() {
	 cout << "# This group of parameters are used to control the common execution of " << endl;
	 cout << "# CENTRE program" << endl;
	 cout << "[control]" << endl;
	 cout << "    # datatype: string. It specifies the path of directory which contains" << endl;
	 cout << "    # input files" << endl;
	 cout << "    #   parameter type      : required" << endl;
	 cout << "    infilepath           = /path/to/directory/of/input/files" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: string. It specifies the path of directory where output " << endl;
	 cout << "    # files will be created" << endl;
	 cout << "    #   parameter type      : required" << endl;
	 cout << "    outfilepath          = /path/to/directory/of/output/files" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: string. It specifies the name of file `outfilepath` " << endl;
	 cout << "    # directory which logs" << endl;
	 cout << "    # information about the progress of an on-going executation of CENTRE" << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : progress.info" << endl;
	 cout << "    infofile             = progress-1.info" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned byte [0, 255] It specifies the total number of " << endl;
	 cout << "    # steps in which entropy is computed for entire trajectory. It must " << endl;
	 cout << "    # be a factor of total_frames." << endl;
	 cout << "    #   parameter type: optional" << endl;
	 cout << "    #   default value : 20" << endl;
	 cout << "    nsteps               = 10" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: boolean true/false. It specifies whether to do " << endl;
	 cout << "    # discretization of BAT trajecory." << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : true" << endl;
	 cout << "    discretize           = true" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: boolean true/false. It specifies whether to do entropy " << endl;
	 cout << "    # calculation from discretized trajectory." << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : true" << endl;
	 cout << "    calcentropy          = true" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: boolean true/false. It specifies whether to create " << endl;
	 cout << "    # convergance report after entropy calculation according to chosen " << endl;
	 cout << "    # method i.e. MIE/MIST/A-MIE/A-MIST " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : true" << endl;
	 cout << "    genconvgdata         = true" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned int. It specifies number of DOFs per CPU to be  " << endl;
	 cout << "    # loaded into memory during discretization stage of execution. A big " << endl;
	 cout << "    # value for this paramter will increase memory required during " << endl;
	 cout << "    # execution of the program. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 5" << endl;
	 cout << "    cachedscrtdimspercpu = 5" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned int. It specifies number of DOFs per CPU to be  " << endl;
	 cout << "    # loaded into memory during entropy calculation stage of execution. " << endl;
	 cout << "    # A big value for this paramter will increase memory required during" << endl;
	 cout << "    # execution of the program. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 10" << endl;
	 cout << "    cacheentcdimspercpu  = 5" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned int. It specifies number how often the progress.info" << endl;
	 cout << "    # has to be updated during discretization stage of execution. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 25" << endl;
	 cout << "    dscrinfofreq         = 25" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned int. It specifies number how often the progress.info" << endl;
	 cout << "    # has to be updated during entropy calculation stage of execution. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 500" << endl;
	 cout << "    entcinfofreq         = 500" << endl;
	 cout << "    " << endl;
	 cout << "[discretization]" << endl;
	 cout << "    # datatype: string. It specifies the name of BAT trajectory file. " << endl;
	 cout << "    #   parameter type      : required" << endl;
	 cout << "    fname           = 3SRI.lcr1-20_3sri_f_dfs_r2-1-5_bat.nc" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned long. It specifies the number of frames in " << endl;
	 cout << "    # BAT trajectory file. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : read-from-BAT-trajectory-file" << endl;
	 cout << "    nframe          = 10000" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned long. It specifies the number of bond DOFs in " << endl;
	 cout << "    # BAT trajectory file. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : read-from-BAT-trajectory-file" << endl;
	 cout << "    nbond           = 364" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned long. It specifies the number of angle DOFs in " << endl;
	 cout << "    # BAT trajectory file. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : read-from-BAT-trajectory-file" << endl;
	 cout << "    nangle          = 363" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned long. It specifies the number of dihedral DOFs in " << endl;
	 cout << "    # BAT trajectory file. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : read-from-BAT-trajectory-file" << endl;
	 cout << "    ndihed          = 362" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: boolean true/false. It specifies whether to shuffle the frames" << endl;
	 cout << "    # of the BAT trajectory during the discretization." << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : false" << endl;
	 cout << "    shuffleframes   = false" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: boolean true/false. It specifies whether to shuffle the DOFs of" << endl;
	 cout << "    # the BAT trajectory independently for every DOF during the discretization." << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : false" << endl;
	 cout << "    shuffledofs     = false" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned long. It specifies the seed for the randomizer" << endl;
	 cout << "    # used to generate randomization sequence to be used for shuffling the DOFs" << endl;
	 cout << "    # during the discretization." << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 0 # i.e. random-seed" << endl;
	 cout << "    randseed        = 1234" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: boolean true/false. It specifies whether to do optimization to" << endl;
	 cout << "    # minimize the width of distribution of dihedrals fluctuations." << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : true" << endl;
	 cout << "    optimizedih     = true" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: choose one from {histogram, vonmiseskde}. " << endl;
	 cout << "    # It specifies the method used for discretization options " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : histogram" << endl;
	 cout << "    pdfmethod       = vonmiseskde" << endl;
	 cout << "" << endl;
	 cout << "[vonmiseskde]" << endl;
	 cout << "    # datatype: boolean true/false. " << endl;
	 cout << "    # It specifies whether to write the frequencies of histogram to files. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : false" << endl;
	 cout << "    writefreq       = true" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: a comma seperated list of one or more of " << endl;
	 cout << "    #           {" << endl;
	 cout << "    #               NOTHING : No histogram" << endl;
	 cout << "    #               D1D  : Dihedral DOFs          1D histogram," << endl;
	 cout << "    #               DD2D : Dihedral-dihedral DOFs 2D histogram, " << endl;
	 cout << "    #           }. " << endl;
	 cout << "    #           where       " << endl;
	 cout << "    #              1D is equivallent to B1D, A1D, D1D" << endl;
	 cout << "    #              2D is equivallent to BB2D, AA2D, DD2D, BA2D, BD2D, AD2D" << endl;
	 cout << "    #   parameter type      : required" << endl;
	 cout << "    writeset        = 2D" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned byte [1, control.nsteps] It specifies the index of " << endl;
	 cout << "    # step from which frequencies are written to hist_*.nc files. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : control.nsteps -1" << endl;
	 cout << "    writestepstart  = 39" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned byte [1, control.nsteps] It specifies the stride of " << endl;
	 cout << "    # step used for writing frequencies to hist_*.nc files. " << endl;
	 cout << "    #      writestepstart to control.nsteps in increment of writestepstride" << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 1" << endl;
	 cout << "    writestepstride = 1" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned byte [1, 255] It specifies the maximum number of " << endl;
	 cout << "    # minimas to search for in the von Mises Distribution of fluctuation" << endl;
	 cout << "    # of dihedral DOFs" << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 5" << endl;
	 cout << "    nmaxconf        = 5" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: float, It specifies the initial guess of the concentration" << endl;
	 cout << "    # parameter for von Mises KDE fitting" << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 1.0" << endl;
	 cout << "    kappa           = 1.0" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned byte [1, 255] It specifies the number of steps in" << endl;
	 cout << "    # Steepest-descent-optimization (sdo)" << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 5" << endl;
	 cout << "    sdosteps        = 5" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned int, It specifies the number of iteration in" << endl;
	 cout << "    # Steepest-descent-optimization (sdo)" << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 1000" << endl;
	 cout << "    sdoiterations   = 1000" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: a small double much less than 1.0, It specifies the tolerance " << endl;
	 cout << "    # to check for convergence during Steepest-descent-optimization (sdo) " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 0.0001" << endl;
	 cout << "    sdoconvlimit    = 0.0001" << endl;
	 cout << "" << endl;
	 cout << "[histogram]" << endl;
	 cout << "    # datatype: boolean true/false. " << endl;
	 cout << "    # It specifies whether to write the frequencies of vmKDE states to files. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : false" << endl;
	 cout << "    writefreq       = true" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: a comma seperated list of one or more of " << endl;
	 cout << "    #           {" << endl;
	 cout << "    #               NOTHING : No histogram" << endl;
	 cout << "    #               B1D  : Bond DOFs              1D histogram, " << endl;
	 cout << "    #               A1D  : Angle DOFs             1D histogram, " << endl;
	 cout << "    #               D1D  : Dihedral DOFs          1D histogram," << endl;
	 cout << "    #               BB2D : Bond-bond DOFs         2D histogram," << endl;
	 cout << "    #               AA2D : Angle-angle DOFs       2D histogram, " << endl;
	 cout << "    #               DD2D : Dihedral-dihedral DOFs 2D histogram, " << endl;
	 cout << "    #               BA2D : Bond-angle DOFs        2D histogram," << endl;
	 cout << "    #               BD2D : Bond-dihedral DOFs     2D histogram," << endl;
	 cout << "    #               AD2D : Angle-dihedral DOFs    2D histogram," << endl;
	 cout << "    #               1D   , " << endl;
	 cout << "    #               2D" << endl;
	 cout << "    #           }. " << endl;
	 cout << "    #           where       " << endl;
	 cout << "    #              1D is equivallent to B1D, A1D, D1D" << endl;
	 cout << "    #              2D is equivallent to 1D, BB2D, AA2D, DD2D, BA2D, BD2D, AD2D" << endl;
	 cout << "    #   parameter type      : required" << endl;
	 cout << "    writeset        = 2D" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned byte [1, control.nsteps] It specifies the index of " << endl;
	 cout << "    # step from which frequencies are written to hist_*.nc files. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : control.nsteps -1" << endl;
	 cout << "    writestepstart  = 39" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned byte [1, control.nsteps] It specifies the stride of " << endl;
	 cout << "    # step used for writing frequencies to hist_*.nc files. " << endl;
	 cout << "    #      writestepstart to control.nsteps in increment of writestepstride" << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 1" << endl;
	 cout << "    writestepstride = 1" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned byte [1, 255] It specifies the number of bins used for" << endl;
	 cout << "    # histogram of fluctuation data in 1D." << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 30" << endl;
	 cout << "    nbins           = 30" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned byte [1, 255] It specifies the number of bins used for" << endl;
	 cout << "    # histogram of fluctuation data in 1D for optimization of dihedral " << endl;
	 cout << "    # distribution width. Usually, 3 * referencenbins bins are used while" << endl;
	 cout << "    # scanning fluctuations-distribtion for optimal width." << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 30" << endl;
	 cout << "    referencenbins  = 30" << endl;
	 cout << "    " << endl;
	 cout << "[entropy]" << endl;
	 cout << "    # datatype: a comma seperated list of one or more of " << endl;
	 cout << "    #           {" << endl;
	 cout << "    #               NOTHING : No entropy calculation" << endl;
	 cout << "    #               B1D  : Bond DOFs              1D entropy, " << endl;
	 cout << "    #               A1D  : Angle DOFs             1D entropy, " << endl;
	 cout << "    #               D1D  : Dihedral DOFs          1D entropy," << endl;
	 cout << "    #               BB2D : Bond-bond DOFs         2D mutual information," << endl;
	 cout << "    #               AA2D : Angle-angle DOFs       2D mutual information, " << endl;
	 cout << "    #               DD2D : Dihedral-dihedral DOFs 2D mutual information, " << endl;
	 cout << "    #               BA2D : Bond-angle DOFs        2D mutual information," << endl;
	 cout << "    #               BD2D : Bond-dihedral DOFs     2D mutual information," << endl;
	 cout << "    #               AD2D : Angle-dihedral DOFs    2D mutual information," << endl;
	 cout << "    #               1D   , " << endl;
	 cout << "    #               2D" << endl;
	 cout << "    #           }. " << endl;
	 cout << "    #           where       " << endl;
	 cout << "    #              1D is equivallent to B1D, A1D, D1D" << endl;
	 cout << "    #              2D is equivallent to 1D, BB2D, AA2D, DD2D, BA2D, BD2D, AD2D" << endl;
	 cout << "    #   parameter type      : required" << endl;
	 cout << "    workset         = 2D" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: boolean true/false. It specifies whether to use only DOFs" << endl;
	 cout << "    # present in the subsetfile. This is intended to filter-out frozen DOFs" << endl;
	 cout << "    # it instructs whether to use only flexible DOFs present in subsetfile." << endl;
	 cout << "    #   parameter type      : required" << endl;
	 cout << "    #   default value       : false" << endl;
	 cout << "    usesubset       = true" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: boolean true/false. It specifies whether to do run a distance" << endl;
	 cout << "    # cutoff based or Neighbor Approximated version of MIE/MIST. The neighbor" << endl;
	 cout << "    # DOFs information would be read from the neighborfile." << endl;
	 cout << "    #   parameter type      : required" << endl;
	 cout << "    #   default value       : false" << endl;
	 cout << "    useneighbor     = true" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: boolean true/false. It specifies whether to do optimization to" << endl;
	 cout << "    # minimize the width of distribution of dihedrals fluctuations." << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : true" << endl;
	 cout << "    jacobian        = true" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: string. It specifies the path of directory which contains " << endl;
	 cout << "    # input files" << endl;
	 cout << "    #   parameter type      : optional if usesubset==false else required" << endl;
	 cout << "    subsetfile      = subset_3sri_f_dfs_r2-1-5.txt" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: string. It specifies the path of directory which contains " << endl;
	 cout << "    # input files" << endl;
	 cout << "    #   parameter type      : optional if useneighbor==false else required" << endl;
	 cout << "    neighborfile    = neigh_3sri_f_dfs_r2-1-5_6.txt" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: choose one from {MIE, MIST}. " << endl;
	 cout << "    # It specifies the method used for entropy scoring.  " << endl;
	 cout << "    # Note: when useneighbor = true the methods MIE/MIST become A-MIE/A-MIST." << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : MIST" << endl;
	 cout << "    scoringmethod   = MIST" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: choose one or more from " << endl;
	 cout << "    #           {" << endl;
	 cout << "    #               ML : Maximum Liklihood," << endl;
	 cout << "    #               MM : Miller Madow, " << endl;
	 cout << "    #               CS : Chao & Shen, " << endl;
	 cout << "    #               JS : James and Stein" << endl;
	 cout << "    #           }. " << endl;
	 cout << "    # It specifies comma seperated list of entropy estimator used" << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : MM" << endl;
	 cout << "    estimator       = ML" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned long. It specifies starting frame number used in   " << endl;
	 cout << "    # entropy calculation. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 0" << endl;
	 cout << "    startframe      = 0" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned long. It specifies every nth frame used in entropy  " << endl;
	 cout << "    # calculation. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : 1" << endl;
	 cout << "    strideframe     = 1" << endl;
	 cout << "    " << endl;
	 cout << "    # datatype: unsigned long. It specifies the number of frames used for " << endl;
	 cout << "    # entropy calculation. " << endl;
	 cout << "    #   parameter type      : optional" << endl;
	 cout << "    #   default value       : read-from-BAT-trajectory-file" << endl;
	 cout << "    nframe          = 10000" << endl;
	 cout << "" << endl;
	 cout << "" << endl;
}

void printFooter() {
	cout << endl << endl
			<< "Kindly Acknowledge the use of CENTRE in your work through citation"
			<< endl
			<< "Author: Shailesh Kr. Panday                                       "
			<< endl
			<< "                               *.*                                "
			<< endl
			<< "                                                                  "
			<< endl;
}