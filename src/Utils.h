/* 
 * File:   Utils.h
 * Author: shailesh
 *
 * Created on 13 July, 2014, 10:52 PM
 */

#ifndef UTILS_H
#define	UTILS_H

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include "Constants.h"

/*! \file Utils.h
 \brief A collection of utility functions

 This class defines important functions which will be used in project,
 in several places.
 */
class Utils {
public:
	/**
	 * The default constructor.
	 */
	Utils();

	/**
	 * A copy constructor.
	 * @param orig
	 */
	Utils(const Utils &orig);

	/**
	 * The destructor of object.
	 */
	virtual ~Utils();

	/**
	 * A function which returns boolean equivallent from the string, valid inputs are:
	 * "true" or "false"
	 * @param str
	 * @return
	 */
	static bool string2bool(std::string str);

	/**
	 * A function which return a string after removing all leading and trailing spaces
	 * from it.
	 * @param str
	 * @return
	 */
	static std::string trim(const std::string &str);

	/**
	 *
	 * @param x
	 * @return
	 */
	static double modifiedBesselI1(const double x);
	/**
	 *
	 * @param x
	 * @return
	 */
	static double modifiedBesselI0(const double x);

	/**
	 *
	 * @param n
	 * @param x
	 * @return
	 */
	static double modifiedBesselI(const int n, const double x);

	/*!
	 Returns a string generated input string \a s, after trimming space and form-feed characters from right.
	 */
	static inline std::string& rtrim(std::string &s, const char *t =
			" \t\n\r\f\v") {
		s.erase(s.find_last_not_of(t) + 1);
		return s;
	}

	/*!
	 Returns a string generated input string \a s, after trimming space and form-feed characters from left.
	 */
	static inline std::string& ltrim(std::string &s, const char *t =
			" \t\n\r\f\v") {
		s.erase(0, s.find_first_not_of(t));
		return s;
	}

	/*!
	 Returns a string generated input string \a s, after trimming space and form-feed characters from both ends.
	 */
	static inline std::string& trim(std::string &s, const char *t =
			" \t\n\r\f\v") {
		return ltrim(rtrim(s, t), t);
	}

	/*!
	 Returns \c true if input string \a name, represents an existing file; otherwise \c false.
	 */
	static inline bool fileExists(const std::string &name) {
		std::ifstream f(name.c_str());
		return f.good();
	}

	/*!
	 Returns \c true if input string \a name, represents an existing directory; otherwise \c false.
	 */
	static inline bool isDirectory(const std::string &name) {
		struct stat st;
		bool res = false;
		if (stat(name.c_str(), &st) == 0)
			if ((st.st_mode & S_IFDIR) != 0)
				res = true;
		return res;
	}

	static void printParsedConfiguration(
			std::pair<std::string, std::vector<std::string>> x);

	static std::map<std::string, std::vector<std::string> > parseConfiguration(
			std::string filename);
};

#endif	/* UTILS_H */

