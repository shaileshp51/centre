/* 
 * File:   Utils.cpp
 * Author: shailesh
 * 
 * Created on 13 July, 2014, 10:52 PM
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm> 
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

#include "Utils.h"
#include "Constants.h"

Utils::Utils() {
}

Utils::Utils(const Utils &orig) {
}

Utils::~Utils() {
}

std::string Utils::trim(const std::string &str) {
	std::string temp(str);
	std::string::size_type pos = temp.find_last_not_of(" \t");
	if (pos != std::string::npos) {
		temp.erase(pos + 1);
		pos = temp.find_first_not_of(" \t");
		if (pos != std::string::npos)
			temp.erase(0, pos);
	} else
		temp.erase(temp.begin(), temp.end());
	return (temp);
}

bool Utils::string2bool(std::string str) {
	if (str.empty())
		return (!str.empty());
	std::transform(str.begin(), str.end(), str.begin(), ::tolower);
	std::istringstream is(str);
	bool b;
	is >> std::boolalpha >> b;
	return b;
}

double Utils::modifiedBesselI(const int n, const double x) {

	int M;
	int IACC = 40;
	double bignumber_positive = 1.e10, BIGNI = 1.e-10;
	double result, TOX, BIM, BI, BIP;

	if (n == 0) {
		result = modifiedBesselI0(x);
		return result;
	}

	if (n == 1) {
		result = modifiedBesselI1(x);
		return result;
	}

	if (x == 0.0) {
		result = 0.0;
		return result;
	}

	TOX = 2.0 / x;
	BIP = 0.0;
	BI = 1.0;
	result = 0.0;
	M = 2 * ((n + (int) (sqrt((float) (IACC * n)))));

	for (int J = M; J >= 1; --J) {
		BIM = BIP + (double) (J) * TOX * BI;
		BIP = BI;
		BI = BIM;
		if (fabs(BI) > bignumber_positive) {
			BI = BI * BIGNI;
			BIP = BIP * BIGNI;
			result = result * BIGNI;
		}
		if (J == n) {
			result = BIP;
		}
	}

	result = result * modifiedBesselI0(x) / BI;

	return result;
}

double Utils::modifiedBesselI0(const double x) {

	double result, y;
	double P[7] = { 1.0e0, 3.5156229e0, 3.0899424e0, 1.2067429e0, 0.2659732e0,
			0.360768e-1, 0.45813e-2 };
	double Q[9] =
			{ 0.39894228e0, 0.1328592e-1, 0.225319e-2, -0.157565e-2,
					0.916281e-2, -0.2057706e-1, 0.2635537e-1, -0.1647633e-1,
					0.392377e-2 };
	double ax, bx;

	if (fabs(x) < 3.750) {
		y = (x / 3.750);
		y = y * y;
		result =
				P[0]
						+ y
								* (P[1]
										+ y
												* (P[2]
														+ y
																* (P[3]
																		+ y
																				* (P[4]
																						+ y
																								* (P[5]
																										+ y
																												* P[6])))));
	} else {
		ax = fabs(x);
		y = 3.750 / ax;
		bx = exp(ax) / sqrt(ax);
		double temp = (Q[5] + y * (Q[6] + y * (Q[7] + y * Q[8])));
		ax = Q[0]
				+ y * (Q[1] + y * (Q[2] + y * (Q[3] + y * (Q[4] + y * temp))));
		result = ax * bx;
	}

	return result;
}

double Utils::modifiedBesselI1(const double x) {
	double result, y;
	double P[7] = { 0.5e0, 0.87890594e0, 0.51498869e0, 0.15084934e0,
			0.2658733e-1, 0.301532e-2, 0.32411e-3 };
	double Q[9] = { 0.39894228e0, -0.3988024e-1, -0.362018e-2, 0.163801e-2,
			-0.1031555e-1, 0.2282967e-1, -0.2895312e-1, 0.1787654e-1,
			-0.420059e-2 };
	double ax, bx;

	if (fabs(x) < 3.75e0) {
		y = (x / 3.75e0);
		y = y * y;
		result =
				x
						* (P[0]
								+ y
										* (P[1]
												+ y
														* (P[2]
																+ y
																		* (P[3]
																				+ y
																						* (P[4]
																								+ y
																										* (P[5]
																												+ y
																														* P[6]))))));
	} else {
		ax = fabs(x);
		y = 3.75e0 / ax;
		bx = exp(ax) / sqrt(ax);
		double temp = (Q[5] + y * (Q[6] + y * (Q[7] + y * Q[8])));
		ax = Q[0]
				+ y * (Q[1] + y * (Q[2] + y * (Q[3] + y * (Q[4] + y * temp))));
		result = ax * bx;
	}
	return result;
}

void printParsedConfiguration(
		std::pair<std::string, std::vector<std::string> > x) {
	std::cout << x.first << ": ";
	if (x.second.size() > 1) {
		std::cout << "[";
		int f = 0;
		for (auto v : x.second) {
			if (f != 0) {
				std::cout << ", " << v;
			} else {
				std::cout << v;
			}
			++f;
		}
		std::cout << "]" << std::endl;
	} else {
		std::cout << x.second[0] << std::endl;
	}
}

std::map<std::string, std::vector<std::string> > Utils::parseConfiguration(
		std::string filename) {
	std::ifstream input(filename.c_str()); //The input stream
	std::map<std::string, std::vector<std::string> > ans; //A map of key-value pairs in the file
	std::string line;
	std::string section;
	while (std::getline(input, line, '\n')) //Keep on going as long as the file stream is good
	{
		std::string key; //The key
		std::string value; // A value
		std::vector<std::string> values; //List of value
		// std::cout << "read: '" << line << "'" << endl;
		line = Utils::trim(line);
		if (line.length() > 0) {
			if (line.at(0) == '#') {
				// std::cout << "COMMENT<<<" << line << ">>>" << std::endl;
				continue;
			} else if (line.length() > 2 && line.at(0) == '['
					&& line.at(line.length() - 1) == ']') {
				section = line.substr(1, line.length() - 2);
			} else {
				std::stringstream sline(line);
				// std::cout << "read-stringstream: '" << sline << "'" << endl;
				std::getline(sline, key, '='); //Read up to the : delimiter into key
				key = Utils::trim(key);
				if (!section.empty()) {
					key = section + '.' + key;
				}
				std::getline(sline, value, '\n'); //Read up to the newline into value
				value = Utils::trim(value);
				std::string::size_type pos1 = value.find_first_of("\""); //Find the first quote in the value
				std::string::size_type pos2 = value.find_last_of("\""); //Find the last quote in the value

				//Check if the found positions are all valid
				if (pos1 != std::string::npos && pos2 != std::string::npos
						&& pos2 > pos1) {
					value = value.substr(pos1 + 1, pos2 - pos1 - 1); //Take a substring of the part between the quotes
					values.push_back(value);

					ans[key] = values; //Store the result in the map

				} else {
					std::string token;
					std::stringstream swords(value);
					int token_no = 0;
					while (swords >> token) {
						if (token.length() > 0 && token.at(0) == '#') {
							// std::cout << "INLINE-COMMENT<<<#" << swords.rdbuf() << ">>>" << std::endl;
							break;
						}
						++token_no;
						token = Utils::trim(token, ",");
						// std::cout << "token(" << token_no << "): |" << token << "| ";
						values.push_back(token);
					}
					ans[key] = values;
				}
			}
		}
	}
	input.close(); //Close the file stream
	return ans; //And return the result
}

