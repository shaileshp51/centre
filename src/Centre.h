/*!\brief This file defines basic class and enum used in project.
 *
 * File:   Centre.h
 * Author: shailesh
 *
 * Created on 5 August, 2016, 10:41 PM
 */

#ifndef INC_CENTRE_H
#define INC_CENTRE_H

#ifndef INC_CENTRE_CONFIG_H
#include "configcentre.h"
#endif

#include <map> // map
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <cstddef>
#include <sstream>
#include <fstream>
#include <string>
#include <limits>

/*!
 An enum CoordSys

 This enum type defines types of Coordinate Systems.
 */
enum CoordSys : u_int8_t {
	CARTESIAN = 0, /*!< The Coordinate System in Use is Cartesian */
	BAT = 1 /*!< The Coordinate System in Use is Bond/Angle/Torsion (BAT) */
};

/*!
 An enum BAT_t

 This enum type defines types of Degrees-of-freedom in BAT Coordinate Systems.
 */
enum BAT_t : u_int8_t {
	NONE = 0, /*!< The BAT Type is unknown */
	BOND = 1, /*!< The BAT Type is Bond */
	ANGLE = 2, /*!< The BAT Type is Angle */
	DIHEDRAL = 3 /*!< The BAT Type is Dihedral/Torsion */
};

/*!
 An enum Entropy4Data

 This enum type describes dataset to be used for entropy estimation.
 */
enum Entropy4Data : u_int8_t {
	ALL = 0, /*!< ALL: Entire data avaliable in input data files will be used for entropy analysis */
	SUBSET = 1 /*!< SUBSET: A subset of entire data avaliable in input data files will be used for entropy analysis */
};

/*!
 An enum PDFMethod

 This enum type describes methods for PDF estimation from input data.
 */
enum PDFMethod : u_int8_t {
	HISTOGRAM = 0, /*!< Histogram method is used */
	vonMisesKDE = 1 /*!< vonMises KDE state is used */
};

/*!
 An enum EntropyMethods

 This enum type describes entropy methods.

 This enum lists all the methods for entropy, However, only MIE/AMIE/MIST implemented.
 */
enum EntropyMethods : u_int8_t {
	QH = 0, /*!< Quasi harmonic */
	MIE = 1, /*!< Mutual Information Expansion */
	AMIE = 2, /*!< Approximate MIE */
	MIST = 4, /*!< Maximum Information Spanning Tree */
	AMIST = 8, /*!< Maximum Local Approximation */
};

/*!
 An enum EntropyEstimators

 This enum type defines entropy estimators implemented.
 */
enum EntropyEstimators : u_int8_t {
	ML = 0, /*!< ML: Maximum Likelihood Estimator */
	MM = 1, /*!< MM: Miller Madow Estimator */
	CS = 2, /*!< CS: Chao Shen */
	JS = 3 /*!< JS: James and Stein Shrinkage Estimator */
};

enum BATSet : u_int16_t {
	NOTHING = 0, /*!< None */
	B1D = 1, /*!< Bonds 1D */
	A1D = 2, /*!< Angles 1D */
	D1D = 4, /*!< Torsions 1D */
	X1D = 7, /*!< Bonds, Angles, Torsions 1D */
	BB2D = 8, /*!< Bonds 2D */
	AA2D = 16, /*!< Angles 2D */
	DD2D = 32, /*!< Torsions 2D */
	XX2D = 63, /*!< Bonds, Angles, Torsions 2D */
	BA2D = 64, /*!< Bonds-Angles 2D */
	BD2D = 128, /*!< Bonds-Torsions 2D */
	AD2D = 256, /*!< Angles-Torsions 2D */
	XY2D = 511 /*!< Bonds, Angles, Torsions 2D, Bonds-Angles 2D, Bonds-Torsions 2D, Angles-Torsions 2D */
};

using namespace std;

//!  A range class.
/*!
 This class represents a range where a range can be generated from stast value till end value
 with elements seperated by step
 */
class Range {
	int start_; /*!< start value for range */
	int end_; /*!< end value for range */
	int step_; /*!< increment between consecutive values in range */

	/*!
	 \fn Range()

	 Default constructor of the class Range
	 */
	Range() {
		start_ = 0;
		end_ = 0;
		step_ = 0;
	}

	/*!
	 \fn Range(int end)

	 A single argument constructor of the class Range with start value is 0, and step is 1
	 \param[in] end the stop value for the range
	 */
	Range(int end) {
		start_ = 0;
		end_ = end;
		step_ = 1;
	}

	/*!
	 \fn Range(int end)

	 A single argument constructor of the class Range with step value 1
	 \param[in] start the start value for the range
	 \param[in] end the stop value for the range
	 */
	Range(int start, int end) {
		start_ = start;
		end_ = end;
		step_ = 1;
	}

	/*!
	 \fn Range(int end)

	 A single argument constructor of the class Range
	 \param[in] start the start value for the range
	 \param[in] end the stop value for the range
	 \param[in] step the step value for the range
	 */
	Range(int start, int end, int step) {
		start_ = start;
		end_ = end;
		step_ = step;
	}

	/*!
	 \fn bool isNone()

	 Return \c true if range is empty; otherwise false
	 */
	bool isNone() {
		bool result = false;
		if (start_ == 0 && end_ == 0 && step_ == 0) {
			result = true;
		}
		return (result);
	}

	/*!
	 \fn void expand(std::vector<int> *vec)

	 Gives a vector of range member values
	 \param[in/out] *vec the vec is modified by pushing range values at the end.
	 */
	void expand(std::vector<int> *vec) {
		int f = 0, t = 0, s = 0;
		if (step_ == 0)
			s = 1;
		if (start_ == -1)
			f = 0;

		vec->clear();
		for (int i = f; i < t; i += s) {
			vec->push_back(i);
		}
	}
};

//!  A Edge class.
/*!
 This class represents a weighted edge in a Graph.
 */
class Edge {
public:
	int source; /*!< start node of the edge */
	int dest; /*!< end node of the edge */
	double weight; /*!< weight of the edge */
	Edge(int s, int d, double w) {
		this->source = s;
		this->dest = d;
		this->weight = w;
	}
};

//!  A Node class.
/*!
 This class represents a node in a Graph.
 */
class Node {
public:
	BAT_t type; /*!< Type of BAT node belongs to e.g. Bond or Angle or Torsion */
	u_int id; /*!< ID of the node */
	Node(BAT_t typ, u_int i) {
		type = typ;
		id = i;
	}
};

//!  A Graph class.
/*!
 This class represents a a connected, undirected
 and weighted graph as a collection of nodes and edges amongst them,
 It is a undirected weighted Graph.
 */
class Graph {
public:
	std::vector<Node> nodes; /*!< A vector of nodes in the Graph */
	std::vector<Edge> edges; /*!< A vector of edges in the Graph */
	Graph() {

	}
};

//!  A Component class.
/*!
 A class to represent a subset for union-find.
 */
class Component {
public:
	u_int parent;
	int rank;
};

/*!
 \fn int findRoot(std::vector<Component>& compnents, int i)

 A utility function to find set of an element i (uses path compression technique).
 Returns \a int the parent node of component \a i in the vector of component \a compnents.
 */
int findRoot(std::vector<Component> &compnents, int i);

/*!
 \fn void unionComponents(std::vector<Component>& compnents, int x, int y)

 A utility function that does union of two componets x and y by attaching smaller rank
 tree under root of high rank tree; Links parents of two components \a x and \a y
 to join two components of the vector \a components.
 */
void unionComponents(std::vector<Component> &compnents, int x, int y);

/*!
 \fn std::vector<Edge> boruvkaMST(Graph& graph)

 Returns a vector of Edges in the Maximum Spanning Tree for the supplied
 Graph \a graph. Where the MST is calculated using boruvka Algorithm
 */
std::vector<Edge> boruvkaMST(Graph &graph);

class nodeMIST {

public:
	BAT_t d1type; /*!< Type of DOF of first DOF in DOF pair */
	BAT_t d2type; /*!< Type of DOF of second DOF in DOF pair */
	u_int d1indx; /*!< Index of first DOF in DOF pair */
	u_int d2indx; /*!< Index of second DOF in DOF pair */

	double MI; /*!< Mutual information for DOF pair */

	nodeMIST();

	nodeMIST(BAT_t d1t, BAT_t d2t, int d1indx, int d2indx, double I);
};

//!  A SubsetData class.
/*!
 A class to store a subset of DOFs to be used for entropy estimation.

 The total number of bonds, angles and dihedrals/torsions in the BAT trajectory
 are stored in fields, The indexes of selected bonds, angles, torsions chosen for
 entropy calculation are stored in bonds angles and torsions vectors.
 */
class SubsetData {
private:
	u_int n_tot_bnd; /*!< Total number of bonds DOFs in the subset */
	u_int n_tot_ang; /*!< Total number of angles DOFs in the subset */
	u_int n_tot_tor; /*!< Total number of torsions DOFs in the subset */
	bool has_bonds_; /*!< Whether bond DOFs is part of the subset */
	bool has_angles_; /*!< Whether angle DOFs is part of the subset */
	bool has_torsions_; /*!< Whether torsion DOFs is part of the subset */
	bool isgood_; /*!< Whether object could be constructed successfully */
	std::vector<u_int> bonds; /*!< A vector of ids of bonds DOFs in the subset */
	std::vector<u_int> angles; /*!< A vector of ids of angles DOFs in the subset */
	std::vector<u_int> torsions; /*!< A vector of ids of torsions DOFs in the subset */
	BATSet entropy_workset;
public:

	/*!
	 \fn SubsetData()

	 Default constructor of the class SubsetData
	 */
	SubsetData() {
	}
	;

	/*!
	 \fn SubsetData(u_int n_bnd, u_int n_ang, u_int n_tor)

	 A constructor of the class SubsetData which accepts number of bonds, angles and torsions
	 \param[in] n_bnd number of bonds
	 \param[in] n_ang number of angles
	 \param[in] n_tor number of torsions
	 */
	SubsetData(u_int n_bnd, u_int n_ang, u_int n_tor, BATSet entropy_workset);

	/*!
	 \fn SubsetData(std::string const &fname)

	 A constructor of the class SubsetData which accepts name of the file (SubsetData file format)
	 from which the subset data is read.
	 \param[in] fname name of the subset file
	 */
	SubsetData(std::string const &fname, BATSet entropy_workset);

	/*!
	 \fn bool isGood() const

	 Return \c true if object constructed successfully; else \c false.
	 */
	bool isGood() const {
		return isgood_;
	}

	/*!
	 \fn bool getHasBonds() const

	 Return \c true if subset has any bonds DOFs; else \c false.
	 */
	bool getHasBonds() const {
		return has_bonds_;
	}

	/*!
	 \fn bool getHasAngles() const

	 Return \c true if subset has any angles DOFs; else \c false.
	 */
	bool getHasAngles() const {
		return has_angles_;
	}

	/*!
	 \fn bool getHasTorsions() const

	 Return \c true if subset has any torsion/dihedral DOFs; else \c false.
	 */
	bool getHasTorsions() const {
		return has_torsions_;
	}

	/*!
	 \fn void setHasBonds(bool v)

	 Set status \a v whether there are bond DOFs in subset.
	 */
	void setHasBonds(bool v) {
		has_bonds_ = v;
	}

	/*!
	 \fn void setHasAngles(bool v)

	 Set status \a v whether there are angle DOFs in subset.
	 */
	void setHasAngles(bool v) {
		has_angles_ = v;
	}

	/*!
	 \fn void setHasTorsions(bool v)

	 Set status \a v whether there are torsion DOFs in subset.
	 */
	void setHasTorsions(bool v) {
		has_torsions_ = v;
	}

	/*!
	 \fn void push_bond(u_int value)

	 Add id of a bond \a value in the vector bonds of subset.
	 */
	void push_bond(u_int value);

	/*!
	 \fn void push_bonds(std::vector<u_int> const& vec)

	 Add ids of a bond \a vec in the vector bonds of subset.
	 */
	void push_bonds(std::vector<u_int> const &vec);

	/*!
	 \fn void push_angle(u_int value)

	 Add id of a angles \a value in the vector bonds of subset.
	 */
	void push_angle(u_int value);

	/*!
	 \fn void push_angles(std::vector<u_int> const& vec)

	 Add ids of a angles \a vec in the vector angles of subset.
	 */
	void push_angles(std::vector<u_int> const &vec);

	/*!
	 \fn void push_torsion(u_int value)

	 Add id of a torsion \a value in the vector torsions of subset.
	 */
	void push_torsion(u_int value);

	/*!
	 \fn void push_torsions(std::vector<u_int> const& vec)

	 Add ids of a torsions \a vec in the vector torsions of subset.
	 */
	void push_torsions(std::vector<u_int> const &vec);

	/*!
	 \fn u_int bondsSize() const

	 Return number of bond DOFs stored in vector bonds of subset.
	 */
	u_int bondsSize() const {
		return bonds.size();
	}

	/*!
	 \fn u_int anglesSize() const

	 Return number of angle DOFs stored in vector angles of subset.
	 */
	u_int anglesSize() const {
		return angles.size();
	}

	/*!
	 \fn u_int torsionsSize() const

	 Return number of torsion DOFs stored in vector torsions of subset.
	 */
	u_int torsionsSize() const {
		return torsions.size();
	}

	/*!
	 \fn const std::vector<u_int>& getAngles() const

	 Return a constant vector reference of angle DOFs stored in vector angles of subset.
	 */
	const std::vector<u_int>& getAngles() const {
		return angles;
	}

	/*!
	 \fn const u_int getAngleId(u_int id) const

	 Return the id of angles stored at i-th index in angles of subset.
	 */
	const u_int getAngleId(u_int id) const {
		return angles[id];
	}

	/*!
	 \fn const u_int getBondId(u_int id) const

	 Return the id of bond stored at i-th index in bonds of subset.
	 */
	const u_int getBondId(u_int id) const {
		return bonds[id];
	}

	/*!
	 \fn const std::vector<u_int>& getBonds() const

	 Return a constant vector reference of bond DOFs stored in vector bonds of subset.
	 */
	const std::vector<u_int>& getBonds() const {
		return bonds;
	}

	/*!
	 \fn const u_int getTorsionId(u_int id) const

	 Return the id of torsion stored at i-th index in torsions of subset.
	 */
	const u_int getTorsionId(u_int id) const {
		return torsions[id];
	}

	const u_int getBATTypeId(BAT_t bat_type, u_int id) const;

	/*!
	 \fn const std::vector<u_int>& getTorsions() const

	 Return a constant vector reference of torsion DOFs stored in vector torsions of subset.
	 */
	const std::vector<u_int>& getTorsions() const {
		return torsions;
	}

	/*!
	 \fn void setAngles(const std::vector<u_int>& angles)

	 Set ids of a angles \a torsions in the vector angles of subset.
	 */
	void setAngles(const std::vector<u_int> &angles) {
		this->angles = angles;
	}

	/*!
	 \fn void setBonds(const std::vector<u_int>& bonds)

	 Set ids of a bonds \a torsions in the vector bonds of subset.
	 */
	void setBonds(const std::vector<u_int> &bonds) {
		this->bonds = bonds;
	}

	/*!
	 \fn void setTorsions(const std::vector<u_int>& torsions)

	 Set ids of a torsions \a torsions in the vector torsions of subset.
	 */
	void setTorsions(const std::vector<u_int> &torsions) {
		this->torsions = torsions;
	}
};

//!  A Neighbors class.
/*!
 A class to store a neighbors of DOFs to be used for entropy estimation.

 The total number of bonds, angles, dihedrals/torsions, and cross-terms
 bond-angle, bond-dihedral and angle-dihedral in the BAT trajectory are
 stored in fields, The indexes of selected bonds, angles, torsions chosen for
 entropy calculation are stored in bonds angles and torsions vectors.
 */
class Neighbors {
private:
	bool has_bonds_;
	bool has_angles_;
	bool has_torsions_;
	bool has_bacross_;
	bool has_bdcross_;
	bool has_adcross_;
	bool isgood_;
	std::map<u_int, std::vector<u_int> > bonds; // ordered: sorted by key
	std::map<u_int, std::vector<u_int> > angles; // ordered: sorted by key
	std::map<u_int, std::vector<u_int> > torsions; // ordered: sorted by key
	std::map<u_int, std::vector<u_int> > bacross; // ordered: sorted by key
	std::map<u_int, std::vector<u_int> > bdcross; // ordered: sorted by key
	std::map<u_int, std::vector<u_int> > adcross; // ordered: sorted by key
	BATSet entropy_workset;
public:
	bool isGood() const {
		return isgood_;
	}

	bool getHasBonds() const {
		return has_bonds_;
	}
	bool getHasAngles() const {
		return has_angles_;
	}
	bool getHasTorsions() const {
		return has_torsions_;
	}

	void setHasBonds(const bool v) {
		has_bonds_ = v;
	}

	void setHasAngles(const bool v) {
		has_angles_ = v;
	}

	void setHasTorsions(const bool v) {
		has_torsions_ = v;
	}

	Neighbors(const u_int n_bnd, const u_int n_ang, const u_int n_tor);

	Neighbors(const SubsetData subset, BATSet entropy_workset);

	Neighbors(std::string const &fname, BATSet entropy_workset);

	Neighbors() {
	}
	;

	~Neighbors() {

	}

	friend std::ostream& operator<<(std::ostream &os, const Neighbors &ht);

	void bondKeys(std::vector<u_int>&) const;

	void angleKeys(std::vector<u_int>&) const;

	void torsionKeys(std::vector<u_int>&) const;

	void bacrossKeys(std::vector<u_int>&) const;

	void bdcrossKeys(std::vector<u_int>&) const;

	void adcrossKeys(std::vector<u_int>&) const;

	const u_int bondsSize() const;

	const u_int anglesSize() const;

	const u_int torsionsSize() const;

	const u_int baCrossSize() const;

	const u_int bdCrossSize() const;

	const u_int adCrossSize() const;

	const u_int baNeighSize(u_int id) const;

	const u_int bdNeighSize(u_int id) const;

	const u_int adNeighSize(u_int id) const;

	const u_int bondNeighSize(u_int id) const;

	const u_int angleNeighSize(u_int id) const;

	const u_int torsionNeighSize(u_int id) const;

	const std::vector<u_int>& getAngle(const u_int id) const;

	void setAngle(const u_int id, const std::vector<u_int> &nei);

	const std::vector<u_int>& getBond(const u_int id) const;

	void setBond(const u_int id, const std::vector<u_int> &nei);

	const std::vector<u_int>& getTorsion(const u_int id) const;

	void setTorsion(const u_int id, const std::vector<u_int> &nei);

	const std::map<u_int, std::vector<u_int> >& getAngles() const;

	const std::map<u_int, std::vector<u_int> >& getBonds() const;

	const std::map<u_int, std::vector<u_int> >& getTorsions() const;

	const std::map<u_int, std::vector<u_int> >& getAdcross() const;

	void setAdcross(const u_int id, const std::vector<u_int> &nei);

	const std::map<u_int, std::vector<u_int> >& getBacross() const;

	const std::vector<u_int>& getBacross(u_int id) const;

	const std::vector<u_int>& getBdcross(u_int id) const;

	const std::vector<u_int>& getAdcross(u_int id) const;

	const std::vector<u_int>& getDimTypeNeighs(const BATSet dimType,
			const u_int id) const;

	void setBacross(const u_int id, const std::vector<u_int> &nei);

	const std::map<u_int, std::vector<u_int> >& getBdcross() const;

	void setBdcross(const u_int id, const std::vector<u_int> &nei);

	bool isHasAdcross() const;

	void setHasAdcross(bool hasAdcross);

	bool isHasBacross() const;

	void setHasBacross(bool hasBacross);

	bool isHasBdcross() const;

	void setHasBdcross(bool hasBdcross);
};

#endif /* INC_CENTREMMD_H */
