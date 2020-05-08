#ifndef INC_COORDINATEINFO_H
#define INC_COORDINATEINFO_H

#include <vector>

class CoordinateInfo {
private:
	int atmStartIdx_;
	int roots_[3];
	int nPseudos_;
	std::vector<long> pseudos_;
	int nBond_; ///< True if coords have associated bonds.
	int nAngle_; ///< True if coords include angles.
	int nDihedral_; ///< True if coords include dihedrals.
	bool hasTime_; ///< True if coords include time info.
	bool hasPseudo_;
public:
	/// CONSTRUCTOR

	CoordinateInfo() :
			nBond_(-1), nAngle_(-1), nDihedral_(-1), hasTime_(false) {
		roots_[0] = 0;
		roots_[1] = 0;
		roots_[2] = 0;
		nPseudos_ = -1;
		hasPseudo_ = false;
		atmStartIdx_ = -1;
		pseudos_ = std::vector<long>();
	}

	bool operator==(const CoordinateInfo &rhs) {
		bool pseudoMatch = pseudos_.size() == rhs.pseudos_.size();
		if (pseudoMatch) {
			for (size_t i = 0; i < pseudos_.size(); ++i) {
				pseudoMatch = pseudoMatch && pseudos_[i] == rhs.pseudos_[i];
			}
		}
		bool isEql = pseudoMatch && roots_[0] == rhs.roots_[0]
				&& roots_[1] == rhs.roots_[1] && roots_[2] == rhs.roots_[2]
				&& nPseudos_ == rhs.nPseudos_ && hasPseudo_ == rhs.hasPseudo_
				&& atmStartIdx_ == rhs.atmStartIdx_;
		isEql = isEql && nBond_ == rhs.nBond_ && nAngle_ == rhs.nAngle_
				&& nDihedral_ == rhs.nDihedral_ && hasTime_ == rhs.hasTime_;
		return isEql;
	}

	bool HasBond() const {
		return (nBond_ > 0);
	}

	bool HasAngle() const {
		return (nAngle_ > 0);
	}

	bool HasDihedral() const {
		return (nDihedral_ > 0);
	}

	bool HasTime() const {
		return hasTime_;
	}

	bool HasPseudo() const {
		return hasPseudo_;
	}

	void SetTime(bool t) {
		hasTime_ = t;
	}

	void SetBond(int nbnd) {
		nBond_ = nbnd;
	}

	int nBonds() const {
		return nBond_;
	}

	int nAngles() const {
		return nAngle_;
	}

	int nDihedrals() const {
		return nDihedral_;
	}

	int nPseudos() const {
		return nPseudos_;
	}

	void rootAtoms(int &r1, int &r2, int &r3) const {
		r1 = roots_[0];
		r2 = roots_[1];
		r3 = roots_[2];
	}

	int atomStartIdx() const {
		return atmStartIdx_;
	}

	void SetAngle(int nang) {
		nAngle_ = nang;
	}

	void SetDihedral(int ndih) {
		nDihedral_ = ndih;
	}

	void SetAtomStartIndex(int s) {
		atmStartIdx_ = s;
	}

	void SetRoots(int r1, int r2, int r3) {
		roots_[0] = r1;
		roots_[1] = r2;
		roots_[2] = r3;
	}

	void SetPseudo(bool b) {
		hasPseudo_ = b;
	}

	void SetPseudo(int num) {
		nPseudos_ = num;
	}

	void SetPseudoBonds(std::vector<long> const &vec) {
		std::vector<long>::const_iterator it;
		for (it = vec.begin(); it != vec.end(); ++it) {
			pseudos_.push_back(*it);
		}
	}

	void getPseudoBonds(std::vector<long> &vec) {
		std::vector<long>::const_iterator it;
		for (it = pseudos_.begin(); it != pseudos_.end(); ++it) {
			vec.push_back(*it);
		}
	}
};
#endif
