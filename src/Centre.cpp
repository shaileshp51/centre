#include "Centre.h"
#include <cstdlib>
#include <limits>

std::vector<Edge> boruvkaMST(Graph &graph) {
	int V = graph.nodes.size();
	int E = graph.edges.size();
	std::vector<Edge> tree;
	std::vector<Edge> edge = graph.edges;

	std::vector<Component> components(V);

	// An array to store index of the cheapest edge of
	// subset.  The stored index for indexing array 'edge[]'
	std::vector<double> cheapest(V, -1.0);

	// Create V subsets with single elements
	for (int v = 0; v < V; ++v) {
		components[v].parent = v;
		components[v].rank = 0;
	}

	int numTrees = V; // Initially there are V different trees.
	int MSTweight = 0;

	// Keep combining components (or sets) until all
	// compnentes are not combined into single MST.
	while (numTrees > 1) {
		// Traverse through all edges and update
		// cheapest of every component
		int iterations_num_cheapest = 0;
		for (int i = 0; i < E; i++) {
			// Find components (or sets) of two corners
			// of current edge
			int set1 = findRoot(components, edge[i].source);
			int set2 = findRoot(components, edge[i].dest);

			// If two corners of current edge belong to
			// same set, ignore current edge
			if (set1 != set2) {
				// Else check if current edge is closer to previous
				// cheapest edges of set1 and set2
				if (cheapest[set1] == -1
						|| edge[cheapest[set1]].weight < edge[i].weight) {
					cheapest[set1] = i;
                    iterations_num_cheapest++;
                }

				if (cheapest[set2] == -1
						|| edge[cheapest[set2]].weight < edge[i].weight) {
					cheapest[set2] = i;
                    iterations_num_cheapest++;
                }
			}
		}
        // The next largets edges-weight equals one already added.
        // Lets break the deadlock by allowing addition of an equal weight edge.
		if (iterations_num_cheapest == 0) {
			for (int i = 0; i < E; i++) {
				// Find components (or sets) of two corners
				// of current edge
				int set1 = findRoot(components, edge[i].source);
				int set2 = findRoot(components, edge[i].dest);
				// Consider an edge only when two ends belog to different components.
				if (set1 != set2) {
					// Else check if current edge is closer to previous
					// cheapest edges of set1 and set2
					if (cheapest[set1] == -1
							|| edge[cheapest[set1]].weight <= edge[i].weight) {
						cheapest[set1] = i;
						iterations_num_cheapest++;
						break;
					}

					if (cheapest[set2] == -1
							|| edge[cheapest[set2]].weight <= edge[i].weight) {
						cheapest[set2] = i;
						iterations_num_cheapest++;
						break;
					}
				}
			}
		}
		if(iterations_num_cheapest==0)
			return tree;
		// Consider the above picked cheapest edges and add them
		// to MST
		for (int i = 0; i < V; i++) {
			// Check if cheapest for current set exists
			if (cheapest[i] != -1) {
				int set1 = findRoot(components, edge[cheapest[i]].source);
				int set2 = findRoot(components, edge[cheapest[i]].dest);
				if (set1 != set2) {
					MSTweight += edge[cheapest[i]].weight;
					tree.push_back(edge[cheapest[i]]);
					// Do a union of set1 and set2 and decrease number
					// of trees
					unionComponents(components, set1, set2);
					numTrees--;
				}
			}
		}
		for (int v = 0; v < V; ++v) {
			cheapest[v] = -1.0;
		}
	}
//   for (int i = 0; i < tree.size(); ++i) {
//      std::cout << graph.nodes[tree[i].source].type << ", " << graph.nodes[tree[i].source].id
//            << ", ";
//      std::cout << graph.nodes[tree[i].dest].type << ", " << graph.nodes[tree[i].dest].id << ", "
//            << tree[i].weight << std::endl;
//   }
	// printf("Weight of MST is %d\n", MSTweight);
	return tree;
}

int findRoot(std::vector<Component> &compnents, int i) {
	// find root and make root as parent of i
	// (path compression)
	if (compnents[i].parent != (u_int) i)
		compnents[i].parent = findRoot(compnents, compnents[i].parent);

	return compnents[i].parent;
}

void unionComponents(std::vector<Component> &comonents, int x, int y) {
	int xroot = findRoot(comonents, x);
	int yroot = findRoot(comonents, y);

	// Attach smaller rank tree under root of high
	// rank tree (Union by Rank)
	if (comonents[xroot].rank < comonents[yroot].rank)
		comonents[xroot].parent = yroot;
	else if (comonents[xroot].rank > comonents[yroot].rank)
		comonents[yroot].parent = xroot;

	// If ranks are same, then make one as root and
	// increment its rank by one
	else {
		comonents[yroot].parent = xroot;
		comonents[xroot].rank++;
	}
}

nodeMIST::nodeMIST() {
	d1type = BAT_t::NONE;
	d2type = BAT_t::NONE;
	d1indx = -1;
	d2indx = -1;
	MI = -std::numeric_limits<float>::max();
}

nodeMIST::nodeMIST(BAT_t d1t, BAT_t d2t, int d1idx, int d2idx, double I) {
	d1type = d1t;
	d2type = d2t;
	d1indx = d1idx;
	d2indx = d2idx;
	MI = I;
}

SubsetData::SubsetData(u_int n_bnd, u_int n_ang, u_int n_tor,
		BATSet entropy_workset) {
	this->entropy_workset = entropy_workset;
	if (entropy_workset & BATSet::B1D) {
		n_tot_bnd = n_bnd;
		bonds.reserve(n_tot_bnd);
		for (u_int i = 0; i < n_tot_bnd; i++) {
			bonds.push_back(i);
		}
		if (n_tot_bnd > 0) {
			setHasBonds(true);
		}
	} else {
		n_tot_bnd = 0;
	}

	if (entropy_workset & BATSet::A1D) {
		n_tot_ang = n_ang;
		angles.reserve(n_tot_ang);
		for (u_int i = 0; i < n_tot_ang; i++) {
			angles.push_back(i);
		}
		if (n_tot_ang > 0) {
			setHasAngles(true);
		}
	} else {
		n_tot_ang = 0;
	}

	if (entropy_workset & BATSet::D1D) {
		n_tot_tor = n_tor;
		torsions.reserve(n_tot_tor);
		for (u_int i = 0; i < n_tot_tor; i++) {
			torsions.push_back(i);

		}
		if (n_tot_tor > 0) {
			setHasTorsions(true);
		}
	} else {
		n_tot_tor = n_tor;
	}
	this->isgood_ = true;
}

SubsetData::SubsetData(std::string const &fname, BATSet entropy_workset) {
	n_tot_bnd = 0;
	n_tot_ang = 0;
	n_tot_tor = 0;
	setHasBonds(false);
	setHasAngles(false);
	setHasTorsions(false);

	std::string str_bond("bond");
	std::string str_angle("angle");
	std::string str_torsion("torsion");
	std::ifstream file(fname.c_str());
	std::string str;

	bool isOK = true;

	while (isOK && (std::getline(file, str, ';'))) {
		if (str[0] == '#') {
			std::cout << "Read a comment line\n" << str << "\n";
		} else {
			std::istringstream stream(str);
			std::string temp;
			u_int total, n_member, v;

			stream >> temp;
			stream >> total;
			stream >> n_member;

			std::vector<u_int> vec;
			for (size_t i = 0; i < n_member; ++i) {
				stream >> v;
				if (!stream.fail()) {
					vec.push_back(v);
				} else {
					std::cout << "SubsetDataReadError:: Expected (" << n_member
							<< ") " << temp << " values found (" << vec.size()
							<< ")\n";
					stream.clear();
					isOK = false;
					break;
				}
			}
			if (isOK && str_bond == temp && (entropy_workset & BATSet::B1D)) {
				n_tot_bnd = total;
				bonds.assign(vec.begin(), vec.end());
				setHasBonds(true);
			}
			if (isOK && str_angle == temp && (entropy_workset & BATSet::A1D)) {
				n_tot_ang = total;
				angles.assign(vec.begin(), vec.end());
				setHasAngles(true);
			}
			if (isOK && str_torsion == temp
					&& (entropy_workset & BATSet::D1D)) {
				n_tot_tor = total;
				torsions.assign(vec.begin(), vec.end());
				setHasTorsions(true);
			}
		}
		file >> std::ws; // read and ignore newline/leading/trailing spaces
	}
	this->isgood_ = isOK;
}

void SubsetData::push_bond(u_int value) {
	bonds.push_back(value);
}

void SubsetData::push_bonds(std::vector<u_int> const &vec) {
	std::vector<u_int>::const_iterator it;
	for (it = vec.begin(); it != vec.end(); ++it) {
		bonds.push_back(*it);
	}
}

void SubsetData::push_angle(u_int value) {
	angles.push_back(value);
}

void SubsetData::push_angles(std::vector<u_int> const &vec) {
	std::vector<u_int>::const_iterator it;
	for (it = vec.begin(); it != vec.end(); ++it) {
		angles.push_back(*it);
	}
}

void SubsetData::push_torsion(u_int value) {
	torsions.push_back(value);
}

void SubsetData::push_torsions(std::vector<u_int> const &vec) {
	std::vector<u_int>::const_iterator it;
	for (it = vec.begin(); it != vec.end(); ++it) {
		torsions.push_back(*it);
	}
}

const u_int SubsetData::getBATTypeId(BAT_t bat_type, u_int id) const {
	u_int found_id = std::numeric_limits<u_int>::max();
	switch (bat_type) {
	case BAT_t::BOND:
		found_id = bonds[id];
		break;
	case BAT_t::ANGLE:
		found_id = angles[id];
		break;
	case BAT_t::DIHEDRAL:
		found_id = torsions[id];
		break;
	default:
		std::cout << "UNKNOWN BAT_t provided, can not continue anymore."
				<< std::endl;
		exit(0);
	}
	return found_id;
}

Neighbors::Neighbors(u_int n_bnd, u_int n_ang, u_int n_tor) {
	if (n_bnd > 0) {
		setHasBonds(true);
	}
	if (n_ang > 0) {
		setHasAngles(true);
	}
	if (n_tor > 0) {
		setHasTorsions(true);
	}
	// TODO: Check if item already exists in list, insert only if not
	for (u_int i = 0; i < n_bnd - 1; i++) {
		std::vector<u_int> bnd;
		for (u_int j = i + 1; j < n_bnd; j++) {
			bnd.push_back(j);
		}
		bonds[i] = bnd;
	}

	for (u_int i = 0; i < n_bnd; i++) {
		std::vector<u_int> ba;
		for (u_int j = 0; j < n_ang; j++) {
			ba.push_back(j);
		}
		bacross[i] = ba;

		std::vector<u_int> bd;
		for (u_int j = 0; j < n_tor; j++) {
			bd.push_back(j);
		}
		bdcross[i] = bd;
	}

	for (u_int i = 0; i < n_ang - 1; i++) {
		std::vector<u_int> ang;
		for (u_int j = i + 1; j < n_ang; j++) {
			ang.push_back(j);
		}
		angles[i] = ang;
	}

	for (u_int i = 0; i < n_ang; i++) {
		std::vector<u_int> ad;
		for (u_int j = 0; j < n_tor; j++) {
			ad.push_back(j);
		}
		adcross[i] = ad;
	}

	for (u_int i = 0; i < n_tor - 1; i++) {
		std::vector<u_int> tor;
		for (u_int j = i + 1; j < n_tor; j++) {
			tor.push_back(j);
		}
		torsions[i] = tor;
	}
	this->isgood_ = true;
}

Neighbors::Neighbors(SubsetData subset, BATSet entropy_workset) {

	if (subset.isGood()) {
		this->entropy_workset = entropy_workset;
		u_int bnd_size = subset.bondsSize();
		u_int ang_size = subset.anglesSize();
		u_int tor_size = subset.torsionsSize();

		setHasAngles(subset.getHasAngles());
		setHasTorsions(subset.getHasTorsions());

		if (entropy_workset & BATSet::B1D) {
			setHasBonds(subset.getHasBonds());
			if (entropy_workset & BATSet::BB2D) {
				for (u_int i = 0; i < bnd_size - 1; i++) {
					std::vector<u_int> b_neigh;
					for (u_int j = i + 1; j < bnd_size; j++) {
						b_neigh.push_back(subset.getBonds()[j]);
					}
					setBond(subset.getBonds()[i], b_neigh);
				}
			}
			if (entropy_workset & BATSet::BA2D) {
				for (u_int i = 0; i < bnd_size; i++) {
					std::vector<u_int> a_neigh;
					for (u_int j = 0; j < ang_size; j++) {
						a_neigh.push_back(subset.getAngles()[j]);
					}
					setBacross(subset.getBonds()[i], a_neigh);
				}
			}
			if (entropy_workset & BATSet::BD2D) {
				for (u_int i = 0; i < bnd_size; i++) {
					std::vector<u_int> d_neigh;
					for (u_int j = 0; j < tor_size; j++) {
						d_neigh.push_back(subset.getTorsions()[j]);
					}
					setBdcross(subset.getBonds()[i], d_neigh);
				}
			}
		}

		if (entropy_workset & BATSet::A1D) {
			if (entropy_workset & BATSet::AA2D) {
				for (u_int i = 0; i < ang_size - 1; i++) {
					std::vector<u_int> a_neigh;
					for (u_int j = i + 1; j < ang_size; j++) {
						a_neigh.push_back(subset.getAngles()[j]);
					}
					setAngle(subset.getAngles()[i], a_neigh);
				}
			}
			if (entropy_workset & BATSet::AD2D) {
				for (u_int i = 0; i < ang_size; i++) {
					std::vector<u_int> d_neigh;
					for (u_int j = 0; j < tor_size; j++) {
						d_neigh.push_back(subset.getTorsions()[j]);
					}
					setAdcross(subset.getAngles()[i], d_neigh);
				}
			}
		}
		if (entropy_workset & BATSet::D1D) {
			if (entropy_workset & BATSet::DD2D) {
				for (u_int i = 0; i < tor_size - 1; i++) {
					std::vector<u_int> t_neigh;
					for (u_int j = i + 1; j < tor_size; j++) {
						t_neigh.push_back(subset.getTorsions()[j]);
					}
					setTorsion(subset.getTorsions()[i], t_neigh);
				}
			}
		}
	}
	this->isgood_ = subset.isGood();
}

Neighbors::Neighbors(std::string const &fname, BATSet entropy_workset) {
	this->entropy_workset = entropy_workset;
 	setHasBonds(false);
	setHasAngles(false);
	setHasTorsions(false);
	std::string str_bond("bond");
	std::string str_angle("angle");
	std::string str_torsion("torsion");
	std::string str_bacross("bacross");
	std::string str_bdcross("bdcross");
	std::string str_adcross("adcross");

	std::ifstream file(fname.c_str());
	std::string str;
	bool isOK = true;

	while (isOK && std::getline(file, str, ';')) {
		if (str[0] == '#') {
			std::cout << "Read a comment line\n" << str << "\n";
		} else {
			std::istringstream stream(str);
			std::string temp;
			u_int _index;
			u_int n_neigh;
			u_int v;
			stream >> temp;
			stream >> _index;
			stream >> n_neigh;

			if (n_neigh == 0) {
				continue;
			}

			std::vector<u_int> vec;
			for (u_int i = 0; i < n_neigh; ++i) {
				stream >> v;
				if (!stream.fail()) {
					vec.push_back(v);
				} else {
					std::cout << "NeighborsReadError:: Expected (" << n_neigh
							<< ") values found (" << vec.size() << ")\n";
					stream.clear();
					isOK = false;
					break;
				}
			}

			if (isOK && str_bond == temp && (entropy_workset & BATSet::BB2D)) {
				setHasBonds(true);
				std::pair<std::map<u_int, std::vector<u_int>>::iterator, bool> ret;
				ret = bonds.insert(
						std::pair<u_int, std::vector<u_int>>(_index, vec));
				if (ret.second == false) {
					std::cout << "element " << _index << " already existed";
					//std::cout << " with a value of " << ret.first->second << '\n';
				}
			} else {
				setHasBonds(false);
			}
			if (isOK && str_angle == temp && (entropy_workset & BATSet::AA2D)) {
				setHasAngles(true);
				std::pair<std::map<u_int, std::vector<u_int>>::iterator, bool> ret;
				ret = angles.insert(
						std::pair<u_int, std::vector<u_int>>(_index, vec));
				if (ret.second == false) {
					std::cout << "element " << _index << " already existed";
					//std::cout << " with a value of " << ret.first->second << '\n';
				}
			} else {
				setHasAngles(false);
			}
			if (isOK && str_torsion == temp
					&& (entropy_workset & BATSet::DD2D)) {
				setHasTorsions(true);

				std::pair<std::map<u_int, std::vector<u_int>>::iterator, bool> ret;
				ret = torsions.insert(
						std::pair<u_int, std::vector<u_int>>(_index, vec));
				if (ret.second == false) {
					std::cout << "element " << _index << " already existed";
					//std::cout << " with a value of " << ret.first->second << '\n';
				}
			} else {
				setHasTorsions(false);
			}

			if (isOK && str_bacross == temp
					&& (entropy_workset & BATSet::BA2D)) {
				setHasBacross(true);

				std::pair<std::map<u_int, std::vector<u_int>>::iterator, bool> ret;
				ret = bacross.insert(
						std::pair<u_int, std::vector<u_int>>(_index, vec));
				if (ret.second == false) {
					std::cout << "element " << _index << " already existed";
					//std::cout << " with a value of " << ret.first->second << '\n';
				}
			} else {
				setHasBacross(false);
			}

			if (isOK && str_bdcross == temp
					&& (entropy_workset & BATSet::BD2D)) {
				setHasBdcross(true);

				std::pair<std::map<u_int, std::vector<u_int>>::iterator, bool> ret;
				ret = bdcross.insert(
						std::pair<u_int, std::vector<u_int>>(_index, vec));
				if (ret.second == false) {
					std::cout << "element " << _index << " already existed";
					//std::cout << " with a value of " << ret.first->second << '\n';
				}
			} else {
				setHasBdcross(false);
			}

			if (isOK && str_adcross == temp
					&& (entropy_workset & BATSet::AD2D)) {
				setHasAdcross(true);

				std::pair<std::map<u_int, std::vector<u_int>>::iterator, bool> ret;
				ret = adcross.insert(
						std::pair<u_int, std::vector<u_int>>(_index, vec));
				if (ret.second == false) {
					std::cout << "element " << _index << " already existed";
					//std::cout << " with a value of " << ret.first->second << '\n';
				}
			} else {
				setHasAdcross(false);
			}
		}
		file >> std::ws; // read and ignore newline/leading/trailing spaces
	}
	this->isgood_ = isOK;
}

void Neighbors::bondKeys(std::vector<u_int> &kints) const {
	for (auto &imap : bonds) {
		kints.push_back(imap.first);
	}
}

void Neighbors::angleKeys(std::vector<u_int> &kints) const {
	for (auto const &imap : angles) {
		kints.push_back(imap.first);
	}
}

void Neighbors::torsionKeys(std::vector<u_int> &kints) const {
	for (auto const &imap : torsions) {
		kints.push_back(imap.first);
	}
}

void Neighbors::bacrossKeys(std::vector<u_int> &kints) const {
	for (auto const &imap : bacross) {
		kints.push_back(imap.first);
	}
}

void Neighbors::bdcrossKeys(std::vector<u_int> &kints) const {
	for (auto const &imap : bdcross) {
		kints.push_back(imap.first);
	}
}

void Neighbors::adcrossKeys(std::vector<u_int> &kints) const {
	for (auto const &imap : adcross) {
		kints.push_back(imap.first);
	}
}

const u_int Neighbors::bondsSize() const {
	return bonds.size();
}

const u_int Neighbors::anglesSize() const {
	return angles.size();
}

const u_int Neighbors::torsionsSize() const {
	return torsions.size();
}

const u_int Neighbors::baCrossSize() const {
	return bacross.size();
}

const u_int Neighbors::bdCrossSize() const {
	return bdcross.size();
}

const u_int Neighbors::adCrossSize() const {
	return adcross.size();
}

const u_int Neighbors::baNeighSize(u_int id) const {
	return (bacross.at(id)).size();
}

const u_int Neighbors::bdNeighSize(u_int id) const {
	return (bdcross.at(id)).size();
}

const u_int Neighbors::adNeighSize(u_int id) const {
	return (adcross.at(id)).size();
}

const u_int Neighbors::bondNeighSize(u_int id) const {
	return (bonds.at(id)).size();
}

const u_int Neighbors::angleNeighSize(u_int id) const {
	return angles.at(id).size();
}

const u_int Neighbors::torsionNeighSize(u_int id) const {
	return torsions.at(id).size();
}

const std::vector<u_int>& Neighbors::getAngle(const u_int id) const {
	return angles.at(id);
}

void Neighbors::setAngle(const u_int id, const std::vector<u_int> &nei) {
	this->angles[id] = nei;
}

const std::vector<u_int>& Neighbors::getBond(const u_int id) const {
	return bonds.at(id);
}

void Neighbors::setBond(const u_int id, const std::vector<u_int> &nei) {
	this->bonds[id] = nei;
}

const std::vector<u_int>& Neighbors::getTorsion(const u_int id) const {
	return torsions.at(id);
}

void Neighbors::setTorsion(const u_int id, const std::vector<u_int> &nei) {
	this->torsions[id] = nei;
}

void Neighbors::setAdcross(const u_int id, const std::vector<u_int> &nei) {
	this->adcross[id] = nei;
}

const std::vector<u_int>& Neighbors::getBacross(u_int id) const {
	return bacross.at(id);
}

const std::vector<u_int>& Neighbors::getBdcross(u_int id) const {
	return bdcross.at(id);
}

const std::vector<u_int>& Neighbors::getAdcross(u_int id) const {
	return adcross.at(id);
}

const std::vector<u_int>& Neighbors::getDimTypeNeighs(const BATSet dimType,
		const u_int id) const {
	const std::vector<u_int> *dtypneigh;
	switch (dimType) {
	case BATSet::BB2D:
		dtypneigh = &bonds.at(id);
		break;
	case BATSet::AA2D:
		dtypneigh = &angles.at(id);
		break;
	case BATSet::DD2D:
		dtypneigh = &torsions.at(id);
		break;
	case BATSet::BA2D:
		dtypneigh = &bacross.at(id);
		break;
	case BATSet::BD2D:
		dtypneigh = &bdcross.at(id);
		break;
	case BATSet::AD2D:
		dtypneigh = &adcross.at(id);
		break;
	default:
		std::cout << "Provided cases is not a valid DimTypePair." << std::endl;
		exit(0);
	}
	return *dtypneigh;
}

void Neighbors::setBacross(const u_int id, const std::vector<u_int> &nei) {
	this->bacross[id] = nei;
}

void Neighbors::setBdcross(const u_int id, const std::vector<u_int> &nei) {
	this->bdcross[id] = nei;
}

bool Neighbors::isHasAdcross() const {
	return has_adcross_;
}

void Neighbors::setHasAdcross(bool hasAdcross) {
	has_adcross_ = hasAdcross;
}

bool Neighbors::isHasBacross() const {
	return has_bacross_;
}

void Neighbors::setHasBacross(bool hasBacross) {
	has_bacross_ = hasBacross;
}

bool Neighbors::isHasBdcross() const {
	return has_bdcross_;
}

void Neighbors::setHasBdcross(bool hasBdcross) {
	has_bdcross_ = hasBdcross;
}

std::ostream& operator<<(std::ostream &os, Neighbors &neig) {

	if (neig.getHasBonds()) {
		// for (const auto& pair : neig.getBonds()) {
		// 	os << "bond " << pair.first << " " << pair.second.size();
		// 	for ( auto& vv :  pair.second) {
		// 	os << " " << vv;
		// 	}
		// 	os << ";" << std::endl;
		// }
	}

	if (neig.getHasAngles()) {
		std::map<u_int, std::vector<u_int>>::const_iterator it;
    //  for (auto it = neig.getAngles().begin(); it != neig.getAngles().end(); ++it) {
    //     auto nei_v = (*it).second;
    //     os << "angle " << (*it).first << " " << nei_v.size();
    //     std::vector<u_int>::iterator it_v;
    //     for (it_v = nei_v.begin(); it_v != nei_v.end(); ++it_v) {
    //        os << " " << *it_v;
    //     }
    //     os << ";" << std::endl;
    //  }
	}

	if (neig.getHasTorsions()) {
		std::map<u_int, std::vector<u_int>>::const_iterator it;
    //  for (auto it = neig.getTorsions().begin(); it != neig.getTorsions().end(); ++it) {
    //     std::vector<u_int> nei_v = (*it).second;
    //     os << "torsion " << (*it).first << " " << nei_v.size();
    //     std::vector<u_int>::iterator it_v;
    //     for (it_v = nei_v.begin(); it_v != nei_v.end(); ++it_v) {
    //        os << " " << *it_v;
    //     }
    //     os << ";" << std::endl;
    //  }
	}
	return os;
}

