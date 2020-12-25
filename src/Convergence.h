/*!\brief This file defines methods which will be used for calculating sum of
 * entropy contributions for first order entropy and mutual information for second order.
 *
 * File:   Convergence.h
 * Author: shailesh
 *
 * Created on 5 August, 2016, 10:41 PM
 */

#ifndef CONVERGENCE_H_
#define CONVERGENCE_H_

#include<vector>

#include "configcentre.h"
#include "Centre.h"

/*!
 \fn void sumContrib1D(const u_int dim_len, const std::vector<double>& dataIn, double& entSum)

 A function for summing 1D entropy contributions
 \param[in] dim_len the number of DOFs of the type in subset
 \param[in] dataIn a constant vector of entropy contribution values for DOFs
 of type bond or angle or torsion
 \param[out] entSum the sum of entropy contributions 1D
 */
void sumContrib1D(const u_int dim_len, const std::vector<double> &dataIn,
		double &entSum);

/*!
 \fn void sumContrib2D(const u_int start_indx, const u_int count,
 const std::vector<double>& dataIn, double& entSum)

 A function for summing 1D entropy contributions
 \param[in] start_indx the index from where summing of data starts
 \param[in] count the number of elements of dataIn to be summed starting from start_indx
 \param[in] dataIn a constant vector of entropy contribution values for DOFs
 of type bond or angle or torsion
 \param[out] entSum the sum of entropy contributions 2D
 */
void sumContrib2D(const u_int start_indx, const u_int count,
		const std::vector<double> &dataIn, double &entSum);

/*!
 \fn void getMaxMI2D(const SubsetData& subset, const Neighbors& neigh,
 const char d1type, const char d2type, const std::vector<double> S1,
 const std::vector<double>& dataIn, std::vector<double>& dataMI)

 A function for getting pair of DOFs of given types with maximum Mutual Information
 \param[in] subset the object of class SubsetData used for guiding calculation
 \param[in] neigh the object of class Neighbors used for guiding calculation
 \param[in] d1type the tyep of first DOF in pair 'B', or 'A' or 'T'
 \param[in] d2type the tyep of second DOF in pair 'B', or 'A' or 'T'
 \param[in] S1 the vector of entropy contributions 1D
 \param[in] dataIn a constant vector of entropy contribution values for DOFs
 of type bond or angle or torsion
 \param[out] mi_data the array of Mutual Information values
 */
void getMaxMI2D(const SubsetData &subset, const Neighbors &neigh,
		const char d1type, const char d2type, const std::vector<double> S1,
		const std::vector<double> &dataIn, std::vector<double> &dataMI);

/*!
 \fn void getMaxMI2D(const SubsetData& subset, const Neighbors& neigh,
 const char d1type, const char d2type, const std::vector<double> S1,
 const std::vector<int>& block, const std::vector<double>& dataIn,
 const u_int dataOffset, std::vector<double>& dataMI)

 A function for getting pair of DOFs of given types with maximum Mutual Information
 \param[in] subset the object of class SubsetData used for guiding calculation
 \param[in] neigh the object of class Neighbors used for guiding calculation
 \param[in] d1type the tyep of first DOF in pair 'B', or 'A' or 'T'
 \param[in] d2type the tyep of second DOF in pair 'B', or 'A' or 'T'
 \param[in] S1 the vector of entropy contributions 1D
 \param[in] block the vector providing the range of dataIn values to be processed in 2D
 \param[in] dataIn a constant vector of entropy contribution values for DOFs
 of type bond or angle or torsion
 \param[in] dataOffset the linear offset for elements of dataIn values to be processed
 \param[out] mi_data the array of Mutual Information values
 */
void getMaxMI2D(const SubsetData &subset, const Neighbors &neigh,
		const char d1type, const char d2type, const std::vector<double> S1,
		const std::vector<int> &block, const std::vector<double> &dataIn,
		const u_int dataOffset, std::vector<double> &dataMI);

void getMaxMI2D(const std::map<u_int, u_int> &id2index_1d1,
		const std::map<u_int, u_int> &id2index_1d2,
		const std::vector<double> &S1, const std::vector<double> &S2,
		const std::vector<u_int> &index2id_2d,
		const std::vector<double> &dataIn, std::vector<double> &dataMI);

/*!
 \fn void getMaxMI2D(const SubsetData& subset, const Neighbors& neigh,
 const char d1type, const char d2type, const std::vector<double> S1,
 const std::vector<int>& block, const std::vector<double>& dataIn,
 const u_int dataOffset, std::vector<double>& dataMI)

 A function for getting pair of DOFs of given types with maximum Mutual Information
 \param[in] subset the object of class SubsetData used for guiding calculation
 \param[in] neigh the object of class Neighbors used for guiding calculation
 \param[in] d1type the tyep of first DOF in pair 'B', or 'A' or 'T'
 \param[in] d2type the tyep of second DOF in pair 'B', or 'A' or 'T'
 \param[in] S1 the vector of entropy contributions 1D for d1type
 \param[in] S2 the vector of entropy contributions 1D for d2type
 \param[in] block the vector providing the range of dataIn values to be processed in 2D
 \param[in] dataIn a constant vector of entropy contribution values for DOFs
 of type bond or angle or torsion
 \param[in] dataOffset the linear offset for elements of dataIn values to be processed
 \param[out] mi_data the array of Mutual Information values
 */
void getMaxCrossMI2D(const SubsetData &subset, const Neighbors &neigh,
		const char d1type, const char d2type, const std::vector<double> S1d1,
		const std::vector<double> S1d2, const std::vector<int> &block,
		const std::vector<double> &dataIn, const u_int dataOffset,
		std::vector<double> &mi_data);

/*!
 \fn void getMaxCrossMI2D(const SubsetData& subset, const Neighbors& neigh,
 const char d1type, const char d2type, const std::vector<double> S1d1,
 const std::vector<double> S1d2, const std::vector<double>& dataIn,
 std::vector<double>& mi_data)

 A function for getting pair of DOFs of given types with maximum Mutual Information
 \param[in] subset the object of class SubsetData used for guiding calculation
 \param[in] neigh the object of class Neighbors used for guiding calculation
 \param[in] d1type the tyep of first DOF in pair 'B', or 'A' or 'T'
 \param[in] d2type the tyep of second DOF in pair 'B', or 'A' or 'T'
 \param[in] S1 the vector of entropy contributions 1D for d1type
 \param[in] S2 the vector of entropy contributions 1D for d2type
 \param[in] dataIn a constant vector of entropy contribution values for DOFs
 of type bond or angle or torsion
 \param[out] mi_data the array of Mutual Information values
 */
void getMaxCrossMI2D(const SubsetData &subset, const Neighbors &neigh,
		const char d1type, const char d2type, const std::vector<double> S1d1,
		const std::vector<double> S1d2, const std::vector<double> &dataIn,
		std::vector<double> &mi_data);

#endif /* CONVERGENCE_H_ */
