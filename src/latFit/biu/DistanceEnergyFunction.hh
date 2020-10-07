// $Id: DistanceEnergyFunction.hh,v 1.2 2016/08/08 12:41:56 mmann Exp $
#ifndef BIU_DISTANCEENERGYFUNCTION_HH_
#define BIU_DISTANCEENERGYFUNCTION_HH_

#include "biu/Alphabet.hh"
#include "biu/LatticeModel.hh"

namespace biu
{

	/*! 
	 * An energy function that evaluates the energy contribution of two monomers
	 * depending on their type and distance.
	 * 
	 * It provides the distance based energy for elements of a specified 
	 * alphabet.
	 * 
	 * @author Martin Mann
	 */
	class DistanceEnergyFunction
	{
	public:
		virtual ~DistanceEnergyFunction();
		
		
			/*! Returns the energy contribution of two elements of the alphabet
			 * that are in a certain distance. 
			 * 
			 * @param seq_i,seq_j have to be elements of the specified alphabet.
			 * @param distance the distance of the two monomers to evaluate
			 * @return the energy contribution
			 */
		virtual
		double getEnergy(	const Alphabet::AlphElem& seq_i, 
							const Alphabet::AlphElem& seq_j,
							const double & distance ) const = 0;
									
			/*! Returns the energy contribution of two elements of the alphabet
			 * with given coordinates. 
			 * 
			 * @param seq_i,seq_j have to be elements of the specified alphabet.
			 * @param cor_i,cor_j the corresponding coordinates of the two 
			 *                    monomers to evaluate
			 * @return the energy contribution
			 */
		virtual
		double getEnergy(	const Alphabet::AlphElem & seq_i, 
							const Alphabet::AlphElem & seq_j,
							const IntPoint & cor_i,
							const IntPoint & cor_j  ) const = 0;
									
			/*! Returns the energy contribution of two elements of the alphabet
			 * with given coordinates. 
			 * 
			 * @param seq_i,seq_j have to be elements of the specified alphabet.
			 * @param cor_i,cor_j the corresponding coordinates of the two 
			 *                    monomers to evaluate
			 * @return the energy contribution
			 */
		virtual
		double getEnergy(	const Alphabet::AlphElem & seq_i, 
							const Alphabet::AlphElem & seq_j,
							const DblPoint & cor_i,
							const DblPoint & cor_j  ) const = 0;
									
			/*! 
			 * Access to the alphabet this energy function is based on.
			 * 
			 * @return the alphabet in use 
			 * */
		virtual
		const Alphabet* const getAlphabet() const = 0;
		
		virtual
		bool operator == ( const DistanceEnergyFunction& toCompare ) const = 0;
		virtual
		bool operator != ( const DistanceEnergyFunction& toCompare ) const = 0;
		
	};

} // biu

#include "biu/Matrix.hh"

namespace biu {

	/*! The contact energy matrix for a ContactEnergyFunction object. */
	typedef biu::Matrix<double> EnergyMatrix;

} // biu


namespace biu
{


		/*! This class implements a contact based energy function. 
		 * 
		 * It provides the contact energy for elements of a specified alphabet.
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	class ContactEnergyFunction : public DistanceEnergyFunction
	{
	private:
	
			//! stores the allowed sequence elements to evaluate
		const Alphabet* const alphabet;
			//! the energy table that contains the energy contributions of 
			//! contacts
		const EnergyMatrix* const energyMat;
			//! the lattice model used to check if two positions are neighbored
		const LatticeModel* const lattice;
			//! the length of the first base vector of the lattice to allow for
			//! a distance based contact evaluation
		const double firstBaseVecLength;

	public:
	
			/*! 
			 * Constructs a new contact energy function for the given alphabet.
			 *  
			 * The alphabet elements of type AlphElem are mapped to
			 * the EnergyMatrix indices via the getIndex(..) of the
			 * Alphabet class.
			 * 
			 * @param alphabet the alphabet the energy function supports
			 * @param energyMat the contact energy matrix used. NOTE: has to 
			 *        have the same dimensions as the alphabet has elements!
			 * @param lattice the lattice model used to check if two positions 
			 *        are neighbored
			 */
		ContactEnergyFunction(	const Alphabet* const alphabet, 
								const EnergyMatrix* const energyMat,
								const LatticeModel* const lattice);
								
		virtual ~ContactEnergyFunction();
			
			/*! Returns the contact energy of two elements of the alphabet. 
			 * 
			 * @param first,second have to be elements of the specified alphabet.*/
		virtual
		double getContactEnergy(	const Alphabet::AlphElem& first, 
									const Alphabet::AlphElem& second) const;
									
			/*! Returns the energy contribution of two elements of the alphabet
			 * that are in a certain distance. 
			 * 
			 * @param seq_i,seq_j have to be elements of the specified alphabet.
			 * @param distance the distance of the two monomers to evaluate
			 * @return the energy contribution
			 */
		virtual
		double getEnergy(	const Alphabet::AlphElem& seq_i, 
							const Alphabet::AlphElem& seq_j,
							const double & distance ) const;
									
			/*! Returns the energy contribution of two elements of the alphabet
			 * with given coordinates. 
			 * 
			 * @param seq_i,seq_j have to be elements of the specified alphabet.
			 * @param cor_i,cor_j the corresponding coordinates of the two 
			 *                    monomers to evaluate
			 * @return the energy contribution
			 */
		virtual
		double getEnergy(	const Alphabet::AlphElem & seq_i, 
							const Alphabet::AlphElem & seq_j,
							const IntPoint & cor_i,
							const IntPoint & cor_j  ) const;
									
			/*! Returns the energy contribution of two elements of the alphabet
			 * with given coordinates. 
			 * 
			 * @param seq_i,seq_j have to be elements of the specified alphabet.
			 * @param cor_i,cor_j the corresponding coordinates of the two 
			 *                    monomers to evaluate
			 * @return the energy contribution
			 */
		virtual
		double getEnergy(	const Alphabet::AlphElem & seq_i, 
							const Alphabet::AlphElem & seq_j,
							const DblPoint & cor_i,
							const DblPoint & cor_j  ) const;
									
			/*! Returns the alphabet this energy function is based on. */
		virtual
		const Alphabet* const getAlphabet() const { 
			return alphabet; 
		}
		
		virtual
		bool operator == (const DistanceEnergyFunction& ef2) const;
		virtual
		bool operator != (const DistanceEnergyFunction& ef2) const;
		
	}; // class

} // namespace biu


namespace biu
{

	/*! 
	 * An distance based energy function that does a discretised energy 
	 * evaluation. For each distance interval a fixed energy table is used.
	 * 
	 * @author Martin Mann
	 */
	class IntervalEnergyFunction : public DistanceEnergyFunction
	{
		
	protected:
		
			//! the alphabet that stores all allowed sequence elements
		const Alphabet* const alphabet;

			//! Contains the energy matrices in use for each interval
		std::vector<const biu::EnergyMatrix*> energyMat;
		
			//! upper distance bounds of the intervals 
		std::vector< double > intervalMax;
		
		
	public:
		
		
			/*!
			 * Creates and initialises a new distance interval energy function.
			 * 
			 * @param alphabet the alphabet the energy function supports
			 */
		IntervalEnergyFunction(	const Alphabet* const alphabet );

			/*!
			 * Copy construction
			 * @param toCopy the object to copy
			 */
		IntervalEnergyFunction(	const IntervalEnergyFunction& toCopy );
		
		virtual ~IntervalEnergyFunction();
		
			
			/*! Returns the energy contribution of two elements of the alphabet
			 * that are in a certain distance. 
			 * 
			 * @param seq_i,seq_j have to be elements of the specified alphabet.
			 * @param distance the distance of the two monomers to evaluate
			 * @return the energy contribution
			 */
		virtual
		double getEnergy(	const Alphabet::AlphElem& seq_i, 
							const Alphabet::AlphElem& seq_j,
							const double & distance ) const;
									
			/*! Returns the energy contribution of two elements of the alphabet
			 * with given coordinates. 
			 * 
			 * @param seq_i,seq_j have to be elements of the specified alphabet.
			 * @param cor_i,cor_j the corresponding coordinates of the two 
			 *                    monomers to evaluate
			 * @return the energy contribution
			 */
		virtual
		double getEnergy(	const Alphabet::AlphElem & seq_i, 
							const Alphabet::AlphElem & seq_j,
							const IntPoint & cor_i,
							const IntPoint & cor_j  ) const;
									
			/*! Returns the energy contribution of two elements of the alphabet
			 * with given coordinates. 
			 * 
			 * @param seq_i,seq_j have to be elements of the specified alphabet.
			 * @param cor_i,cor_j the corresponding coordinates of the two 
			 *                    monomers to evaluate
			 * @return the energy contribution
			 */
		virtual
		double getEnergy(	const Alphabet::AlphElem & seq_i, 
							const Alphabet::AlphElem & seq_j,
							const DblPoint & cor_i,
							const DblPoint & cor_j  ) const;
									
			/*! 
			 * Access to the alphabet this energy function is based on.
			 * 
			 * @return the alphabet in use 
			 * */
		virtual
		const Alphabet* const getAlphabet() const;
		
		bool operator == (const DistanceEnergyFunction& cef2) const;
		bool operator != (const DistanceEnergyFunction& cef2) const;
		

		
			/**
			 * Adds a contact energy matrix to the energy function that is used
			 * for distance values that are greater than the upper bound of the
			 * last added interval and smaller or equal to the given upperBound. 
			 * 
			 * @param energies the contact energy matrices to use for this
			 *                 interval 
			 *                 NOTE: has to have the same 
			 *                 dimensions as the used alphabet has elements!
			 * @param upperBound the upper distance bound of this interval. 
			 *        NOTE: Has to be at higher than last upper bound so far!
			 * @return the index of this interval, i.e. the last interval index 
			 */
		virtual 
		size_t
		addInterval(	const biu::EnergyMatrix& energies,
						const double upperBound );
		
			/**
			 * Access to the number of intervals in use.
			 * @return the number of intervals
			 */
		virtual
		size_t
		getIntervalNum(void) const;
		
		
			/**
			 * Access to the upper interval bound.
			 * @param index the index of the interval of intrest (has to be less
			 *              than getIntervalNum())
			 * @return the upper bound of the specified interval
			 */
		virtual
		double
		getIntervalMax(const size_t index) const;
		
		
			/**
			 * Access to the energy table used for the given interval.
			 * @param index the index of the interval of intrest (has to be less
			 *              than getIntervalNum())
			 * @return the energy table of the specified interval
			 */
		virtual
		const biu::EnergyMatrix* const
		getIntervalMatrix(const size_t index) const;
		
		
			/**
			 * Access to the interval index of a given distance.
			 * @param distance the distance to check the interval for
			 * @return the interval index the distance falls into OR UINT_MAX if
			 *         the distance is higher than the highest distance handled
			 */
		virtual
		size_t
		getInterval( double distance) const;
	};

} // biu



#endif /*DISTANCEENERGYFUNCTION_HH_*/
