// $Id: BioMolecule.hh,v 1.2 2016/08/08 12:41:56 mmann Exp $
#ifndef BIU_BIOMOLECULE_H_
#define BIU_BIOMOLECULE_H_

#include "biu/Alphabet.hh"

namespace biu
{
		/*! A type for sequence representation */	
	typedef Alphabet::Sequence Sequence;
	
		/*! A type for structure representation */	
	typedef Alphabet::Sequence Structure;
	
		/*! A BioMolecule represents a state in a biological
		 *  energy landscape. 
		 * 
		 * It consists of its sequence and structure information.
		 * 
		 * @author Martin Mann, Sebastian Will, Andreas Richter
		 */
	class BioMolecule
	{
	public:
		BioMolecule()
		{}
		
		virtual ~BioMolecule()
		{}
		
		virtual Sequence	getSequence() const = 0;
		
		virtual Structure	getStructure() const = 0;
		
			//! Returns the length of the Biomolecule, i.e. the number of 
			//! monomers.
		virtual size_t		getLength() const = 0;
		
			/*! Returns the specific energy of the BioMolecule. */
		virtual double	getEnergy() const = 0;
		
			/*! Returns whether or not the BioMolecule is valid
			 *  concerning some internal criterias
			 *  (e.g. selfavoidingness of lattice proteins).
			 */
		virtual bool	isValid() const = 0;
		
			//! Returns a specific std::string representation of this
			//!  BioMolecule. 
		virtual std::string	getStringRepresentation() const = 0;	
	};

} // namespace biu

#endif /*BIOMOLECULE_H_*/
