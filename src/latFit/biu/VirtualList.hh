// $Id: VirtualList.hh,v 1.2 2016/08/08 12:42:01 mmann Exp $
/************************************************************
 * 
 * Implementation of virtual lists that can be iterated STL-ish
 *
 * @author Sebastian Will <will@@informatik.uni-freiburg.de>
 *
 * Copyright Sebastian Will, 2006
 *
 ************************************************************/

#ifndef BIU_VIRTUAL_LIST_H
#define BIU_VIRTUAL_LIST_H

namespace biu {

    //------------------------------------------------------------
    // CLASS VirtualList
    
	//! Implementation of virtual lists that can be iterated STL-ish
    //!
    //! An example of usage is provided in testVirtualList.cc
    //!
    //! The elements of virtual lists are constructed on demand and exist
	//! only as long as an iterator points to them
    //!
    //! The class VirtualList<T> is an abstract class that serves
    //! as base class for concrete virtual lists of elements T.
    //!
    //! Child classes have to implement the abstract methods
    //!
    //!   virtual T* next(ItState *itstate, T* elem) const = 0;
    //!   virtual T* first(ItState **itstate) const = 0;
    //!
    //! and usually refine the iterator state (class ItState) by
    //! deriving from the nested class ItState.
    //!
	//! class T must implement a clone() function
    //!
	//! @author Sebastian Will <will@@informatik.uni-freiburg.de>
	//!
	template<class T> class VirtualList {

	public:
		//------------------------------------------------------------
		// Nested CLASS ItState
		
		//! State of an iterator. Usually the information in state is
		//! used for computing the next list entry 
		class ItState { 
			friend class VirtualList<T>;
		protected:
			ItState() {}
		public:
			ItState(const ItState &) {}
			virtual ~ItState();
			
			//! iterator states must be clonable
			virtual ItState * clone() const {return new ItState(*this);}
		};
  
		//------------------------------------------------------------
		// Nested CLASS Iterator
		
		//! the iterator for traversing a list
		//!
		//! the iterator can be used in an STL-ish way:
		//!  e.g. for(Iterator it=list->begin(); list->end()!=it; ++it)
		//!               {...*it...}
		//!
		//! since the iterator holds a pointer to its associated list and uses
		//! virtual methods of VirtualList<T>, it behaves "polymorph" due to
		//! the polymorphism of its associated list
		//!
		class Iterator {
			friend class VirtualList<T>;
		public:
			const VirtualList<T> *parent; //!< virtual list of iterator 
			ItState *itstate;             //!< information on iterator state
			T *elem;                      //!< element referenced by iterator
		private:
			//! Constructor that is used in begin() and end()
			//!
			//! Since only VirtualList<T> should use this, this constructor is private
			Iterator( const VirtualList<T> *parent_,
					  ItState *itstate_,
					  T* elem_ )
				: parent(parent_), itstate(itstate_), elem(elem_) {}

		public:
			//! copying iterators. Since element and iterator state
			//! are cloned we get an independent iterator
			Iterator(const Iterator &it)
				: parent(it.parent),
				  itstate(it.itstate->clone()),
				  elem(it.elem->clone())
				{}
			
			virtual ~Iterator();
			
			//! switch iterator to point to next list element
			const Iterator &operator ++() {
				T* nextElem = parent->next(itstate,elem);
				if (0==nextElem && 0!=elem) delete elem;
				elem=nextElem;
				return *this;
			}
        
			//! compare iterators, only two end-iterators are equal
			//! this behavior allows for loop-termination,
			//! but differs from the standard behavior of STL-iterators
			bool operator != (const Iterator &it) const {
				return this->elem!=0 || it.elem!=0; 
			}
			
			//! assigning iterators. Since element and iterator state
			//! are cloned we get an independent iterator
			Iterator &operator =(const Iterator &it) {
				parent=it.parent;
			
				if (itstate) delete itstate;
				itstate=(it.itstate)?it.itstate->clone():0;
			
				if (elem) delete elem;
				elem=(it.elem)?it.elem->clone():0;
			
				return *this;
			}

			//! returns the element referenced by the iterator
			const T& operator *() const {
				return *elem;
			}

			//! accesses the element referenced by the iterator
			const T* operator ->() const {
				return elem;
			}

		};
	  
	  VirtualList() {}
	  virtual ~VirtualList() {}

		//! @returns iterator pointing at first element of list 
		Iterator begin() const {
			ItState *itstate = 0;
		
			T *f = first(&itstate);
		
			return Iterator(this,itstate,f);
		}
  
		//! @returns iterator pointing at entry after last element of list
		//! end() can be compared to an iterator by !=
		Iterator end() const {
			return Iterator(this,0,0);
		}

	protected:
		// ------------------------------------------------------------
		// abstract methods that define the iteration method
		//

		//! returns pointer to first element of the virtual list
		//! and a pointer to a new, initialized itstate
		//!
		//!
		//! first() must delete (or overwrite) itstate if itstate!=0
		//!  @param itstate is an in-out parameter
		virtual T* first(ItState **itstate) const = 0;
  
		// --------------------
		//! Changes elem to be next element and returns the elem pointer.
		//! Position specific information should
		//! be drawn from itstate and/or elem.
		//! next may overwrite *elem (or delete old elem and create new one)
		//! but NEVER delete elem completely. This is done internally.
		//! @return NULL if no next element is available or the changed elem
		//!         pointer
		virtual T* next(ItState *itstate, T* const elem) const = 0;
    
	};

	template<class T>
	VirtualList<T>::ItState::~ItState() {}

	template<class T>
	VirtualList<T>::Iterator::~Iterator() {
		if (elem) delete elem;
		if (itstate) delete itstate;
	}

} // namespace biu

#endif // VIRTUAL_LIST_H
