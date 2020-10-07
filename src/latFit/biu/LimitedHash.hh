// $Id: LimitedHash.hh,v 1.2 2016/08/08 12:42:00 mmann Exp $
#ifndef BIU_LIMITEDHASH_HH_
#define BIU_LIMITEDHASH_HH_

#include <utility>


// get information on available hashes or maps
#include "biu/HashMap.hh"

// include best available hash or map
#if HAVE_UNORDERED_MAP == 1
	#include <unordered_map>
#else
#if HAVE_TR1_UNORDERED_MAP == 1
	#include <tr1/unordered_map>
#else
#if HAVE_GNU_HASH_MAP == 1
	#include <ext/hash_map>
#else
	#include <map>
namespace biu {
	//! dummy hashing object for default use when a std::map has to be used
	template<class T> 
	struct hash_dummy
	{
		//! @param s the value to hash
		//! @return 0 all the time
		unsigned int operator()(const T& s) const
		{ return 0; }
	};
} // namespace biu
#endif
#endif
#endif



namespace biu
{


	/*!
	 * An unordered hash implementation with limited capacity. It handles up to
	 * a given maximum number (maxSize) of elements and stores at least the last
	 * maxSize/2 elements inserted. It utilizes two unordered_map containers
	 * which are filled up to a maximum size of maxSize/2. If one container is
	 * full the other is filled. If this is full as well, the first container
	 * is cleared and filled afterwards. That way, the last maxSize/2 inserted 
	 * elements are maintained for sure. 
	 * 
	 * @param Key the key type to use that is hashed
	 * @param T the mapped type to store in the hash
	 * @param Hash the hash function object to apply for hashing
	 * @param Pred the comparison function object to apply to check for equality
	 * @param Alloc the memory allocator to use to allocate an object of 
	 *              value_type
	 * 
	 * @author Martin Mann http://www.bioinf.uni-freiburg.de/~mmann/
	 * 
	 */
	template <	class Key,
			    class T,
#if HAVE_UNORDERED_MAP == 1 
			    class Hash = std::hash<Key>,
			    class Pred = std::equal_to<Key>,
#else
#if HAVE_TR1_UNORDERED_MAP == 1
			    class Hash = std::tr1::hash<Key>,
			    class Pred = std::equal_to<Key>,
#else
#if HAVE_GNU_HASH_MAP == 1
			    class Hash  = __gnu_cxx::hash<Key>,
			    class Pred = std::equal_to<Key>,
#else
			    class Hash = biu::hash_dummy<Key>,
			    class Pred = std::less<Key>,
#endif
#endif
#endif
			    class Alloc = std::allocator< typename std::pair<const Key, T> > >
	class LimitedHash
	{

	public:
		
		
		// set typedef for best available hash or map
		#if HAVE_UNORDERED_MAP == 1
			//! the type of hash container internally used
			typedef std::unordered_map< Key, T, Hash, Pred, Alloc > hash_type;
		#else
		#if HAVE_TR1_UNORDERED_MAP == 1
			//! the type of hash container internally used
			typedef std::tr1::unordered_map< Key, T, Hash, Pred, Alloc > hash_type;
		#else
		#if HAVE_GNU_HASH_MAP == 1
			//! the type of hash container internally used
			typedef __gnu_cxx::hash_map< Key, T, Hash, Pred, Alloc > hash_type;
		#else
			//! the type of hash container internally used
			typedef std::map< Key, T, Pred, Alloc > hash_type;
		#endif
		#endif
		#endif
		
		//! iterator on the internal hash_type
		typedef typename hash_type::iterator internal_iterator;
		//! const_iterator on the internal hash_type
		typedef typename hash_type::iterator const_internal_iterator;
		
		/*!
		 * Generic iterator class that is able to iterator on all elements 
		 * stored in both internal hashes to represent them as one container.
		 * 
		 * @param iterator_type the internal iterator type masked by this object
		 * 
		 * @author Martin Mann
		 */
	    template< class iterator_type >
	    class iteratorT {
	    	  //! to enable access to non-public constructors
	    	friend class LimitedHash < Key, T, Hash, Pred, Alloc >;
	    protected:
	    	  //! the internal iterator that gives the current element
	    	iterator_type cur;
	    	  //! the internal iterator that marks the end of the first internal
	    	  //! hash container
	    	iterator_type firstEnd;
	    	  //! the internal iterator that marks the beginning of the second
	    	  //! internal hash container
	    	iterator_type secondStart;
	    	  //! a flag the is TRUE if cur is still in the first container or
	    	  //! FALSE if it is already iterating on the second container
	    	bool firstNotEnded;
	    	
	    	  /*! Construction
	    	   * 
	    	   * @param cur_ the beginning of the range to iterate
	    	   * @param firstEnd_ the end of the first range to iterate
	    	   * @param secondStart_ the begin of the second range to iterate
	    	   * @param firstNotEnded_ if TRUE we assume cur will reach firstEnd
	    	   *        such that we will have to switch to secondStart, 
	    	   *        if FALSE we assume cur is equal to secondStart or higher
	    	   */
	    	explicit
	    	iteratorT(	iterator_type cur_
	    				, iterator_type firstEnd_
	    				, iterator_type secondStart_
	    				, bool firstNotEnded_ = true );
	    	
	    public:
	    	
	    	//! Default constructor
	    	iteratorT();
	    	//! Destructor
	    	~iteratorT();
	    	
	    	//! Copy construction
	    	iteratorT(const iteratorT< iterator_type >& toCopy);
	    	
	        //! iterator dereferencing
	    	//! @return a reference of a value_type object this iterator points to
	    	typename iterator_type::reference
	        operator*() const;

	        //! iterator pointer dereferencing
	    	//! @return a pointer to the value_type object this iterator points to
	        typename iterator_type::pointer
	        operator->() const;

	        //! pre-incrementation (first increment, than return changed state)
	        //! e.g. (++it)
	        //! @return a reference to (*this) AFTER incrementation
	        iteratorT<iterator_type>&
	        operator++();

	        //! post-incrementation (first return current state, than increment)
	        //! e.g. (it++)
	        //! @return a copy of (*this) BEFORE incrementation
	        iteratorT<iterator_type>
	        operator++(int);

	        //! equality comparision 
	        //! @param toCompare the iterator to compare this to
	        //! @return true if all internal iterators are equal, false otherwise
	        bool
	        operator==( const iteratorT< iterator_type >& toCompare );
	    	
	        //! inequality comparision 
	        //! @param toCompare the iterator to compare this to
	        //! @return returns (!(*this == toCompare))
	        bool
	        operator!=( const iteratorT< iterator_type >& toCompare );
	    	
	    	
	    };

	    //! type of key values
	    typedef typename hash_type::key_type		key_type;
	    //! type of the stored objects : pair<key_type,mapped_type>
	    typedef typename hash_type::value_type		value_type;
	    //! type of mapped values
	    typedef typename hash_type::mapped_type		mapped_type;
	    //! type of the hashing function object for key values
	    typedef typename hash_type::hasher			hasher;
	    //! type of the comparison function object for equality checks
	    typedef typename hash_type::key_equal		key_equal;
	    //! type of the memory allocation function object for value_types
	    typedef typename hash_type::allocator_type	allocator_type;
	    
	    //! type of a value_type object pointer
	    typedef typename allocator_type::pointer         	pointer;             
	    //! type of a value_type object constant pointer
	    typedef typename allocator_type::const_pointer   	const_pointer;   
	    //! type of a value_type object reference
	    typedef typename allocator_type::reference			reference;
	    //! type of a value_type object constant reference
	    typedef typename allocator_type::const_reference	const_reference;
	    //! type of size describing values
	    typedef typename hash_type::size_type				size_type;
	    
	    //! type of forward iterator to traverse the container content
	    typedef iteratorT< internal_iterator >			iterator;
	    //! type of constant forward iterator to traverse the container content
	    typedef iteratorT< const_internal_iterator >	const_iterator;

	protected:
		
		//! the maximal number of elements this hash should contain
		size_type maxSize;
		//! the maximal number of elements each of the internal hashs should
		//! contain
		size_type maxSizeHalf;
		
		//! internal hash 1
		hash_type hashObject1;
		//! internal hash 2
		hash_type hashObject2;
		
		//! access to the internal hash currently to fill
		hash_type * hash2Fill;
		//! access to the internal hash that is full and used for lookups only
		hash_type * hash2Check;
		
	public:
		
	    // construct/destroy/copy
		/*! Construction of a new hash with limited size.
		 * @param maxSize the maximal number of elements allowed
		 * @param n the initial size
		 * @param hf the hash function object to use
		 * @param eql the comparison object to use
		 * @param a the allocator to use
		 */
	    explicit LimitedHash(	const size_type maxSize, 
	    						const size_type n = 3,
	                           const hasher& hf = hasher(),
	                           const key_equal& eql = key_equal(),
	                           const allocator_type& a = allocator_type());

	    /*! Copy construction
	     * @param toCopy the LimitedHash to copy
	     */
	    LimitedHash(const LimitedHash& toCopy);
	    
	    //! Destruction
	    ~LimitedHash();
	    
	    /*! Assignment operator
	     * @param toCopy the Limitedhash to copy
	     * @return access to this object which is now a copy of the given hash
	     */
	    LimitedHash& 
	    operator=(const LimitedHash& toCopy);
	    
	    
	    
	    /*! Access to the allocator used.
	     * @return the allocator in use
	     */
	    Alloc 
	    get_allocator() const;

	    /*! Whether or not an element is present or not.
	     * @return true if no element is stored in the hash, false otherwise.
	     */
	    bool 
	    empty() const;
	    
	    /*! Access to the current number of elements stored.
	     * @return the number of elements in the hash.
	     */
	    size_type 
	    size() const;
	    
	    /*! Access to the maximal number of elements allowed in the hash.
	     * @return the maximal number elements allowed.
	     */
	    size_type 
	    max_size() const;

	    /*! Sets and returns the maximal number of elements allowed in the hash.
	     * @param newMax the new maximal number of elements allowed
	     * @return the maximal number elements allowed.
	     */
	    size_type 
	    max_size( const size_type newMax );

	    // iterators
	    
	    //! Iterator that points to the begin of container
	    //! @return iterator to the containers begin
	    iterator       
	    begin();
	    
	    //! Constant iterator that points to the begin of container
	    //! @return constant iterator to the containers begin
	    const_iterator 
	    begin() const;
	    
	    //! Iterator that points to the end of container, i.e. AFTER
	    //! the last element accessible
	    //! @return iterator to the containers end
	    iterator       
	    end();
	    
	    //! Constant iterator that points to the end of container, i.e. AFTER
	    //! the last element accessible
	    //! @return constant iterator to the containers end
	    const_iterator 
	    end() const;

	    // modifiers
	    
	    /*! Inserts the given value object into the hash. If its key is already
	     * present, only the mapped value is updated. 
	     * 
	     * @param obj the value object to insert (key,val)
	     * @return first is the iterator to the inserted object, second is TRUE
	     * if the object was inserted, or FALSE if it was updated only
	     */
	    std::pair<iterator, bool> 
	    insert(const value_type& obj);
	    
//	    iterator insert(const_iterator hint, const value_type& obj);
	    
	    /*! Inserts successively a list of value objects into the hash utilizing
	     * the insert(value_type) method.
	     * 
	     * @param first the begin of the list of values to insert
	     * @param last the end of the list of values to insert
	     */
	    template <class InputIterator>
	    void 
	    insert(InputIterator first, InputIterator last);

	    //! Erases the value object the iterator points to
	    //! @param position the value object to remove
	    void 
	    erase(const_iterator position);
	    
	    //! Erases all value objects that equal the key.
	    //! @param k the key to be equal
	    //! @return the number of elements removed
	    size_type 
	    erase(const key_type& k);
	    
	    /*! Erases successively a list of value objects from the hash utilizing
	     * the erase(const_iterator) method.
	     * 
	     * @param first the begin of the list of values to erase
	     * @param last the end of the list of values to erase
	     */
	    void 
	    erase(const_iterator first, const_iterator last);
	    
	    //! Removes all elements from the hash.
	    void 
	    clear();

	    //! Swaps the content with the given hash.
	    //! @param toSwap the hash to swap content with
	    void 
	    swap(LimitedHash& toSwap);

	    // observers
	    
	    //! Access to the hashing function object used to hash the key values
	    //! @return the hashing object
	    hasher 
	    hash_function() const;
	    
	    
	    //! Access to the comparison function object used to check the key 
	    //! values upon equality
	    //! @return the equality comparison object
	    key_equal 
	    key_eq() const;

	    // lookup
	    
	    /*! Access to the value object that matches the given key. If none is
	     * present this->end() is returned. 
	     * 
	     * @param k the key to find
	     * @return an iterator pointing to the element in the hash matching the
	     *         key or this->end() if none exists
	     */
	    iterator       
	    find(const key_type& k);
	    
	    /*! Constant access to the value object that matches the given key. 
	     * If none is present this->end() is returned. 
	     * 
	     * @param k the key to find
	     * @return a constant iterator pointing to the element in the hash 
	     *         matching the key or this->end() if none exists
	     */
	    const_iterator 
	    find(const key_type& k) const;
	    
	    /*! The number of value objects that match the given key. 
	     * 
	     * @param k the key to find
	     * @return number of element matching the key value
	     */
	    size_type 
	    count(const key_type& k) const;
	    
	    /*! Access to all elements that match the given key value. 
	     * 
	     * @param k the key to find
	     * @return an iterator range (first,last) of elements that equal the 
	     * given key, or (end(),end()) if none exists
	     */
	    std::pair<iterator, iterator>
	    equal_range(const key_type& k);
	    
	    /*! Constant access to all elements that match the given key value. 
	     * 
	     * @param k the key to find
	     * @return a constant iterator range (first,last) of elements that  
	     * equal the given key, or (end(),end()) if none exists
	     */
	    std::pair<const_iterator, const_iterator>
	    equal_range(const key_type& k) const;

	    /*! The member function determines the iterator it as the return value 
	     * of insert( value_type(k, T()) ). (It inserts an element with the
	     * specified key if no such element exists.) It then returns a 
	     * reference to (*it).second.
	     */
	    mapped_type& 
	    operator[](const key_type& k);

	    /*! The member function rehashes the unordered map, ensuring that it 
	     * contains at least n buckets.
	     * @param n the minimum number of buckets to use
	     */
	    void 
	    rehash(size_type n);
		
	};

} // namespace biu

// include template implementations
#include "biu/LimitedHash.icc"

#endif /*BIU_LIMITEDHASH_HH_*/
