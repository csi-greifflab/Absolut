// $Id: Graph_UD.hh,v 1.2 2016/08/08 12:42:01 mmann Exp $
#ifndef BIU_GRAPH_UD_HH_
#define BIU_GRAPH_UD_HH_

#include <vector>
#include <set>
#include <string>
#include <ostream>

namespace biu
{


/*!
 * An undirected graph implementation. Nodes are labeled using integers starting
 * from 0. The numbering has to be succesively.
 * @author Martin Mann
 */
class Graph_UD
{
	
protected:

	  //! the internal adjacency list representation 
	std::vector< std::set<size_t> > adjList;  
	

public:
	
	  // component label representation
	typedef std::vector<size_t> CompLabel;
	
	  //! iterator on the adjacency list of a node 
	typedef std::set<size_t>::const_iterator AdjacencyIterator;
	  //! pair of begin/end adjacency iterator
	typedef std::pair<AdjacencyIterator,AdjacencyIterator> AdjItPair;

	  //! constructs an undirected graph and initializes the node list with
	  //! nodeNumber nodes.
	Graph_UD(size_t nodeNumber = 0);
	
	  //! destruction
	virtual ~Graph_UD();
	
	  //! updates the node list to indices 0..(n-1) and removes all edges 
	  //! previously connected to nodes with index >= n
	  //! @param n the number of nodes after the changes
	void setNodeNumber(const size_t n);
	
	  //! Adds an edge between node n1 and n2.
	void addEdge(const size_t n1, const size_t n2);
	
	  //! Removes an edge from the graph between node n1 and n2.
	void remEdge(const size_t n1, const size_t n2);
	
	  //! Tests whether or not an edge in the graph between node n1 and n2.
	  //! @return true if the edge exists, false otherwise
	bool isEdge(const size_t n1, const size_t n2);
	
	  //! Returns the number of nodes in the graph.
	  //! @return number of nodes
	size_t nodes() const;
	
	  //! Returns the number of edges in the graph.
	  //! @return number of edges
	size_t edges() const;

	  //! Returns the node degree in the graph.
	  //! @param n node index
	  //! @return node degree of the node
	size_t degree(const size_t n) const;

	  //! calculates the number of connected components and writes the nodeID
	  //! componentID labeling into the given label vector.
	  //! @return the number of independent connected components
	size_t connectedComponents(CompLabel & compID) const;
	
	  //! prints the graph in DOT-format to the given output stream.
	void printDOT( std::ostream & out ) const;
	
	  //! Gives the begin/end iterator of the adjacency list of the specified
	  //! node.
	  //! @param n node index
	  //! @return the begin and end iterator of the adjacency list of node
	  //!         indices
	AdjItPair
	adjNodes( const size_t n ) const;

protected:

	  //! Labels all uncolored direct or indirect to curNode adjacent nodes and 
	  //! therefore the corresponding conntected component. The compID of 
	  //! uncolored nodes is UINT_MAX.
	  //! The method is internally called by method 'connectedComponents(..)'. 
	  //! @param curNode the initial node to explore the connected component
	  //! @param compID the labels set so far
	  //! @param label the label of this connected component to set
	void labelAdjacentNodes( const size_t curNode, CompLabel & compID, const size_t label) const;

};

}

#endif /*GRAPH_UD_HH_*/
