// $Id: Graph_UD.cc,v 1.2 2016/08/08 12:41:57 mmann Exp $

#include "biu/Graph_UD.hh"

#include <stack>
#include "biu/assertbiu.hh"
#include <limits.h>


namespace biu
{

	Graph_UD::Graph_UD( size_t nodeNumber) 
	  :	adjList(nodeNumber)
	{
	}
	
	Graph_UD::~Graph_UD()
	{
	}
	
	void
	Graph_UD::addEdge(const size_t n1, const size_t n2) {
		assertbiu(std::max(n1,n2) < adjList.size(), "node index out of range");
		adjList[n1].insert(n2);
		adjList[n2].insert(n1);
	}

	void 
	Graph_UD::remEdge(const size_t n1, const size_t n2) {
		assertbiu(std::max(n1,n2) < adjList.size(), "node index out of range");
		adjList[n1].erase(n2);
		adjList[n2].erase(n1);
	}
	
	bool
	Graph_UD::isEdge(const size_t n1, const size_t n2) {
		assertbiu(std::max(n1,n2) < adjList.size(), "node index out of range");
		return adjList[n1].find(n2) != adjList[n1].end();
	}
	
	size_t
	Graph_UD::nodes(void) const {
		return adjList.size();
	}

	size_t
	Graph_UD::edges(void) const {
		size_t num = 0;
		for(size_t i = 0; i< adjList.size(); i++) {
			num += adjList[i].size();
		}
		return num/2; // each edges was counted twice
	}

	size_t 
	Graph_UD::degree(const size_t n) const {
		return adjList[n].size();
	}

	size_t 
	Graph_UD::connectedComponents( Graph_UD::CompLabel & compID ) const {
		  // initializing the labels
		compID.resize(adjList.size()); // resizing to neccessary size
		  // check if anything to do
		if (adjList.size() == 0)
			return 0;
		  // set initial labels
		for (Graph_UD::CompLabel::iterator it1 = compID.begin(); it1!=compID.end(); it1++) {
			*it1 = UINT_MAX;
		}
		  // temp data
		size_t curLabel = 0;
		size_t lastInitID = 0;
		
		// recursion for complete graph exploration
		do {
			  // color the connected component adjacent to node lastInitID
			labelAdjacentNodes(lastInitID, compID, curLabel);
			  // shift lastInitID to next unseen node and therefore next component
			for (lastInitID++; lastInitID < compID.size() && compID[lastInitID] != UINT_MAX; lastInitID++);
			  // increase connected component label
			curLabel++;
		} while (lastInitID < adjList.size());
		
		return curLabel;
	}
	
	// NON-RECURSIVE VERSION
	void 
	Graph_UD::labelAdjacentNodes( const size_t curNode, Graph_UD::CompLabel & compID, const size_t label) const {
		assertbiu( curNode < compID.size() , "initial node out of node index range");
		if (compID[curNode] != UINT_MAX)	// nothing to do
			return;
		  // temporary stack and initialization
		typedef std::set<size_t>::const_iterator adjIT;
		std::stack<adjIT> path;
		compID[curNode] = label;
		
		if (adjList[curNode].empty()) // abortion if no adjacent nodes
			return;
		else
			path.push( adjList[curNode].begin() );	// add first adjacent node
		
		while ( !path.empty() ) {
			if ( compID[*path.top()] == UINT_MAX) {
				compID[*path.top()] = label;	// set label
				if ( ! adjList[*path.top()].empty() )	// add adjacency iterator
					path.push( adjList[*path.top()].begin() );
			} else {
				adjIT t = path.top(); // get head element
				path.pop();
				t++; // shift to next in list
				// if a valid neighbor add to stack again
				if ( (path.empty() && t != adjList[curNode].end()) || 
					(!path.empty() && t != adjList[*path.top()].end()) ) { 
					path.push(t);
				}
			}
		}
	}
/*  // RECURSIVE VERSION
	void 
	Graph_UD::labelAdjacentNodes( const size_t curNode, Graph_UD::CompLabel & compID, const size_t label) const {
		assertbiu( curNode < compID.size() , "initial node out of node index range");
		  // recursion abortion
		if (compID[curNode] != UINT_MAX)
			return;
		  // label current node
		compID[curNode] = label;
		  // recursive labeling of all adjacent nodes if not already labeled.
		for (std::set<size_t>::const_iterator it = adjList[curNode].begin();
			it != adjList[curNode].end(); it++) 
		{
			labelAdjacentNodes(*it, compID, label);
		}
	}
*/
	void 
	Graph_UD::printDOT( std::ostream & out ) const {
		size_t i = 0;
		  // initial header
		out <<"graph G {\n";
		  // print node indices
		for (i=0; i<nodes(); i++)
			out <<" "<<i <<";";
		out <<"\n";
		  // print edges
		for (i=0; i<nodes(); i++)
			for (std::set<size_t>::const_iterator it = adjList[i].begin(); it != adjList[i].end(); it++) 
				if (i <= (*it))
					out <<" " <<i <<"--" <<(*it) <<";\n";
		  // closing
		out <<"}";
	}
	
	void
	Graph_UD::setNodeNumber(const size_t max) {
		  // resize adjacency list
		adjList.resize(max);
		  // remove all nodes with index >= max
		for (size_t i=0; i< adjList.size(); i++) {
			adjList[i].erase(adjList[i].find(max), adjList[i].end());
		}
	}
	
	Graph_UD::AdjItPair
	Graph_UD::adjNodes( const size_t n ) const
	{
		return AdjItPair(adjList[n].begin(), adjList[n].end());
	}


}
