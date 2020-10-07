#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>			// << setprecision() <<
#include <cmath>			// min(x,y)
#include <queue>			// pour dijkstra
#include <fstream>			// i/o fichiers par << et >>
#include <cstdlib>			// for linux, atoi ...
using namespace std;

typedef long double valeur;

/// A - First, definition of internal classes (to store a graph in memory).
///    		a graph is a set of nodes (also called vertices) and edges (weighted arrows).
///            In memory,
///    		1:nodes - a node is an integer : let say a graph is of size S, we will call the nodes 1, 2, 3 ..., S.
///    		2:edges - first, an edge is defined as:
///                a destination node,
//////                a weight.
///                (the starting node is omitted because each starting node will be given its edges from it).

struct edge {
    int dest;
    valeur weight;
    edge(int d, valeur w): dest(d), weight(w) {}            /// Constructor (int destination_node, valeur weight)
};


///	3:graph - the graph is stored in a 'adjacence list', meaning : for each node, we associate the list of edges starting from it.
typedef vector<edge> list;
struct adjacency_list {
    adjacency_list();
    adjacency_list(int _taille);       /// constructor of a new graph with _taille nodes

    vector<list *> table;              /// adjacency list (vector of list of edges)
    int nNodes;                         /// number of nodes in the graph
    int size();                         /// get the number of nodes

    void add(int i, int j, valeur w);   /// to add an edge, from node i to node j with weight w
    list & operator [](int i);         /// to access the list of edges starting from node i

    void resize(int _taille);           /// rebuild the graph with an other number of nodes
    ~adjacency_list();
};


///	4:paths - In order to manipulate optimal pathes in the microcanonical problem, a path (chemin) is defined as a table of the nodes it goes through

struct path {
    vector<int> content;                     /// a path = a list of nodes ID (int) in a graph

    void copy(const path &ch2);       /// to erase and copy from an other path
    int operator [] (int i);
    path(vector<int> ch2);
    path();
    int size();
    void push_back(int z);              /// add a new node to the end of the path
    void show();
};

/// 5:set of pathes - In order to manipulate the set of all optimal paths and to display it,

/// ths structure stores all the pathes from the source (or dhe sink), to the node 'destination'.
struct ensemble_paths {

    int destination;        /// common destination of all the paths
    vector<path> table;   /// table of all paths in the set
    int size;

    void add(path &ch);    /// add a new path in the set
    path operator[](int i);
    ensemble_paths();
    void set_destination(int _destination);    /// set the destination destination to _destination

    void clear();
    void add(struct ensemble_paths &ens);         /// add all the paths from an other set to me
    void addElementPath();                      ///
    void show();
};


///	6:grid-graphs -	Since we will manipulate graph that look like 2D grids, here is a structure to convert
/// from (i,j) positions into a grid to the corresponding ID (number) of the node.

struct grille {
    int width;                /// with of the 2D grid
    grille(int width);        /// constructor for a new grid of width largeur (nb of columns).
    int ID(int i, int j);       /// get the node ID (number) sitting at position (i, j) on the grid (i=row, j=column)
    int pos_x(int ID);          /// get the x position (column) on the grid of node ID
    int pos_y(int ID);          /// get the y position (row) on the grid of node ID

    void showpos(int ID);       /// displays  cout << "$ID ($x,$y) >/t"
    int source;                 /// ID of the special node that will hold the source (ex : ID(N+1, 1))
    int arrivee;                /// ID of the special node that will hold the sink  (ex : ID(N+1, 2))

    void convert_to_positions_x(vector<int> &IDs, vector<int> &pos_to_fill);            /// convert a vector of IDs to corresponding x positions (filling an already defined vector)

    /// trick (or cheating ...), to say 'I want to display the node id, but with a static function'.
    static void showpostatic(int id);
};


/// 7:Big General Class of graphs that include topolocigal sorting and dijkstra from the source or from the node.

class Graph {

private:
    int size;									/// Number of edges
public:
    adjacency_list t;							/// the graph
    adjacency_list reverse;					/// the reverted graph (all arrows opposite direction)
public:
    Graph(int size);							/// to clear t and reverse
    void add_to_graph(int i, int j, valeur w);	/// to add i--(w)-->j to t and j--(w)-->i to the reverse
    int nb_sommets();
private:
    struct vert{ int num; valeur dist; vert(int n, valeur d): num(n), dist(d) { } };
    struct comp { bool operator() (vert &v1, vert &v2){ return (v1.dist > v2.dist); } };
    void DFS(int start);
    void DFS_reverse(int start);



private:
    vector<int> Dpred;							/// PrŽcŽcesseur de chaque sommet (aprs dijkstra)
    vector<valeur> Ddist;						/// Distance de chaque sommet ˆ la source (aprs dijkstra)
    vector<bool> Dvisited;						/// Table des sommets visitŽs par dijkstra au cours de son exŽcution
    bool dijkstra_done;
    int dijkstra_source;
    int dijkstra_dest;
public:
    bool dijkstra(int start, int dest);			/// Calcule dist et pred ˆ partir du sommet start. Ds que toutes les manires optimales de relier destination ont Ã©tÃ© trouvÃ©es, l'algorithme s'arrÃªte.
    int pred(int sommet);
    valeur dist(int sommet);



private:
    vector<int> Dpred_reverse;					/// PrŽcŽcesseur de chaque sommet (aprs dijkstra)
    vector<valeur> Ddist_reverse;				/// Distance de chaque sommet ˆ la source (aprs dijkstra)
    vector<bool> Dvisited_reverse;				/// Table des sommets visitŽs par dijkstra au cours de son exŽcution
    bool dijkstra_reverse_done;
    int dijkstra_reverse_source;
    int dijkstra_reverse_dest;
public:
    bool dijkstra_reverse(int start, int dest);	/// Calcule dist et pred ˆ partir du sommet start. Ds que toutes les manires optimales de relier destination ont ŽtŽ trouvŽes, l'algorithme s'arrte.
    int pred_reverse(int sommet);
    valeur dist_reverse(int sommet);



private:
    vector<int> Tpred_topo;
    int current_index_topo;
    vector<int> Tordre;
    bool tri_topologique_done;
    int tri_topologique_source;
public:
    void tri_topologique(int start);
    int pred_topo(int);
    int nb_sommets_ordonnes();
    int ordre(int nb);
    int number(int sommet);



private:
    vector<int> Tpred_topo_reverse;
    int current_index_topo_reverse;
    vector<int> Tordre_reverse;
    bool tri_topologique_reverse_done;
    int tri_topologique_reverse_source;
public:
    void tri_topologique_reverse(int start);
    int pred_topo_reverse(int);
    int nb_sommets_ordonnes_reverse();
    int ordre_reverse(int nb);
    int number_reverse(int sommet);

    void show_graph(void (* conv)(int ID) = NULL);


private:
    int cpt_interne;

private:
    typedef vector<int> solution;
    void trace_all_paths(int ending_node, int source);
    vector<ensemble_paths> big_table;
    vector<solution> ensemble_solutions;
    bool all_paths_computed;


    // careful, dijkstra and topological have to be done from the same source and destination
    vector<double> density_graph;
    void compute_density_from_nb_paths(int source, int destination);
    vector<int> trace_one_optimal_path(int dest);
    void show_one_best_solution(int dest);
public:

};

#endif
