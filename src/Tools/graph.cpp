// Created by Robert Philippe on 18/03/11.

#include "graph.h"

#define my_precision 0.000001
bool egal(valeur a, valeur b){
    return (abs(a-b) < my_precision);
}

#define DBG false

/* --- class edge is already defined in the .h file --- */



/* --- class liste d'adjacence --- */

// liste := vector<edge>
int adjacency_list::size(){
    return nNodes;
}
adjacency_list::adjacency_list() {nNodes = 0; table.clear();}
adjacency_list::adjacency_list(int _taille) {resize(_taille);}
void adjacency_list::resize(int _taille){
    nNodes = _taille;
    table.resize(nNodes);
    for(int i = 0; i < nNodes; ++i){
        table[i] = new vector<edge>();
    }
}
list & adjacency_list::operator [](int i){		//passe la référence
    if((i < 0) || (i >= nNodes)) cerr << "Unauthorized call to adjacency table (index " << i << ")\n";
    return *(table[i]);
}
void adjacency_list::add(int i, int j, valeur w){
    table[i]->push_back(edge(j, w));
}
adjacency_list::~adjacency_list(){
    for(int i = 0; i < nNodes; ++i){
        table[i]->clear();
        delete table[i];
    }
}



/* ---- class Graph ---- */

Graph::Graph(int _size){
    size = _size;

    t.resize(size);
    reverse.resize(size);

    Dpred.resize(size);
    Ddist.resize(size);
    Dvisited.resize(size);
    dijkstra_done = false;
    dijkstra_source = -1;
    dijkstra_dest = -1;

    Dpred_reverse.resize(size);
    Ddist_reverse.resize(size);
    Dvisited_reverse.resize(size);
    dijkstra_reverse_done = false;
    dijkstra_reverse_source = -1;
    dijkstra_reverse_dest = -1;

    Tpred_topo.resize(size);
    current_index_topo = 0;
    Tordre.resize(size);					// only the beginning used (0..current_index_topo - 1)
    tri_topologique_done = false;
    tri_topologique_source = -1;

    Tpred_topo_reverse.resize(size);
    current_index_topo_reverse = 0;
    Tordre_reverse.resize(size);			// only the beginning used (0..current_index_topo_reverse - 1)
    tri_topologique_reverse_done = false;
    tri_topologique_reverse_source = -1;

    cpt_interne = 0;
}

void Graph::add_to_graph(int i, int j, valeur w){
    t[i].push_back(edge(j, w));
    reverse[j].push_back(edge(i,w));
    dijkstra_done = false;					// we could also reset the sources/destinations, but not necessary
    dijkstra_reverse_done = false;
    tri_topologique_done = false;
    tri_topologique_reverse_done = false;
}

int Graph::nb_sommets(){
    return size;
}


void Graph::tri_topologique(int start){
    if(tri_topologique_done && (start == tri_topologique_source)){
        cerr << "Tri Topologique already done from the same source. Not doing it again.\n";
        return;
    }
    if((start < 0) || (start >= size)) {
        cerr << "tri_topologique : out of range for source vertex number. Abort." << endl;
        return;
    }
    tri_topologique_done = false;

    for(int i = 0; i < size; ++i){
        Tordre[i] = 0;
        Tpred_topo[i] = -1;
    }
    Tpred_topo[start] = start;
    current_index_topo = 0;

    cpt_interne = 0;
    DFS(start);
        //cerr << "    -> DFS finished ; " << current_index_topo << " vertices seen\n";
    //if(current_index_topo < size) cerr << "WARNING : tri_topologique didn't manage to join all vertices of the graph. It may be not connex.\n";

    vector<int> cp = Tordre;
    for(int i = 0; i < current_index_topo; ++i){
        Tordre[i] = cp[current_index_topo - i - 1];
    }
    tri_topologique_done = true;
    tri_topologique_source = start;
}

void Graph::DFS(int start){
    if(size > 10000) {
        if((cpt_interne % 1000) == 0) cerr << "... " << cpt_interne << " / " << size << " sommets visites\n";
    }
    list::iterator it;
    for(it = t[start].begin(); it != t[start].end(); ++it){
        int d = it->dest;
        if(Tpred_topo[d] < 0){
            Tpred_topo[d] = start;
            cpt_interne++;
            DFS(d);
        }
    }
    Tordre[current_index_topo] = start;
    current_index_topo++;
}

void Graph::tri_topologique_reverse(int start){
    if(tri_topologique_reverse_done && (start == tri_topologique_reverse_source)){
        cerr << "Tri Topologique (reverse) already done from the same source. Not doing it again.\n";
        return;
    }
    if((start < 0) || (start >= size)) {
        cerr << "tri_topologique_reverse : out of range for source vertex number. Abort." << endl;
        return;
    }
    tri_topologique_reverse_done = false;
    for(int i = 0; i < size; ++i){
        Tordre_reverse[i] = 0;
        Tpred_topo_reverse[i] = -1;
    }
    Tpred_topo_reverse[start] = start;
    current_index_topo_reverse = 0;

    cpt_interne = 0;
    DFS_reverse(start);
        ///cerr << "    -> DFS_reverse finished ; " << current_index_topo_reverse << " vertices seen\n";
    //if(current_index_topo_reverse < size) cerr << "WARNING : tri_topologique_reverse didn't manage to join all vertices of the graph. It may be not connex.\n";

    vector<int> cp = Tordre_reverse;
    for(int i = 0; i < current_index_topo_reverse; ++i){
        Tordre_reverse[i] = cp[current_index_topo_reverse - i - 1];
    }
    tri_topologique_reverse_done = true;
    tri_topologique_reverse_source = start;
}

void Graph::DFS_reverse(int start){
    if(size > 10000) {
        if((cpt_interne % 1000) == 0) cerr << "... " << cpt_interne << " / " << size << " sommets visites\n";
    }
    list::iterator it;
    for(it = reverse[start].begin(); it != reverse[start].end(); ++it){
        int d = it->dest;
        if(Tpred_topo_reverse[d] < 0){
            Tpred_topo_reverse[d] = start;
            cpt_interne++;
            DFS_reverse(d);
        }
    }
    Tordre_reverse[current_index_topo_reverse] = start;
    current_index_topo_reverse++;
}





bool Graph::dijkstra(int start, int destination){
    if(dijkstra_done && (start == dijkstra_source) && (destination == dijkstra_dest)){
        cerr << "dijkstra already done from the same source and the same destination. Not doing it again." << endl;
    }
    if((start < 0) || (start >= size) || (destination < 0) || (destination >= size) ) {
        cerr << "dijkstra : out of range for source or destination vertex number. Abort." << endl;
        return false;
    }
    dijkstra_done = false;

    vector<bool> added(size, false);	// To remember who is in the queue without looking in all the queue
    valeur best_dist_found = 1e25;		// to stop when all the vertices of distance dist[destination] has been found (useless to go on then)
    cpt_interne = 0;					// to display a message all 1000 iterations if the graph is big

    int u;
    priority_queue <vert, vector<vert>, comp> q;	// in this queue, vertices will be sorted by distance from the source
    list::iterator it;
    for(int i = 0; i < size; ++i){
        Ddist[i] = -1;
        Dpred[i] = -1;
        Dvisited[i] = false;
    }
    Ddist[start] = 0;
    q.push(vert(start,0));


    while(!q.empty()){

        /*  ---- Takes the first vertex out of the list. This is the closest unseen vertex from the source. --- */
        u = q.top().num;
        q.pop();

        /*  ---- to display something every 1000 iterations ---- */
        cpt_interne++;
        if(size > 10000){
            if((cpt_interne % 1000) == 0) cerr << "... " << cpt_interne << " / " << size << " sommets visites " << q.size() << "\n";
        }
        /*  ---- to stop when all vertices of distance <= dist[destination] has been found  --- */
        if(Ddist[u] > best_dist_found) {
            dijkstra_done = true;
            dijkstra_source = start;
            dijkstra_dest = destination;
            return true;
        }
        if(u == destination) best_dist_found = Ddist[u];


        /* ---- the big loop ---- */
        if(!Dvisited[u]){
            Dvisited[u] = true;
            for(it = t[u].begin() ; it != t[u].end(); ++it){
                int v = it->dest;
                valeur d2 = Ddist[u] + it->weight;		// does the distance of v is better going through source -> u -> v ?
                if((Ddist[v] < 0) || (Ddist[v] > d2)){	// if yes, (two cases : either v is unseen, either it is already seen (in the list), but we can improve its distance
                    Ddist[v] = d2;
                    Dpred[v] = u;
                    if(!Dvisited[v]) {
                        q.push(vert(v, Ddist[v]));
                        added[v] = true;
                    }
                }
                if(!Dvisited[v] && (!added[v])) {		// In the case v is already in the list (added[]), but is not improved, we don't repeat it in the list, in order to keep the list small.
                    q.push(vert(v, Ddist[v]));
                    added[v] = true;
                }
            }
        }
    }

    /* ---- if we arrive here, we have not managed to reach the destination ---- */
    dijkstra_done = true;
    dijkstra_source = start;
    dijkstra_dest = destination;
    return false;
}


bool Graph::dijkstra_reverse(int start, int destination){
    if(dijkstra_reverse_done && (start == dijkstra_reverse_source) && (destination == dijkstra_reverse_dest)){
        cerr << "dijkstra_reverse already done from the same source and the same destination. Not doing it again." << endl;
    }
    if((start < 0) || (start >= size) || (destination < 0) || (destination >= size) ) {
        cerr << "dijkstra_reverse : out of range for source or destination vertex number. Abort." << endl;
        return false;
    }
    dijkstra_reverse_done = false;

    vector<bool> added(size, false);	// To remember who is in the queue without looking in all the queue
    valeur best_dist_found = 1e25;		// to stop when all the vertices of distance dist[destination] has been found (useless to go on then)
    cpt_interne = 0;					// to display a message all 1000 iterations if the graph is big

    int u;
    priority_queue <vert, vector<vert>, comp> q;	// in this queue, vertices will be sorted by distance from the source
    list::iterator it;
    for(int i = 0; i < size; ++i){
        Ddist_reverse[i] = -1;
        Dpred_reverse[i] = -1;
        Dvisited_reverse[i] = false;
    }
    Ddist_reverse[start] = 0;
    q.push(vert(start,0));


    while(!q.empty()){

        /*  ---- Takes the first vertex out of the list. This is the closest unseen vertex from the source. --- */
        u = q.top().num;
        q.pop();

        /*  ---- to display something every 1000 iterations ---- */
        cpt_interne++;
        if(size > 10000){
            if((cpt_interne % 1000) == 0) cerr << "... " << cpt_interne << " / " << size << " sommets visites " << q.size() << "\n";
        }
        /*  ---- to stop when all vertices of distance <= dist[destination] has been found  --- */
        if(Ddist_reverse[u] > best_dist_found) {
            dijkstra_reverse_done = true;
            dijkstra_reverse_source = start;
            dijkstra_reverse_dest = destination;
            return true;
        }
        if(u == destination) best_dist_found = Ddist_reverse[u];


        /* ---- the big loop ---- */
        if(!Dvisited_reverse[u]){
            Dvisited_reverse[u] = true;
            for(it = reverse[u].begin() ; it != reverse[u].end(); ++it){
                int v = it->dest;
                valeur d2 = Ddist_reverse[u] + it->weight;		// does the distance of v is better going through source -> u -> v ?
                if((Ddist_reverse[v] < 0) || (Ddist_reverse[v] > d2)){	// if yes, (two cases : either v is unseen, either it is already seen (in the list), but we can improve its distance
                    Ddist_reverse[v] = d2;
                    Dpred_reverse[v] = u;
                    if(!Dvisited_reverse[v]) {
                        q.push(vert(v, Ddist_reverse[v]));
                        added[v] = true;
                    }
                }
                if(!Dvisited_reverse[v] && (!added[v])) {		// In the case v is already in the list (added[]), but is not improved, we don't repeat it in the list, in order to keep the list small.
                    q.push(vert(v, Ddist_reverse[v]));
                    added[v] = true;
                }
            }
        }
    }

    /* ---- if we arrive here, we have not managed to reach the destination ---- */
    dijkstra_reverse_done = true;
    dijkstra_reverse_source = start;
    dijkstra_reverse_dest = destination;
    return false;
}


int Graph::pred(int sommet){
    if(dijkstra_done){
        if((sommet >= 0) && (sommet < size)){
            return Dpred[sommet];
        } else cerr << "Out of bounds, pred(" << sommet << ")" << endl;
    } else cerr << "Can't access predecessors before having done dijkstra." << endl;
    return -1;
}

valeur Graph::dist(int sommet){
    if(dijkstra_done){
        if((sommet >= 0) && (sommet < size)){
            return Ddist[sommet];
        } else cerr << "Out of bounds, dist(" << sommet << ")" << endl;
    } else cerr << "Can't access distances before having done dijkstra." << endl;
    return -1;
}

int Graph::pred_reverse(int sommet){
    if(dijkstra_reverse_done){
        if((sommet >= 0) && (sommet < size)){
            return Dpred_reverse[sommet];
        } else cerr << "Out of bounds, pred_reverse(" << sommet << ")" << endl;
    } else cerr << "Can't access predecessors (reverse) before having done dijkstra." << endl;
    return -1;
}

valeur Graph::dist_reverse(int sommet){
    if(dijkstra_reverse_done){
        if((sommet >= 0) && (sommet < size)){
            return Ddist_reverse[sommet];
        } else cerr << "Out of bounds, dist_reverse(" << sommet << ")" << endl;
    } else cerr << "Can't access distances (reverse) before having done dijkstra." << endl;
    return -1;
}


int Graph::pred_topo(int sommet){
    if(tri_topologique_done){
        if((sommet >= 0) && (sommet < size)){
            return Tpred_topo[sommet];
        } else cerr << "Out of bounds, pred_topo(" << sommet << ")" << endl;
    } else cerr << "Can't access predecessors (topological) before having done dijkstra." << endl;
    return -1;
}

int Graph::nb_sommets_ordonnes(){
    if(tri_topologique_done){
        return current_index_topo;
    } else cerr << "Topological sorting not done. Can't call nb_sommets_ordonnes yet." << endl;
    return 0;
}

int Graph::ordre(int nb){
    if(tri_topologique_done){
        if((nb >= 0) && (nb < current_index_topo)){
            return Tordre[nb];
        } else cerr << "Can't access the order of vertices whereas tri_topologique has not been done" << endl;
    } else cerr << "Topological sorting not done. Can't call nb_sommets_ordonnes yet." << endl;
    return 0;
}

int Graph::number(int sommet){
    if((sommet < 0) || (sommet > size)) {
        cerr << "number() : Out of range in vertex index" << endl;
        return -1;
    }
    for(int i = 0; i < nb_sommets_ordonnes() ; ++i){
        if(ordre(i) == sommet) return i;
    }
    return -2; //case the vertex has not been seen by topological sorting
}

int Graph::pred_topo_reverse(int sommet){
    if(tri_topologique_reverse_done){
        if((sommet >= 0) && (sommet < size)){
            return Tpred_topo_reverse[sommet];
        } else cerr << "Out of bounds, pred_topo_reverse(" << sommet << ")" << endl;
    } else cerr << "Can't access predecessors (topological, reverse) before having done dijkstra." << endl;
    return -1;
}

int Graph::nb_sommets_ordonnes_reverse(){
    if(tri_topologique_reverse_done){
        return current_index_topo_reverse;
    } else cerr << "Reverse topological sorting not done. Can't call nb_sommets_ordonnes yet." << endl;
    return 0;
}

int Graph::ordre_reverse(int nb){
    if(tri_topologique_reverse_done){
        if((nb >= 0) && (nb < current_index_topo_reverse)){
            return Tordre_reverse[nb];
        } else cerr << "Can't access the order of vertices whereas tri_topologique has not been done" << endl;
    } else cerr << "Topological sorting not done. Can't call nb_sommets_ordonnes yet." << endl;
    return 0;
}

int Graph::number_reverse(int sommet){
    if((sommet < 0) || (sommet > size)) {
        cerr << "number_reverse() : Out of range in vertex index" << endl;
        return -1;
    }
    for(int i = 0; i < nb_sommets_ordonnes_reverse() ; ++i){
        if(ordre_reverse(i) == sommet) return i;
    }
    return -2; //case the vertex has not been seen by topological sorting
}


void Graph::show_graph(void (* conv)(int ID)){
    cout << " ==============  Structure du graphe généré :  ============== \n";
    cout << "   - nb sommets : " << size << endl;
    if(dijkstra_done) cout << "    -> Dijkstra done from " << dijkstra_source << " to " << dijkstra_dest << endl;
    else cout << "    -> Dijkstra not done." << endl;
    if(dijkstra_reverse_done) cout << "    -> Dijkstra (in reverse) done from " << dijkstra_reverse_source << " to " << dijkstra_reverse_dest << endl;
    else cout << "    -> Dijkstra (reverse) not done." << endl;
    if(tri_topologique_done) cout << "    -> Topological sorting done from " << tri_topologique_source << endl;
    else cout << "    -> Topological sorting not done." << endl;
    if(tri_topologique_reverse_done) cout << "    -> Topological sorting (in reverse) done from " << tri_topologique_reverse_source << endl;
    else cout << "    -> Topological sorting (in reverse) not done." << endl;

    cout << "   - Displayed :\n";
    cout << "\tVertex_Number\n";
    if(dijkstra_done)					cout << "\tDdist := Dijkstra_distance_from " << dijkstra_source << "\tDpred := Dijkstra_predecessor\n";
    if(dijkstra_reverse_done)			cout << "\tDRdist:= Dijkstra_distance_IN_REVERSE_GRAPH_from" << dijkstra_reverse_source << "\tDpred := Dijkstra_predecessor_IN_REVERSE\n";
    if(tri_topologique_done)			cout << "\tTordre := number_of_this_vertex_in_topological_order_from" << tri_topologique_source<< "\n";
    if(tri_topologique_reverse_done)	cout << "\tTRordre:= number_of_this_vertex_in_REVERSE_topological_order_from " << tri_topologique_reverse_source << "\n";
    cout << "\tSons := neighbours in (oriented) graph(weight)\n";
    cout << "\tPredecessors := neighbours in REVERSE graph(weight)\n";
    cout << " ==============............................... ============== \n";

    for(int i = 0; i < size; ++i){
        if(conv) {conv(i); cout << endl;}
        else cout << "S" << i	<< "\n";
        if(dijkstra_done)					cout << "\tDdist  : " << dist(i) << "\tDpred : " << pred(i) << "\n";
        if(dijkstra_reverse_done)			cout << "\tDRdist : " << dist_reverse(i) << "\t DRpred:" << pred_reverse(i) << "\n";
        if(tri_topologique_done)			cout << "\tTordre : " << number(i) << "\n";
        if(tri_topologique_reverse_done)	cout << "\tTRordre: " << number_reverse(i) << "\n";
        cout << "\tSons : ";
                for(int j = 0; j < (int) t[i].size(); ++j){
            if(conv) {conv(t[i][j].dest); cout << "(" << t[i][j].weight << "),   ";}
            else cout << "S" << t[i][j].dest << "(" << t[i][j].weight << "),   ";
        }
        cout << endl;
        cout << "\tPredecessors : ";
                for(int j = 0; j < (int) reverse[i].size(); ++j){
            if(conv) {conv(reverse[i][j].dest); cout << "(" << reverse[i][j].weight << "),   ";}
            else cout << reverse[i][j].dest << "(" << reverse[i][j].weight << "),   ";
        }
        cout << endl;
    }
};




/* void Graphe::show_order(){
    cout << "Order of the " << current_index_topo << " vertices seen from the source (";
    showpos(source);
    cout << "): " << endl;
    for(int i = 0; i < current_index_topo; ++i){
        cout << i+1 << "\t";
        ordre[i];
        cout << endl;
    }
} */



/* ---- class grille 2D ---- */


grille::grille(int _largeur) : width(_largeur){
    source = -1;
    arrivee = -1;
        //showpostatic(_largeur);
};

int grille::ID(int i, int j){
    return (i - 1) * width + (j - 1);
}
int grille::pos_x(int ID){
    return (ID % width) + 1;
}
int grille::pos_y(int ID){
    return (ID / width) + 1;
}
void grille::showpos(int ID){
    cerr << ID << " (" << pos_y(ID) << "," << pos_x(ID) << ") ";
}

void grille::convert_to_positions_x(vector<int> &IDs, vector<int> &pos_to_fill){
    int lg = IDs.size();
    pos_to_fill.resize(lg);
    for(int i = 0; i < lg; ++i){
        pos_to_fill[i] = pos_x(IDs[i]);
    }
}

void grille::showpostatic(int ID){
    ID++;
        /*static int largeur_static = -2;
    if(largeur_static == -2) largeur_static = ID;
    else {
               cout << ID << "[" << (ID / largeur_static) + 1 << "," << (ID % largeur_static) + 1 << "]";
        }*/
}





void Graph::compute_density_from_nb_paths(int source, int destination){
    int size = nb_sommets();
    vector<int> nb_paths(size, 0);
    vector<int> nb_paths_reverse(size, 0);
    density_graph.resize(size);

    nb_paths[source] = 1;
    nb_paths_reverse[destination] = 1;

    vector<bool> test_visited(size, false);			// attention, si on teste que les prédécesseurs ont été déjà vus,
                                // ça ne marche pas parce qu'ileexiste des prédécesseurs non accessibles depuis la source.

    for(int i = 0; i < nb_sommets_ordonnes(); ++i){
        list::iterator it;
        int u = ordre(i);
        test_visited[u] = true;
        for(it = reverse[u].begin(); it != reverse[u].end(); ++it){
            int prev = it->dest;
            if((test_visited[prev] == false) && (pred_topo(prev) >= 0)){		// this should never occur ...
                cerr << "ERR : topological order problem : vertex " << u << " has an unseen predecessor :" << prev << " Please check the graph is acyclic\n";
            }
            if(egal(dist(prev) + it->weight, dist(u))){
                nb_paths[u] += nb_paths[prev];
            }
        }
    }
    test_visited.clear();
    test_visited.resize(size, false);

    for(int i = 0; i < nb_sommets_ordonnes_reverse(); ++i){
        list::iterator it;
        int u = ordre_reverse(i);
        test_visited[u] = true;
        for(it = t[u].begin(); it != t[u].end(); ++it){
            int prev = it->dest;
            if((test_visited[prev] == false) && (pred_topo_reverse(prev) >= 0)) { //ie if also visited during topological algo
                cerr << "ERR : topological REVERSE order problem : vertex " << u << " has an unseen predecessor :" << prev <<
                        " Please check the graph is acyclic\n";
            }
            if(egal(dist(prev) - it->weight, dist(u))){
                nb_paths_reverse[u] += nb_paths_reverse[prev];
            }
        }
    }

    int nb_paths_src_dest = nb_paths[destination];
    if(nb_paths_reverse[source] != nb_paths_src_dest) {cerr << "ERR : finds more paths in one way than in another\n";}

    if(nb_paths_src_dest == 0) cerr << "ERR compute_density_from_recursion() : destination not reached.\n";
    for(int i = 0; i < size; ++i){
        density_graph[i] = ((valeur) nb_paths[i] * (valeur) nb_paths_reverse[i] ) / ((valeur) nb_paths_src_dest);
    }

}



vector<int> Graph::trace_one_optimal_path(int dest){
    vector<int> une_solution_optimale;
    int i = dest;
    int cpt = 0;
    while ((pred(i) != -1) && (cpt < 100)){
        cpt++;
        if(cpt > 2) {
            une_solution_optimale.push_back(i);
        }
        i = pred(i);
    }
    return une_solution_optimale;
}

void Graph::show_one_best_solution(int dest){

    vector<int> une_solution_optimale = trace_one_optimal_path(dest);
    //cerr << "    => Coût optimal :" << (dist(dest)) - CstAntiNeg << endl;
    cerr << "    => One optimal path : (" << une_solution_optimale.size() << " particles) :\t";
        for(int i = 0; i < (int) une_solution_optimale.size(); ++i){
        cerr << une_solution_optimale[i] << "\t";
    }
    cerr << endl;
    int j = une_solution_optimale.size() - 1;
    cerr << "\t\t";
    int N = une_solution_optimale.size();
    for(int i = 0; i < N; ++i){
        if((j >= 0) && (une_solution_optimale[j] == i+1)){
            cerr << "N";
            j--;
        } else {
            cerr << ".";
        }
    }
    cerr << endl;
}








/* -------------------- Structure for storing a path inside a graph --------------- */

void path::copy(const path &ch2){
    content = ch2.content;
/*    int N = ch2.size();
    ch.clear();
    ch.resize(N);
    for(int i = 0; i < N; ++i){
        ch[i] = ch2.ch[i];
    }*/
}
int path::operator [] (int i) {
    return content[i];
}
path::path(vector<int> ch2){
    content.clear();
    content.resize(ch2.size());
        for(int i = 0; i < (int) ch2.size(); ++i){
        content[i] = ch2[i];
    }
}
path::path(){
    content.clear();
}
int path::size(){
    return content.size();
}
void path::push_back(int z){
    content.push_back(z);
}
void path::show(){
    for(int i = 0; i < size(); ++i){
        cerr << content[i] << ",   ";
    }
}



/* -------------------- Structure for manipulating a list of pathes inside a graph, with common operations on this listt ---------------- */

void ensemble_paths::add(path &ch){
    size++;
    table.resize(size);
    table[size-1].copy(ch);
}
path ensemble_paths::operator[](int i){
    return table[i];
}
ensemble_paths::ensemble_paths() {table.clear(); destination = 0; size = 0;}
void ensemble_paths::set_destination(int _destination){
    destination = _destination;
}
void ensemble_paths::clear(){
    size = 0;
    table.clear();
}
void ensemble_paths::add(struct ensemble_paths &ens){
    int to_add = ens.size;
    size = size + to_add;
    table.resize(size);
    for(int i = 0; i < to_add; ++i){
        table[size - to_add + i].copy(ens[i]);
        table[size - to_add + i].push_back(destination);
    }
}
void ensemble_paths::addElementPath(){
    size = size + 1;
    table.resize(size);
    vector<int> start (1,destination);
    table[size - 1].copy(start);
}
void ensemble_paths::show(){
    for(int i = 0; i < size; ++i){
        table[i].show();
        cerr << endl;
    }
}


/* ------- Exemple d'utilisation de la structure ens_chemins ---------
  l<-----	4 <----- 2 <----- 1
 l		^				  l
 6 <-l		l------- 3 <------l
 l                l
 l<----- 5 <------l

 ens_chemins big_table[MAX];
 int N = 10;
 for(int i = 0; i < N; ++i){
    big_table[i].set_destination(i);
 }
 big_table[1].addElementPath();
 big_table[2].add(big_table[1]);
 big_table[3].add(big_table[1]);
 big_table[4].add(big_table[2]);
 big_table[4].add(big_table[3]);
 big_table[5].add(big_table[3]);
 big_table[6].add(big_table[4]);
 big_table[6].add(big_table[5]);
 for(int i = 1; i <= 6; ++i){
    cout << "Chemins to " << i << " :\n";
    big_table[i].show();
    cout << "\n";
 }
 ----------------------------------------------------------------------- */



// all the paths with optimal distance from the source that can reach 'ending node'. Usually called with (destination, source)
void Graph::trace_all_paths(int ending_node, int source){

    int size = nb_sommets();

    // prepare the table of paths from the source to the
    big_table.clear();
    big_table.resize(size);
    for(int i = 0; i < size; ++i){
        big_table[i].set_destination(i);
    }
    big_table[ending_node].addElementPath();// the path which is just the destination

    vector<bool> test_visited(size, false);
    for(int i = 0; i < nb_sommets_ordonnes_reverse(); ++i){
        int u = ordre_reverse(i);
        test_visited[u] = true;

        list::iterator it;
        for(it = t[u].begin(); it != t[u].end(); ++it){
            int prev = it->dest;
            if((test_visited[prev] == false) && (pred_topo_reverse(prev) >= 0)){
                cerr << "ERR : topological reverse order problem : vertex " << u;
                cerr << " has an unseen predecessor :" << prev << " Please check the graph is acyclic\n";
            }
            // only takes optimal moves
            if(egal(dist_reverse(prev) + it->weight, dist_reverse(u))){
                big_table[u].add(big_table[prev]);
            }
        }
    }

    if(DBG) big_table[source].show();

    ensemble_solutions.clear();
    int TOT = big_table[source].size;
    ensemble_solutions.resize(TOT);

    for(int k = 0; k < TOT; ++k){
        int ch_size = big_table[source][k].size() - 2;	// remove arrivee and the state just before
        //cerr << "    .  solution optimale n " << k << " : " <<  ch_size << " particules) :\t";
        for(int i = 1; i < ch_size; ++i){		// remove source (first one).
            ensemble_solutions[k].push_back(big_table[source][k][ch_size - i + 1]);
        }
    }
}


/*void graphe::show_all_best_solutions(){
    int nb_solutions = ensemble_solutions.size();
    cerr << "    => Ensemble des solutions optimales (" << nb_solutions << "(" << endl;
    for(int k = 0; k < nb_solutions; ++k){
        int j =0;
        cerr << "\t\t";
        for(int i = 0; i < N; ++i){
            if((j < (int) ensemble_solutions[k].size()) && (ensemble_solutions[k][j] == i+1)){
                cerr << "N";
                j++;
            } else {
                cerr << ".";
            }
        }
        cerr << "\t(" << ensemble_solutions[k].size() << " part.)\t";
        for(int i = 0; i < (int) ensemble_solutions[k].size(); ++i){
            cerr << ensemble_solutions[k][i] << ", ";
        }
        cerr << endl;
    }
}*/


