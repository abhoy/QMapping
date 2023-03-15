/*--------------------------------------------------------------------------------
Copyright (C) 2023, Abhoy Kole, Kamalika Datta, Indranil Sengupta, Rolf Drechsler

Quantum Mapping Tool (QMT) for mapping Quantum circuits to IBM architecture
released under the MIT licence.

All right reserved.
----------------------------------------------------------------------------------*/

#include <iomanip>
#include <map>
#include <vector>

#include "circuit.hpp"

#ifndef GRAPH_HPP
#define GRAPH_HPP

#define FCOST 0    // FORWARD COUPLING COST
#define RCOST 0.6  // REVERSE COUPLING COST

namespace grph {
enum { OLD = 0,
       NEW = 1 };
enum class window_type { fixed,
                         variable };
typedef long dimension_t;  // data type for representing width and height of a rectangular regular architecture
typedef unsigned nnc_t;    // nearest neighbor cost type: OLD (gate count based) or NEW (priority based)
typedef unsigned node_t;
typedef int map_t;
typedef int indeg_t;
typedef int outdeg_t;
typedef float weight_t;
typedef double distance_t;
typedef bool reverse_t;
typedef unsigned window_t;
typedef int deg_t;  // degree type to keep track of number of adjacent nodes
typedef unsigned base_t;
const int INVALID = -20;

int roundNum(double num);

class Coordinate {
   public:
    int x, y, z;
    Coordinate();
    Coordinate(int x, int y, int z);
    bool isIdentical(Coordinate& pos);
};
typedef struct {
    node_t node;
    indeg_t indeg;
    outdeg_t outdeg;
} Degree_t;

class Graph {
   public:
    Graph();
    Graph(const Graph& g);
    Graph(rev::Circuit&);
    Graph(rev::Circuit&, window_t w);                        // constant size window
    Graph(rev::Circuit&, deg_t maxdeg, window_type type);    // variable size window
    Graph(rev::Circuit&, nnc_t type, window_t w, base_t x);  // for coupling cost
    Graph(std::vector<rev::line_t>& lines, std::vector<rev::Gate>& gates, int index, nnc_t type, window_t w, base_t x);
    Graph(std::vector<rev::line_t>& lines, std::vector<rev::Gate>& gates, int index, window_t w);                      // Constant size window
    void init();
    void addNode(node_t nd);
    void addEdge(node_t nd1, node_t nd2);                                // For Qubit Coupling Graph
    void addDirEdge(node_t nd1, node_t nd2);                             // Directed graph
    void addEdge2(node_t nd1, node_t nd2, base_t x, window_t y, int i);  // For Qubit Interaction Graph
    void addEdge2(node_t nd1, node_t nd2, int weight);                   // For Qubit Interaction Graph
    void addEdge3(node_t nd1, node_t nd2, int weight);                   // For Qubit Interaction Graph
    weight_t getWeight(node_t nd1, node_t nd2);
    void setWeight(node_t nd1, node_t nd2, weight_t w);
    void displayGraph();
    void displayDistance();
    void displayPath();
    std::vector<node_t>& getNodes();
    std::map<node_t, std::vector<node_t> >& getEdges();
    std::map<node_t, std::vector<node_t> >& getRedges();
    std::map<node_t, std::map<node_t, weight_t> >& getEdgeWeights();
    node_t getNodeInPath(node_t p1, node_t p2);
    node_t getNodeInPath2(node_t p1, node_t p2);
    distance_t getDistance(node_t p1, node_t p2);
    void readLayout(string file);                      // Directed Graph Layout
    void readLayout2(std::istream& is);                // Undirected Graph Layout
    void computeDistance();                            // For Directed Graph Layout
    void computeDistance2();                           // For Undirected Graph Layout
    void computePath();                                // For Directed Graph Layout
    void computePath2();                               // For Undirected Graph Layout
    void mapLayout(Graph, std::map<node_t, node_t>&);  // map circuit qubits based on their degree of association
    void mapLayout2(Graph, std::map<node_t, node_t>&);
    void mapQubits(Graph& pg, std::map<node_t, node_t>& phmap, std::vector<node_t>& pqubits);  // map circuit qubits to specified qubits
    void removeCNOTEdge(rev::Gate& g, window_t w, base_t x);
    bool isEmpty();
    void operator=(const Graph& g);
    void updateEdgeWeight(node_t n, node_t m, double weight);

   private:
    std::vector<node_t> nodes;                      // all nodes
    std::map<node_t, std::vector<node_t> > edges;   // node to edges mapping
    std::map<node_t, std::vector<node_t> > redges;  // node to edges mapping

    std::map<node_t, std::map<node_t, weight_t> > weight;      // edge weight
    std::map<node_t, std::map<node_t, distance_t> > distance;  // distance between nodes
    std::map<node_t, std::map<node_t, node_t> > path;          // path
    std::map<node_t, std::map<node_t, node_t> > path2;         // path
    std::map<node_t, std::map<node_t, reverse_t> > reverse;    // reverse directed edge
    std::map<node_t, deg_t> degs;                              // degree of each node for undirected graph
};

}  // namespace grph
#endif /* GRAPH_HPP */
