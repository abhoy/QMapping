/*--------------------------------------------------------------------------------
Copyright (C) 2023, Abhoy Kole, Kamalika Datta, Indranil Sengupta, Rolf Drechsler

Quantum Mapping Tool (QMT) for mapping Quantum circuits to IBM architecture
released under the MIT licence.

All right reserved.
----------------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include "circuit.hpp"
#include "graph.hpp"

#ifndef OPTIMIZE_HPP
#define OPTIMIZE_HPP
enum class Metric_Type { RCNOT_COST,
                         COUPLING_COST };
enum class Approach_Type { A1,
                           A2,
                           A3 };
typedef struct {
    grph::node_t control;
    grph::node_t target;
    double swap_cost;
    double swaps;
} mincost_t;

void displayPath(std::vector<grph::node_t>& path);
bool isInvalidPath(std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& plmap);
void removeDuplicatePath(std::vector<std::vector<grph::node_t> >& paths);
void addInitNode(std::vector<std::vector<grph::node_t> >& paths, grph::node_t nd);
void displayAllPaths(std::vector<std::vector<grph::node_t> >& paths);
void getAllPaths(grph::node_t p1, grph::node_t p2, std::vector<std::vector<grph::node_t> >& paths, grph::Graph& pg, int i);

void getAPath(grph::node_t p1, grph::node_t p2, std::vector<grph::node_t>& path,
              std::pair<int, int> dim);
void getPath(grph::node_t p1, grph::node_t p2, std::vector<grph::node_t>& path, grph::Graph& pg, int i);

mincost_t minCostPathIndex(std::vector<grph::node_t>& path, std::map<grph::node_t, grph::node_t>& temp_phmap, grph::Graph& pg, grph::Graph& temp_g);

void updateMapping(std::vector<grph::node_t>& path, std::map<grph::node_t, grph::node_t>& temp_phmap, grph::Graph& pg,
                   std::vector<grph::node_t>& nodes, grph::node_t control, grph::node_t target);

// minimal estimation of coupling cost exploring all possible paths
double swapCost(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& cg, int base, unsigned depth);

// local ordering: estimation using coupling cost metric
double mappingCost(std::map<grph::node_t, grph::node_t>& phmap, std::map<rev::line_t, rev::line_t>& lgmap, grph::Graph& pg, grph::Graph& cg);  // IBM QX

// global ordering: estimation using coupling cost metric
double mappingCost(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& cg);  // ibm qx

// local ordering: estimation using remote CNOT cost metric
double remoteCNOTCost(std::map<grph::node_t, grph::node_t>& phmap, std::map<rev::line_t, rev::line_t>& lgmap, grph::Graph& pg, grph::Graph& cg);  // IBM QX

// global ordering: estimation using remote CNOT cost metric
double remoteCNOTCost(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& cg);  // IBM QX

// estimation using WCOST metric
double getWCost(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& cg);  // IBM QX

double computeDepth(rev::Circuit& ckt);
long computeCNOTGates(rev::Circuit& ckt);
double move_cost(pair<grph::node_t, grph::node_t>& swapped_nodes, grph::Graph& p);
void generateRemoteCNOT(std::vector<rev::Gate>& gates, rev::Gate& g, std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg);

#endif /* OPTIMIZE_HPP */
