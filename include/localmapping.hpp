/*--------------------------------------------------------------------------------
Copyright (C) 2023, Abhoy Kole, Kamalika Datta, Indranil Sengupta, Rolf Drechsler

Quantum Mapping Tool (QMT) for mapping Quantum circuits to IBM architecture
released under the MIT licence.

All right reserved.
----------------------------------------------------------------------------------*/

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <vector>

#include "circuit.hpp"
#include "costmetrics.hpp"
#include "graph.hpp"

#ifndef LOCAL_HPP
#define LOCAL_HPP


// IBM QX
mincost_t minCostPathIndex(std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& lgmap,
                           std::map<rev::line_t, rev::line_t>& glmap, std::map<rev::line_t, rev::line_t>& plmap,
                           std::map<grph::node_t, grph::node_t>& phmap, std::vector<rev::Gate>& gates,
                           std::vector<rev::line_t>& lines, grph::Graph& pg, grph::Graph& lg, unsigned pos, int base);
// IBM QX
void localorder(rev::Circuit& ckt, std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& lg, int base);

// IBM QX
mincost_t minCostPathIndexv2(std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& lgmap,
                             std::map<rev::line_t, rev::line_t>& glmap, std::map<rev::line_t, rev::line_t>& plmap,
                             std::map<grph::node_t, grph::node_t>& phmap, std::vector<rev::Gate>& gates,
                             std::vector<rev::line_t>& lines, grph::Graph& pg, grph::Graph& lg, unsigned pos, int base);

// IBM QX using only SWAP operation
void localorderv2(rev::Circuit& ckt, std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& lg, int base, rev::window_t window);

// insert Swap as cascade of CNOT gates
void insertSWAP(std::vector<rev::Gate>& swap_gates, std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& lgmap,
                std::map<rev::line_t, rev::line_t>& glmap, std::map<rev::line_t, rev::line_t>& plmap,
                rev::line_t control, rev::line_t target);
// insert Swap as an abstract gate
void insertSWAPv2(std::vector<rev::Gate>& swap_gates, std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& lgmap,
                  std::map<rev::line_t, rev::line_t>& glmap, std::map<rev::line_t, rev::line_t>& plmap,
                  rev::line_t control, rev::line_t target);

#endif /* LOCAL_HPP */
