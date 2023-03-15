/*--------------------------------------------------------------------------------
Copyright (C) 2023, Abhoy Kole, Kamalika Datta, Indranil Sengupta, Rolf Drechsler

Quantum Mapping Tool (QMT) for mapping Quantum circuits to IBM architecture
released under the MIT licence.

All right reserved.
----------------------------------------------------------------------------------*/

#include "costmetrics.hpp"

void getPath(grph::node_t p1, grph::node_t p2, std::vector<grph::node_t>& path, grph::Graph& pg) {
    if (pg.getNodeInPath(p1, p2) != p2) {
        getPath(p1, pg.getNodeInPath(p1, p2), path, pg);
        getPath(pg.getNodeInPath(p1, p2), p2, path, pg);
    } else
        path.push_back(p2);
}

// ibm qx architecture
void getAllPaths(grph::node_t p1, grph::node_t p2, std::vector<std::vector<grph::node_t> >& paths, grph::Graph& pg, int i) {
    grph::node_t px = pg.getNodeInPath(p1, p2);
    if (px != p2) {  // pg.getNodeInPath(p1, p2)
        paths[i].push_back(p1);
        getAllPaths(px, p2, paths, pg, i);

        grph::node_t py = pg.getNodeInPath2(p1, p2);
        if (pg.getDistance(px, p2) == pg.getDistance(py, p2)) {
            bool isParent = false;
            for (int i = 0; i < paths.size(); i++) {
                for (int j = 0; j < paths[i].size(); j++) {
                    if (paths[i][j] == py) isParent = true;
                    break;
                }
            }
            if (py != -1 && isParent == false) {
                paths.push_back(std::vector<grph::node_t>());
                int k = paths.size() - 1;

                for (std::vector<grph::node_t>::iterator p = paths[i].begin(); p != paths[i].end(); ++p) {
                    paths[k].push_back(*p);
                    if (*p == p1) break;
                }

                getAllPaths(py, p2, paths, pg, k);
            }
        }
    } else {
        paths[i].push_back(p1);
        paths[i].push_back(p2);
    }
}

void removeDuplicatePath(std::vector<std::vector<grph::node_t> >& paths) {
    for (std::vector<std::vector<grph::node_t> >::iterator p1 = paths.begin(); p1 != paths.end() - 1; ++p1) {
        unsigned pos1 = std::distance(paths.begin(), p1);
        for (std::vector<std::vector<grph::node_t> >::iterator p2 = p1 + 1; p2 != paths.end(); ++p2) {
            bool identical = true;
            for (std::vector<grph::node_t>::iterator q1 = (*p1).begin(), q2 = (*p2).begin(); q1 != (*p1).end(); ++q1, ++q2) {
                if (*q1 != *q2) identical = false;
            }
            if (identical) {
                unsigned pos2 = std::distance(paths.begin(), p2);
                paths.erase(p2);
                p2 = paths.begin() + pos2 - 1;
                if (p2 == paths.end()) break;
            }
        }
        p1 = paths.begin() + pos1;
        if (p1 + 1 == paths.end()) break;
    }
}

void addInitNode(std::vector<std::vector<grph::node_t> >& paths, grph::node_t nd) {
    for (int i = 0; i < paths.size(); i++) paths[i].insert(paths[i].begin(), nd);
}

bool isInvalidPath(std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& plmap) {
    bool invalidPath = false;
    for (const auto& p : path) {
        if (plmap.find(p) == plmap.end()) {
            invalidPath = true;
            break;
        }
    }
    return invalidPath;
}

void displayAllPaths(std::vector<std::vector<grph::node_t> >& paths) {
    for (std::vector<std::vector<grph::node_t> >::iterator p = paths.begin(); p != paths.end(); ++p) {
        std::cout << "Path: ";
        for (std::vector<grph::node_t>::iterator q = (*p).begin(); q != (*p).end(); ++q)
            std::cout << *q << " ";
        std::cout << std::endl;
    }
}

void displayPath(std::vector<grph::node_t>& path) {
    std::cout << "Path: ";
    for (std::vector<grph::node_t>::iterator p = path.begin(); p != path.end(); ++p)
        std::cout << *p << " ";
    std::cout << std::endl;
}

// IBM QX Mapping
mincost_t minCostPathIndex(std::vector<grph::node_t>& path, std::map<grph::node_t, grph::node_t>& temp_phmap, grph::Graph& pg, grph::Graph& temp_g) {
    std::vector<grph::node_t> nodes = temp_g.getNodes();
    mincost_t mincost;
    bool empty = true;
    for (std::vector<grph::node_t>::iterator p = path.begin(); p != path.end() - 1; ++p) {
        grph::node_t control = *p;
        grph::node_t target = *(p + 1);

        std::map<grph::node_t, grph::node_t> temp_phmap2 = std::map<grph::node_t, grph::node_t>();
        std::map<grph::node_t, grph::node_t> temp_phmap3 = std::map<grph::node_t, grph::node_t>();
        for (std::vector<grph::node_t>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
            temp_phmap2[*n] = temp_phmap[*n];
            temp_phmap3[temp_phmap[*n]] = *n;
        }

        double cur_swap_cost = 0;
        // forward move in path
        for (std::vector<grph::node_t>::iterator pf = path.begin(); pf != p; ++pf) {
            // perform swaping in temporary location
            cur_swap_cost += 1;

            grph::node_t temp = temp_phmap2[temp_phmap3[*pf]];
            temp_phmap2[temp_phmap3[*pf]] = temp_phmap2[temp_phmap3[*(pf + 1)]];
            temp_phmap2[temp_phmap3[*(pf + 1)]] = temp;

            temp = temp_phmap3[*pf];
            temp_phmap3[*pf] = temp_phmap3[*(pf + 1)];
            temp_phmap3[*(pf + 1)] = temp;
        }
        // backward move in path
        for (std::vector<grph::node_t>::reverse_iterator pb = path.rbegin(), end(p + 2); pb != end; ++pb) {
            // perform swaping in temporary location
            cur_swap_cost += 1;

            grph::node_t temp = temp_phmap2[temp_phmap3[*pb]];
            temp_phmap2[temp_phmap3[*pb]] = temp_phmap2[temp_phmap3[*(pb + 1)]];
            temp_phmap2[temp_phmap3[*(pb + 1)]] = temp;

            temp = temp_phmap3[*pb];
            temp_phmap3[*pb] = temp_phmap3[*(pb + 1)];
            temp_phmap3[*(pb + 1)] = temp;
        }

        // current control line is not valid for physical coupling map, i.e *(p+1)->*p
        if (pg.getDistance(*p, *(p + 1)) == RCOST) {
            // perform swaping in temporary location
            cur_swap_cost += RCOST;

            grph::node_t temp = temp_phmap2[temp_phmap3[*p]];
            temp_phmap2[temp_phmap3[*p]] = temp_phmap2[temp_phmap3[*(p + 1)]];
            temp_phmap2[temp_phmap3[*(p + 1)]] = temp;

            temp = temp_phmap3[*p];
            temp_phmap3[*p] = temp_phmap3[*(p + 1)];
            temp_phmap3[*(p + 1)] = temp;
        }
        double coupling_cost = mappingCost(temp_phmap2, pg, temp_g);
        if (empty || coupling_cost < mincost.swap_cost) {
            mincost.control = control;
            mincost.target = target;
            mincost.swaps = cur_swap_cost;
            mincost.swap_cost = coupling_cost;
            empty = false;
        }
    }
    return mincost;
}

// IBM QX
void updateMapping(std::vector<grph::node_t>& path, std::map<grph::node_t, grph::node_t>& temp_phmap, grph::Graph& pg,
                   std::vector<grph::node_t>& nodes, grph::node_t control, grph::node_t target) {
    std::map<grph::node_t, grph::node_t> temp_phmap2 = std::map<grph::node_t, grph::node_t>();
    for (std::vector<grph::node_t>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        temp_phmap2[temp_phmap[*n]] = *n;
    }
    // forward move in path
    for (std::vector<grph::node_t>::iterator pf = path.begin(); *pf != control; ++pf) {
        grph::node_t temp = temp_phmap[temp_phmap2[*pf]];
        temp_phmap[temp_phmap2[*pf]] = temp_phmap[temp_phmap2[*(pf + 1)]];
        temp_phmap[temp_phmap2[*(pf + 1)]] = temp;
        temp = temp_phmap2[*pf];
        temp_phmap2[*pf] = temp_phmap2[*(pf + 1)];
        temp_phmap2[*(pf + 1)] = temp;
    }
    // backward move in path
    for (std::vector<grph::node_t>::reverse_iterator pb = path.rbegin(); *pb != target; ++pb) {
        grph::node_t temp = temp_phmap[temp_phmap2[*pb]];
        temp_phmap[temp_phmap2[*pb]] = temp_phmap[temp_phmap2[*(pb + 1)]];
        temp_phmap[temp_phmap2[*(pb + 1)]] = temp;

        temp = temp_phmap2[*pb];
        temp_phmap2[*pb] = temp_phmap2[*(pb + 1)];
        temp_phmap2[*(pb + 1)] = temp;
    }
    // current control line is not valid for physical coupling map, i.e *(p+1)->*p
    if (pg.getDistance(control, target) == RCOST) {
        grph::node_t temp = temp_phmap[temp_phmap2[control]];
        temp_phmap[temp_phmap2[control]] = temp_phmap[temp_phmap2[target]];
        temp_phmap[temp_phmap2[target]] = temp;

        temp = temp_phmap2[control];
        temp_phmap2[control] = temp_phmap2[target];
        temp_phmap2[target] = temp;
    }
}

// IBM QX: computing coupling cost
double swapCost(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& cg, int base, unsigned depth) {
    double swap_cost = 0;
    grph::Graph temp_g = cg;
    std::vector<grph::node_t> nodes = temp_g.getNodes();
    std::map<grph::node_t, std::vector<grph::node_t> > edge_list = temp_g.getEdges();
    std::map<grph::node_t, bool> rev_map = std::map<grph::node_t, bool>();
    std::map<grph::node_t, grph::node_t> temp_phmap = std::map<grph::node_t, grph::node_t>();
    std::map<grph::node_t, grph::node_t> temp_phmap2 = std::map<grph::node_t, grph::node_t>();

    for (std::vector<grph::node_t>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
        temp_phmap[*n] = phmap[*n];
        temp_phmap2[temp_phmap[*n]] = *n;
        rev_map[phmap[*n]] = false;
    }
    for (unsigned d = depth - 1; d >= 0; d--) {
        double cur_weight = std::pow(base, d);
        bool edge_present = false;
        for (std::vector<grph::node_t>::iterator n = nodes.begin(); n != nodes.end(); ++n) {
            for (std::vector<grph::node_t>::iterator m = edge_list[*n].begin(); m != edge_list[*n].end(); ++m) {
                if (temp_g.getWeight(*n, *m) > 0) edge_present = true;
                if (temp_g.getWeight(*n, *m) >= cur_weight) {
                    temp_g.updateEdgeWeight(*n, *m, cur_weight);
                    if (pg.getDistance(temp_phmap[*n], temp_phmap[*m]) > 0) {
                        std::vector<std::vector<grph::node_t> > paths = std::vector<std::vector<grph::node_t> >();
                        paths.push_back(std::vector<grph::node_t>());
                        getAllPaths(temp_phmap[*n], temp_phmap[*m], paths, pg, 0);
                        removeDuplicatePath(paths);
                        mincost_t mincost;
                        bool empty = true;
                        int index;
                        for (int i = 0; i < paths.size(); i++) {
                            if (isInvalidPath(paths[i], temp_phmap2)) continue;
                            mincost_t tempcost = minCostPathIndex(paths[i], temp_phmap, pg, temp_g);
                            if (empty || mincost.swap_cost > tempcost.swap_cost || (mincost.swap_cost == tempcost.swap_cost && mincost.swaps > tempcost.swaps)) {
                                mincost.control = tempcost.control;
                                mincost.target = tempcost.target;
                                mincost.swaps = tempcost.swaps;
                                mincost.swap_cost = tempcost.swap_cost;
                                empty = false;
                                index = i;
                            }
                        }
                        updateMapping(paths[index], temp_phmap, pg, nodes, mincost.control, mincost.target);
                        swap_cost += mincost.swaps;
                        if (pg.getDistance(mincost.control, mincost.target) == RCOST) {
                            if (rev_map[mincost.control]) swap_cost -= RCOST / 2;
                            if (rev_map[mincost.target]) swap_cost -= RCOST / 2;
                            rev_map[mincost.control] = true;
                            rev_map[mincost.target] = true;
                        } else {
                            rev_map[mincost.control] = false;
                            rev_map[mincost.target] = false;
                        }
                    }
                }
            }
        }
        if (!edge_present) break;
    }
    return swap_cost;
}

// IBM QX : estimation using coupling cost metric for global ordering
double mappingCost(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& cg) {
    double cost = 0;
    std::vector<grph::node_t> nodes = cg.getNodes();
    std::map<grph::node_t, std::vector<grph::node_t> > edge_list = cg.getEdges();
    for (const auto& n : nodes) {
        for (const auto& m : edge_list[n]) {
            cost += pg.getDistance(phmap[n], phmap[m]) * cg.getWeight(n, m);
        }
    }

    return cost;
}

// IBM QX : estimation using coupling cost metric for local ordering
double mappingCost(std::map<grph::node_t, grph::node_t>& phmap, std::map<rev::line_t, rev::line_t>& lgmap,
                   grph::Graph& pg, grph::Graph& cg) {
    double cost = 0;
    std::vector<grph::node_t> nodes = cg.getNodes();
    std::map<grph::node_t, std::vector<grph::node_t> > edge_list = cg.getEdges();
    for (const auto& n : nodes) {
        for (const auto& m : edge_list[n])
            cost += pg.getDistance(phmap[lgmap[n]], phmap[lgmap[m]]) * cg.getWeight(n, m);
    }

    return cost;
}

// IBM QX : estimation using remote CNOT cost metric for global ordering
double remoteCNOTCost(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& cg) {
    double cost = 0;
    std::vector<grph::node_t> nodes = cg.getNodes();
    std::map<grph::node_t, std::vector<grph::node_t> > edge_list = cg.getEdges();
    for (const auto& n : nodes) {
        for (const auto& m : edge_list[n]) {
            cost += 4 * pg.getDistance(phmap[n], phmap[m]) * cg.getWeight(n, m);
        }
    }

    return cost;
}

// IBM QX : estimation using remote CNOT cost metric for global ordering
double getWCost(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& cg) {
    double cost = 0;
    std::vector<grph::node_t> nodes = cg.getNodes();
    std::map<grph::node_t, std::vector<grph::node_t> > edge_list = cg.getEdges();
    for (const auto& n : nodes) {
        for (const auto& m : edge_list[n]) {
            grph::distance_t dis = pg.getDistance(phmap[n], phmap[m]);
            if (dis == 0)
                cost += cg.getWeight(n, m);
            else
                cost += 4 * dis * cg.getWeight(n, m);
        }
    }

    return cost;
}

// IBM QX : estimation using remote CNOT cost metric for local ordering
double remoteCNOTCost(std::map<grph::node_t, grph::node_t>& phmap, std::map<rev::line_t, rev::line_t>& lgmap,
                      grph::Graph& pg, grph::Graph& cg) {
    double cost = 0;
    std::vector<grph::node_t> nodes = cg.getNodes();
    std::map<grph::node_t, std::vector<grph::node_t> > edge_list = cg.getEdges();
    for (const auto& n : nodes) {
        for (const auto& m : edge_list[n]) {
            cost += 4 * pg.getDistance(phmap[lgmap[n]], phmap[lgmap[m]]) * cg.getWeight(n, m);
        }
    }
    return cost;
}

double move_cost(pair<grph::node_t, grph::node_t>& swapped_nodes, grph::Graph& pg) {
    return (pg.getDistance(swapped_nodes.first, swapped_nodes.second) + 1) * 3;  // Number of CNOT operations for SWAPPing state
}

double computeDepth(rev::Circuit& ckt) {
    std::map<grph::node_t, double> depths = std::map<grph::node_t, double>();
    std::vector<rev::Gate> gates = ckt.getGates();
    std::vector<rev::line_t> lines = ckt.getLines();
    for (std::vector<rev::line_t>::iterator l = lines.begin(); l != lines.end(); ++l) {
        depths[*l] = 0;
    }
    for (std::vector<rev::Gate>::iterator g = gates.begin(); g != gates.end(); ++g) {
        if (g->getType() == rev::CX) {
            double cdepth = depths[g->getControls()[0]] > depths[g->getTargets()[0]] ? depths[g->getControls()[0]] + 1 : depths[g->getTargets()[0]] + 1;
            depths[g->getControls()[0]] = cdepth;
            depths[g->getTargets()[0]] = cdepth;
        } else {
            depths[g->getTargets()[0]] = depths[g->getTargets()[0]] + 1;
        }
    }
    double maxDepth = 0;
    for (std::vector<rev::line_t>::iterator l = lines.begin(); l != lines.end(); ++l) {
        if (maxDepth < depths[*l]) maxDepth = depths[*l];
    }
    return maxDepth;
}

long computeCNOTGates(rev::Circuit& ckt) {
    long CNOT_Gates = 0;
    std::vector<rev::Gate> gates = ckt.getGates();
    for (auto& g : gates) {
        if (g.getType() == rev::CX) CNOT_Gates += 1;
    }
    return CNOT_Gates;
}

void generateRemoteCNOT(std::vector<rev::Gate>& gates, rev::Gate& g, std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg) {
    rev::line_t control = phmap[g.getControls()[0]];  // mapped physical control qubit
    rev::line_t target = phmap[g.getTargets()[0]];    // mapped physical target qubit
    grph::distance_t d = pg.getDistance(control, target);
    while (d > 0) {  // if distance between control and target qubit is greater than 0
        for (const auto& i : phmap) {
            grph::distance_t d_new = pg.getDistance(i.second, target);
            // find a mapped qubit adjacent to control qubit and closer to target qubit
            if (control != i.second && target != i.second && (d - d_new) == 1) {
                gates.push_back(rev::Gate(rev::CX, control, i.second));
                control = i.second;
                d = d_new;
                break;
            }
        }
    }
    std::vector<rev::Gate> gates_seq = std::vector<rev::Gate>();
    for (auto g = gates.rbegin(); g != gates.rend(); ++g)
        gates_seq.push_back(rev::Gate(rev::CX, g->getControls()[0], g->getTargets()[0]));

    gates.push_back(rev::Gate(rev::CX, control, target));
    gates.insert(gates.end(), gates_seq.begin(), gates_seq.end());
    std::vector<rev::Gate> gates_seq2 = std::vector<rev::Gate>();  //(gates.begin() + 1, gates.end());
    for (auto g = gates.begin() + 1; g + 1 != gates.end(); ++g)
        gates_seq2.push_back(rev::Gate(rev::CX, g->getControls()[0], g->getTargets()[0]));
    gates.insert(gates.end(), gates_seq2.begin(), gates_seq2.end());
}
