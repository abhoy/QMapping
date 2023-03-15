/*--------------------------------------------------------------------------------
Copyright (C) 2023, Abhoy Kole, Kamalika Datta, Indranil Sengupta, Rolf Drechsler

Quantum Mapping Tool (QMT) for mapping Quantum circuits to IBM architecture
released under the MIT licence.

All right reserved.
----------------------------------------------------------------------------------*/

#include "localmapping.hpp"

// IBM QX
mincost_t minCostPathIndex(std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& lgmap,
                           std::map<rev::line_t, rev::line_t>& glmap, std::map<rev::line_t, rev::line_t>& plmap,
                           std::map<grph::node_t, grph::node_t>& phmap, std::vector<rev::Gate>& gates,
                           std::vector<rev::line_t>& lines, grph::Graph& pg, grph::Graph& lg, unsigned pos, int base) {
    mincost_t mincost;
    mincost.swap_cost = 0;
    mincost.swaps = 0;
    for (std::vector<grph::node_t>::iterator p = path.begin(); p != path.end() - 1; ++p) {
        rev::line_t cur_control = *p;
        rev::line_t cur_target = *(p + 1);
        std::map<rev::line_t, rev::line_t> lgmap_tmp = std::map<rev::line_t, rev::line_t>();  // temporary local to global line map
        std::map<rev::line_t, rev::line_t> glmap_tmp = std::map<rev::line_t, rev::line_t>();  // temporary global to local line map
        for (std::vector<rev::line_t>::iterator l = lines.begin(); l != lines.end(); ++l) {   // copy the current mapping to temporary location
            lgmap_tmp[*l] = lgmap[*l];
            glmap_tmp[*l] = glmap[*l];
        }
        double cur_swaps = 0;
        // forward move in path
        for (std::vector<grph::node_t>::iterator pf = path.begin(); pf != p; ++pf) {
            // perform swaping in temporary location
            cur_swaps += 1;
            // updating local to global map
            lgmap_tmp[glmap_tmp[plmap[*pf]]] = plmap[*(pf + 1)];
            lgmap_tmp[glmap_tmp[plmap[*(pf + 1)]]] = plmap[*pf];

            // updating global to local map
            rev::line_t l = glmap_tmp[plmap[*pf]];
            glmap_tmp[plmap[*pf]] = glmap_tmp[plmap[*(pf + 1)]];
            glmap_tmp[plmap[*(pf + 1)]] = l;
        }
        // backward move in path
        for (std::vector<grph::node_t>::reverse_iterator pb = path.rbegin(), end(p + 2); pb != end; ++pb) {
            // perform swaping in temporary location
            cur_swaps += 1;
            // updating local to global map
            lgmap_tmp[glmap_tmp[plmap[*pb]]] = plmap[*(pb + 1)];
            lgmap_tmp[glmap_tmp[plmap[*(pb + 1)]]] = plmap[*pb];

            // updating global to local map
            rev::line_t l = glmap_tmp[plmap[*pb]];
            glmap_tmp[plmap[*pb]] = glmap_tmp[plmap[*(pb + 1)]];
            glmap_tmp[plmap[*(pb + 1)]] = l;
        }

        // current control line is not valid for physical coupling map, i.e *(p+1)->*p
        if (pg.getDistance(*p, *(p + 1)) == RCOST) {
            // perform swaping in temporary location
            cur_swaps += RCOST;
            // updating local to global map
            lgmap_tmp[glmap_tmp[plmap[*p]]] = plmap[*(p + 1)];
            lgmap_tmp[glmap_tmp[plmap[*(p + 1)]]] = plmap[*p];

            // updating global to local map
            rev::line_t l = glmap_tmp[plmap[*p]];
            glmap_tmp[plmap[*p]] = glmap_tmp[plmap[*(p + 1)]];
            glmap_tmp[plmap[*(p + 1)]] = l;
        }
        double cur_swap_cost;
        if (lg.isEmpty() == false)
            cur_swap_cost = mappingCost(phmap, lgmap_tmp, pg, lg);
        else
            cur_swap_cost = cur_swaps;
        if (p == path.begin() || mincost.swap_cost > cur_swap_cost || (mincost.swap_cost == cur_swap_cost && mincost.swaps > cur_swaps)) {
            mincost.swap_cost = cur_swap_cost;
            mincost.swaps = cur_swaps;
            mincost.control = cur_control;
            mincost.target = cur_target;
        }
    }
    return mincost;
}

// IBM QX
mincost_t minCostPathIndexv2(std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& lgmap,
                             std::map<rev::line_t, rev::line_t>& glmap, std::map<rev::line_t, rev::line_t>& plmap,
                             std::map<grph::node_t, grph::node_t>& phmap, std::vector<rev::Gate>& gates,
                             std::vector<rev::line_t>& lines, grph::Graph& pg, grph::Graph& lg, unsigned pos, int base) {
    mincost_t mincost;
    mincost.swap_cost = 0;
    mincost.swaps = 0;
    // computing minimum cost insertion path index
    for (std::vector<grph::node_t>::iterator p = path.begin(); p != path.end() - 1; ++p) {
        rev::line_t cur_control = *p;
        rev::line_t cur_target = *(p + 1);
        std::map<rev::line_t, rev::line_t> lgmap_tmp = std::map<rev::line_t, rev::line_t>();  // temporary local to global line map
        std::map<rev::line_t, rev::line_t> glmap_tmp = std::map<rev::line_t, rev::line_t>();  // temporary global to local line map
        for (std::vector<rev::line_t>::iterator l = lines.begin(); l != lines.end(); ++l) {   // copy the current mapping to temporary location
            lgmap_tmp[*l] = lgmap[*l];
            glmap_tmp[*l] = glmap[*l];
        }
        double cur_swaps = 0;
        // forward move in path
        for (std::vector<grph::node_t>::iterator pf = path.begin(); pf != p; ++pf) {
            // perform swaping in temporary location
            cur_swaps += 1;
            // updating local to global map
            lgmap_tmp[glmap_tmp[plmap[*pf]]] = plmap[*(pf + 1)];
            lgmap_tmp[glmap_tmp[plmap[*(pf + 1)]]] = plmap[*pf];

            // updating global to local map
            rev::line_t l = glmap_tmp[plmap[*pf]];
            glmap_tmp[plmap[*pf]] = glmap_tmp[plmap[*(pf + 1)]];
            glmap_tmp[plmap[*(pf + 1)]] = l;
        }
        // backward move in path
        for (std::vector<grph::node_t>::reverse_iterator pb = path.rbegin(), end(p + 2); pb != end; ++pb) {
            // perform swaping in temporary location
            cur_swaps += 1;
            // updating local to global map
            lgmap_tmp[glmap_tmp[plmap[*pb]]] = plmap[*(pb + 1)];
            lgmap_tmp[glmap_tmp[plmap[*(pb + 1)]]] = plmap[*pb];

            // updating global to local map
            rev::line_t l = glmap_tmp[plmap[*pb]];
            glmap_tmp[plmap[*pb]] = glmap_tmp[plmap[*(pb + 1)]];
            glmap_tmp[plmap[*(pb + 1)]] = l;
        }

    
        double cur_swap_cost;
        if (lg.isEmpty() == false)
            cur_swap_cost = remoteCNOTCost(phmap, lgmap_tmp, pg, lg);  
        else
            cur_swap_cost = cur_swaps;
        if (p == path.begin() || mincost.swap_cost > cur_swap_cost || (mincost.swap_cost == cur_swap_cost && mincost.swaps > cur_swaps)) {
            mincost.swap_cost = cur_swap_cost;
            mincost.swaps = cur_swaps;
            mincost.control = cur_control;
            mincost.target = cur_target;
        }
    }
    return mincost;
}

// IBM QX
void localorder(rev::Circuit& ckt, std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& lg, int base) {
    std::vector<rev::Gate> gates = ckt.getGates();
    std::vector<rev::line_t> lines = ckt.getLines();
    unsigned size = gates.size();
    std::map<rev::line_t, rev::line_t> lgmap = std::map<rev::line_t, rev::line_t>();  // local to global line map
    std::map<rev::line_t, rev::line_t> glmap = std::map<rev::line_t, rev::line_t>();  // global to local line map
    std::map<rev::line_t, rev::line_t> plmap = std::map<rev::line_t, rev::line_t>();  // physical to logical line map
    std::vector<rev::line_t> new_lines = std::vector<rev::line_t>();

    for (std::vector<rev::line_t>::iterator l = lines.begin(); l != lines.end(); ++l) {
        lgmap[*l] = *l;
        glmap[*l] = *l;
        plmap[phmap[*l]] = *l;
        new_lines.push_back(phmap[*l]);
    }

    unsigned gcount = 0;
    unsigned tot_gate = gates.size();
    for (std::vector<rev::Gate>::iterator g = gates.begin(); g != gates.end(); ++g) {
        if (1 == g->getControls().size()) {
            gcount++;
            if (pg.getDistance(phmap[lgmap[g->getControls()[0]]], phmap[lgmap[g->getTargets()[0]]]) == 0) {  // adjacent gate
                g->updateLines(phmap[lgmap[g->getControls()[0]]], phmap[lgmap[g->getTargets()[0]]]);
            } else {  // qubits are not adjacent

                std::vector<std::vector<grph::node_t> > paths = std::vector<std::vector<grph::node_t> >();
                paths.push_back(std::vector<grph::node_t>());
                getAllPaths(phmap[lgmap[g->getControls()[0]]], phmap[lgmap[g->getTargets()[0]]], paths, pg, 0);
                removeDuplicatePath(paths);
                unsigned window;
                if (base == 2)
                    window = 20;
                else if (base == 3)
                    window = 80;
                else if (base == 4)
                    window = 60;
                unsigned pos = std::distance(gates.begin(), g);
                grph::Graph gr;
                if (g + 1 != gates.end()) {
                    gr = grph::Graph(lines, gates, pos + 1, grph::NEW, window, base);
                }
                int k = 0;  // path index
                mincost_t mincost;
                bool empty = true;
                for (int i = 0; i < paths.size(); i++) {
                    mincost_t tempcost;
                    if (isInvalidPath(paths[i], plmap)) {
                        continue;
                    }
                    tempcost = minCostPathIndex(paths[i], lgmap, glmap, plmap, phmap, gates, lines, pg, gr, pos, base);
                    if (empty || mincost.swap_cost > tempcost.swap_cost || (mincost.swap_cost == tempcost.swap_cost && mincost.swaps > tempcost.swaps)) {
                        mincost.swap_cost = tempcost.swap_cost;
                        mincost.swaps = tempcost.swaps;
                        mincost.control = tempcost.control;
                        mincost.target = tempcost.target;
                        k = i;
                        empty = false;
                    }
                }
                std::vector<rev::Gate> swap_gates = std::vector<rev::Gate>();
                // Actual swap gate insertion in forward path
                for (std::vector<grph::node_t>::iterator pf = paths[k].begin(); *pf != mincost.control; ++pf) {
                    // updating local to global map
                    lgmap[glmap[plmap[*pf]]] = plmap[*(pf + 1)];
                    lgmap[glmap[plmap[*(pf + 1)]]] = plmap[*pf];

                    // updating global to local map
                    rev::line_t l = glmap[plmap[*pf]];
                    glmap[plmap[*pf]] = glmap[plmap[*(pf + 1)]];
                    glmap[plmap[*(pf + 1)]] = l;

                    if (pg.getDistance(*pf, *(pf + 1)) == 0) {
                        swap_gates.push_back(rev::Gate(rev::CX, *pf, *(pf + 1)));
                        swap_gates.push_back(rev::Gate(rev::H, *pf));
                        swap_gates.push_back(rev::Gate(rev::H, *(pf + 1)));
                        swap_gates.push_back(rev::Gate(rev::CX, *pf, *(pf + 1)));
                        swap_gates.push_back(rev::Gate(rev::H, *pf));
                        swap_gates.push_back(rev::Gate(rev::H, *(pf + 1)));
                        swap_gates.push_back(rev::Gate(rev::CX, *pf, *(pf + 1)));
                    } else {
                        swap_gates.push_back(rev::Gate(rev::CX, *(pf + 1), *pf));
                        swap_gates.push_back(rev::Gate(rev::H, *(pf + 1)));
                        swap_gates.push_back(rev::Gate(rev::H, *pf));
                        swap_gates.push_back(rev::Gate(rev::CX, *(pf + 1), *pf));
                        swap_gates.push_back(rev::Gate(rev::H, *(pf + 1)));
                        swap_gates.push_back(rev::Gate(rev::H, *pf));
                        swap_gates.push_back(rev::Gate(rev::CX, *(pf + 1), *pf));
                    }
                }
                // Actual swap gate insertion in backward path
                for (std::vector<grph::node_t>::reverse_iterator pb = paths[k].rbegin(); *pb != mincost.target; ++pb) {
                    // updating local to global map
                    lgmap[glmap[plmap[*pb]]] = plmap[*(pb + 1)];
                    lgmap[glmap[plmap[*(pb + 1)]]] = plmap[*pb];

                    // updating global to local map
                    rev::line_t l = glmap[plmap[*pb]];
                    glmap[plmap[*pb]] = glmap[plmap[*(pb + 1)]];
                    glmap[plmap[*(pb + 1)]] = l;

                    if (pg.getDistance(*pb, *(pb + 1)) == 0) {
                        swap_gates.push_back(rev::Gate(rev::CX, *pb, *(pb + 1)));
                        swap_gates.push_back(rev::Gate(rev::H, *pb));
                        swap_gates.push_back(rev::Gate(rev::H, *(pb + 1)));
                        swap_gates.push_back(rev::Gate(rev::CX, *pb, *(pb + 1)));
                        swap_gates.push_back(rev::Gate(rev::H, *pb));
                        swap_gates.push_back(rev::Gate(rev::H, *(pb + 1)));
                        swap_gates.push_back(rev::Gate(rev::CX, *pb, *(pb + 1)));
                    } else {
                        swap_gates.push_back(rev::Gate(rev::CX, *(pb + 1), *pb));
                        swap_gates.push_back(rev::Gate(rev::H, *(pb + 1)));
                        swap_gates.push_back(rev::Gate(rev::H, *pb));
                        swap_gates.push_back(rev::Gate(rev::CX, *(pb + 1), *pb));
                        swap_gates.push_back(rev::Gate(rev::H, *(pb + 1)));
                        swap_gates.push_back(rev::Gate(rev::H, *pb));
                        swap_gates.push_back(rev::Gate(rev::CX, *(pb + 1), *pb));
                    }
                }
                // current control line is not valid for physical coupling map, i.e target->control
                if (pg.getDistance(mincost.control, mincost.target) > 0) {
                    swap_gates.push_back(rev::Gate(rev::H, mincost.control));
                    swap_gates.push_back(rev::Gate(rev::H, mincost.target));
                    swap_gates.push_back(rev::Gate(rev::CX, mincost.target, mincost.control));
                    swap_gates.push_back(rev::Gate(rev::H, mincost.control));
                    swap_gates.push_back(rev::Gate(rev::H, mincost.target));

                } else
                    swap_gates.push_back(rev::Gate(rev::CX, mincost.control, mincost.target));

                gates.insert(g, swap_gates.begin(), swap_gates.end());
                g = gates.begin() + pos + swap_gates.size();
                gates.erase(g);
                g = gates.begin() + pos + swap_gates.size() - 1;
                if (g == gates.end()) break;
            }
        } else {
            g->updateLines(phmap[lgmap[g->getTargets()[0]]]);
        }
    }
    ckt.setGates(gates);
    ckt.setLines(new_lines);
}

// IBM QX
void localorderv2(rev::Circuit& ckt, std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& lg, int base, rev::window_t window) {
    std::vector<rev::Gate> gates = ckt.getGates();
    std::vector<rev::line_t> lines = ckt.getLines();
    std::map<rev::line_t, rev::line_t> measures = ckt.getMeasures();
    unsigned size = gates.size();
    std::map<rev::line_t, rev::line_t> lgmap = std::map<rev::line_t, rev::line_t>();  // local to global line map
    std::map<rev::line_t, rev::line_t> glmap = std::map<rev::line_t, rev::line_t>();  // global to local line map
    std::map<rev::line_t, rev::line_t> plmap = std::map<rev::line_t, rev::line_t>();  // physical to logical line map
    std::vector<rev::line_t> new_lines = std::vector<rev::line_t>();
    std::map<rev::line_t, rev::line_t> new_measures = std::map<rev::line_t, rev::line_t>();
    std::vector<rev::Gate> new_gates = std::vector<rev::Gate>();
    bool isInitSwap;
    do {
        isInitSwap = false;
        for (auto l = lines.begin(); l != lines.end(); ++l) {
            lgmap[*l] = *l;
            glmap[*l] = *l;
            plmap[phmap[*l]] = *l;
            new_lines.push_back(phmap[*l]);
        }
        unsigned gcount = 0;
        unsigned tot_gate = gates.size();
        bool isFirstGate = true;
        for (auto g = gates.begin(); g != gates.end(); ++g) {
            if (g->getType() == rev::CX) {
                gcount++;
                // when gate is properly mapped
                if (pg.getDistance(phmap[lgmap[g->getControls()[0]]], phmap[lgmap[g->getTargets()[0]]]) == 0) {  // adjacent gate
                    rev::Gate ng = *g;
                    ng.updateLines(phmap[lgmap[g->getControls()[0]]], phmap[lgmap[g->getTargets()[0]]]);
                    new_gates.push_back(ng);
                    isFirstGate = false;
                } else {  // qubits are not adjacent

                    std::vector<std::vector<grph::node_t> > paths = std::vector<std::vector<grph::node_t> >();
                    paths.push_back(std::vector<grph::node_t>());
                    getAllPaths(phmap[lgmap[g->getControls()[0]]], phmap[lgmap[g->getTargets()[0]]], paths, pg, 0);

                    unsigned pos = std::distance(gates.begin(), g);
                    grph::Graph gr;
                    if (g + 1 != gates.end()) {
                        gr = grph::Graph(lines, gates, pos + 1, grph::NEW, window, base);
                    }
                    int k = 0;  // path index
                    mincost_t mincost;
                    bool empty = true;
                    for (int i = 0; i < paths.size(); i++) {
                        mincost_t tempcost;
                        if (isInvalidPath(paths[i], plmap) && i != paths.size() - 1) {
                            continue;
                        }
                        if (gr.isEmpty() == false) {
                            tempcost = minCostPathIndexv2(paths[i], lgmap, glmap, plmap, phmap, gates, lines, pg, gr, pos, base);
                        } else {  // reverse traverse prosessed gates based on current control and target
                            bool c_flag = false, t_flag = false;
                            unsigned long path_0_c_count = 0, path_0_t_count = 0;
                            unsigned long path_1_c_count = 0, path_1_t_count = 0;
                            unsigned long path_n_c_count = 0, path_n_t_count = 0;
                            unsigned long path_n1_c_count = 0, path_n1_t_count = 0;
                            unsigned long path_len = paths[i].size() - 1;
                            auto c_gate = new_gates.rend();
                            auto t_gate = new_gates.rend();
                            for (auto ng = new_gates.rbegin(); ng != new_gates.rend(); ng++) {
                                if (ng->getType() == rev::CX) {
                                    if (ng->getControls()[0] == paths[i][0]) path_0_c_count += 1;
                                    if (ng->getTargets()[0] == paths[i][0]) path_0_t_count += 1;
                                    if (ng->getControls()[0] == paths[i][1]) path_1_c_count += 1;
                                    if (ng->getTargets()[0] == paths[i][1]) path_1_t_count += 1;

                                    if (ng->getControls()[0] == paths[i][path_len]) path_n_c_count += 1;
                                    if (ng->getTargets()[0] == paths[i][path_len]) path_n_t_count += 1;
                                    if (ng->getControls()[0] == paths[i][path_len - 1]) path_n1_c_count += 1;
                                    if (ng->getTargets()[0] == paths[i][path_len - 1]) path_n1_t_count += 1;

                                    if ((ng->getControls()[0] == paths[i][0] || ng->getTargets()[0] == paths[i][0]) &&
                                        (ng->getControls()[0] == paths[i][1] || ng->getTargets()[0] == paths[i][1])) {
                                        c_gate = ng;
                                        c_flag = true;
                                    }

                                    if ((ng->getControls()[0] == paths[i][path_len] || ng->getTargets()[0] == paths[i][path_len]) &&
                                        (ng->getControls()[0] == paths[i][path_len - 1] || ng->getTargets()[0] == paths[i][path_len - 1])) {
                                        t_gate = ng;
                                        t_flag = true;
                                    }

                                } else {
                                    if (ng->getTargets()[0] == paths[i][0]) path_0_t_count += 1;
                                    if (ng->getTargets()[0] == paths[i][1]) path_1_t_count += 1;

                                    if (ng->getTargets()[0] == paths[i][path_len]) path_n_t_count += 1;
                                    if (ng->getTargets()[0] == paths[i][path_len - 1]) path_n1_t_count += 1;
                                }
                                if (c_flag && t_flag) break;
                            }
                            tempcost.swaps = path_len - 1;
                            if (c_flag && !t_flag) {
                                if (path_0_c_count == 0 && path_1_c_count == 0)
                                    tempcost.swap_cost = 0;
                                else
                                    tempcost.swap_cost = std::max(path_0_c_count + path_0_t_count, path_1_c_count + path_1_t_count);
                                tempcost.control = paths[i][path_len - 1];
                                tempcost.target = paths[i][path_len];
                            } else if (!c_flag && t_flag) {
                                if (path_n_c_count == 0 && path_n1_c_count == 0)
                                    tempcost.swap_cost = 0;
                                else
                                    tempcost.swap_cost = std::max(path_n_c_count + path_n_t_count, path_n1_c_count + path_n1_t_count);
                                tempcost.control = paths[i][path_len - 1];
                                tempcost.target = paths[i][path_len];
                            } else if (c_flag && t_flag) {
                                double control_cost, target_cost;
                                if (c_gate->getControls()[0] == paths[i][0])
                                    control_cost = std::max(path_0_t_count, path_1_c_count * new_gates.size());

                                if (c_gate->getControls()[0] == paths[i][1])
                                    control_cost = std::max(path_1_t_count, path_0_c_count * new_gates.size());

                                if (t_gate->getControls()[0] == paths[i][path_len])
                                    target_cost = std::max(path_n_t_count, path_n1_c_count * new_gates.size());

                                if (t_gate->getControls()[0] == paths[i][path_len - 1])
                                    target_cost = std::max(path_n1_t_count, path_n_c_count * new_gates.size());

                                if (control_cost < target_cost) {
                                    tempcost.swap_cost = control_cost;
                                    tempcost.control = paths[i][path_len - 1];
                                    tempcost.target = paths[i][path_len];
                                } else {
                                    tempcost.swap_cost = target_cost;
                                    tempcost.control = paths[i][0];
                                    tempcost.target = paths[i][1];
                                }
                            } else {
                                tempcost.swap_cost = new_gates.size();
                                tempcost.control = paths[i][path_len - 1];
                                tempcost.target = paths[i][path_len];
                            }
                        }
                        if (empty || mincost.swap_cost > tempcost.swap_cost || (mincost.swap_cost == tempcost.swap_cost && mincost.swaps > tempcost.swaps)) {
                            mincost.swap_cost = tempcost.swap_cost;
                            mincost.swaps = tempcost.swaps;
                            mincost.control = tempcost.control;
                            mincost.target = tempcost.target;
                            k = i;
                            empty = false;
                        }
                    }

                    std::vector<rev::Gate> swap_gates = std::vector<rev::Gate>();

                    insertSWAPv2(swap_gates, paths[k], lgmap, glmap, plmap, mincost.control, mincost.target);

                    if (isFirstGate == true) {
                        std::map<grph::node_t, grph::node_t> new_phmap;
                        for (auto elem : glmap)
                            new_phmap[elem.second] = phmap[elem.first];
                        phmap = new_phmap;

                        lgmap.clear();
                        glmap.clear();
                        plmap.clear();
                        new_lines.clear();
                        new_measures.clear();
                        new_gates.clear();
                        isInitSwap = true;
                        break;
                    }
                    rev::Gate ng = *g;
                    ng.updateLines(mincost.control, mincost.target);
                    new_gates.insert(new_gates.end(), swap_gates.begin(), swap_gates.end());
                    new_gates.push_back(ng);
                }
            } else {
                rev::Gate ng = *g;
                ng.updateLines(phmap[lgmap[g->getTargets()[0]]]);
                new_gates.push_back(ng);
            }
        }
    } while (isInitSwap == true);
    for (std::map<rev::line_t, rev::line_t>::iterator l = measures.begin(); l != measures.end(); ++l) {
        new_measures[phmap[lgmap[l->first]]] = l->second;
    }
    ckt.setGates(new_gates);
    ckt.setLines(new_lines);
    ckt.setMeasures(new_measures);
}

void insertRCNOT(std::vector<rev::Gate>& rcnot_gates, std::vector<grph::node_t>& path,
                 rev::line_t control, rev::line_t target) {
    bool flag = false;
    for (int i = 0; i < path.size() - 1; i++) {
        if (path[i] == control) flag = true;
        if (flag == true) {
            rcnot_gates.push_back(rev::Gate(rev::CX, path[i], path[i + 1]));
        }
        if (path[i + 1] == target) break;
    }
    for (int i = rcnot_gates.size() - 2; i >= 0; i--) {
        rcnot_gates.push_back(rev::Gate(rev::CX, rcnot_gates[i].getControls()[0], rcnot_gates[i].getTargets()[0]));
    }
    for (int i = rcnot_gates.size() - 2; i >= 1; i--) {
        rcnot_gates.push_back(rev::Gate(rev::CX, rcnot_gates[i].getControls()[0], rcnot_gates[i].getTargets()[0]));
    }
}

// insertion of Swap as a cascade of CNOT gates
void insertSWAP(std::vector<rev::Gate>& swap_gates, std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& lgmap,
                std::map<rev::line_t, rev::line_t>& glmap, std::map<rev::line_t, rev::line_t>& plmap,
                rev::line_t control, rev::line_t target) {
    for (int i = 0; path[i] != control; i++) {
        // updating local to global map
        lgmap[glmap[plmap[path[i]]]] = plmap[path[i + 1]];
        lgmap[glmap[plmap[path[i + 1]]]] = plmap[path[i]];

        // updating global to local map
        rev::line_t l = glmap[plmap[path[i]]];
        glmap[plmap[path[i]]] = glmap[plmap[path[i + 1]]];
        glmap[plmap[path[i + 1]]] = l;

        swap_gates.push_back(rev::Gate(rev::CX, path[i], path[i + 1]));
        swap_gates.push_back(rev::Gate(rev::CX, path[i + 1], path[i]));
        swap_gates.push_back(rev::Gate(rev::CX, path[i], path[i + 1]));
    }
    // Actual swap gate insertion in backward path
    for (int j = path.size() - 1; path[j] != target; j--) {
        // updating local to global map
        lgmap[glmap[plmap[path[j]]]] = plmap[path[j - 1]];
        lgmap[glmap[plmap[path[j - 1]]]] = plmap[path[j]];

        // updating global to local map
        rev::line_t l = glmap[plmap[path[j]]];
        glmap[plmap[path[j]]] = glmap[plmap[path[j - 1]]];
        glmap[plmap[path[j - 1]]] = l;

        swap_gates.push_back(rev::Gate(rev::CX, path[j], path[j - 1]));
        swap_gates.push_back(rev::Gate(rev::CX, path[j - 1], path[j]));
        swap_gates.push_back(rev::Gate(rev::CX, path[j], path[j - 1]));
    }
}

// insertion of Swap as an abstract gate
void insertSWAPv2(std::vector<rev::Gate>& swap_gates, std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& lgmap,
                  std::map<rev::line_t, rev::line_t>& glmap, std::map<rev::line_t, rev::line_t>& plmap,
                  rev::line_t control, rev::line_t target) {
    for (int i = 0; path[i] != control; i++) {
        // updating local to global map
        lgmap[glmap[plmap[path[i]]]] = plmap[path[i + 1]];
        lgmap[glmap[plmap[path[i + 1]]]] = plmap[path[i]];

        // updating global to local map
        rev::line_t l = glmap[plmap[path[i]]];
        glmap[plmap[path[i]]] = glmap[plmap[path[i + 1]]];
        glmap[plmap[path[i + 1]]] = l;

        swap_gates.push_back(rev::Gate(rev::SWAP, path[i], path[i + 1]));
    }

    // Actual swap gate insertion in backward path
    for (int j = path.size() - 1; path[j] != target; j--) {
        // updating local to global map
        lgmap[glmap[plmap[path[j]]]] = plmap[path[j - 1]];
        lgmap[glmap[plmap[path[j - 1]]]] = plmap[path[j]];

        // updating global to local map
        rev::line_t l = glmap[plmap[path[j]]];
        glmap[plmap[path[j]]] = glmap[plmap[path[j - 1]]];
        glmap[plmap[path[j - 1]]] = l;

        swap_gates.push_back(rev::Gate(rev::SWAP, path[j], path[j - 1]));
    }
}