/*--------------------------------------------------------------------------------
Copyright (C) 2023, Abhoy Kole, Kamalika Datta, Indranil Sengupta, Rolf Drechsler

Quantum Mapping Tool (QMT) for mapping Quantum circuits to IBM architecture
released under the MIT licence.

All right reserved.
----------------------------------------------------------------------------------*/

#include "initialmapping.hpp"

#include <cmath>
#include <cstdlib>

// Exhaustive search for a minimal cost global mapping of logical qubits to physical qubits.
double ex_search(std::map<grph::node_t, grph::node_t>& phmap, std::map<grph::node_t, grph::node_t>& phmap_min,
                 grph::Graph& pg, grph::Graph& lg, double mcost, int pos) {
    if (pos == phmap.size() - 1) {
        double cost = mappingCost(phmap, pg, lg);
        if (cost < mcost) {
            mcost = cost;
            phmap_min.clear();

            for (std::map<grph::node_t, grph::node_t>::iterator l = phmap.begin(); l != phmap.end(); ++l)
                phmap_min[l->first] = l->second;
        }
        return mcost;
    } else {
        std::map<grph::node_t, grph::node_t>::iterator i = phmap.begin();
        for (std::advance(i, pos); i != phmap.end(); ++i) {
            grph::node_t tmp = phmap[pos];
            phmap[pos] = i->second;
            phmap[i->first] = tmp;
            mcost = ex_search(phmap, phmap_min, pg, lg, mcost, pos + 1);
            tmp = phmap[pos];
            phmap[pos] = i->second;
            phmap[i->first] = tmp;
        }
    }
    return mcost;
}

// Generate the initial population, each solution being of size equal to the size of any element from phmaps.
void gen_initial_pop(std::vector<std::map<grph::node_t, grph::node_t> >& phmaps, std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& lg) {
    for (std::vector<std::map<grph::node_t, grph::node_t> >::iterator i = phmaps.begin(); i != phmaps.end(); ++i) {
        physicalMap(*i, phmap, lg);
    }
}

// Mapping logical qubits to physial qubits
void physicalMap(std::map<grph::node_t, grph::node_t>& phmap_new, std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& lg) {
    int r;
    int x = lg.getNodes().size();
    for (std::map<grph::node_t, grph::node_t>::iterator n = phmap.begin(); n != phmap.end(); ++n) {
        do {
            r = rand() % x;
        } while (phmap_new.find(r) != phmap_new.end());
        phmap_new[r] = n->second;
    }
}

// IBM QX mapping: Calculate the fitness values of the solutions in "phmaps"
void calculate_fitness(std::vector<std::map<grph::node_t, grph::node_t> >& phmaps, std::vector<double>& fitness,
                       grph::Graph& pg, grph::Graph& lg, int base, unsigned depth) {
    fitness.clear();
    for (auto& m : phmaps) {
        fitness.push_back(remoteCNOTCost(m, pg, lg));
    }
}

// Calculate percentage fitness of the solutions
void calculate_percent_fitness(std::vector<double>& fitness, std::vector<double>& percent_fitness) {
    double total = 0;
    for (std::vector<double>::iterator i = fitness.begin(); i != fitness.end(); ++i) {
        total += *i;
    }
    percent_fitness.clear();
    for (std::vector<double>::iterator i = fitness.begin(); i != fitness.end(); ++i) {
        if (total == 0.)
            percent_fitness.push_back((*i));
        else
            percent_fitness.push_back((*i) / total);
    }
}

// Sort the solutions in "phmaps" in ascending order of fitness values
void sort_gener_cur(std::vector<std::map<grph::node_t, grph::node_t> >& phmaps, std::vector<double>& fitness) {
    for (std::vector<double>::iterator i = fitness.begin(); i != fitness.end() - 1; ++i) {
        for (std::vector<double>::iterator j = i + 1; j != fitness.end(); ++j) {
            if (*i > *j) {
                double t = *i;
                *i = *j;
                *j = t;
                for (std::map<grph::node_t, grph::node_t>::iterator m = phmaps[std::distance(fitness.begin(), i)].begin(),
                                                                    n = phmaps[std::distance(fitness.begin(), j)].begin();
                     m != phmaps[std::distance(fitness.begin(), i)].end(); ++m, ++n) {
                    grph::node_t tmp = m->second;
                    m->second = n->second;
                    n->second = tmp;
                }
            }
        }
    }
}

// Select a solution from "phmaps" using roulette wheel method
int selection(std::vector<double>& percent_fitness) {
    int i;
    long double r, cumul_sum;

    r = rand() / (double)RAND_MAX;
    cumul_sum = 0.0;

    for (std::vector<double>::iterator i = percent_fitness.begin(); i != percent_fitness.end(); ++i) {
        cumul_sum += *i;
        if (r < cumul_sum) return std::distance(percent_fitness.begin(), i);
    }
    return 0;
}

/*Perform crossover of two solutions of "phmaps" at indices "pos1" and "pos2", and
store the new solutions in "phmaps_next" at the end.*/
void crossover_xy(std::vector<std::map<grph::node_t, grph::node_t> >& phmaps,
                  std::vector<std::map<grph::node_t, grph::node_t> >& phmaps_next, int pos1, int pos2, int x) {
    int crosspoint;
    do {
        crosspoint = rand() % (x - 1);
    } while (x == crosspoint + 2);  // Crossover requires at least two items
    std::map<grph::node_t, grph::node_t> phmap1, phmap2;

    for (auto i = phmaps[pos1].begin(), j = phmaps[pos2].begin(); i != phmaps[pos1].end() &&
                                                                  std::distance(phmaps[pos1].begin(), i) <= crosspoint;
         ++i, ++j) {  // Copy all elements upto the crossover point
        phmap1[i->first] = i->second;
        phmap2[j->first] = j->second;
    }
    int index1, index2;
    std::map<int, bool> snodes1, snodes2;
    std::map<grph::node_t, bool> crossover1, crossover2;
    int items = x - crosspoint - 1;   // number of items to be cross over
    for (auto i = 0; i < items; i++)  // initially no crossover
        snodes1[i] = false, snodes2[i] = false;
    while (items > 0) {  // perform crossover
        do {
            index1 = rand() % (x - crosspoint - 1);  // select a mapping that is not crossover
        } while (snodes1[index1]);
        do {
            index2 = rand() % (x - crosspoint - 1);  // select another mapping that is not crossover
        } while (snodes1[index2] || index1 == index2);
        snodes1[index1] = true, snodes1[index2] = true;
        auto i = phmaps[pos1].begin(), j = phmaps[pos1].begin();
        std::advance(i, crosspoint + 1 + index1), std::advance(j, crosspoint + 1 + index2);
        phmap1[i->first] = j->second, phmap1[j->first] = i->second;
        crossover1[i->second] = true, crossover1[j->second] = true;
        do {
            index1 = rand() % (x - crosspoint - 1);  // select a mapping that is not crossover
        } while (snodes2[index1]);
        do {
            index2 = rand() % (x - crosspoint - 1);  // select another mapping that is not crossover
        } while (snodes2[index2] || index1 == index2);
        snodes2[index1] = true, snodes2[index2] = true;
        i = phmaps[pos2].begin(), j = phmaps[pos2].begin();
        std::advance(i, crosspoint + 1 + index1), std::advance(j, crosspoint + 1 + index2);
        phmap2[i->first] = j->second, phmap2[j->first] = i->second;
        crossover2[i->second] = true, crossover2[j->second] = true;
        items -= 2;        // two items crossovered
        if (items == 1) {  // can not perform crossover
            i = phmaps[pos1].begin(), j = phmaps[pos2].begin();
            for (std::advance(i, crosspoint + 1); i != phmaps[pos1].end(); i++) {  // copy items that are not crossover
                if (crossover1.find(i->second) == crossover1.end()) {
                    phmap1[i->first] = i->second;  // copy the item
                    break;
                }
            }
            for (std::advance(j, crosspoint + 1); j != phmaps[pos2].end(); j++) {
                if (crossover2.find(j->second) == crossover2.end()) {
                    phmap2[j->first] = j->second;
                    break;
                }
            }
            break;
        }
    }

    phmaps_next.push_back(phmap1);
    phmaps_next.push_back(phmap2);
}

/*Generate next generation "phmaps_next" from the current generation "phmaps"
Copy the top NBEST solutions as it is; rest are generated through crossover using probability CROSSPROB */
void generation_next(std::vector<std::map<grph::node_t, grph::node_t> >& phmaps,
                     std::vector<std::map<grph::node_t, grph::node_t> >& phmaps_next,
                     std::vector<double>& percent_fitness, int x) {
    phmaps_next.clear();
    for (std::vector<std::map<grph::node_t, grph::node_t> >::iterator i = phmaps.begin(); i != phmaps.end() &&
                                                                                          std::distance(phmaps.begin(), i) < NBEST;
         ++i) {
        std::map<grph::node_t, grph::node_t> phmap;
        for (std::map<grph::node_t, grph::node_t>::iterator j = (*i).begin(); j != (*i).end(); ++j) phmap[j->first] = j->second;
        phmaps_next.push_back(phmap);
    }
    int num_crossover = (POPSIZE - NBEST) / 2;  // Number of crossovers to perform POPSIZE & NBIST should be even

    for (int k = 0; k < num_crossover; k++) {
        int pos1 = selection(percent_fitness);
        int pos2 = selection(percent_fitness);
        if ((rand() / (float)RAND_MAX) < CROSSPROB) {
            crossover_xy(phmaps, phmaps_next, pos1, pos2, x);  // Perform crossover with probability
        } else {
            std::map<grph::node_t, grph::node_t> phmap1, phmap2;
            for (std::map<grph::node_t, grph::node_t>::iterator i = phmaps[pos1].begin(), j = phmaps[pos2].begin(); i != phmaps[pos1].end(); ++i, ++j) {
                phmap1[i->first] = i->second;
                phmap2[j->first] = j->second;
            }
            phmaps_next.push_back(phmap1);  // Else copy as it is
            phmaps_next.push_back(phmap2);
        }
    }
}

// Perform mutation in the solution at "index" in "gener_1d_next"
void mutation(std::vector<std::map<grph::node_t, grph::node_t> >& phmaps_next, int index, int x) {
    grph::node_t first, second;
    int type = rand() % 6;  // a number between 0 and 5
    switch (type) {
        case 0:  // swap two adjacent qubits along x direction
        case 1: {
            first = rand() % x;
            if (first + 1 < x)
                second = first + 1;
            else {
                second = first;
                first = second - 1;
            }
            exchange(phmaps_next[index], first, second);
        } break;

        case 2:
        case 3: {  // rotate solution right
            auto n = phmaps_next[index].begin();
            auto m = n;
            auto tmp = m->second;
            for (++m; m != phmaps_next[index].end(); ++n, ++m)
                n->second = m->second;
            n->second = tmp;
        } break;

        case 4: {  // rotate solution left
            auto n = phmaps_next[index].rbegin();
            auto m = n;
            auto tmp = m->second;
            for (++m; m != phmaps_next[index].rend(); ++n, ++m)
                n->second = m->second;
            n->second = tmp;
        } break;

        case 5: {  // swap two arbitrarily randomly chosen entries
            do {
                first = rand() % x;
                second = rand() % x;
            } while (first == second);
            exchange(phmaps_next[index], first, second);
        } break;
    }
}

// Exange two qubits in the solution of "phmap" at positions "a" and "b".
void exchange(std::map<grph::node_t, grph::node_t>& phmap, grph::node_t a, grph::node_t b) {
    grph::node_t n = phmap[a];
    phmap[a] = phmap[b];
    phmap[b] = n;
}

// Perform mutation with probability MUTPROB, and copy back  "phmaps_next" into "phmaps"
void mutation_and_copy(std::vector<std::map<grph::node_t, grph::node_t> >& phmaps, std::vector<std::map<grph::node_t, grph::node_t> >& phmaps_next, int x) {
    for (int i = NBEST; i < POPSIZE; i++) {
        if ((rand() / (float)RAND_MAX) < MUTPROB)
            mutation(phmaps_next, i, x);  // Mutate i-th solution in "gener_nxt"

        for (std::map<grph::node_t, grph::node_t>::iterator j = phmaps_next[i].begin(); j != phmaps_next[i].end(); ++j)
            phmaps[i][j->first] = j->second;
    }
}

// IBM QX mapping : GA based ordering of the qubits so as to optimize the cost
void ga_search(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& lg, int base, unsigned depth) {
    std::vector<std::map<grph::node_t, grph::node_t> > phmaps(POPSIZE, std::map<grph::node_t, grph::node_t>());
    std::vector<std::map<grph::node_t, grph::node_t> > phmaps_next(POPSIZE, std::map<grph::node_t, grph::node_t>());
    std::vector<double> fitness(POPSIZE, 0);
    std::vector<double> percent_fitness(POPSIZE, 0);
    srand(time(NULL));
    gen_initial_pop(phmaps, phmap, lg);

    for (int numgen = 1; numgen <= MAXGEN; numgen++) {  // MAXGEN
        calculate_fitness(phmaps, fitness, pg, lg, base, depth);
        sort_gener_cur(phmaps, fitness);
        if (fitness[0] == 0) break;  // minimal solution
        calculate_percent_fitness(fitness, percent_fitness);
        generation_next(phmaps, phmaps_next, percent_fitness, lg.getNodes().size());
        mutation_and_copy(phmaps, phmaps_next, lg.getNodes().size());
    }

    phmap.clear();
    if (fitness[0] != 0) {
        calculate_fitness(phmaps, fitness, pg, lg, base, depth);
        sort_gener_cur(phmaps, fitness);
    }
    for (std::map<grph::node_t, grph::node_t>::iterator i = phmaps[0].begin(); i != phmaps[0].end(); ++i)
        phmap[i->first] = i->second;
}

// Display the population and their corresponding cost.
void display_all_pop(std::vector<std::map<grph::node_t, grph::node_t> >& phmaps, grph::Graph& pg,
                     grph::Graph& lg, Metric_Type metric) {
    for (std::vector<std::map<grph::node_t, grph::node_t> >::iterator i = phmaps.begin(); i != phmaps.end(); ++i) {
        displayMapping(*i, lg);
        if (metric == Metric_Type::COUPLING_COST)
            std::cout << "Cost: " << mappingCost(*i, pg, lg) << std::endl;
        else
            std::cout << "Cost: " << remoteCNOTCost(*i, pg, lg) << std::endl;
    }
}

// Display placement details in 1D grid.
void displayMapping(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& lg) {
    std::cout << "#Qubit mapping details: (logical) = Physical" << std::endl;
    for (std::vector<grph::node_t>::iterator i = lg.getNodes().begin(); i != lg.getNodes().end(); ++i)
        std::cout << (*i) << "=" << phmap[*i] << std::endl;
}
