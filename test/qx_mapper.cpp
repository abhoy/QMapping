/*--------------------------------------------------------------------------------
Copyright (C) 2023, Abhoy Kole, Kamalika Datta, Indranil Sengupta, Rolf Drechsler

Quantum Mapping Tool (QMT) for mapping Quantum circuits to IBM architecture
released under the MIT licence.

All right reserved.
----------------------------------------------------------------------------------*/

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>

#include "circuit.hpp"
#include "costmetrics.hpp"
#include "gate.hpp"
#include "graph.hpp"
#include "initialmapping.hpp"
#include "localmapping.hpp"

#define BASE 2  // edge weight of IG is power(BASE, depth - pos);


int main(int argc, char** argv) {
    std::string layout;
    std::string circuit;
    std::string out_circuit;
    int runs = 1;
    if (argc != 4 && argc != 5) {
        std::cout << "Usage:\n\tlayoutname[*.lut] \t circuitname[*.qasm] \t outcircuitname[*.qasm] \t runs " << std::endl;
        std::exit(EXIT_FAILURE);
    } else {
        layout = argv[1];
        circuit = argv[2];
        out_circuit = argv[3];
        if (argc == 5) {
            runs = std::stoi(argv[4]);
        }
    }

    rev::Circuit ckt = rev::Circuit();

    std::size_t dot = circuit.find_last_of('.');

    std::string extension = circuit.substr(dot + 1);

    std::transform(extension.cbegin(), extension.cend(), extension.begin(), [](unsigned char ch) { return ::tolower(ch); });

    if (extension == "qasm") {
        auto ifs = std::ifstream(circuit);

        if (!ifs.good()) {
            std::cout << "[readCircuit] unable to open file " + circuit << std::endl;
            exit(EXIT_FAILURE);
        }

        ckt.readQASM(ifs);

        ifs.close();

    } else {
        std::cout << "[readCircuit] extension " + extension + " not recognized" << std::endl;
        exit(EXIT_FAILURE);
    }
    grph::Graph pg = grph::Graph();

    dot = layout.find_last_of('.');

    extension = layout.substr(dot + 1);

    std::transform(extension.cbegin(), extension.cend(), extension.begin(), [](unsigned char ch) { return ::tolower(ch); });

    if (extension == "lut") {
        auto ifs = std::ifstream(layout);

        if (!ifs.good()) {
            std::cout << "[readLayout] unable to open file " + layout << std::endl;
            exit(EXIT_FAILURE);
        }

        pg.readLayout2(ifs);

        ifs.close();

    } else {
        std::cout << "[readLayout] extension " + extension + " not recognized" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (ckt.getLines().size() > pg.getNodes().size()) {
        std::cout << "Unable to map " << circuit.substr(circuit.find_last_of("/") + 1, circuit.find_last_of(".") - circuit.find_last_of("/") - 1) << " to large (qubits=" << ckt.getLines().size() << ") for the current architecture" << std::endl;
        exit(EXIT_FAILURE);
    }
    unsigned CNOTs = ckt.twoQubitGateCount();

    std::cout << "#Circuit:\t" << circuit.substr(circuit.find_last_of("/") + 1, circuit.find_last_of(".") - circuit.find_last_of("/") - 1) << "\t#Lines:\t" << ckt.getLines().size() << "\t#Gates:\t" << ckt.getGates().size() << "\t#CNOTs:\t" << CNOTs;

    rev::window_t window = 0;

    window = ckt.getCNOTDepth();

    if (window > 20) window = 20;

    clock_t start, end;
    double avgtime = 0.;
    rev::Circuit min_ckt{};

    std::map<grph::node_t, grph::node_t> phmap_min = std::map<grph::node_t, grph::node_t>();

    for (int k = 0; k < runs; k++) {
        rev::Circuit tmp_ckt = ckt;
        start = clock();
        grph::Graph lg = grph::Graph(tmp_ckt, grph::NEW, window, BASE);

        std::map<grph::node_t, grph::node_t> phmap = std::map<grph::node_t, grph::node_t>();
        lg.mapLayout2(pg, phmap);
        for (std::map<grph::node_t, grph::node_t>::iterator l = phmap.begin(); l != phmap.end(); ++l)
            phmap_min[l->first] = l->second;

        ga_search(phmap_min, pg, lg, BASE, window);

        localorderv2(tmp_ckt, phmap_min, pg, lg, BASE, window);

        end = clock();
        if (k == 0 || tmp_ckt.getGates().size() < min_ckt.getGates().size()) {
            min_ckt = tmp_ckt;
        }
        avgtime += ((end - start) / (double)CLOCKS_PER_SEC);
    }
    avgtime = avgtime / runs;

    auto ofs = std::ofstream(out_circuit);

    if (!ofs.good()) {
        std::cout << "[writeCircuit] unable to open file " + out_circuit << std::endl;
        exit(EXIT_FAILURE);
    }

    min_ckt.displayQASM(ofs, pg.getNodes().size());

    unsigned SWAPs = min_ckt.twoQubitGateCount() - CNOTs;

    std::cout << "\t#SWAPs:\t" << SWAPs << "\t#Depth:\t" << computeDepth(min_ckt) << "\t#Time:\t" << avgtime << std::endl;
}
