/*--------------------------------------------------------------------------------
Copyright (C) 2023, Abhoy Kole, Kamalika Datta, Indranil Sengupta, Rolf Drechsler

Quantum Mapping Tool (QMT) for mapping Quantum circuits to IBM architecture
released under the MIT licence.

All right reserved.
----------------------------------------------------------------------------------*/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdlib>
#include <nlohmann/json.hpp>
#include <pybind11_json/pybind11_json.hpp>

#include "circuit.hpp"
#include "costmetrics.hpp"
#include "gate.hpp"
#include "graph.hpp"
#include "initialmapping.hpp"
#include "localmapping.hpp"

#define BASE 2

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;
namespace nl = nlohmann;
using namespace pybind11::literals;


// Read both layout and circuit from file and write to file
void qx_mapper_v1(std::string& layout, std::string& circuit, std::string& out_circuit, int runs = 1) {
    rev::Circuit ckt = rev::Circuit();
    grph::Graph pg = grph::Graph();
    
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
        std::cout << "Unable to map, imput circuit is too large (qubits=" << ckt.getLines().size() << ") for the current architecture" << std::endl;
        exit(EXIT_FAILURE);
    }

    rev::window_t window = 0;
    window = ckt.getCNOTDepth();
    if (window > 20) window = 20;

    rev::Circuit min_ckt;
    std::map<grph::node_t, grph::node_t> phmap_min = std::map<grph::node_t, grph::node_t>();

    for (int k = 0; k < runs; k++) {
        rev::Circuit tmp_ckt = ckt;
        grph::Graph lg = grph::Graph(tmp_ckt, grph::NEW, window, BASE);

        std::map<grph::node_t, grph::node_t> phmap = std::map<grph::node_t, grph::node_t>();

        lg.mapLayout2(pg, phmap);
        for (auto elem : phmap) {
            phmap_min[elem.first] = elem.second;
        }

        ga_search(phmap_min, pg, lg, BASE, window);
        localorderv2(tmp_ckt, phmap_min, pg, lg, BASE, window);
        if (k == 0 || tmp_ckt.getGates().size() < min_ckt.getGates().size())
            min_ckt = tmp_ckt;
    }
    auto ofs = std::ofstream(out_circuit);

    if (!ofs.good()) {
        std::cout << "[writeCircuit] unable to open file " + out_circuit << std::endl;
        exit(EXIT_FAILURE);
    }

    min_ckt.displayQASM(ofs, pg.getNodes().size());

    ofs.close();
}

// Read both layout and circuit as parameters and return circuit as qasm string
nl::json qx_mapper_v2(const nl::json& layout, const nl::json& circuit, int runs = 1) {
    grph::Graph pg = grph::Graph();
    auto is1 = std::istringstream(layout.get<std::string>());
    pg.readLayout2(is1);

    rev::Circuit ckt = rev::Circuit();
    auto is2 = std::istringstream(circuit.get<std::string>());
    ckt.readQASM(is2);

    if (ckt.getLines().size() > pg.getNodes().size()) {
        std::cout << "Unable to map, imput circuit is too large (qubits=" << ckt.getLines().size() << ") for the current architecture" << std::endl;
        exit(EXIT_FAILURE);
    }
    rev::window_t window = 0;
    window = ckt.getCNOTDepth();
    if (window > 20) window = 20;

    rev::Circuit min_ckt;
    std::map<grph::node_t, grph::node_t> phmap_min = std::map<grph::node_t, grph::node_t>();

    for (int k = 0; k < runs; k++) {
        rev::Circuit tmp_ckt = ckt;
        grph::Graph lg = grph::Graph(tmp_ckt, grph::NEW, window, BASE);

        std::map<grph::node_t, grph::node_t> phmap = std::map<grph::node_t, grph::node_t>();

        lg.mapLayout2(pg, phmap);
        for (auto elem : phmap) {
            phmap_min[elem.first] = elem.second;
        }

        ga_search(phmap_min, pg, lg, BASE, window);

        localorderv2(tmp_ckt, phmap_min, pg, lg, BASE, window);

        if (k == 0 || tmp_ckt.getGates().size() < min_ckt.getGates().size())
            min_ckt = tmp_ckt;
    }

    auto os = std::ostringstream();
    min_ckt.displayQASM(os, pg.getNodes().size());

    nl::json out_circuit = os.str();
    return out_circuit;
}

PYBIND11_MODULE(PyQMapping, m) {
    m.doc() = R"pbdoc(
        Python interface for the QMapping Layout Mapping library
        ---------------------------------------------------
        .. autosummary::
           :toctree: _generate
           qx_mapper_v1
           qx_mapper_v2
    )pbdoc";

    m.def("qx_mapper_file", &qx_mapper_v1, "map a quantum circuit to a lyout by reading circuit and layout from files and writing circuit to file");
    m.def("qx_mapper", &qx_mapper_v2, "map a quantum circuit to a lyout by reading circuit and layout from string and returning circuit as string");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}