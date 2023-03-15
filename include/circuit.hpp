/*--------------------------------------------------------------------------------
Copyright (C) 2023, Abhoy Kole, Kamalika Datta, Indranil Sengupta, Rolf Drechsler

Quantum Mapping Tool (QMT) for mapping Quantum circuits to IBM architecture
released under the MIT licence.

All right reserved.
----------------------------------------------------------------------------------*/

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

#include "gate.hpp"

#ifndef CIRCUIT_HPP
#define CIRCUIT_HPP
using namespace std;
namespace rev {
typedef unsigned window_t;

int convertToInteger(string);
void exchange(vector<line_t>& lines, int i, int j);
class Circuit {
   private:
    int numvars;

    vector<string> variables;      // backward compitable
    vector<line_t> lines;          // quantum registers
    vector<line_t> clines;         // classical registers
    map<line_t, line_t> measures;  // quantum to classical register mapping for measurement

    vector<Gate> gates;
    map<line_t, string> var_names;  // qubit names
    map<string, line_t> var_indices;

    map<line_t, string> var_cnames;  // cbit names
    map<string, line_t> var_cindices;

    map<line_t, double> CNOTDepth;

   public:
    Circuit();
    Circuit(const Circuit& ckt);
    ~Circuit();
    void init();
    void readQASM(std::istream& myfile);
    void readReal(string);
    void writeQASM(string, unsigned);
    void writeReal(string file);
    void displayQASM(ostream& os, unsigned size);
    void displayQisKit(ostream& os, unsigned size);
    vector<Gate> getGates();
    vector<line_t> getLines();
    map<line_t, line_t> getMeasures();
    void setGates(vector<Gate>& gates);
    void setLines(vector<line_t>& lines);
    void setMeasures(map<line_t, line_t>& measures);
    unsigned twoQubitGateCount();
    unsigned getCNOTDepth();
    unsigned maxCNOTDepth();
    void reduceCNOTDepth(Gate& g);
    unsigned twoQubitGateCount(unsigned pos);
    void addQubit(string);
    void addCbit(string);
    void addGate(Gate&);
    void addGate(Gate);
    void clear();
    void operator=(const Circuit& ckt);
    friend void exchange(vector<line_t>& lines, int i, int j);
    friend int convertToInteger(string);
};
}  // namespace rev
#endif /* CIRCUIT_HPP */
