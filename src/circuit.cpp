/*--------------------------------------------------------------------------------
Copyright (C) 2023, Abhoy Kole, Kamalika Datta, Indranil Sengupta, Rolf Drechsler

Quantum Mapping Tool (QMT) for mapping Quantum circuits to IBM architecture
released under the MIT licence.

All right reserved.
----------------------------------------------------------------------------------*/

#include "circuit.hpp"

using namespace rev;

/* Interchange the position of i-th qubit with j-th qubit*/
void rev::exchange(vector<line_t>& lines, int i, int j) {
    line_t line = lines[i];
    lines[i] = lines[j];
    lines[j] = line;
}
/* Display the Quantum circuit in QASM format*/
void Circuit::displayQASM(ostream& os, unsigned size) {
    os << "OPENQASM 2.0;" << std::endl;
    os << "include \"qelib1.inc\";" << std::endl;
    os << "qreg q[" << size << "];" << std::endl;
    os << "creg meas[" << size << "];" << std::endl;
    for (vector<Gate>::iterator i = this->gates.begin(); i != this->gates.end(); ++i)
        (*i).display(os);
    os << "barrier q[";
    for (std::vector<rev::line_t>::iterator l = this->lines.begin(); l != this->lines.end(); ++l)
        if ((l + 1) != this->lines.end())
            os << *l << "],q[";
        else
            os << *l << "];" << std::endl;
    for (std::map<rev::line_t, rev::line_t>::iterator l = this->measures.begin(); l != this->measures.end(); ++l) {
        os << "measure q[" << l->first << "] -> meas[" << l->second << "];" << std::endl;
    }
}
/* Display the Quantum circuit in IBM QisKit format*/
void Circuit::displayQisKit(ostream& os, unsigned size) {
    os << "#QISKIT Format" << std::endl;
    os << "from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, Aer" << std::endl;
    os << "from qiskit.compiler import transpile" << std::endl;
    os << "qc = QuantumCircuit(" << size << ", " << size << ")" << std::endl;
    for (vector<Gate>::iterator i = this->gates.begin(); i != this->gates.end(); ++i)
        (*i).display2(os);
    os << "qc.qasm(formatted=True)" << std::endl;
}
void Circuit::readQASM(std::istream& myfile) {
    string line;
    bool readGate = false;
    if (myfile.good()) {
        while (getline(myfile, line)) {
            if (0 == line.size() || 1 == line.size()) continue;  // skip empty lines
            if ('#' == line[0]) continue;                        // skip comments
            istringstream iss(line);
            string command;
            iss >> command;
            if ("OPENQASM" == command)
                continue;  // Version
            else if ("include" == command)
                continue;  // header file
            else if ("qreg" == command)
                continue;
            else if ("creg" == command)
                continue;  // skip classical registers
            else if ("barrier" == command)
                continue;                     // skip
            else if ("measure" == command) {  // continue;      //skip
                iss >> command;
                string qubit = command;
                addQubit(qubit);
                iss >> command;  // skip string ->
                iss >> command;
                string cbit = command;
                addCbit(cbit);
                this->measures[this->var_indices[qubit]] = this->var_cindices[cbit];

            } else {
                iss >> command;
                while (iss) {
                    auto pos1 = string::npos, pos2 = command.find(",");
                    // int i =0;
                    for (; pos2 != string::npos;) {
                        string qubit = pos1 != string::npos ? command.substr(pos1 + 1, pos2 - pos1 - 1) : command.substr(0, pos2);
                        addQubit(qubit);
                        pos1 = pos2, pos2 = command.find(",", pos1 + 1);
                        if (pos2 == string::npos) {
                            if ((pos2 = command.find(";")) != string::npos) break;
                            iss >> command;
                            pos1 = string::npos, pos2 = command.find(",");
                            if (pos2 == string::npos && (pos2 = command.find(";")) != string::npos) break;
                        }
                    }
                    if (pos2 == string::npos) pos2 = command.find(";");
                    if (pos2 != string::npos) {
                        string qubit = pos1 != string::npos ? command.substr(pos1 + 1, pos2 - pos1 - 1) : command.substr(0, pos2);

                        addQubit(qubit);
                    }
                    iss >> command;
                }

                this->gates.push_back(Gate(line, this->var_indices));
            }
        }
    } else
        cout << "File can not be opened." << endl;
}

void Circuit::addGate(Gate& g) {
    this->gates.push_back(g);
}

void Circuit::addGate(Gate g) {
    this->gates.push_back(g);
}

void Circuit::readReal(string file) {
    string line;
    ifstream myfile(file.c_str());
    bool readGate = false;
    if (myfile.is_open()) {
        int qubits = 0;
        while (getline(myfile, line)) {
            if (0 == line.size() || 1 == line.size()) continue;  // skip empty lines
            if ('#' == line[0]) continue;                        // skip comments
            istringstream iss(line);
            string command;
            iss >> command;
            if (command == ".version")
                continue;  // Version
            else if (command == ".mode")
                continue;
            else if (command == ".numvars") {
                iss >> command;
                qubits = std::stoi(command);
                continue;
            } else if (command == ".variables") {
                while (iss >> command) addQubit(command);
                continue;
            } else if (command == ".inputs")
                continue;
            else if (command == ".outputs")
                continue;
            else if (command == ".constants")
                continue;
            else if (command == ".garbage")
                continue;
            else if (command == ".inputbus")
                continue;
            else if (command == ".outputbus")
                continue;
            else if (command == ".state")
                continue;
            else if (command == ".module")
                continue;
            else if (command == ".begin")
                continue;
            else if (command == ".end")
                continue;
            else if (command == ".define") {
                do {
                    getline(myfile, line);
                    istringstream iss(line);
                    iss >> command;
                } while (command != ".enddefine");
                continue;
            } else {
                if (command[0] == 't') {  // Toffoli Gate
                    int count = std::stoi(command.substr(1));
                    if (count == 1) {
                        iss >> command;
                        this->gates.push_back(Gate(rev::X, this->var_indices[command]));
                    } else if (count == 2) {
                        string control, target;
                        iss >> control;
                        iss >> target;
                        this->gates.push_back(Gate(rev::CX, this->var_indices[control], this->var_indices[target]));
                    } else {
                        std::vector<line_t> controls;
                        string control, target;
                        for (auto i = 0; i < count - 1; i++) {
                            iss >> control;
                            controls.push_back(this->var_indices[control]);
                        }
                        iss >> target;
                        rev::gate_t type;
                        if (count == 3)
                            type = rev::CCX;
                        else
                            type = rev::MCX;
                        this->gates.push_back(Gate(type, controls, this->var_indices[target]));
                    }
                } else if (command[0] == 'f') {  // Fredkin Gate
                    int count = std::stoi(command.substr(1));
                    if (count == 2) {
                        string target1, target2;
                        iss >> target1;
                        iss >> target2;
                        this->gates.push_back(Gate(rev::SWAP, this->var_indices[target1], this->var_indices[target2]));
                    } else {
                        std::vector<line_t> controls;
                        string control, target;
                        for (auto i = 0; i < count - 1; i++) {
                            iss >> control;
                            controls.push_back(this->var_indices[control]);
                        }
                        iss >> target;
                        rev::gate_t type;
                        if (count == 3)
                            type = rev::CCX;
                        else
                            type = rev::MCX;
                        this->gates.push_back(Gate(rev::CX, this->var_indices[target], controls[count - 2]));
                        this->gates.push_back(Gate(type, controls, this->var_indices[target]));
                        this->gates.push_back(Gate(rev::CX, this->var_indices[target], controls[count - 2]));
                    }
                } else if (command == "p") {  // Peres Gate
                    std::vector<line_t> controls;
                    string control, target;
                    for (auto i = 0; i < 2; i++) {
                        iss >> control;
                        controls.push_back(this->var_indices[control]);
                    }
                    iss >> target;
                    this->gates.push_back(Gate(rev::CCX, controls, this->var_indices[target]));
                    this->gates.push_back(Gate(rev::CX, controls[0], controls[1]));
                } else if (command == "pi") {  // Peres Gate
                    std::vector<line_t> controls;
                    string control, target;
                    for (auto i = 0; i < 2; i++) {
                        iss >> control;
                        controls.push_back(this->var_indices[control]);
                    }
                    iss >> target;
                    this->gates.push_back(Gate(rev::CX, controls[0], controls[1]));
                    this->gates.push_back(Gate(rev::CCX, controls, this->var_indices[target]));
                } else if (command == "v") {  // V Gate
                    string control, target;
                    iss >> control;
                    iss >> target;
                    this->gates.push_back(Gate(rev::H, this->var_indices[target]));
                    this->gates.push_back(Gate(rev::T, this->var_indices[control]));
                    this->gates.push_back(Gate(rev::CX, this->var_indices[target], this->var_indices[control]));
                    this->gates.push_back(Gate(rev::T, this->var_indices[target]));
                    this->gates.push_back(Gate(rev::TDG, this->var_indices[control]));
                    this->gates.push_back(Gate(rev::CX, this->var_indices[target], this->var_indices[control]));
                    this->gates.push_back(Gate(rev::H, this->var_indices[target]));
                } else if (command == "v+") {  // V Gate
                    string control, target;
                    iss >> control;
                    iss >> target;
                    this->gates.push_back(Gate(rev::H, this->var_indices[target]));
                    this->gates.push_back(Gate(rev::TDG, this->var_indices[control]));
                    this->gates.push_back(Gate(rev::CX, this->var_indices[target], this->var_indices[control]));
                    this->gates.push_back(Gate(rev::TDG, this->var_indices[target]));
                    this->gates.push_back(Gate(rev::T, this->var_indices[control]));
                    this->gates.push_back(Gate(rev::CX, this->var_indices[target], this->var_indices[control]));
                    this->gates.push_back(Gate(rev::H, this->var_indices[target]));
                } else
                    cout << "Invalid gate: " << command << endl;
            }
        }
    } else
        cout << "File can not be opened." << endl;
}

void Circuit::addQubit(string qubit) {
    if (this->var_indices.find(qubit) == this->var_indices.end()) {
        this->var_names[this->lines.size()] = qubit;
        this->var_indices[qubit] = this->lines.size();
        this->lines.push_back(this->lines.size());
    }
}
void Circuit::addCbit(string cbit) {
    if (this->var_cindices.find(cbit) == this->var_cindices.end()) {
        this->var_cnames[this->clines.size()] = cbit;
        this->var_cindices[cbit] = this->clines.size();
        this->clines.push_back(this->clines.size());
    }
}
unsigned Circuit::twoQubitGateCount() {
    unsigned count = 0;
    for (std::vector<rev::Gate>::iterator i = this->gates.begin(); i != this->gates.end(); ++i)
        if (CX == i->getType() || SWAP == i->getType()) count++;
    return count;
}
unsigned Circuit::getCNOTDepth() {
    this->CNOTDepth = std::map<rev::line_t, double>();
    for (std::vector<rev::line_t>::iterator l = this->lines.begin(); l != this->lines.end(); ++l) {
        this->CNOTDepth[*l] = 0;
    }
    for (std::vector<rev::Gate>::iterator g = this->gates.begin(); g != this->gates.end(); ++g) {
        if (CX == g->getType()) {
            double cdepth = this->CNOTDepth[g->getControls()[0]] > this->CNOTDepth[g->getTargets()[0]] ? this->CNOTDepth[g->getControls()[0]] + 1 : this->CNOTDepth[g->getTargets()[0]] + 1;
            this->CNOTDepth[g->getControls()[0]] = cdepth;
            this->CNOTDepth[g->getTargets()[0]] = cdepth;
        }
    }
    return maxCNOTDepth();
}

unsigned Circuit::maxCNOTDepth() {
    double maxDepth = 0;
    for (std::vector<rev::line_t>::iterator l = this->lines.begin(); l != this->lines.end(); ++l) {
        if (maxDepth < this->CNOTDepth[*l]) maxDepth = CNOTDepth[*l];
    }
    return maxDepth;
}
void Circuit::reduceCNOTDepth(Gate& g) {
    double cdepth = this->CNOTDepth[g.getControls()[0]] < this->CNOTDepth[g.getTargets()[0]] ? this->CNOTDepth[g.getControls()[0]] - 1 : this->CNOTDepth[g.getTargets()[0]] - 1;
    this->CNOTDepth[g.getControls()[0]] = cdepth;
    this->CNOTDepth[g.getTargets()[0]] = cdepth;
    std::cout << "Cnot Depth";
    for (std::vector<rev::line_t>::iterator l = this->lines.begin(); l != this->lines.end(); ++l) {
        std::cout << this->CNOTDepth[*l] << " ";
    }
    std::cout << std::endl;
}
unsigned Circuit::twoQubitGateCount(unsigned pos) {
    unsigned count = 0;
    for (std::vector<rev::Gate>::iterator i = this->gates.begin() + pos; i != this->gates.end(); ++i)
        if (CX == i->getType()) count++;
    return count;
}

void Circuit::writeQASM(string file, unsigned size) {
    char fileName[50];
    int i, ln = file.length();
    for (i = 0; i < ln; i++)
        fileName[i] = file[i];
    fileName[i] = '\0';
    ofstream myfile(fileName);
    std::map<rev::line_t, string> vars;

    if (myfile.is_open()) {
        myfile << "OPENQASM 2.0;" << std::endl;
        myfile << "include \"qelib1.inc\";" << std::endl;
        myfile << "qreg q[" << size << "];" << std::endl;
        myfile << "creg q[" << size << "];" << std::endl;
        for (vector<Gate>::iterator i = this->gates.begin(); i != this->gates.end(); ++i)
            (*i).display(myfile);
        myfile.close();
    } else
        cout << "File can not be opened." << endl;
}

void Circuit::writeReal(string file) {
    char fileName[50];
    int i, ln = file.length();
    for (i = 0; i < ln; i++)
        fileName[i] = file[i];
    fileName[i] = '\0';
    ofstream myfile(fileName);
    std::map<rev::line_t, string> vars;
    variables.clear();
    int vcount = 0, gcount = 0;
    for (vector<line_t>::iterator i = lines.begin(); i != lines.end(); ++i) {
        std::ostringstream vrname;
        vrname << "x" << vcount++;
        vars[*i] = vrname.str();
    }

    if (myfile.is_open()) {
        myfile << ".version 1.0" << endl;
        myfile << ".numvars " << lines.size() << endl;
        myfile << ".variables";
        for (vector<line_t>::iterator i = lines.begin(); i != lines.end(); ++i)
            myfile << " " << vars[(*i)];
        myfile << endl;
        myfile << ".begin" << endl;
        for (vector<Gate>::iterator i = gates.begin(); i != gates.end(); ++i) {
            vector<line_t> controls = (*i).getControls();
            vector<line_t> targets = (*i).getTargets();
            if (X == (*i).getType())
                myfile << "x";
            else if (Y == (*i).getType())
                myfile << "y";
            else if (Z == (*i).getType())
                myfile << "z";
            else if (H == (*i).getType())
                myfile << "h";
            else if (S == (*i).getType())
                myfile << "s";
            else if (SDG == (*i).getType())
                myfile << "s+";
            else if (T == (*i).getType())
                myfile << "p";
            else if (TDG == (*i).getType())
                myfile << "p+";
            else if (RX == (*i).getType())
                myfile << "rx";
            else if (RY == (*i).getType())
                myfile << "ry";
            else if (RZ == (*i).getType())
                myfile << "rz";
            else if (U1 == (*i).getType())
                myfile << "u1";
            else if (U2 == (*i).getType())
                myfile << "u2";
            else if (U3 == (*i).getType())
                myfile << "u3";
            else if (CX == (*i).getType() || CCX == (*i).getType() || MCX == (*i).getType())
                myfile << "t";
            else if (SWAP == (*i).getType())
                myfile << "f";
            else {
                myfile << "unknown";
                return;
            }
            if (CX == (*i).getType() || CCX == (*i).getType() || MCX == (*i).getType() || SWAP == (*i).getType())
                myfile << controls.size() + targets.size();
            for (vector<line_t>::iterator i = controls.begin(); i != controls.end(); ++i)
                myfile << " " << vars[(*i)];
            for (vector<line_t>::iterator i = targets.begin(); i != targets.end(); ++i)
                myfile << " " << vars[(*i)];
            myfile << endl;
        }
        // std::cout<<"Hello3"<<std::endl;
        myfile << ".end" << endl;
        myfile.close();
    } else
        cout << "File can not be opened." << endl;
}

vector<Gate> Circuit::getGates() {
    return gates;
}
vector<line_t> Circuit::getLines() {
    return lines;
}

map<line_t, line_t> Circuit::getMeasures() {
    return measures;
}

void Circuit::clear() {
    this->gates.clear();
    this->lines.clear();
}

void Circuit::setGates(vector<Gate>& ngates) {
    this->gates.clear();
    for (vector<Gate>::iterator g = ngates.begin(); g != ngates.end(); ++g)
        this->gates.push_back(*g);
}
void Circuit::setLines(vector<line_t>& nlines) {
    this->lines.clear();
    for (vector<line_t>::iterator l = nlines.begin(); l != nlines.end(); ++l)
        this->lines.push_back(*l);
}

void Circuit::setMeasures(map<line_t, line_t>& nmeasures) {
    this->measures.clear();
    for (map<line_t, line_t>::iterator l = nmeasures.begin(); l != nmeasures.end(); ++l)
        this->measures[l->first] = l->second;
}

void Circuit::init() {
    this->variables = vector<string>();
    this->lines = vector<line_t>();
    this->clines = vector<line_t>();
    this->measures = map<line_t, line_t>();
    this->gates = vector<Gate>();
    this->var_names = map<line_t, string>();
    this->var_indices = map<string, line_t>();
    this->var_cnames = map<line_t, string>();
    this->var_cindices = map<string, line_t>();
}
Circuit::Circuit() {
    init();
}
Circuit::Circuit(const Circuit& ckt) {
    init();
    for (vector<string>::const_iterator l = ckt.variables.begin(); l != ckt.variables.end(); ++l)
        this->variables.push_back(*l);
    for (vector<line_t>::const_iterator l = ckt.lines.begin(); l != ckt.lines.end(); ++l)
        this->lines.push_back(*l);
    for (vector<line_t>::const_iterator l = ckt.clines.begin(); l != ckt.clines.end(); ++l)
        this->clines.push_back(*l);
    for (map<line_t, line_t>::const_iterator vn = ckt.measures.begin(); vn != ckt.measures.end(); ++vn)
        this->measures[vn->first] = vn->second;
    for (vector<Gate>::const_iterator g = ckt.gates.begin(); g != ckt.gates.end(); ++g)
        this->gates.push_back(*g);
    for (map<line_t, string>::const_iterator vn = ckt.var_names.begin(); vn != ckt.var_names.end(); ++vn)
        this->var_names[vn->first] = vn->second;
    for (map<string, line_t>::const_iterator vi = ckt.var_indices.begin(); vi != ckt.var_indices.end(); ++vi)
        this->var_indices[vi->first] = vi->second;
    for (map<line_t, string>::const_iterator vn = ckt.var_cnames.begin(); vn != ckt.var_cnames.end(); ++vn)
        this->var_cnames[vn->first] = vn->second;
    for (map<string, line_t>::const_iterator vi = ckt.var_cindices.begin(); vi != ckt.var_cindices.end(); ++vi)
        this->var_cindices[vi->first] = vi->second;
    for (map<rev::line_t, double>::const_iterator cd = ckt.CNOTDepth.begin(); cd != ckt.CNOTDepth.end(); ++cd)
        this->CNOTDepth[cd->first] = cd->second;
}
void Circuit::operator=(const Circuit& ckt) {
    init();
    for (vector<string>::const_iterator l = ckt.variables.begin(); l != ckt.variables.end(); ++l)
        this->variables.push_back(*l);
    for (vector<line_t>::const_iterator l = ckt.lines.begin(); l != ckt.lines.end(); ++l)
        this->lines.push_back(*l);
    for (vector<line_t>::const_iterator l = ckt.clines.begin(); l != ckt.clines.end(); ++l)
        this->clines.push_back(*l);
    for (map<line_t, line_t>::const_iterator vn = ckt.measures.begin(); vn != ckt.measures.end(); ++vn)
        this->measures[vn->first] = vn->second;
    /*for (map<line_t, line_t>::const_iterator vn = ckt.flgmap.begin(); vn != ckt.flgmap.end(); ++vn)
      this->flgmap[vn->first] = vn->second;*/
    for (vector<Gate>::const_iterator g = ckt.gates.begin(); g != ckt.gates.end(); ++g)
        this->gates.push_back(*g);
    for (map<line_t, string>::const_iterator vn = ckt.var_names.begin(); vn != ckt.var_names.end(); ++vn)
        this->var_names[vn->first] = vn->second;
    for (map<string, line_t>::const_iterator vi = ckt.var_indices.begin(); vi != ckt.var_indices.end(); ++vi)
        this->var_indices[vi->first] = vi->second;
    for (map<line_t, string>::const_iterator vn = ckt.var_cnames.begin(); vn != ckt.var_cnames.end(); ++vn)
        this->var_cnames[vn->first] = vn->second;
    for (map<string, line_t>::const_iterator vi = ckt.var_cindices.begin(); vi != ckt.var_cindices.end(); ++vi)
        this->var_cindices[vi->first] = vi->second;
    for (map<rev::line_t, double>::const_iterator cd = ckt.CNOTDepth.begin(); cd != ckt.CNOTDepth.end(); ++cd)
        this->CNOTDepth[cd->first] = cd->second;
}
Circuit::~Circuit() {
}
int rev::convertToInteger(string value) {
    int i = 0, result = 0, len = value.length();
    while (' ' == value[i]) i++;  // skipping initial spaces
    for (; i < len; i++) {
        if (48 <= value[i] && 57 >= value[i])
            result = 10 * result + value[i] - 48;
    }
    return result;
}
