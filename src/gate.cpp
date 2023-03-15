/*--------------------------------------------------------------------------------
Copyright (C) 2023, Abhoy Kole, Kamalika Datta, Indranil Sengupta, Rolf Drechsler

Quantum Mapping Tool (QMT) for mapping Quantum circuits to IBM architecture
released under the MIT licence.

All right reserved.
----------------------------------------------------------------------------------*/

#include "gate.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

using namespace rev;

int Gate::gate_count = 0;

string Gate::removeZero(string str) {
    int i = str.length();
    while (str[i] == '0')
        i--;

    str.erase(i);
    return str;
}

double Gate::round_sig(double num, int sig) {
    double result = static_cast<double>(static_cast<int>(num * pow(10, sig))) / pow(10, sig);
    return result;
}

bool GateCompare::operator()(const Gate& lhs, const Gate& rhs) {
    return lhs.getId() < rhs.getId();
}

void Gate::updateLines(line_t control, line_t target) {
    controls[0] = control;
    targets[0] = target;
}

void Gate::updateLines(line_t target) {
    targets[0] = target;
}

vector<line_t> Gate::getControls() {
    return controls;
}

vector<line_t> Gate::getTargets() {
    return targets;
}

gate_t Gate::getType() {
    return type;
}

int Gate::getId() const {
    return gate_id;
}

Gate::Gate(gate_t type, line_t target) {
    gate_id = ++gate_count;
    targets.push_back(target);
    this->type = type;
}

void Gate::makeEquivalent() {
    switch (this->type) {
        case U1:
        case U2:
        case U3:
        case RZ:
            break;
        case X: {
            this->type = U3;
            this->theta = "PI";
            this->phi = "0";
            this->lamda = "PI";
            break;
        }
        case Y: {
            this->type = U3;
            this->theta = "PI";
            this->phi = "PI/2";
            this->lamda = "PI/2";
            break;
        }
        case Z: {
            this->type = U1;
            this->theta = "0";
            this->phi = "0";
            this->lamda = "PI";
            break;
        }
        case H: {
            this->type = U2;
            this->theta = "PI/2";
            this->phi = "0";
            this->lamda = "PI";
            break;
        }
        case S: {
            this->type = U1;
            this->theta = "0";
            this->phi = "0";
            this->lamda = "PI/2";
            break;
        }
        case SDG: {
            this->type = U1;
            this->theta = "0";
            this->phi = "0";
            this->lamda = "-PI/2";
            break;
        }
        case T: {
            this->type = U1;
            this->theta = "0";
            this->phi = "0";
            this->lamda = "PI/4";
            break;
        }
        case TDG: {
            this->type = U1;
            this->theta = "0";
            this->phi = "0";
            this->lamda = "-PI/4";
            break;
        }
        default: {
            std::cout << "INVALID GATE TYPE" << std::endl;
        }
    }
}

Gate& Gate::operator=(const Gate& g) {
    this->gate_id = g.gate_id;
    this->type = g.type;
    this->theta = g.theta;
    this->phi = g.phi;
    this->lamda = g.lamda;
    for (auto control : g.controls)
        this->controls.push_back(control);
    for (auto target : g.targets)
        this->targets.push_back(target);

    return *this;
}

Gate::Gate(const Gate& g) {
    this->gate_id = g.gate_id;
    this->type = g.type;
    this->theta = g.theta;
    this->phi = g.phi;
    this->lamda = g.lamda;
    for (auto control : g.controls)
        this->controls.push_back(control);
    for (auto target : g.targets)
        this->targets.push_back(target);
}

Gate::Gate(gate_t type, line_t control, line_t target) {
    gate_id = ++gate_count;
    if (rev::SWAP == type)
        targets.push_back(control);
    else
        controls.push_back(control);
    targets.push_back(target);
    this->type = type;
}

Gate::Gate(std::string gate, std::map<string, line_t>& var_indices) {
    gate_id = ++gate_count;
    gate = gate.substr(gate.find_first_not_of(" \n\r\t\f\v"));

    stringstream check(gate);
    string gate_type;
    check >> gate_type;

    if (gate_type == "x")
        type = X;
    else if (gate_type == "y")
        type = Y;
    else if (gate_type == "z")
        type = Z;
    else if (gate_type == "h")
        type = H;
    else if (gate_type == "s")
        type = S;
    else if (gate_type == "sdg")
        type = SDG;
    else if (gate_type == "t")
        type = T;
    else if (gate_type == "tdg")
        type = TDG;
    else if (gate_type == "cx")
        type = CX;
    else if (gate_type == "ccx")
        type = CCX;
    else if (gate_type == "mcx")
        type = MCX;
    else if (gate_type == "swap")
        type = SWAP;
    else if (gate_type.substr(0, 2) == "rz") {
        type = U1;  // RZ;
        theta = "0.";
        phi = "0.";
        lamda = gate_type.substr(3, gate_type.length() - 4);
    } else if (gate_type.substr(0, 2) == "u1") {
        type = U1;
        theta = "0.";
        phi = "0.";
        lamda = gate_type.substr(3, gate_type.length() - 4);
    } else if (gate_type.substr(0, 2) == "u2") {
        type = U2;
        theta = "PI/2";
        string angles = gate_type.substr(3, gate_type.length() - 4);
        size_t start, end = 0;
        start = angles.find_first_not_of(',', end);
        end = angles.find(',', start);
        phi = angles.substr(start, end - start);
        start = angles.find_first_not_of(',', end);
        end = angles.find(',', start);
        lamda = angles.substr(start, end - start);
    } else if (gate_type.substr(0, 2) == "u3") {
        type = U3;
        string angles = gate_type.substr(3, gate_type.length() - 4);
        size_t start, end = 0;
        start = angles.find_first_not_of(',', end);
        end = angles.find(',', start);
        theta = angles.substr(start, end - start);
        start = angles.find_first_not_of(',', end);
        end = angles.find(',', start);
        phi = angles.substr(start, end - start);
        start = angles.find_first_not_of(',', end);
        end = angles.find(',', start);
        lamda = angles.substr(start, end - start);
    } else
        cout << "Invalid gate: " << gate_type << endl;

    string qubits;
    check >> qubits;
    if (CX == type || CCX == type || MCX == type || SWAP == type) {
        for (auto pos1 = string::npos, pos2 = qubits.find(","); pos2 != string::npos;) {
            string qubit = pos1 != string::npos ? qubits.substr(pos1 + 1, pos2 - pos1 - 1) : qubits.substr(0, pos2);
            if (SWAP == type)
                targets.push_back(var_indices[qubit]);
            else
                controls.push_back(var_indices[qubit]);
            pos1 = pos2, pos2 = qubits.find(",", pos1 + 1);
            if (pos2 == string::npos) {
                if ((pos2 = qubits.find(";")) == string::npos) {
                    check >> qubits;
                    pos1 = string::npos, pos2 = qubits.find(",");
                    if (pos2 == string::npos) {
                        qubit = qubits.substr(0, qubits.find(";"));
                        targets.push_back(var_indices[qubit]);
                        break;
                    }
                } else {
                    qubit = qubits.substr(pos1 + 1, pos2 - pos1 - 1);
                    targets.push_back(var_indices[qubit]);
                    break;
                }
            }
        }
    } else {  // arbitrary single qubit gate
        string qubit = qubits.substr(0, qubits.find(";"));
        targets.push_back(var_indices[qubit]);
    }
}
// Real type
Gate::Gate(gate_t type, std::vector<line_t>& controls, line_t target) {
    gate_id = ++gate_count;
    for (auto c : controls)
        this->controls.push_back(c);
    targets.push_back(target);
    this->type = type;
}
void Gate::display(ostream& os) {
    switch (type) {
        case X:
            os << "x ";
            break;
        case CX:
            os << "cx ";
            break;
        case CCX:
            os << "ccx ";
            break;
        case MCX:
            os << "mcx ";
            break;
        case Y:
            os << "y ";
            break;
        case Z:
            os << "z ";
            break;
        case H:
            os << "h ";
            break;
        case S:
            os << "s ";
            break;
        case SDG:
            os << "sdg ";
            break;
        case T:
            os << "t ";
            break;
        case TDG:
            os << "tdg ";
            break;
        case RZ: {
            std::ostringstream strs;
            strs << lamda;
            std::string str = strs.str();
            os << "rz(" << removeZero(str) << ") ";
            break;
        }
        case U1: {
            std::ostringstream strs;
            strs << lamda;
            std::string str = strs.str();
            os << "u1(" << removeZero(str) << ") ";
            break;
        }
        case U2: {
            std::ostringstream strs2;
            strs2 << phi;
            std::string str2 = strs2.str();
            std::ostringstream strs3;
            strs3 << lamda;
            std::string str3 = strs3.str();
            os << "u2(" << removeZero(str2) << "," << removeZero(str3) << ") ";
            break;
        }
        case U3: {
            std::ostringstream strs1;
            strs1 << theta;
            std::string str1 = strs1.str();
            std::ostringstream strs2;
            strs2 << phi;
            std::string str2 = strs2.str();
            std::ostringstream strs3;
            strs3 << lamda;
            std::string str3 = strs3.str();
            os << "u3(" << removeZero(str1) << "," << removeZero(str2) << "," << removeZero(str3) << ") ";
            break;
        }
        case SWAP:
            os << "swap ";
            break;
        default:
            os << "Invalid choice" << endl;
    }
    for (vector<line_t>::iterator i = controls.begin(); i != controls.end(); ++i)
        os << "q[" << *i << "],";
    for (vector<line_t>::iterator i = targets.begin(); i != targets.end(); ++i)
        os << "q[" << *i << "]" << (targets.size() > 1 && (i + 1) != targets.end() ? "," : ";");
    os << endl;
}
void Gate::display2(ostream& os) {
    switch (type) {
        case X:
            os << "qc.x(";
            break;
        case CX:
            os << "qc.cx(";
            break;
        case CCX:
            os << "qc.ccx(";
            break;
        case MCX:
            os << "qc.mct([";
            break;
        case Y:
            os << "qc.y(";
            break;
        case Z:
            os << "qc.z(";
            break;
        case H:
            os << "qc.h(";
            break;
        case S:
            os << "qc.s(";
            break;
        case SDG:
            os << "qc.sdg(";
            break;
        case T:
            os << "qc.t(";
            break;
        case TDG:
            os << "qc.tdg(";
            break;
        case RZ: {
            std::ostringstream strs;
            strs << lamda;
            std::string str = strs.str();
            os << "qc.u1(" << removeZero(str) << ", ";
            break;
        }
        case U1: {
            std::ostringstream strs;
            strs << lamda;
            std::string str = strs.str();
            os << "qc.u1(" << removeZero(str) << ", ";
            break;
        }
        case U2: {
            std::ostringstream strs2;
            strs2 << phi;
            std::string str2 = strs2.str();
            std::ostringstream strs3;
            strs3 << lamda;
            std::string str3 = strs3.str();
            os << "qc.u2(" << removeZero(str2) << "," << removeZero(str3) << ", ";
            break;
        }
        case U3: {
            std::ostringstream strs1;
            strs1 << theta;
            std::string str1 = strs1.str();
            std::ostringstream strs2;
            strs2 << phi;
            std::string str2 = strs2.str();
            std::ostringstream strs3;
            strs3 << lamda;
            std::string str3 = strs3.str();
            os << "qc.u3(" << removeZero(str1) << "," << removeZero(str2) << "," << removeZero(str3) << ", ";
            break;
        }
        case SWAP:
            cout << "qc.swap(";
            break;
        default:
            os << "Invalid choice" << endl;
    }
    if (controls.size() > 0) {
        for (vector<line_t>::iterator i = controls.begin(); i != controls.end() - 1; ++i)
            os << *i << ",";
        if (type == MCX)
            os << controls[controls.size() - 1] << "],";
        else
            os << controls[controls.size() - 1] << ",";
    }
    for (vector<line_t>::iterator i = targets.begin(); i != targets.end(); ++i)
        os << *i << (targets.size() > 1 && (i + 1) != targets.end() ? ", " : ")");
    os << endl;
}
