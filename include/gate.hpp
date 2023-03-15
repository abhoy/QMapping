/*--------------------------------------------------------------------------------
Copyright (C) 2023, Abhoy Kole, Kamalika Datta, Indranil Sengupta, Rolf Drechsler

Quantum Mapping Tool (QMT) for mapping Quantum circuits to IBM architecture
released under the MIT licence.

All right reserved.
----------------------------------------------------------------------------------*/

#include <cmath>
#include <map>
#include <vector>

#define PI 3.142857143

#ifndef GATE_HPP
#define GATE_HPP
using namespace std;
namespace rev {

enum { X,
       Y,
       Z,
       H,
       S,
       SDG,
       T,
       TDG,
       RX,
       RY,
       RZ,
       U1,
       U2,
       U3,
       CX,
       CCX,
       MCX,
       SWAP };

typedef unsigned gate_t;
typedef unsigned line_t;
typedef string angle_t;
class Gate;
struct GateCompare {
    bool operator()(const Gate& lhs, const Gate& rhs);
};
class Gate {
   private:
    static int gate_count;
    int gate_id;
    vector<line_t> controls;
    vector<line_t> targets;
    gate_t type;
    angle_t theta;
    angle_t phi;
    angle_t lamda;

   public:
    Gate(const Gate& g);
    Gate(gate_t type, line_t);

    Gate(gate_t, line_t, line_t);
    Gate(std::string gate, std::map<string, line_t>& var_indices);    // QASM type
    Gate(gate_t type, std::vector<line_t>& controls, line_t target);  // Real type
    Gate& operator=(const Gate& g);
    void display(ostream& os);
    void display2(ostream& os);
    gate_t getType();
    int getId() const;
    vector<line_t> getControls();
    vector<line_t> getTargets();
    void updateLines(line_t control, line_t target);
    void updateLines(line_t target);
    string removeZero(string str);
    double round_sig(double num, int sig);
    void makeEquivalent();  // Need to redefine due angle type changed from double to string
};
}  // namespace rev
#endif /* GATE_HPP */
