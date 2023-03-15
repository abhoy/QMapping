# QMapping

QMapping: A tool for architectural mapping of quantum circuits

### Installation:
    Download (clone) the source to a local directory and run the following commands
    for source build:

    1. Create a python virtual environment
        python3 -m venv <env name> 

    2. Activate the environment
        source <env name>/bin/activate

    3. Install the following packages 
        pip install cmake ninja setuptools setuptools_scm 

    4. Install the QMapping tool
        python3 setup.py build 

### Example
    
####    1. Using Python Interface:

        import PyQMapping
        
        #Layout as undirected graph
        layout = '''
        #L Q5
        Qubit: Q0 Q1 Q2 Q3 Q4 Q5 
        Coupling:
        Q0 Q1
        Q1 Q4
        Q1 Q2
        Q2 Q3
        Q3 Q5
        end
        '''

        input = '''
        OPENQASM 2.0;
        include "qelib1.inc";
        qreg q[27];
        creg meas[27];
        z q[0];
        x q[1];
        u3(pi/2,0,pi) q[2];
        cx q[0],q[2];
        h q[2];
        cx q[0],q[2];
        cx q[1],q[2];
        '''
        print(PyQMapping.qx_mapper(layout, input, 1))

####    2. Using Direct Executable

        ./build/qx_mapper <path to layout file> <path to qasm file> <path to output file> 
