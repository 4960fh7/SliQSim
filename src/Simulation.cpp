#include "Simulator.h"
#include "util_sim.h"


/**Function*************************************************************

  Synopsis    [Initailize simulator]

  Description [This function will set #qubits n, construct initial state, and enable dynamic reordering]

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::init_simulator(int nQubits)
{
    n = nQubits; // set the number n here
    manager = Cudd_Init(n, n, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
    int *constants = new int[n];
    for (int i = 0; i < n; i++)
        constants[i] = 0; // TODO: costom initial state
    measured_qubits_to_clbits = std::vector<std::vector<int>>(n, std::vector<int>(0));
    init_state(constants);
    delete[] constants;
    if (isReorder) Cudd_AutodynEnable(manager, CUDD_REORDER_SYMM_SIFT);
}


/**Function*************************************************************

  Synopsis    [parse and simulate the qasm file]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::sim_qasm_pass(std::string qasm, int forced_val, int mid_measure_qubit)
{
    this->p0_mid = 0.0;
    this->p1_mid = 0.0;
    std::string inStr;
    std::stringstream inFile_ss(qasm);
    while (getline(inFile_ss, inStr))
    {
        inStr = inStr.substr(0, inStr.find("//"));
        
        inStr.erase(0, inStr.find_first_not_of(" \t\r\n"));
        if (inStr.empty()) continue;

        if (inStr.find("if") != std::string::npos) {
            size_t start = inStr.find("==");
            size_t end = inStr.find(")");
            if (start != std::string::npos && end != std::string::npos && start < end) {
                try {
                    int condition_val = std::stoi(inStr.substr(start + 2, end - start - 2));
                    if (forced_val != condition_val) continue; 
                } catch (...) { continue; } 
                
                inStr = inStr.substr(end + 1); 
                inStr.erase(0, inStr.find_first_not_of(" \t")); 
            }
        }
        // ------------------------------------------------

        if (inStr.find_first_not_of("\t\n\r ") != std::string::npos)
        {
            std::stringstream inStr_ss(inStr);
            getline(inStr_ss, inStr, ' ');
            if (inStr == "qreg")
            {
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                init_simulator(stoi(inStr));
                
                if (sim_type == 1 && n > 50) 
                {
                    std::cerr << "[Error]: The all-amplitude mode will print the whole state vector, which is too long for large qubit number. Please consider using the sampling mode." << std::endl;
                    assert(sim_type != 1 || n <= 50); 
                }
            }
            else if (inStr == "creg")
            {
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                nClbits = stoi(inStr);                
            }
            else if (inStr == "OPENQASM"){;}
            else if (inStr == "include"){;}
            else if (inStr == "measure")
            {
                isMeasure = 1;
                int qIndex = 0, cIndex = 0;
                try {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    qIndex = std::stoi(inStr);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    cIndex = std::stoi(inStr);
                } catch (...) { continue; }

                if (qIndex == mid_measure_qubit) {
                    // Force the mid-circuit measurement to the target pass value
                    execute_mid_circuit_measure(qIndex, forced_val);
                } else {
                    measure(qIndex, cIndex); 
                }
            }
            else
            {
                if (inStr == "x")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    PauliX(stoi(inStr));
                }
                else if (inStr == "y")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    PauliY(stoi(inStr));
                }
                else if (inStr == "z")
                {
                    std::vector<int> iqubit(1);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    iqubit[0] = stoi(inStr);
                    PauliZ(iqubit);
                    iqubit.clear();
                }
                else if (inStr == "h")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    Hadamard(stoi(inStr));
                }
                else if (inStr == "s")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift(2, iqubit);
                }
                else if (inStr == "sdg")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift_dagger(-2, iqubit);
                }
                else if (inStr == "t")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift(4, iqubit);
                }
                else if (inStr == "tdg")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    int iqubit = stoi(inStr);
                    Phase_shift_dagger(-4, iqubit);
                }
                else if (inStr == "rx(pi/2)")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    rx_pi_2(stoi(inStr));
                }
                else if (inStr == "ry(pi/2)")
                {
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    ry_pi_2(stoi(inStr));
                }
                else if (inStr == "cx")
                {
                    std::vector<int> cont(1);
                    std::vector<int> ncont(0);
                    int targ;
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    cont[0] = stoi(inStr);
                    getline(inStr_ss, inStr, '[');
                    getline(inStr_ss, inStr, ']');
                    targ = stoi(inStr);
                    Toffoli(targ, cont, ncont);
                    cont.clear();
                    ncont.clear();
                }
                else if (inStr == "cz")
                {
                    std::vector<int> iqubit(2);
                    for (int i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        iqubit[i] = stoi(inStr);
                    }
                    PauliZ(iqubit);
                    iqubit.clear();
                }
                else if (inStr == "swap")
                {
                    int swapA, swapB;
                    std::vector<int> cont(0);
                    for (int i = 0; i < 2; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        if (i == 0)
                            swapA = stoi(inStr);
                        else
                            swapB = stoi(inStr);
                    }
                    Fredkin(swapA, swapB, cont);
                    cont.clear();
                }
                else if (inStr == "cswap")
                {
                    int swapA, swapB;
                    std::vector<int> cont(1);
                    for (int i = 0; i < 3; i++)
                    {
                        getline(inStr_ss, inStr, '[');
                        getline(inStr_ss, inStr, ']');
                        if (i == 0)
                            cont[i] = stoi(inStr);
                        else if (i == 1)
                            swapA = stoi(inStr);
                        else
                            swapB = stoi(inStr);
                    }
                    Fredkin(swapA, swapB, cont);
                    cont.clear();
                }
                else if (inStr == "ccx" || inStr == "mcx")
                {
                    std::vector<int> cont(0);
                    std::vector<int> ncont(0);
                    int targ;
                    getline(inStr_ss, inStr, '[');
                    while(getline(inStr_ss, inStr, ']'))
                    {
                        cont.push_back(stoi(inStr));
                        getline(inStr_ss, inStr, '[');
                    }
                    targ = cont.back();
                    cont.pop_back();
                    Toffoli(targ, cont, ncont);
                    cont.clear();
                    ncont.clear();
                }
                else
                {
                    std::cerr << std::endl
                            // << "[warning]: Gate \'" << inStr << "\' is not supported in this simulator. The gate is ignored ..." << std::endl;
                            << "[warning]: Syntax \'" << inStr << "\' is not supported in this simulator. The line is ignored ..." << std::endl;
                }
            }
        }
    }
    if (isReorder) Cudd_AutodynDisable(manager);
}

/**Function*************************************************************

  Synopsis    [simulate the circuit described by a qasm file]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::sim_qasm(std::string qasm)
{
    if (qasm.find("while") != std::string::npos || qasm.find("reset") != std::string::npos) {
        sim_dynamic(qasm);
        return;
    }

    int mid_measure_qubit = -1;
    bool is_dynamic = false;
    std::stringstream scan_ss(qasm);
    std::string line;
    
    while (getline(scan_ss, line)) {
        line = line.substr(0, line.find("//")); 
        if (line.find_first_not_of("\t\n\r ") == std::string::npos) continue;

        if (line.find("if") != std::string::npos) {
            is_dynamic = true;
            break; 
        }
        if (line.find("measure") != std::string::npos) {
            size_t start = line.find("[");
            size_t end = line.find("]");
            if (start != std::string::npos && end != std::string::npos && start < end) {
                try {
                    mid_measure_qubit = std::stoi(line.substr(start + 1, end - start - 1));
                } catch (...) { }
            }
        }
    }
    
    if (!is_dynamic || mid_measure_qubit == -1) {
        sim_qasm_pass(qasm, 0, -1);
        if (sim_type == 0) {
            create_bigBDD();
            if (this->shots > 0) measurement();
        } else if (sim_type == 1) {
            getStatevector();
        }
        print_results();
        return;
    }

    int total_shots = shots;
    std::unordered_map<std::string, int> final_counts;

    // PASS 0: Force mid-circuit measurement to 0
    state_count.clear(); 
    sim_qasm_pass(qasm, 0, mid_measure_qubit); 
    
    if (sim_type == 0) {
        create_bigBDD();
        this->shots = std::round(total_shots * p0_mid); 
        if (this->shots > 0) measurement();
        for (auto const& [key, val] : state_count) final_counts[key] += val;
    }
    
    double saved_p0_mid = p0_mid; 
    clear(); 
    state_count.clear(); 
    
    // PASS 1: Force mid-circuit measurement to 1
    sim_qasm_pass(qasm, 1, mid_measure_qubit); 
    
    if (sim_type == 0) {
        create_bigBDD();
        this->shots = total_shots - std::round(total_shots * saved_p0_mid); 
        if (this->shots > 0) measurement();
        for (auto const& [key, val] : state_count) final_counts[key] += val;
    }

    state_count = final_counts;
    this->shots = total_shots;
    print_results();
}

/**Function*************************************************************

  Synopsis    [print state vector and distribution of sampled outcomes]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void Simulator::print_results()
{
    std::unordered_map<std::string, int>::iterator it;
    
    run_output = "{";
    if (state_count.begin() != state_count.end()){
        run_output += "\"counts\": { ";
        for (it = state_count.begin(); it != state_count.end(); it++)
        {
            if (std::next(it) == state_count.end())
                run_output = run_output + "\"" + it->first + "\": " + std::to_string(it->second);
            else
                run_output = run_output + "\"" + it->first + "\": " + std::to_string(it->second) + ", ";
        }
        run_output += " }";
        run_output += (statevector != "null") ? ", " : ""; 
    }    

    run_output += (statevector != "null") ? "\"statevector\": " + statevector + " }" : " }";
    //return;
    std::cout << run_output << std::endl;
}

/**Function*************************************************************
  Synopsis    [Active qubit reset: Measure and apply conditional X]
***********************************************************************/
void Simulator::reset(int qIndex)
{
    int outcome = measure_and_collapse(qIndex);
    if (outcome == 1) {
        PauliX(qIndex);
    }
}

/**Function*************************************************************
  Synopsis    [Evaluate classical condition against live register]
***********************************************************************/
int Simulator::evaluate_condition(std::string condition_str)
{
    condition_str.erase(remove_if(condition_str.begin(), condition_str.end(), ::isspace), condition_str.end());
    size_t eq_pos = condition_str.find("==");
    if (eq_pos != std::string::npos) {
        int target_val = std::stoi(condition_str.substr(eq_pos + 2));
        int current_val = 0;
        for (int i = 0; i < live_creg.size(); i++) {
            current_val += live_creg[i] * (1 << i);
        }
        return (current_val == target_val) ? 1 : 0;
    }
    return 0;
}

/**Function*************************************************************
  Synopsis    [Recursive QASM Block Interpreter for Loops]
***********************************************************************/
void Simulator::execute_qasm_block(std::string block)
{
    std::string inStr;
    std::stringstream inFile_ss(block);
    
    while (getline(inFile_ss, inStr))
    {
        inStr = inStr.substr(0, inStr.find("//"));
        inStr.erase(0, inStr.find_first_not_of(" \t\r\n"));
        inStr.erase(inStr.find_last_not_of(" \t\r\n") + 1);
        if (inStr.empty()) continue;

        if (inStr.find("while") == 0) { 
            size_t start_cond = inStr.find("(");
            size_t end_cond = inStr.find(")");
            std::string condition = inStr.substr(start_cond + 1, end_cond - start_cond - 1);
            
            std::string loop_body = "";
            int bracket_count = 0;
            if (inStr.find("{") != std::string::npos) {
                bracket_count = 1;
            }
            std::string loop_line;
            while (getline(inFile_ss, loop_line)) {
                if (loop_line.find("{") != std::string::npos) bracket_count++;
                if (loop_line.find("}") != std::string::npos) bracket_count--;
                if (bracket_count == 0) break;
                loop_body += loop_line + "\n";
            }
            
            while (evaluate_condition(condition)) {
                execute_qasm_block(loop_body);
            }
            continue; 
        }

        if (inStr.find("if") == 0) {
            size_t start = inStr.find("(");
            size_t end = inStr.find(")");
            std::string condition = inStr.substr(start + 1, end - start - 1);
            
            if (!evaluate_condition(condition)) continue; 
            
            inStr = inStr.substr(end + 1); 
            inStr.erase(0, inStr.find_first_not_of(" \t\r\n")); 
        }

        std::stringstream inStr_ss(inStr);
        
        getline(inStr_ss, inStr, ' ');

        if (inStr == "OPENQASM" || inStr == "include" || inStr == "qreg" || inStr == "creg" || inStr == "}") continue;
        if (inStr == "measure") {
            int qIndex = 0, cIndex = 0;
            getline(inStr_ss, inStr, '['); getline(inStr_ss, inStr, ']'); qIndex = std::stoi(inStr);
            getline(inStr_ss, inStr, '['); getline(inStr_ss, inStr, ']'); cIndex = std::stoi(inStr);
            
            int outcome = measure_and_collapse(qIndex);
            if (cIndex >= live_creg.size()) live_creg.resize(cIndex + 1, 0);
            live_creg[cIndex] = outcome;
        }
        else if (inStr == "reset") {
            getline(inStr_ss, inStr, '['); getline(inStr_ss, inStr, ']');
            reset(std::stoi(inStr));
        }
        else if (inStr == "x") {
            getline(inStr_ss, inStr, '[');
            getline(inStr_ss, inStr, ']');
            PauliX(stoi(inStr));
        }
        else if (inStr == "y")
        {
            getline(inStr_ss, inStr, '[');
            getline(inStr_ss, inStr, ']');
            PauliY(stoi(inStr));
        }
        else if (inStr == "z")
        {
            std::vector<int> iqubit(1);
            getline(inStr_ss, inStr, '[');
            getline(inStr_ss, inStr, ']');
            iqubit[0] = stoi(inStr);
            PauliZ(iqubit);
            iqubit.clear();
        }
        else if (inStr == "h")
        {
            getline(inStr_ss, inStr, '[');
            getline(inStr_ss, inStr, ']');
            Hadamard(stoi(inStr));
        }
        else if (inStr == "s")
        {
            getline(inStr_ss, inStr, '[');
            getline(inStr_ss, inStr, ']');
            int iqubit = stoi(inStr);
            Phase_shift(2, iqubit);
        }
        else if (inStr == "sdg")
        {
            getline(inStr_ss, inStr, '[');
            getline(inStr_ss, inStr, ']');
            int iqubit = stoi(inStr);
            Phase_shift_dagger(-2, iqubit);
        }
        else if (inStr == "t")
        {
            getline(inStr_ss, inStr, '[');
            getline(inStr_ss, inStr, ']');
            int iqubit = stoi(inStr);
            Phase_shift(4, iqubit);
        }
        else if (inStr == "tdg")
        {
            getline(inStr_ss, inStr, '[');
            getline(inStr_ss, inStr, ']');
            int iqubit = stoi(inStr);
            Phase_shift_dagger(-4, iqubit);
        }
        else if (inStr == "rx(pi/2)")
        {
            getline(inStr_ss, inStr, '[');
            getline(inStr_ss, inStr, ']');
            rx_pi_2(stoi(inStr));
        }
        else if (inStr == "ry(pi/2)")
        {
            getline(inStr_ss, inStr, '[');
            getline(inStr_ss, inStr, ']');
            ry_pi_2(stoi(inStr));
        }
        else if (inStr == "cx")
        {
            std::vector<int> cont(1);
            std::vector<int> ncont(0);
            int targ;
            getline(inStr_ss, inStr, '[');
            getline(inStr_ss, inStr, ']');
            cont[0] = stoi(inStr);
            getline(inStr_ss, inStr, '[');
            getline(inStr_ss, inStr, ']');
            targ = stoi(inStr);
            Toffoli(targ, cont, ncont);
            cont.clear();
            ncont.clear();
        }
        else if (inStr == "cz")
        {
            std::vector<int> iqubit(2);
            for (int i = 0; i < 2; i++)
            {
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                iqubit[i] = stoi(inStr);
            }
            PauliZ(iqubit);
            iqubit.clear();
        }
        else if (inStr == "swap")
        {
            int swapA, swapB;
            std::vector<int> cont(0);
            for (int i = 0; i < 2; i++)
            {
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                if (i == 0)
                    swapA = stoi(inStr);
                else
                    swapB = stoi(inStr);
            }
            Fredkin(swapA, swapB, cont);
            cont.clear();
        }
        else if (inStr == "cswap")
        {
            int swapA, swapB;
            std::vector<int> cont(1);
            for (int i = 0; i < 3; i++)
            {
                getline(inStr_ss, inStr, '[');
                getline(inStr_ss, inStr, ']');
                if (i == 0)
                    cont[i] = stoi(inStr);
                else if (i == 1)
                    swapA = stoi(inStr);
                else
                    swapB = stoi(inStr);
            }
            Fredkin(swapA, swapB, cont);
            cont.clear();
        }
        else if (inStr == "ccx" || inStr == "mcx")
        {
            std::vector<int> cont(0);
            std::vector<int> ncont(0);
            int targ;
            getline(inStr_ss, inStr, '[');
            while(getline(inStr_ss, inStr, ']'))
            {
                cont.push_back(stoi(inStr));
                getline(inStr_ss, inStr, '[');
            }
            targ = cont.back();
            cont.pop_back();
            Toffoli(targ, cont, ncont);
            cont.clear();
            ncont.clear();
        }
    }
}

/**Function*************************************************************
  Synopsis    [Monte Carlo Execution Entry Point]
***********************************************************************/
void Simulator::sim_dynamic(std::string qasm)
{
    std::unordered_map<std::string, int> final_counts;

    int total_qubits = 0;
    int total_clbits = 0;
    std::stringstream scan_ss(qasm);
    std::string line;
    while(getline(scan_ss, line)) {
        if (line.find("qreg") != std::string::npos) {
            size_t start = line.find("[");
            size_t end = line.find("]");
            if(start != std::string::npos && end != std::string::npos)
                total_qubits = std::stoi(line.substr(start + 1, end - start - 1));
        }
        if (line.find("creg") != std::string::npos) {
            size_t start = line.find("[");
            size_t end = line.find("]");
            if(start != std::string::npos && end != std::string::npos)
                total_clbits = std::stoi(line.substr(start + 1, end - start - 1));
        }
    }
    
    this->n = total_qubits;
    this->nClbits = total_clbits;
    
    init_simulator(total_qubits); 

    for (int s = 0; s < shots; s++) {
        reset_shot(total_qubits);
        this->nClbits = total_clbits;
        live_creg.assign(nClbits, 0);
        normalize_factor = 1.0;

        execute_qasm_block(qasm);

        std::string final_outcome = "";
        for (int i = 0; i < nClbits; i++) {
            final_outcome = std::to_string(live_creg[i]) + final_outcome; 
        }
        final_counts[final_outcome]++;
    }

    this->state_count = final_counts;
    print_results();
    clear();
}