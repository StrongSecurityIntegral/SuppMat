#include <vector>
#include <iostream>
#include <cstdint>
#include <string>
#include <array>
#include <memory>

#include "gurobi_c++.h"

using namespace std;

string exec(const char* cmd){
//Exec cmd and grab the stdout
    array<char, 128> buffer;
    string result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result += buffer.data();
    }
    return result;
}

int main(){

	uint const rMax = 13;
	uint const blockSize = 64;
	int rMax_minus_1 = (int)(rMax) - 1;

	//Generate the model file
	exec("sage modelPresent.sage");

	//Input and output patterns
	vector<uint8_t> input(64,1);
	input[0] = 0;
	vector<uint8_t> output(64,0);
	output[0] = 1;

	//Key patterns (first one fron PRESENT_13r_keyPatterns.txt") 
	vector<string> kp_str({
	"1000000010000000000000000000000000000000000000000000000000000000",
	"0000000000000000000000000000000000000000000000000000000000001000",
	"0000000000000000000000000000000000000000000000000000000000000000",
	"0000000000000000000000000000000000000000000000000000000000000000",
	"0001000100000000000000000000000000010000000100000000000000000000",
	"1111111000011100101010000100101010110000100100001110110000001101",
	"0100000000000001110000000000000011010000000011000001000100000000",
	"0000000010010000000000000000000000000001101000000000000000000000",
	"0000000000000000000000000000000000000000000000000000000000000000",
	"0000000000000000000000000000000000000000000000000000000000000000",
	"0000000000000000001000000000000000000000000000000000000000000000",
	"0010000000000000000000000000000000000000000000000000000000000000",});
	vector<vector<uint8_t>> kp(rMax_minus_1, vector<uint8_t>(blockSize));
	for(int i = 0; i < rMax_minus_1; i++){
		for(uint j = 0; j < blockSize; j++){
			if(kp_str[i][j] == '0') kp[i][j] = 0;
			else kp[i][j] = 1;
		}
	}

	//Read the model
	GRBEnv gurobiEnv;
	gurobiEnv.set(GRB_IntParam_OutputFlag,0);
	GRBModel m(gurobiEnv, "PRESENT_13r.mps");
	

	//Fix the value of input/output/ke patterns
	for(uint i = 0; i < blockSize; i++){
		m.addConstr(m.getVarByName("x_0_"+to_string(i)) == input[i]);
		m.addConstr(m.getVarByName("y_"+to_string(rMax_minus_1)+"_"+to_string(i)) == output[i]);
	}
	for(int r = 0; r < rMax_minus_1; r++){
		for(uint i = 0; i < blockSize; i++){
			m.addConstr(m.getVarByName("k_"+to_string(r)+"_"+to_string(i)) == kp[r][i]);
		}
	}

	//Count trails
	m.set(GRB_IntParam_PoolSearchMode, 2);
	m.set(GRB_IntParam_PoolSolutions, 2000000000); //maximum number of solution allowed by Gurobi

	//Dummy objective (sometime help speeding up the search)
	GRBLinExpr objExpr;
	for(uint i = 0; i < blockSize; i++){
		objExpr += m.getVarByName("x_"+to_string(rMax/2)+"_"+to_string(i));
	}
	m.setObjective(objExpr);
	m.update();

	m.optimize();

	int solcount = m.get(GRB_IntAttr_SolCount);
	cout << "Number of trails : " << solcount << endl;
}