#include"gurobi_c++.h"
#include<vector>

using namespace std;

// linear inequalities for CRAFT (equivalent) Sbox with logic condition.
// let (x0,x1,x2,x3) --> (y0,y1,y2,y3)
// the first one means  "x1 + x2 + (1-x3) + y0 + y1 + y2 >= 1"  
vector<vector<int>> table_sbox = {
	{1,2,0,0,0,2,0,0}, {0,0,1,0,0,2,2,0}, {0,0,1,1,0,2,2,1}, {2,1,2,1,0,1,1,0},
	{2,1,1,1,0,1,1,2}, {2,2,2,1,0,0,0,0}, {0,0,0,0,2,2,2,1}, {2,1,2,0,0,0,2,0},
	{1,2,2,1,2,0,0,0}, {2,2,1,0,2,0,1,0}, {2,1,1,0,0,0,2,2}, {1,2,1,2,0,0,0,2},
	{2,2,0,0,2,1,1,1}, {1,1,2,2,2,0,0,1}, {1,1,1,2,0,2,0,2}, {1,2,2,1,0,0,1,2},
	{2,1,1,0,1,2,2,0}, {2,1,1,1,0,2,2,0}, {2,1,2,1,2,0,0,1}, {2,2,1,1,2,1,1,0},
	{1,1,1,1,2,0,2,2}, {0,0,0,2,2,1,2,2}, {2,0,0,0,2,1,2,2}, {2,0,2,0,2,1,2,1},
	{0,2,0,2,2,2,1,1}, {1,1,2,1,2,2,0,2}, {2,0,2,2,1,2,1,1}, {0,2,0,2,1,2,2,2},
	{0,0,2,2,2,2,1,2}, {0,2,2,0,2,2,1,2}, {2,0,0,2,1,2,2,2}, {0,2,2,2,1,1,2,2},
	{2,2,0,2,1,2,2,1}, {2,0,2,2,1,1,2,2}
};

// linear inequalities for CRAFT MixColumns with logic condition.
// let (x0,x1,x2,x3) --> (y0,y1,y2,y3)
// the first one means  "x1 + x2 + (1-x3) + y0 + y1 + y2 >= 1"  
vector<vector<int>> table_lbox = {
	{2,2,1,1,2,0,0,0}, {0,0,2,2,1,1,1,2}, {2,1,1,1,2,2,0,0}, {0,2,0,0,1,2,2,2},
	{2,2,2,1,0,0,2,0}, {2,1,2,1,0,2,2,0}, {1,2,2,1,2,0,2,0}, {0,2,2,2,1,2,1,1},
	{1,1,2,1,2,2,2,0}, {2,2,1,2,0,2,0,2}, {2,0,2,0,2,1,2,2}, {1,2,1,2,2,2,0,2},
	{2,2,2,0,2,2,2,1}, {2,2,0,2,2,2,1,2}, {1,2,2,2,0,2,2,2}, {2,1,2,2,2,0,2,2}
};

// shiftrows
vector<int> perm = { 15,10,9,4,3,6,5,8,7,2,1,12,11,14,13,0 };

// add linear constraints for Sbox to model
void sboxConstr(GRBModel& model, vector<vector<int>> table, vector<GRBVar> inVar, vector<GRBVar> outVar) {
	for (int i = 0; i < table.size(); i++) {
		GRBLinExpr tmp = 0;
		for (int j = 0; j < 4; j++) {
			if (table[i][j] == 0)
				tmp += inVar[j];
			else if (table[i][j] == 1)
				tmp += (1 - inVar[j]);
		}
		for (int j = 0; j < 4; j++) {
			if (table[i][4 + j] == 0)
				tmp += outVar[j];
			else if (table[i][4 + j] == 1)
				tmp += (1 - outVar[j]);
		}
		model.addConstr(tmp >= 1);
	}
}

// set pattern
void setPattern(GRBModel& model, vector<vector<GRBVar>> var, unsigned long long int v) {
	for (int i = 0; i < 64; i++) {
		if (((v >> i) & 1) == 1) {
			model.addConstr(var[15 - i / 4][i % 4] == 1);
		}
		else {
			model.addConstr(var[15 - i / 4][i % 4] == 0);
		}
	}
}

// main function
int main(void) {

	// rounds that we show the lower bound
	int rounds = 14;

	// input/key/output pattern
	unsigned long long int in = 0xEFFFFFFFFFFFFFFF;
	vector<unsigned long long int> key = {
		0x6000000000000000,
		0x0000100000000000,
		0x0000000000000000,
		0x0000000000000000,
		0x0000000000C80008,
		0x0008000C00000080,
		0x489C080CC8C80C88,
		0x88CC0C0CCCC88CC8,
		0x0400C08884000400,
		0x00C8000000000000,
		0x0000000000000000,
		0x0000000000000000,
		0x0000000000000006
	};
	unsigned long long int out = 0x1000000000000000;

	// start gurobi 
	try {
		// create the environment
		GRBEnv env = GRBEnv();

		// enumurate solutions up to 2000000000
		env.set(GRB_IntParam_PoolSearchMode, 2);
		env.set(GRB_IntParam_PoolSolutions, 2000000000);
		env.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);

		// create the model
		GRBModel model = GRBModel(env);

		// create variables
		vector<vector<vector<GRBVar>>> X(rounds, vector<vector<GRBVar>>(16, vector<GRBVar>(4)));
		vector<vector<vector<GRBVar>>> Y(rounds, vector<vector<GRBVar>>(16, vector<GRBVar>(4)));
		vector<vector<vector<GRBVar>>> Z(rounds - 1, vector<vector<GRBVar>>(16, vector<GRBVar>(4)));
		vector<vector<vector<GRBVar>>> W(rounds - 1, vector<vector<GRBVar>>(16, vector<GRBVar>(4)));
		vector<vector<vector<GRBVar>>> K(rounds - 1, vector<vector<GRBVar>>(16, vector<GRBVar>(4)));
		for (int r = 0; r < rounds; r++) {
			for (int i = 0; i < 16; i++) {
				for (int j = 0; j < 4; j++) {
					if (r == 0) {
						X[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
					}
					Y[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
					if (r < rounds - 1) {
						Z[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
						W[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
						K[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
					}
				}
			}
		}

		// set pattern
		setPattern(model, X[0], in);
		setPattern(model, Y[rounds - 1], out);
		for (int r = 0; r < rounds - 1; r++)
			setPattern(model, K[r], key[r]);

		// create constraints
		for (int r = 0; r < rounds - 1; r++) {

			// sbox
			for (int i = 0; i < 16; i++)
				sboxConstr(model, table_sbox, X[r][i], Y[r][i]);

			//MixColumns
			for (int col = 0; col < 4; col++) {
				int i = 4 * col;
				for (int j = 0; j < 4; j++) {
					vector<GRBVar> inVar = { Y[r][i][j], Y[r][i + 1][j], Y[r][i + 2][j], Y[r][i + 3][j] };
					vector<GRBVar> outVar = { Z[r][i][j], Z[r][i + 1][j], Z[r][i + 2][j], Z[r][i + 3][j] };
					sboxConstr(model, table_lbox, inVar, outVar);
				}
			}

			// addkey
			for (int i = 0; i < 16; i++) {
				for (int j = 0; j < 4; j++) {
					model.addConstr(Z[r][i][j] + K[r][i][j] == W[r][i][j]);
				}
			}

			// ShiftRows
			for (int i = 0; i < 16; i++) {
				X[r + 1][i] = W[r][perm[i]];
			}

		}

		// last sbox
		for (int i = 0; i < 16; i++)
			sboxConstr(model, table_sbox, X[rounds - 1][i], Y[rounds - 1][i]);

		// solve this model 
		model.optimize();

		int solcount = model.get(GRB_IntAttr_SolCount);
		cout << "the number of trails : " << solcount << endl;

	}
	catch (GRBException e) {
		cerr << "Error code = " << e.getErrorCode() << endl;
		cerr << e.getMessage() << endl;
	}
	catch (...) {
		cerr << "Exception during optimization" << endl;
	}

	return 0;
}

