#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
#include <ctime>

#include <cryptominisat5/cryptominisat.h>
//#include <boost/program_options.hpp>
//namespace po = boost::program_options;

using KeyPattern = std::vector<std::vector<uint64_t >>;
using KeyPatternList = std::vector<KeyPattern>;


// Parse input file with the key pattern
uint64_t parseFile(const std::string &fileName, KeyPatternList &keyPatternList,
                   uint64_t n, uint64_t rounds){

    std::ifstream file(fileName, std::ios::in);

    if(file.is_open()){
        std::string line;

        while(std::getline(file, line)){
            if(line.empty()){
                file.close();
                return 0;
            }
            else{

                KeyPattern keyPattern;

                for(uint64_t i = 0; i < rounds; i++){
                    std::vector<uint64_t> roundKeyPattern;

                    for(uint64_t k = 0; k < n; k++){
                        roundKeyPattern.push_back(std::stoull(std::string(1, line[k])));
                    }
                    std::getline(file, line);
                    keyPattern.push_back(roundKeyPattern);
                }

                keyPatternList.push_back(keyPattern);
                if(line != "------"){
                    std::cerr << "read file error" << std::endl;
                    file.close();
                    return -1;
                }
            }
        }

    }
    else {
        std::cerr << "Error opening file!" << std::endl;
        return -1;
    }
    return 0;
}


// Model COPY for (a) -> (b_1, b_2)
inline std::string copySAT(uint64_t a, uint64_t b1, uint64_t b2){
    std::string satStr = "-" + std::to_string(a) + " " + std::to_string(b1) + " " + std::to_string(b2) + " 0\n";
    satStr += std::to_string(a) + " -" + std::to_string(b1) + " 0\n";
    satStr += std::to_string(a) + " -" + std::to_string(b2) + " 0\n";

    return satStr;
}


// Model XOR for (a_1, a_2) -> (b)
inline std::string xorSAT(uint64_t a1, uint64_t a2, uint64_t b){
    std::string satStr = std::to_string(a1) + " " + std::to_string(a2) + " -" + std::to_string(b) + " 0\n";
    satStr += std::to_string(b) + " -" + std::to_string(a1) + " 0\n";
    satStr += std::to_string(b) + " -" + std::to_string(a2) + " 0\n";
    satStr += "-" + std::to_string(a1) + " -" + std::to_string(a2) + " 0\n";

    return satStr;
}


// Model AND for (a_1, a_2) -> (b)
inline std::string andSAT(uint64_t a1, uint64_t a2, uint64_t b){
    std::string satStr = std::to_string(b) + " -" + std::to_string(a1) + " 0\n";
    satStr += std::to_string(b) + " -" + std::to_string(a2) + " 0\n";
    satStr += std::to_string(a1) + " -" + std::to_string(b) + " 0\n";
    satStr += std::to_string(a2) + " -" + std::to_string(b) + " 0\n";

    return satStr;
}


// Generate CNF string for a given simon instance with a fixed key pattern
std::string generateSimonSatFormula(uint64_t n, uint64_t rounds, uint64_t shift1, uint64_t shift2, uint64_t shift3,
                                    const KeyPattern &keyPattern,
                                    const std::vector<uint64_t> &inputMonomial,
                                    uint64_t outputBit){

    uint64_t nVariables = 11*n*rounds + 2*n;

    // Clauses for set input, output and key
    uint64_t nClauses = 4*n + n*rounds;
    // Clauses for COPY
    nClauses += 3 * 3 * n * rounds;
    // Clauses for AND
    nClauses += 4 * n * rounds;
    // Clauses for XOR
    nClauses += 3 * 4 * n * rounds;

    // Round offset
    const uint64_t R_OFF = 11 * n;

    std::string cnfString = "p cnf " + std::to_string(nVariables) + " " + std::to_string(nClauses) + "\n";

    // Set input monomial
    for(uint64_t i = 0; i < 2*n; i++){
        if(inputMonomial[i] == 0){
            cnfString += "-" + std::to_string(i+1) + " 0\n";
        }
        else{
            assert(inputMonomial[i] == 1);
            cnfString += std::to_string(i+1) + " 0\n";
        }
    }


    // Propagate the rounds
    for(uint64_t r = 0; r < rounds; r++){
        // First copy
        for(uint64_t i = 1; i <= n; i++){
            uint64_t shiftedPos = r*R_OFF + 2*n + ( (i - shift1 - 1) % n) + 1;
            cnfString += copySAT(r*R_OFF + i, shiftedPos, r*R_OFF + 3*n + i);
        }
        // Second copy
        for(uint64_t i = 1; i <= n; i++){
            uint64_t shiftedPos = r*R_OFF + 4*n + ( (i - shift2 - 1) % n) + 1;
            cnfString += copySAT(r*R_OFF + 3*n + i, shiftedPos, r*R_OFF + 5*n + i);
        }
        // And
        for(uint64_t i = 1; i <= n; i++){
            cnfString += andSAT(r*R_OFF + 2*n + i, r*R_OFF + 4*n + i, r*R_OFF + 6*n + i);
        }
        // First Xor
        for(uint64_t i = 1; i <= n; i++){
            cnfString += xorSAT(r*R_OFF + n + i, r*R_OFF + 6*n + i, r*R_OFF + 7*n + i);
        }
        // Third copy
        for(uint64_t i = 1; i <= n; i++){
            uint64_t shiftedPos = r*R_OFF + 8*n + ( (i - shift3 - 1) % n) + 1;
            cnfString += copySAT(r*R_OFF + 5*n + i, shiftedPos, (r+1)*R_OFF + n + i);
        }
        // Second Xor
        for(uint64_t i = 1; i<= n; i++){
            cnfString += xorSAT(r*R_OFF + 7*n + i, r*R_OFF + 8*n + i, r*R_OFF + 9*n + i);
        }
        // Set Key Pattern
        for(uint64_t i = 0; i < n; i++){
            if(keyPattern[r][i] == 0){
                cnfString += "-" + std::to_string(r*R_OFF + 10*n + i + 1) + " 0\n";
            }
            else{
                assert(keyPattern[r][i] == 1);
                cnfString += std::to_string(r*R_OFF + 10*n + i + 1) + " 0\n";
            }
        }
        // Third Xor
        for(uint64_t i = 1; i<= n; i++){
            cnfString += xorSAT(r*R_OFF + 9*n + i, r*R_OFF + 10*n + i, (r+1)*R_OFF + i);
        }
    }

    // Set output bit
    assert(outputBit < 2*n);
    for(uint64_t i = 0; i < 2*n; i++){
        if(i == outputBit){
            cnfString += std::to_string(rounds * R_OFF + i + 1) + " 0\n";
        }
        else{
            cnfString += "-" + std::to_string(rounds * R_OFF + i + 1) + " 0\n";
        }
    }

    return cnfString;
}


// Convert a CNF string to a cryptominisat (SAT solver) instance
CMSat::SATSolver cnfToSatInstance(const std::string &cnfStr){
    CMSat::SATSolver solver;

    std::string line;
    std::stringstream sstream(cnfStr);

    //std::getline(sstream, line, '\n');
    std::getline(sstream, line, '\n');

    line = line.substr(6);
    std::size_t pos = line.find(' ');

    uint64_t noClauses = std::stoll(line.substr(pos+1));
    uint64_t noLiterals = std::stoll(line.substr(0, pos));

    uint64_t iterations = 0;

    solver.new_vars(noLiterals);

    while(std::getline(sstream, line, '\n')){

        std::istringstream istream(line);
        int64_t literal;

        std::vector<CMSat::Lit> clause;
        while( istream >> literal ) {
            if(literal != 0){

                if(literal > 0){
                    clause.emplace_back(literal-1, false);
                }
                else{
                    clause.emplace_back(-literal-1, true);
                }
            }
        }

        solver.add_clause(clause);
        iterations++;
    }

    if(iterations != noClauses){
        std::cerr << "SAT conversion error" << std::endl;
        exit(-1);
    }

    return solver;
}


// Count the number of solutions for a given cryptominisat instance
// If the time limit is reached, output -1
// If reduced = true, then output the number of solutions modulo 2
int64_t count_solutions(CMSat::SATSolver &sat_solver, bool reduced=false, int timeLimit=0){
    uint64_t solutions = 0;
    CMSat::lbool ret = sat_solver.solve();
    std::vector<CMSat::Lit> ban_solution;

    std::time_t startTime = std::time(nullptr);

    while(ret == CMSat::l_True){
        ban_solution.clear();

        if(std::time(nullptr) - startTime > timeLimit){
            return -1;
        }

        for (uint64_t var = 0; var < sat_solver.nVars(); var++) {
            if (sat_solver.get_model()[var] != CMSat::l_Undef) {
                ban_solution.emplace_back(var, (sat_solver.get_model()[var] == CMSat::l_True )) ;
            }
            else{
                std::cout << "Error " << std::endl;
            }
        }

        sat_solver.add_clause(ban_solution);
        ret = sat_solver.solve();
        solutions++;
    }
    if(reduced) {
        return solutions % 2;
    }
    else{
        return solutions;
    }
}


// Count the number of trails for a given Simon instance with cryptominisat
// If reduced = true, then output the number of solutions modulo 2
int64_t countTrailsSat(uint64_t n, uint64_t shift1, uint64_t shift2, uint64_t shift3, uint64_t rounds,
                       uint64_t outputBit, const std::vector<uint64_t> &inputMonomial,
                       const KeyPattern &keyPattern, bool reduced, int timeLimit){

    auto cnfStr = generateSimonSatFormula(n, rounds, shift1, shift2, shift3, keyPattern, inputMonomial, outputBit);
    auto solver = cnfToSatInstance(cnfStr);
    return count_solutions(solver, reduced, timeLimit);
}


inline uint64_t rotatePosition(uint64_t position, uint64_t n, uint64_t shiftAmount){
    assert(position < 2*n);

    if(position >= n){
        return ((position + shiftAmount) % n) + n;
    }
    else{
        return (position + shiftAmount) % n;
    }
}


std::vector<std::vector<int64_t>> rotateTrailVector(const std::vector<std::vector<int64_t>> &trailVector,
                                                    uint64_t n, uint64_t shiftAmount){

    std::vector<std::vector<int64_t>> rotatedTrailVector(2*n, std::vector<int64_t>(2*n, -5));

    for(uint64_t inputBit = 0; inputBit < 2*n; inputBit++){
        for(uint64_t outputBit = 0; outputBit < 2*n; outputBit++){

            uint64_t newInputPosition = rotatePosition(inputBit, n, shiftAmount);
            uint64_t newOutputPosition = rotatePosition(outputBit, n, shiftAmount);

            rotatedTrailVector[newInputPosition][newOutputPosition] = trailVector[inputBit][outputBit];
        }
    }

    return rotatedTrailVector;
}


void printTrailVector(const std::vector<std::vector<int64_t>> &trailVector, std::ostream &stream){
    stream << "[";

    for(auto const &row : trailVector){
        for(auto const elem: row){
            stream << elem << ", ";
        }
    }
    stream << "]," << std::endl;
}


int main(int argc, char **argv){

    if(argc != 7){
        std::cerr << "Usage: " << argv[0] << " input n rounds threads timelimit cipher" << std::endl;
        return -1;
    }

    std::string inputFile = std::string(argv[1]);
    const uint64_t n = atoll(argv[2]);
    const uint64_t rounds = atoll(argv[3]);
    const uint64_t numberOfThreads = atoll(argv[4]);
    const uint64_t timeLimit = atoll(argv[5]);

    std::string outputFileName = "matrix_" + std::to_string(n) + "_" + std::to_string(rounds)
            + "_" + std::string(argv[6]) + ".py";

    uint64_t shift1, shift2, shift3;

    if(std::string(argv[6]) == "simon"){
        shift1 = 1;
        shift2 = 8;
        shift3 = 2;
    }
    else if(std::string(argv[6]) == "simeck"){
        shift1 = 0;
        shift2 = 5;
        shift3 = 1;
    }
    else{
        std::cerr << "Invalid cipher" << std::endl;
        return -1;
    }

    std::cout << "Verification program" << std::endl;
    std::cout << "n: " << n << std::endl << "rounds: " << rounds << std::endl;
    std::cout << "Number of threads: " << numberOfThreads << std::endl;
    std::cout << "Time limit: " << timeLimit << " seconds" << std::endl;
    std::cout << "Output file: " << outputFileName << std::endl;


    KeyPatternList keyPatternList;

    if(parseFile(inputFile, keyPatternList, n, rounds) == 0){

        std::ofstream outputFile(outputFileName);

        outputFile << "arr = [";

        // Compute the big matrix
        #pragma omp parallel for num_threads(numberOfThreads) schedule(dynamic)\
            shared(n, rounds, shift1, shift2, shift3, keyPatternList, outputFile, timeLimit) default(none)

        for (auto const &keyPattern : keyPatternList){

            // Compute matrix row
            std::vector<std::vector<int64_t>> trailVector(2*n, std::vector<int64_t>(2*n, -5));

            for(uint64_t inputBit = 0; inputBit < 2*n; inputBit++){
                std::vector<uint64_t> inputMonomial(2*n, 1);
                inputMonomial[inputBit] = 0;

                for(uint64_t outputBit = 0; outputBit < 2*n; outputBit++){

                    int64_t trails = countTrailsSat(n, shift1, shift2, shift3, rounds, outputBit, inputMonomial,
                                                    keyPattern, true, timeLimit);

                    trailVector[inputBit][outputBit] = trails;
                }

            }

            // Shift row because of key invariance property
            for(uint64_t shiftAmount = 0; shiftAmount < n; shiftAmount++){
                auto rotatedTrailVector = rotateTrailVector(trailVector, n, shiftAmount);
                #pragma omp critical
                {
                    printTrailVector(rotatedTrailVector, outputFile);
                }
            }
        }

        outputFile << "]" << std::endl;
    }

    return 0;
}