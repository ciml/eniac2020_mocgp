#include <bdd.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>

class Circuit{
    public:
        int num_inputs;                                 /**< Number of inputs */
        int num_outputs;                                /**< Number of outputs */
        int num_gates;                                  /**< Number of gates on seeding population */
        std::ifstream file;                             /**< Circuits input file */
        std::vector<bdd> outputs;                       /**< Vector of outputs from ESPRESSO's circuit used as comparison */
        std::vector<std::string> boolean_expression;    /**< Vector of boolean expression from ESPRESSO's input file */
        
        /**
		* @brief Constructor
		* @param filename - The circuit definition file
		* @return none
		*/
        Circuit(const char *filename);

        /**
		* @brief Destructor
		* @return none
		*/
        ~Circuit();

        /**
		* @brief Check if the input variable is complemented or not  and return the 
        * respective bdd variable
		* @param buffer - Boolean expression containing an input
		* @return The respective bdd variable complemented or not
		*/
        bdd parse_variable(std::string buffer);

        /**
		* @brief Breaks the buffer in the '*' character, separate it on strings, 
        * send each string to parse_variable, create the bdd according to the 
        * expression and return the resultant bdd.
		* @param buffer - Boolean expression containing only AND functions
		* @return The bdd created according to the boolean expression received as input
		*/
        bdd parse_and(std::string buffer);

        /**
		* @brief Breaks the buffer in the '+' character, separate it on strings, 
        * send each string to parse_and function, create the bdd according to the expression 
        * and return the resultant bdd.
		* @param buffer - Boolean expression containing the complete boolean expression
		* @return The bdd created according to the boolean expression received as input
		*/
        bdd parse_or(std::string buffer);

        /**
		* @brief Read the ESPRESSO boolean expressions from the input file, send it to 
        * parse_or function and stores the bdd returned.
		* @return none
		*/
        void make_circuit_bdd();

        void update_base_circuit(std::vector<bdd> new_outputs);

        /**
		* @brief Print circuit information.
		* @param filename - The circuit definition file
		* @return none
		*/
        void print_infos(const char *filename);
};

Circuit::Circuit(const char *filename){
    std::string buffer;
    this->file.open(filename);

    if(!this->file.is_open()){
        std::cout << "File " << filename << " cannot be opened!" << std::endl;
        exit(0);
    }

    while(std::getline(this->file, buffer)){
        if(buffer.find(".p") != std::string::npos){
            sscanf(buffer.c_str(), ".p %d", &this->num_gates);
        }
        else if(buffer.find(".i") != std::string::npos){
            sscanf(buffer.c_str(), ".i %d", &this->num_inputs);
        }
        else if(buffer.find(".o") != std::string::npos){
            sscanf(buffer.c_str(), ".o %d", &this->num_outputs);
            break;
        }
    }

    /**< Setup Binary Decision Diagram>*/
    bdd_setvarnum(this->num_inputs);
    this->outputs = std::vector<bdd>(this->num_outputs);

    print_infos(filename);
}

Circuit::~Circuit(){
    std::vector<bdd>().swap(this->outputs);
    std::vector<std::string>().swap(this->boolean_expression);
}

bdd Circuit::parse_variable(std::string buffer){
    int var = -1;
    if(buffer.find("~") != std::string::npos){
        sscanf(buffer.c_str(), "~i%d", &var);
        return bdd_not(bdd_ithvar(var));
    }
    else
    {
        sscanf(buffer.c_str(), "i%d", &var);
        return bdd_ithvar(var);
    }
}

bdd Circuit::parse_and(std::string buffer){
    std::vector<std::string> splited_expression;
    std::string temp;
    bdd and_exp;

    /*< Replace + with space >*/
    std::replace(buffer.begin(), buffer.end(), '*', ' ');
    std::stringstream str(buffer);

    /*< Split string in space character and save on vector >*/
    while (str >> temp) splited_expression.push_back(temp);

    and_exp = parse_variable(splited_expression[0]);
    for (int i = 1; i < (int)splited_expression.size(); i++)
    {
        and_exp = bdd_and(and_exp, parse_variable(splited_expression[i]));
    }
    return and_exp;
}

bdd Circuit::parse_or(std::string buffer){
    std::vector<std::string> splited_expression;
    std::string temp;
    bdd or_exp;

    /*< Replace * with space >*/
    std::replace(buffer.begin(), buffer.end(), '+', ' ');
    std::stringstream str(buffer);

    /*< Split string in space character and save on vector >*/
    while (str >> temp) splited_expression.push_back(temp);

    or_exp = parse_and(splited_expression[0]);
    for(int i = 1; i < (int)splited_expression.size(); i++){
        or_exp = bdd_or(or_exp, parse_and(splited_expression[i]));
    }
    return or_exp;
}

void Circuit::make_circuit_bdd(){
    int counter = 0;
    std::string buffer;
    
    while (std::getline(this->file, buffer)){
        this->boolean_expression.push_back(buffer);
        this->outputs[counter] = parse_or(buffer);
        counter++;
    }
    file.close();
}

void Circuit::update_base_circuit(std::vector<bdd> new_outputs){
    for(int i = 0 ; i < this->num_outputs; i++){
        this->outputs[i]= new_outputs[i];
    }
}

void Circuit::print_infos(const char *filename){
    std::stringstream file(filename);
    std::string temp;
    std::vector<std::string> list;
    while(std::getline(file, temp, '/')){
        list.push_back(temp);
    }
    std::cout << "Filename: " << list[list.size() - 1] << std::endl;
    std::cout << "Number of Inputs: " << this->num_inputs << "\tNumber of Outputs: " 
    << this->num_outputs << std::endl;
}