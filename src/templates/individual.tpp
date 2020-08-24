#include <time.h>
#include <bdd.h>
#include <cstring>
#include <math.h>
#include <iterator>
#include "node.tpp"
#include "circuit.tpp"

#define MR 0.05

#define AND 1
#define OR 2
#define NOT 3
#define NAND 4
#define NOR 5
#define XOR 6
#define XNOR 7

int NCOL;
int LB;
std::vector<int> functions = {AND, OR, XOR, NOT, NAND, NOR, XNOR};
int num_functions = 7;

#ifdef MO
///Propagation delay(ns)
#define tdAND 1.7
#define tdNAND 1.8
#define tdOR 1.7
#define tdNOR 1.7
#define tdXOR 1.9
#define tdXNOR 72.5
#define tdNOT 4.5

///Operational parameters
#define freq 50000000
#define Capa_Load 0.00000000005
#define Vcc 5
#endif


template<typename S, typename T>
class Individual{
    public:
        std::vector<Node<T>> genotype;      /**< Vector of nodes */
        std::vector<int> outputs;           /**< Vector of outputs */
        std::vector<S> outputs_error;       /**< Vector of error for each output */
        int num_outputs;                    /**< Number of outputs */
        int transistors;                    /**< Number of transistors */
        S error;                            /**< Total error */
        #ifdef MO
            double power;                   /**< Total pwoer */
            double delay;                   /**< Maximum delay */
            int rank;                       /**< Pareto's Front Rank */
            double mre;
        #endif

        /**
		* @brief Constructor
		* @param num_outputs - The circuit's number of outputs
		* @return None
		*/
        Individual(int num_outputs);
        
        /**
		* @brief Destructor
		* @return None
		*/
        ~Individual();

        /**
		* @brief Initialize individual genotype
		* @param circuit - The circuit
		* @return None
		*/
        void initialize(Circuit *circuit);
        
        /**
		* @brief Print individual information
		* @param circuit - The circuit
		* @return none
		*/
        void print(Circuit *circuit);

        /**
		* @brief Set the active nodes at an output path
		* @param circuit - The circuit
		* @param node - The actual node
		* @return none
		*/
        void set_active_nodes_output(Circuit *circuit, int node);

        /**
		* @brief Set the active nodes at all outputs path
		* @param circuit - The circuit
		* @return none
		*/
        void set_active_nodes(Circuit *circuit);

        /**
		* @brief Set all nodes as inactive
		* @return none
		*/
        void clear_active_nodes();

        /**
		* @brief Mutate an output gene
		* @param circuit - The circuit
		* @param sorted_node - The sorted output node to be mutated
		* @return none
		*/
        void mutate_output(Circuit *circuit, int sorted_node);

        /**
		* @brief Mutate a regular node
		* @param circuit - The circuit
		* @param sorted_node - The sorted node to be mutated
		* @return none
		*/
        void mutate_node(Circuit *circuit, int sorted_node);

        /**
		* @brief Mutate a node
		* @param circuit - The circuit
		* @return true if the node is active and false if not
		*/
        bool mutate_individual(Circuit *circuit);

        /**
		* @brief Apply the Point Mutation
		* @param circuit - The circuit
		* @return none
		*/
        void apply_PM(Circuit *circuit);

        /**
		* @brief Apply the Single Active Mutation
		* @param circuit - The circuit
		* @return none
		*/
        void apply_SAM(Circuit *circuit);

        /**
		* @brief Apply the Guided Active Mutation
		* @param circuit - The circuit
		* @return none
		*/
        void apply_GAM(Circuit *circuit);

        /**
		* @brief Return the number of transistors according to the logic gate
		* @param function - The logic gate
		* @return The number of transistors according to the logic gate
		*/
        int get_num_transistors(int function);

        /**
		* @brief Count the number of transistors
		* @return none
		*/
        void count_transistors();

        void rand_nodes_not_used(Circuit *circuit, int start);
        int get_next_available_position(Circuit *circuit, int start);
        
        void set_nots(Circuit *circuit);
        int parse_pla_input(const char *str, Circuit *circuit);
        int parse_only_ands(const char *str, Circuit *circuit);
        int parse_only_ors(const char *str, Circuit *circuit);
        int parse_ands_ors(const char *str, Circuit *circuit);        
        void seed_pla(Circuit *circuit);
        int parse_evoapprox8b_input(Circuit *circuit, std::string buffer, std::vector<std::string> temp);
        void parse_evoapprox8b_output(Circuit *circuit, std::string buffer, std::vector<std::string> temp);
        void seed_evoapprox8b(Circuit *circuit, const char *filename);
        void seed(Circuit *circuit, std::string filename); 

        /**
		* @brief Create and return the bdd according to the logic function of the node
        and its inputs
		* @param left - The bdd of the left input
		* @param function - The logic gate
		* @param circuit - The bdd of the right input
		* @return The resulting bdd
		*/
        bdd get_function_output(bdd left, int function, bdd right);

        /**
		* @brief Create the bdd of each active node according to the logic function of the node
        and its inputs
		* @param circuit - The circuit
		* @param node - The actual node 
		* @return The node output bdd
		*/
        bdd make_bdd_per_output(Circuit *circuit, int node);

        /**
		* @brief Calculate the maximum delay at each active node
		* @param circuit - The circuit
		* @param node - The actual node 
		* @return The maximum delay at that node 
		*/
        double eval_delay(Circuit *circuit, int node);

        /**
		* @brief Calculate the power at each active node
		* @param circuit - The circuit
		* @param node - The actual node 
		* @return The power at that node 
		*/
        double eval_power(Circuit *circuit, int node);

        /**
		* @brief Calculate individual's sat count
		* @param circuit - The circuit
		* @return none
		*/
        void evaluate_sat_count(Circuit *circuit);
        #ifdef MO
        /**
		* @brief Calculate individual's maximum delay
		* @param circuit - The circuit
		* @return none
		*/
        void evaluate_delay(Circuit *circuit);

        /**
		* @brief Calculate individual's total power
		* @param circuit - The circuit
		* @return none
		*/
        void evaluate_power(Circuit *circuit);

        /**
		* @brief Print individual's information
		* @return none 
		*/
        void print_MO();
        #endif
};

template<typename S, typename T>
Individual<S, T>::Individual(int num_outputs){
    this->num_outputs = num_outputs;
    this->error = -1;
    this->transistors = 0;

    #ifdef MO
    this->power = 1000.0;
    this->delay = 1000.0;
    this->rank = 10;
    this->mre = 100.0;
    #endif

    this->genotype = std::vector<Node<T>>(NCOL, Node<T>());
    this->outputs = std::vector<int>(this->num_outputs, -1);
    this->outputs_error = std::vector<S>(this->num_outputs, -1);
}

template<typename S, typename T>
Individual<S, T>::~Individual(){
    // for(int i = 0; i < NCOL; i++)this->genotype[i].~Node();
    std::vector<Node<T>>().swap(this->genotype);
    std::vector<int>().swap(this->outputs);
    std::vector<S>().swap(this->outputs_error);
}

template<typename S, typename T>
void Individual<S, T>::initialize(Circuit *circuit){
    /*define individual genotype*/
    int temp;
    for(int i = 0; i < NCOL; i++){
        for(int j = 0; j < INPUTS; j++){
            while (1){
                temp = rand() % (circuit->num_inputs + i);
                if (temp < circuit->num_inputs || temp >= (i + circuit->num_inputs - LB)){
                    this->genotype[i].inputs[j] = temp;
                    break;
                }
            }
        }
        temp = rand()%num_functions;
        this->genotype[i].function = functions[temp];
    }

    /*define individual outputs*/
    for(int i = 0; i < this->num_outputs; i++){
        this->outputs[i] = rand()%(circuit->num_inputs + NCOL);
    }
}

template<typename S, typename T>
void Individual<S, T>::print(Circuit *circuit){
    set_active_nodes(circuit);
    std::cout << "GENOTYPE" << std::endl;
    for(int i = 0; i < NCOL; i+= 1){
        std::cout << i+circuit->num_inputs <<": ["; 
        for(int j = 0; j < INPUTS; j++){
            std::cout << this->genotype[i].inputs[j] << ", ";
        }
        std::cout << this->genotype[i].function << "] " << this->genotype[i].active << std::endl;
    }

    std::cout << "OUTPUTS:" << std::endl;
    for(int i = 0; i < this->num_outputs; i++){
        std::cout << this->outputs[i] << "\t";
    }
    std::cout << std::endl << std::endl;
    clear_active_nodes();
}

template<typename S, typename T>
void Individual<S, T>::set_active_nodes_output(Circuit *circuit, int node){
    if(node >= circuit->num_inputs){
        int pos = node - circuit->num_inputs;
        for(int i = 0; i < INPUTS; i++){
            set_active_nodes_output(circuit, this->genotype[pos].inputs[i]);
        }
        this->genotype[pos].active = true;
    }
}

template<typename S, typename T>
void Individual<S, T>::set_active_nodes(Circuit *circuit){
    for(int i = 0; i < this->num_outputs; i++){
        set_active_nodes_output(circuit, this->outputs[i]);
    }
}

template<typename S, typename T>
void Individual<S, T>::clear_active_nodes(){
    for(int i = 0; i < NCOL; i++){
        this->genotype[i].active = false;
    }
}

template<typename S, typename T>
void Individual<S, T>::mutate_output(Circuit *circuit, int sorted_node){
    int temp = 0;

    sorted_node = sorted_node - NCOL;

    if(sorted_node < 0 || sorted_node >= (int)this->num_outputs){
        std::cout << "Invalid value sorted as output!" << std::endl;
        exit(0);
    }
    
    while (1)
    {
        temp = rand()%(circuit->num_inputs + NCOL);
        if(temp != this->outputs[sorted_node]){
            this->outputs[sorted_node] = temp;
            break;
        }
    }
}

template<typename S, typename T>
void Individual<S, T>::mutate_node(Circuit *circuit, int sorted_node){
    int temp = rand() % (INPUTS + 1);

    if(temp < INPUTS){
        int input = temp;
        while (1){
            temp = rand() % (circuit->num_inputs + sorted_node);

            if ((temp < circuit->num_inputs || temp >= (sorted_node + circuit->num_inputs - LB)) && (temp != this->genotype[sorted_node].inputs[input]))
            {
                this->genotype[sorted_node].inputs[input] = temp;
                break;
            }
        }
    }
    else if (temp == INPUTS){
        while (1)
        {
            temp = rand()%num_functions;
            if (functions[temp] != this->genotype[sorted_node].function){
                this->genotype[sorted_node].function = functions[temp];
                break;
            }
        }
        
    }
    else{
        std::cout << "Sorted value when mutating node isn't valid!" << std::endl;
        exit(0);
    }
}

template<typename S, typename T>
bool Individual<S, T>::mutate_individual(Circuit *circuit){
    int sorted_node = rand() % (NCOL + this->num_outputs);
    
    if(sorted_node >= 0 && sorted_node < NCOL){
        mutate_node(circuit, sorted_node);
    }
    else if(sorted_node >= NCOL && sorted_node < (NCOL + this->num_outputs)){
        mutate_output(circuit, sorted_node);
        return true;
    }
    else{
        std::cout << "Sorted value isn't valid!" << std::endl; 
    }

    if(this->genotype[sorted_node].active) return true;

    return false;
}

template<typename S, typename T>
void Individual<S, T>::apply_PM(Circuit *circuit){
    for(int i = 0; i < (int)(MR*NCOL); i++){
        mutate_individual(circuit);
    }
}

template<typename S, typename T>
void Individual<S, T>::apply_SAM(Circuit *circuit){
    set_active_nodes(circuit);
    while (!mutate_individual(circuit));
    clear_active_nodes();
}

template<typename S, typename T>
void Individual<S, T>::apply_GAM(Circuit *circuit){
    std::vector<int> temparray;
    int worst_output = 0;
    temparray.push_back(0);

    for(int i = 1; i < this->num_outputs; i++){
        if(this->outputs_error[i] > this->outputs_error[worst_output]){
            worst_output = i;
            temparray.clear();
            temparray.push_back(i);
        }
        else if(this->outputs_error[i] == this->outputs_error[worst_output]){
            temparray.push_back(i);
        }
    }

    if((int)temparray.size() == 1) worst_output = temparray[0];
    else if((int)temparray.size() > 1){
        worst_output = temparray[rand() % (int)temparray.size()];
    }
    else{
        std::cout << "Array size is wrong!" << std::endl;
    }

    set_active_nodes_output(circuit, this->outputs[worst_output]);
    int temp;
    while (1){
        temp = rand() % (NCOL + this->num_outputs);
        
        if(temp >= NCOL){
            mutate_output(circuit, temp);
            break;
        }
        else if(this->genotype[temp].active == true){
            mutate_node(circuit, temp);
            break;
        }
    }
    clear_active_nodes();
}

template<typename S, typename T>
int Individual<S, T>::get_num_transistors(int function){
    switch(function){
    case 1: //AND
        return 2;
        break;
    case 2: //OR
        return 2;
        break;
    case 3: //NOT
        return 1;
        break;
    case 4: //NAND
        return 2;
        break;
    case 5: //NOR
        return 1;
        break;
    case 6: //XOR
        return 3;
        break;
    case 7: //XNOR
        return 4;
        break;
    default:
        std::cout << "Gate code unknow!\n";
        exit(0);
        break;
    }
}

template<typename S, typename T>
void Individual<S, T>::count_transistors(){
    this->transistors = 0;
    for(int i = 0; i < NCOL; i++){
        if(this->genotype[i].active) 
            this->transistors += get_num_transistors(this->genotype[i].function);
    }
}

template<typename S, typename T>
void Individual<S, T>::rand_nodes_not_used(Circuit *circuit, int start){
    int nodes_not_used = NCOL - circuit->num_gates - start;
    int temp;

    if(nodes_not_used < 0){
        std::cout << "Number of columns is lower than needed!" << std::endl;
        exit(0);
    }

    int count = 0;
    while (1){
        temp = start + rand() % (NCOL - start);

        if(this->genotype[temp].active == false){
            this->genotype[temp].active = true;
            count += 1;
        }

        if(count == nodes_not_used) break;
    }
    
}

template<typename S, typename T>
void Individual<S, T>::set_nots(Circuit *circuit){
    for(int i = 0; i < circuit->num_inputs; i++){
        this->genotype[i].function = NOT;
        this->genotype[i].inputs[0] = i;
        this->genotype[i].active = true;
    }
}

template<typename S, typename T>
int Individual<S, T>::get_next_available_position(Circuit *circuit, int start){
    for(int i = start; i < NCOL; i++){
        if(this->genotype[i].active == 0) return i;
    }

    std::cout << "Couldn't find any free space on genotype!" << std::endl;
    exit(0);

    return 0;
}

template<typename S, typename T>
int Individual<S, T>::parse_pla_input(const char *str, Circuit *circuit){
    int pos = 0;
    
    if(strstr(str, "~") != NULL){
        sscanf(str, "~i%d", &pos);
        return pos + circuit->num_inputs;
    }
    else{
        sscanf(str, "i%d", &pos);
        return pos;
    }
}

template<typename S, typename T>
int Individual<S, T>::parse_only_ands(const char *str, Circuit *circuit){
    char *token;
    char *temp;
    char *saveptr;
    int flag = 0;
    int input0 = 0;
    int input1 = 0;
    int pos = 0;

    if(strstr(str, "*") == NULL) return parse_pla_input(str, circuit);

    for (temp = (char *)str;; temp = NULL){
        token = strtok_r(temp, "*", &saveptr);
        if (token == NULL) break;
        if (flag == 0){
            input0 = parse_pla_input(token, circuit);
            flag = 1;
        }
        else if(flag == 1){
            input1 = parse_pla_input(token, circuit);
            pos = get_next_available_position(circuit, circuit->num_inputs);

            this->genotype[pos].function = AND;
            this->genotype[pos].inputs[0] = input0;
            this->genotype[pos].inputs[1] = input1;
            this->genotype[pos].active = true;

            input1 = pos + circuit->num_inputs;
            flag = 2;
        }
        else if(flag == 2){
            input0 = parse_pla_input(token, circuit);
            pos = get_next_available_position(circuit, circuit->num_inputs);

            this->genotype[pos].function = AND;
            this->genotype[pos].inputs[0] = input0;
            this->genotype[pos].inputs[1] = input1;
            this->genotype[pos].active = true;

            input1 = pos + circuit->num_inputs;
        }       
    }
    return pos + circuit->num_inputs;
}

template<typename S, typename T>
int Individual<S, T>::parse_only_ors(const char *str, Circuit *circuit){
    char *token;
    char *temp;
    char *saveptr;
    int flag = 0;
    int input0 = 0;
    int input1 = 0;
    int pos = 0;

    if (strstr(str, "+") == NULL) return parse_pla_input(str, circuit);

    for (temp = (char *)str;; temp = NULL){
        token = strtok_r(temp, "+", &saveptr);
        if (token == NULL) break;
        if (flag == 0){
            input0 = parse_pla_input(token, circuit);
            flag = 1;
        }
        else if (flag == 1){
            input1 = parse_pla_input(token, circuit);
            pos = get_next_available_position(circuit, circuit->num_inputs);

            this->genotype[pos].function = OR;
            this->genotype[pos].inputs[0] = input0;
            this->genotype[pos].inputs[1] = input1;
            this->genotype[pos].active = true;

            input1 = pos + circuit->num_inputs;
            flag = 2;
        }
        else if (flag == 2){
            input0 = parse_pla_input(token, circuit);
            pos = get_next_available_position(circuit, circuit->num_inputs);

            this->genotype[pos].function = OR;
            this->genotype[pos].inputs[0] = input0;
            this->genotype[pos].inputs[1] = input1;
            this->genotype[pos].active = true;

            input1 = pos + circuit->num_inputs;
        }
    }
    return pos + circuit->num_inputs;
}

template<typename S, typename T>
int Individual<S, T>::parse_ands_ors(const char *str, Circuit *circuit){
    char *token;
    char *temp;
    char *saveptr;
    int flag = 0;
    int input0 = 0;
    int input1 = 0;
    int pos = 0;

    for (temp = (char *)str;; temp = NULL){
        token = strtok_r(temp, "+", &saveptr);
        if (token == NULL) break;
        if (flag == 0){
            input0 = parse_only_ands(token, circuit);
            flag = 1;
        }
        else if (flag == 1){
            input1 = parse_only_ands(token, circuit);
            pos = get_next_available_position(circuit, circuit->num_inputs);

            this->genotype[pos].function = OR;
            this->genotype[pos].inputs[0] = input0;
            this->genotype[pos].inputs[1] = input1;
            this->genotype[pos].active = true;

            input1 = pos + circuit->num_inputs;
            flag = 2;
        }
        else if (flag == 2){
            input0 = parse_only_ands(token, circuit);
            pos = get_next_available_position(circuit, circuit->num_inputs);

            this->genotype[pos].function = OR;
            this->genotype[pos].inputs[0] = input0;
            this->genotype[pos].inputs[1] = input1;
            this->genotype[pos].active = true;

            input1 = pos + circuit->num_inputs;
        }
    }
    return pos + circuit->num_inputs;
}

template<typename S, typename T>
void Individual<S, T>::seed_pla(Circuit *circuit){

    rand_nodes_not_used(circuit, circuit->num_inputs);
    set_nots(circuit);


    for (int i = 0; i < circuit->num_outputs; i++)
    {
        if (strstr(circuit->boolean_expression[i].c_str(), "*") == NULL && strstr(circuit->boolean_expression[i].c_str(), "+") == NULL)
        {
            this->outputs[i] = parse_pla_input(circuit->boolean_expression[i].c_str(), circuit);
        }
        else if(strstr(circuit->boolean_expression[i].c_str(), "+") == NULL)
        {
            this->outputs[i] = parse_only_ands(circuit->boolean_expression[i].c_str(), circuit);
        }
        else if(strstr(circuit->boolean_expression[i].c_str(), "*") == NULL)
        {
            this->outputs[i] = parse_only_ors(circuit->boolean_expression[i].c_str(), circuit);
        }
        else if (strstr(circuit->boolean_expression[i].c_str(), "*") != NULL && strstr(circuit->boolean_expression[i].c_str(), "+") != NULL)
        {
            this->outputs[i] = parse_ands_ors(circuit->boolean_expression[i].c_str(), circuit);
        }
    }
}

template<typename S, typename T>
int Individual<S, T>::parse_evoapprox8b_input(Circuit *circuit, std::string buffer, std::vector<std::string> temp){
    std::vector<std::string>::iterator it = std::find(temp.begin(), temp.end(), buffer);
    if (it != temp.end()) return std::distance(temp.begin(), it);
    else{
        std::cout << "EvoApprox8b input not recognized: " << buffer << std::endl;
        exit(0);
    }
}

template<typename S, typename T>
void Individual<S, T>::parse_evoapprox8b_output(Circuit *circuit, std::string buffer, std::vector<std::string> temp){
    int node = 0, output = 0;
    
    buffer.erase(std::remove(buffer.begin(), buffer.end(), ';'), buffer.end());
    buffer.erase(std::remove(buffer.begin(), buffer.end(), '('), buffer.end());
    buffer.erase(std::remove(buffer.begin(), buffer.end(), ')'), buffer.end());
    std::istringstream str(buffer);
    std::vector<std::string> splited_buffer{std::istream_iterator<std::string>(str), {}};

    std::vector<std::string>::iterator it = std::find(temp.begin(), temp.end(), splited_buffer[2]);
    if (it != temp.end()) node = std::distance(temp.begin(), it);
    else{
        std::cout << "EvoApprox8b node not recognized: " << buffer << std::endl;
        exit(0);
    }
    output = circuit->num_outputs - 1 - std::stoi(splited_buffer[6]);
    this->outputs[output] = node;
}

template<typename S, typename T>
void Individual<S, T>::seed_evoapprox8b(Circuit *circuit, const char *filename){
    std::ifstream file;
    std::string buffer;
    std::vector<std::string> temp(NCOL + circuit->num_inputs, "");
    int xor1 = 0;
    
    file.open(filename);
    if(!file.is_open()){
        std::cout << "File cannot be opened!" << std::endl;
        exit(0);
    }

    while(std::getline(file, buffer)){
        if(buffer.find("Nodes") != std::string::npos){
            sscanf(buffer.c_str(), "///  Nodes = %d", &circuit->num_gates);
            rand_nodes_not_used(circuit, 0);
        }
        else if(buffer.find("c |=") != std::string::npos){
            parse_evoapprox8b_output(circuit, buffer, temp);
        }
        else if(buffer.find(">>") != std::string::npos){
            int input = 0;

            buffer.erase(std::remove(buffer.begin(), buffer.end(), '('), buffer.end());
            buffer.erase(std::remove(buffer.begin(), buffer.end(), ')'), buffer.end());
            
            std::istringstream str(buffer);
            std::vector<std::string> splited_buffer{std::istream_iterator<std::string>(str), {}};

            if(splited_buffer[3].find("a") != std::string::npos){
                input = circuit->num_inputs/2 - 1;
                input -= atoi(splited_buffer[5].c_str());
            }
            else if(splited_buffer[3].find("b") != std::string::npos){
                input = circuit->num_inputs - 1;
                input -= atoi(splited_buffer[5].c_str());
            }
            temp[input] = splited_buffer[1];
        }
        else if((buffer.find("uint8_t") == std::string::npos) && 
        (buffer.find("&") != std::string::npos || buffer.find("|") != std::string::npos 
        || buffer.find("^") != std::string::npos || buffer.find("~") != std::string::npos)){
            int num_and = 0, num_or = 0, num_xor = 0, num_not = 0;

            num_and = std::count(buffer.begin(), buffer.end(), '&');
            num_or  = std::count(buffer.begin(), buffer.end(), '|');
            num_xor = std::count(buffer.begin(), buffer.end(), '^');
            num_not = std::count(buffer.begin(), buffer.end(), '~');

            buffer.erase(std::remove(buffer.begin(), buffer.end(), ';'), buffer.end());
            buffer.erase(std::remove(buffer.begin(), buffer.end(), '('), buffer.end());
            buffer.erase(std::remove(buffer.begin(), buffer.end(), ')'), buffer.end());
            buffer.erase(std::remove(buffer.begin(), buffer.end(), '~'), buffer.end());

            std::istringstream str(buffer);
            std::vector<std::string> splited_buffer{std::istream_iterator<std::string>(str), {}};

            if(num_and == 2 && num_or == 1 && num_not == 1){
                std::vector<int> input;
                std::vector<int> node;

                for(int i = 2; i < 9; i+=2) input.push_back(parse_evoapprox8b_input(circuit, splited_buffer[i], temp));
                for(int i = 0; i < 4; i++){
                    node.push_back(get_next_available_position(circuit, 0));
                    this->genotype[node[i]].active = true;
                }

                this->genotype[node[0]].inputs[0] = input[0];
                this->genotype[node[0]].inputs[1] = input[1];
                this->genotype[node[0]].function = AND;
                
                this->genotype[node[1]].inputs[0] = input[2];
                this->genotype[node[1]].function = NOT;

                this->genotype[node[2]].inputs[0] = node[1] + circuit->num_inputs;
                this->genotype[node[2]].inputs[1] = input[3];
                this->genotype[node[2]].function = AND;

                this->genotype[node[3]].inputs[0] = node[0] + circuit->num_inputs;
                this->genotype[node[3]].inputs[1] = node[2] + circuit->num_inputs;
                this->genotype[node[3]].function = OR;
                temp[node[3] + circuit->num_inputs] = splited_buffer[0];
            }
            else if(num_and == 3 && num_or == 2){
                int input1 = 0, input2 = 0, node = 0;
                std::vector<int> nodes;

                for(int i = 2; i < 13; i+=4){
                    input1 = parse_evoapprox8b_input(circuit, splited_buffer[i], temp);
                    input2 = parse_evoapprox8b_input(circuit, splited_buffer[i+2], temp);
                    node = get_next_available_position(circuit, 0);
                    this->genotype[node].inputs[0] = input1;
                    this->genotype[node].inputs[1] = input2;
                    this->genotype[node].function = AND;
                    this->genotype[node].active = true;
                    nodes.push_back(node);
                }

                for(int i = 0; i < 4; i+=2){
                    input1 = nodes[i] + circuit->num_inputs;
                    input2 = nodes[i+1] + circuit->num_inputs;
                    node = get_next_available_position(circuit, 0);
                    this->genotype[node].inputs[0] = input1;
                    this->genotype[node].inputs[1] = input2;
                    this->genotype[node].function = OR;
                    this->genotype[node].active = true;
                    nodes.push_back(node);
                }
                temp[nodes[4] + circuit->num_inputs] = splited_buffer[0];
            }
            else if(num_and == 2 && num_xor == 2){
                std::vector<int> input;
                std::vector<int> node;

                input.push_back(parse_evoapprox8b_input(circuit, splited_buffer[2], temp));
                input.push_back(parse_evoapprox8b_input(circuit, splited_buffer[4], temp));
                input.push_back(parse_evoapprox8b_input(circuit, splited_buffer[10], temp));
            
                for(int i = 0; i < 3; i++){
                    node.push_back(get_next_available_position(circuit, 0));
                    this->genotype[node[i]].active = true;
                }

                this->genotype[node[0]].inputs[0] = input[0];
                this->genotype[node[0]].inputs[1] = input[1];
                this->genotype[node[0]].function = AND;

                this->genotype[node[1]].inputs[0] = input[2];
                this->genotype[node[1]].inputs[1] = xor1 + circuit->num_inputs;
                this->genotype[node[1]].function = AND;

                this->genotype[node[2]].inputs[0] = node[0] + circuit->num_inputs;
                this->genotype[node[2]].inputs[1] = node[1] + circuit->num_inputs;
                this->genotype[node[2]].function = XOR;
                temp[node[2] + circuit->num_inputs] = splited_buffer[0];
            
            }
            else if(num_and == 2 && num_xor == 1 && num_or == 1){
                std::vector<int> input;
                std::vector<int> node;

                input.push_back(parse_evoapprox8b_input(circuit, splited_buffer[2], temp));
                input.push_back(parse_evoapprox8b_input(circuit, splited_buffer[4], temp));
                input.push_back(parse_evoapprox8b_input(circuit, splited_buffer[10], temp));
            
                for(int i = 0; i < 3; i++){
                    node.push_back(get_next_available_position(circuit, 0));
                    this->genotype[node[i]].active = true;
                }

                this->genotype[node[0]].inputs[0] = input[0];
                this->genotype[node[0]].inputs[1] = input[1];
                this->genotype[node[0]].function = AND;

                this->genotype[node[1]].inputs[0] = input[2];
                this->genotype[node[1]].inputs[1] = xor1 + circuit->num_inputs;
                this->genotype[node[1]].function = AND;

                this->genotype[node[2]].inputs[0] = node[0] + circuit->num_inputs;
                this->genotype[node[2]].inputs[1] = node[1] + circuit->num_inputs;
                this->genotype[node[2]].function = OR;
                temp[node[2] + circuit->num_inputs] = splited_buffer[0];
            
            }
            else if(num_xor == 2){
                int input1 = 0, input2 = 0, node = 0;
                
                input1 = parse_evoapprox8b_input(circuit, splited_buffer[2], temp);
                input2 = parse_evoapprox8b_input(circuit, splited_buffer[4], temp);
                node = get_next_available_position(circuit, 0);
                this->genotype[node].inputs[0] = input1;
                this->genotype[node].inputs[1] = input2;
                this->genotype[node].function = XOR;
                this->genotype[node].active = true;
                xor1 = node;

                input1 = node + circuit->num_inputs;
                input2 = parse_evoapprox8b_input(circuit, splited_buffer[6], temp);
                node = get_next_available_position(circuit, 0);
                this->genotype[node].inputs[0] = input1;
                this->genotype[node].inputs[1] = input2;
                this->genotype[node].function = XOR;
                this->genotype[node].active = true;
                temp[node + circuit->num_inputs] = splited_buffer[0];
            }
            else if(num_and == 1 || num_xor == 1 || num_or == 1){
                int input1 = 0, input2 = 0, node = 0;

                input1 = parse_evoapprox8b_input(circuit, splited_buffer[2], temp);
                input2 = parse_evoapprox8b_input(circuit, splited_buffer[4], temp);
                node = get_next_available_position(circuit, 0);
                this->genotype[node].inputs[0] = input1;
                this->genotype[node].inputs[1] = input2;
                if(num_and) this->genotype[node].function = AND;
                if(num_or) this->genotype[node].function = OR;
                if(num_xor) this->genotype[node].function = XOR;
                this->genotype[node].active = true;
                temp[node + circuit->num_inputs] = splited_buffer[0];
            }
            else if(num_not == 1){
                int input1 = 0, node = 0;

                input1 = parse_evoapprox8b_input(circuit, splited_buffer[2], temp);
                node = get_next_available_position(circuit, 0);
                this->genotype[node].inputs[0] = input1;
                this->genotype[node].function = NOT;
                this->genotype[node].active = true;
                temp[node + circuit->num_inputs] = splited_buffer[0];
            }
            else{
                std::cout << "Invalid operation (" << buffer <<") option on seeding individual with EvoApprox8b circuit!" << std::endl;
                exit(0);
            }
        }
    }
}

template<typename S, typename T>
void Individual<S, T>::seed(Circuit *circuit, std::string filename){
    if(filename.find(".ep") != std::string::npos){
        seed_pla(circuit);
    }
    else if(filename.find(".c") != std::string::npos){
        seed_evoapprox8b(circuit, filename.c_str());
    }
}

template<typename S, typename T>
bdd Individual<S, T>::get_function_output(bdd left, int function, bdd right){
    switch (function)
    {
    case 1: //AND
        return bdd_and(left, right);
        break;
    case 2: //OR
        return bdd_or(left, right);
        break;
    case 3: //NOT
        return bdd_not(left);
        break;
    case 4: //NAND
        return bdd_not(bdd_and(left, right));
        break;
    case 5: //NOR
        return bdd_not(bdd_or(left, right));
        break;
    case 6: //XOR
        return bdd_xor(left, right);
        break;
    case 7: //XNOR
        return bdd_not(bdd_xor(left, right));
        break;
    default:
        std::cout << "Gate code unknow!" << std::endl;
        exit(0);
        break;
    }
}

template<typename S, typename T>
bdd Individual<S, T>::make_bdd_per_output(Circuit *circuit, int node){
    if(node < 0 || node >= (NCOL + circuit->num_inputs)){
        std::cout << "Output value isn't a valid one!" << std::endl;
        exit(0);
    }
    else if(node >= 0 && node < circuit->num_inputs){
        return bdd_ithvar(node);
    }
    
    int pos = node - circuit->num_inputs;

    if(this->genotype[pos].active){
        return this->genotype[pos].output;
    }
    if(this->genotype[pos].function == NOT){
        this->genotype[pos].output = bdd_not(make_bdd_per_output(circuit, this->genotype[pos].inputs[0]));
    }
    else{
        bdd left, right;
        left = make_bdd_per_output(circuit, this->genotype[pos].inputs[0]);
        right = make_bdd_per_output(circuit, this->genotype[pos].inputs[1]);

        this->genotype[pos].output = get_function_output(left, this->genotype[pos].function, right);
    }
    this->genotype[pos].active = true;
    return this->genotype[pos].output;
}

template<typename S, typename T>
void Individual<S, T>::evaluate_sat_count(Circuit *circuit){
    bdd temp;
    this->error = 0;
    this->mre = 0.0;
    clear_active_nodes();
    for(int i = 0; i < this->num_outputs; i++){
        temp = make_bdd_per_output(circuit, this->outputs[i]);
        temp = bdd_xor(temp, circuit->outputs[i]);
        this->outputs_error[i] = bdd_satcount(temp);
        this->error += this->outputs_error[i];
        this->mre += this->outputs_error[i]/pow(2,circuit->num_inputs);
    }
    this->mre = this->mre/circuit->num_outputs;
    count_transistors();
    clear_active_nodes();
}

#ifdef MO
template<typename S, typename T>
double Individual<S, T>::eval_delay(Circuit *circuit, int node){
    if(node >= 0 && node < circuit->num_inputs){
        return 0.0;
    }
    else if(node < 0){
        std::cout << "Output value isn't a valid one! (delay)" << std::endl;
        exit(0);
    }

    int pos = node - circuit->num_inputs;

    if(this->genotype[pos].active){
        return this->genotype[pos].delay;
    }
    else if(this->genotype[pos].function == NOT){
        this->genotype[pos].delay = eval_delay(circuit, this->genotype[pos].inputs[0]) + tdNOT;
    }
    else{
        double left = 0 , right = 0 , temp = 0;
        left = eval_delay(circuit, this->genotype[pos].inputs[0]);
        right = eval_delay(circuit, this->genotype[pos].inputs[1]);

        if(right > left) left = right;
        
        if(this->genotype[pos].function == AND) temp = tdAND;
        if(this->genotype[pos].function == NAND) temp = tdNAND;
        if(this->genotype[pos].function == OR) temp = tdOR;
        if(this->genotype[pos].function == NOR) temp = tdNOR;
        if(this->genotype[pos].function == XOR) temp = tdXOR;
        if(this->genotype[pos].function == XNOR) temp = tdXNOR;
        this->genotype[pos].delay = left + temp;
    }
    this->genotype[pos].active = true;
    return this->genotype[pos].delay;
}

template<typename S, typename T>
double Individual<S, T>::eval_power(Circuit *circuit, int node){
    if(node>= 0 && node < circuit->num_inputs){
        return 0.5;
    }
    else if(node < 0){
        std::cout << "Output value isn't a valid one! (power)" << std::endl;
        exit(0);
    }

    int pos = node - circuit->num_inputs;
    
    if(this->genotype[pos].active){
        return this->genotype[pos].power;
    }
    else if(this->genotype[pos].function == NOT){
        double left = eval_power(circuit, this->genotype[pos].inputs[0]);
        this->genotype[pos].power = 1.0 - left;
    }
    else{
        double left = 0, right = 0;
        left = eval_power(circuit, this->genotype[pos].inputs[0]);
        right = eval_power(circuit, this->genotype[pos].inputs[1]);

        if(this->genotype[pos].function == AND) 
            this->genotype[pos].power = left * right;
        if(this->genotype[pos].function == NAND) 
            this->genotype[pos].power = 1.0 - (left * right);
        if(this->genotype[pos].function == OR) 
            this->genotype[pos].power = 1.0 - (1.0 - left)*(1.0 - right);
        if(this->genotype[pos].function == NOR) 
            this->genotype[pos].power = (1.0 - left)*(1.0 - right);
        if(this->genotype[pos].function == XOR)
            this->genotype[pos].power = 1.0 - ((1.0 - left)*(1.0 - right) + (left * right));
        if(this->genotype[pos].function == XNOR) 
            this->genotype[pos].power = (1.0 - left)*(1.0 - right) + (left * right);   
    }
    this->genotype[pos].active = true;
    return this->genotype[pos].power;
}

template<typename S, typename T>
void Individual<S, T>::evaluate_delay(Circuit *circuit){
    double delayt = 0.0;
    this->delay = 0.0;
    clear_active_nodes();
    for(int i = 0; i < this->num_outputs; i++){
        delayt = eval_delay(circuit, this->outputs[i]);
        if(delayt > this->delay) this->delay = delayt;
        // std::cout << this->outputs[i] << " " << this->delay << " " << delayt << std::endl;
    }
    clear_active_nodes();
}

template<typename S, typename T>
void Individual<S, T>::evaluate_power(Circuit *circuit){
    this->power = 0.0;
    clear_active_nodes();
    for(int i = 0; i < this->num_outputs; i++){
        eval_power(circuit, this->outputs[i]);
    }

    for(int i = 0; i < NCOL; i++){
        if(this->genotype[i].active == true){
            this->genotype[i].power = this->genotype[i].power * (1.0 - this->genotype[i].power);
            this->power += this->genotype[i].power;
        }
    }
    this->power = this->power * Capa_Load * freq * Vcc * Vcc;
}

template <typename S, typename T>
void Individual<S, T>::print_MO(){
    std::cout << "MRE: " << this->mre << " Error: " << this->error << " Delay: " << 
    this->delay << "ns Power: " << this->power << 
        "mW Transistors: " << this->transistors << std::endl;
}
#endif