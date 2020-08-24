#include "population.tpp"

int setup_mutation;
int eval_flag;

template <typename S = long int, typename T = bdd>
class Cgp{
    public:
        Population<S, T> *population;                                               /**< Population */
        Circuit *circuit;                                                           /**< Circuit */
        int evaluations;                                                            /**< Number of evaluations */
        void (Population<S, T>::*mutation)(Circuit *circuit, int first_invidual);   /**< Mutation algorithm used */
        void (Population<S, T>::*select)();                                         /**< Selection Algorithm used */
        
        #ifdef MO
        bool seed_individual;                                                       /**< Seed individual flag */
        std::string seed_file;                                                      /**< Seed individual filename */
        #endif  

        /**
		* @brief Constructor
		* @param filename - The circuit definition file
		* @param setupfile - The CGP parameters setup file
		* @return none
		*/
        Cgp(int argc, char const *argv[]);

        /**
		* @brief Destructor
		* @return none
		*/
        ~Cgp();

        /**
		* @brief Initializes the population according to the parameter SeedIndividual, 
        * if SeedIndividual is defined as true the population is initialized based on the 
        * ESPRESSO circuit and them returns true. If SeedIndividual is definid as false, the 
        * population is initialized randomly and the proccess of search for a feasible 
        * solution starts. 
		* @return True if a feasible circuit was found and false otherwise.
		*/
        bool evolve();

        /**
		* @brief Starts with a feasible solution and search for the most optimized circuit
        * regarding the number of transistors.
		* @return none
		*/
        void optimize();

        /**
		* @brief Starts with a feasible solution and search for the mosts optimized circuits
        * regarding the delay, the error and the power.
		* @return none
		*/
        void optimize_MO();
};

template <typename S, typename T>
Cgp<S, T>::Cgp(int argc, char const *argv[]) {    
    
    #ifdef MO
    seed_individual = 0;
    #endif

    this->circuit = new Circuit(argv[2]);
    this->circuit->make_circuit_bdd();

    std::cout << "INPUTS: " << INPUTS << std::endl;
    NCOL = atoi(argv[7]);
    std::cout << "NCOL: " << NCOL << std::endl;
    LB = atoi(argv[7]);
    std::cout << "LB: " << LB << std::endl;


    this->evaluations = atoi(argv[3]);
    std::cout << "NumEvaluations: " << this->evaluations << std::endl;

    int num_individuals = atoi(argv[4]);
    std::cout << "PopSize: " << num_individuals << std::endl;
    this->population = new Population<S, T>(num_individuals, this->circuit->num_outputs);

    std::cout << argv[5] << std::endl;
    std::cout << "Mutation: ";
    if(!strcmp(argv[5], "SAM")){
        this->mutation = &Population<S, T>::apply_SAM;
        std::cout << "SAM" << std::endl;
    }
    else if(!strcmp(argv[5], "PM")){
        this->mutation = &Population<S, T>::apply_PM;
        std::cout << "PM" << std::endl;    
    }
    else if(!strcmp(argv[5], "SG")){
        this->mutation = &Population<S, T>::apply_SG;
        std::cout << "SG" << std::endl;
    }
    else{
        std::cout << "Mutation Method not recognized!" << std::endl;
        exit(0);
    }

    std::cout << "Selection: ";
    if(!strcmp(argv[6], "APS")){
        this->select = &Population<S, T>::select_aps;
        std::cout << "APS" << std::endl;
    }
    else if(!strcmp(argv[6], "NSGA2")){
        std::cout << "NSGA2" << std::endl;
        this->select = &Population<S, T>::select_nsga2;
    }
    else{
        std::cout << "Selection Method not recognized!" << std::endl;
        exit(0);
    }
    

    if(argc == 9){
        this->seed_file = argv[8];
        this->seed_individual = true;
        std::cout << "Seeding Population: " << this->seed_file << std::endl;
    }

    std::cout << "Cgp setup finished successfully!" << std::endl;
    std::cout << "*~*~*~* *~*~*~* *~*~*~* *~*~*~* *~*~*~* *~*~*~* *~*~*~*" << std::endl;
}

template <typename S, typename T>
Cgp<S, T>::~Cgp(){
    this->population->~Population();
    this->circuit->~Circuit();
}

template <typename S, typename T>
bool Cgp<S, T>::evolve(){

    int best_indv = 0;
    eval_flag = this->evaluations - 100000;
    int first_individual = 1;

    this->population->initialize(circuit);

    #ifdef MO
    if(this->seed_individual){
        this->population->individuals[0].seed(circuit, this->seed_file);
        this->population->individuals[0].clear_active_nodes();
        this->population->individuals[0].evaluate_sat_count(circuit);
        this->population->individuals[0].evaluate_delay(circuit);
        this->population->individuals[0].evaluate_power(circuit);   

        if(this->population->individuals[0].error != 0){
            std::cout << "Seeding population didn't work! Error = " << this->population->individuals[0].error << std::endl;
            
            this->population->individuals[0].print(circuit);
            for(int i = 0; i < circuit->num_outputs; i++)
            {
                std::cout << this->population->individuals[0].outputs_error[i] << "  ";
            }
            exit(0);
        }
        else{
            std::vector<bdd> new_outputs;
            for(int i = 0; i < circuit->num_outputs; i++){
                new_outputs.push_back(this->population->individuals[0].make_bdd_per_output(circuit, this->population->individuals[0].outputs[i]));
            }
            circuit->update_base_circuit(new_outputs);

            std::cout << "Seeding population worked successfully!" << std::endl;
            std::cout << "Evolution Final Solution:" << std::endl;
            this->population->individuals[0].print_MO();
            std::cout << "*~*~*~* *~*~*~* *~*~*~* *~*~*~* *~*~*~* *~*~*~* *~*~*~*" << std::endl;
            return true;
        }
    }
    else{
        std::cout << "Initializing population worked successfully!" << std::endl;
        std::cout << "Evolution Final Solution:" << std::endl;
        this->population->individuals[0].evaluate_sat_count(circuit);
        this->population->individuals[0].evaluate_power(circuit);
        this->population->individuals[0].evaluate_delay(circuit);
        this->population->individuals[0].print_MO();
        std::cout << "*~*~*~* *~*~*~* *~*~*~* *~*~*~* *~*~*~* *~*~*~* *~*~*~*" << std::endl;
        return true;
    }
    #endif

    this->population->individuals[0].evaluate_sat_count(circuit);
    this->population->individuals[0].evaluate_power(circuit);
    this->population->individuals[0].evaluate_delay(circuit);
    std::cout << "Evaluations: " << this->evaluations << " Individual: " << 0 << " ";
    this->population->individuals[0].print_MO();

    while (1){
        this->population->clone_best_individual(best_indv);
        (this->population->*mutation)(circuit, first_individual);
        this->population->evaluate(circuit, first_individual);
        evaluations -= this->population->size - 1;
        best_indv = this->population->get_best_sat_count();

        if(this->population->individuals[best_indv].error == 0){
            this->population->clone_best_individual(best_indv);
            
            std::vector<bdd> new_outputs;
            for(int i = 0; i < circuit->num_outputs; i++){
                new_outputs.push_back(this->population->individuals[0].make_bdd_per_output(circuit, this->population->individuals[0].outputs[i]));
            }
            circuit->update_base_circuit(new_outputs);     

            this->population->individuals[0].evaluate_power(circuit);
            this->population->individuals[0].evaluate_delay(circuit);
            std::cout << "Evolution Final Solution:" << std::endl;
            #ifdef MO
            this->population->individuals[0].print_MO();
            #endif
            #ifndef MO
            this->population->individuals[0].print(circuit);
            #endif
            std::cout << "*~*~*~* *~*~*~* *~*~*~* *~*~*~* *~*~*~* *~*~*~* *~*~*~*" << std::endl;

            return true;
        }
        if(this->evaluations - (this->population->size - 1) <= 0) break;
        if(this->evaluations <= eval_flag){
            this->population->individuals[0].evaluate_power(circuit);
            this->population->individuals[0].evaluate_delay(circuit);
            std::cout << "Evaluations: " << this->evaluations << " Individual: " << 
            best_indv << " ";
            this->population->individuals[0].print_MO();
            eval_flag = this->evaluations - 100000;
            fflush(stdout);
        }
        if(bdd_getnodenum() >= (int) (0.8 * bdd_getallocnum())) bdd_gbc();
    }
    for(int i = 0; i < this->population->size; i++){
        this->population->individuals[i].evaluate_power(circuit);
        this->population->individuals[i].evaluate_delay(circuit);
        std::cout << "Individual: " << i << " ";
        this->population->individuals[i].print_MO();
    }

    return false;
}

template <typename S, typename T>
void Cgp<S, T>::optimize(){
    int best_indv = -1;
    int count = 0;
    if(setup_mutation == 3) this->mutation = &Population<S, T>::apply_SAM;
 
    while (1){
        (this->population->*mutation)(circuit, this->population->size/2);
        this->population->evaluate(circuit);
        this->evaluations -= this->population->size - 1;
        best_indv = this->population->get_optimized();
        this->population->clone_best_individual(best_indv);

        if(best_indv != 0) count = 0;
        if(count >= 100000) break;
        if(this->evaluations - (this->population->size - 1) <= 0) break;
        if(this->evaluations <= eval_flag){
            std::cout << "Evaluations: " << this->evaluations << " Error: " << 
            this->population->individuals[0].error << " Transistors: " << 
            this->population->individuals[0].transistors << std::endl;
            fflush(stdout);
        }
        if(bdd_getnodenum() >= (int) (0.8 * bdd_getallocnum())) bdd_gbc();

        count++;
    }
    std::cout << "Optimization Final Solution:" << std::endl;
    this->population->individuals[0].print(circuit);
}

template <typename S, typename T>
void Cgp<S, T>::optimize_MO(){
    int first_individual = 1;
    if(setup_mutation == 3){
        this->mutation = &Population<S, T>::apply_SAM;
        std::cout << "Mutation changed to SAM!" << std::endl;
    }
    while (1){
        
        this->population->clone_MO();
        (this->population->*mutation)(circuit, first_individual);
        this->population->evaluate_MO(circuit, first_individual);
        this->evaluations -= this->population->size/2;
        (this->population->*select)();

        if(this->evaluations - this->population->size <= 0) break;
        if(this->evaluations <= eval_flag){
            std::cout << "Evaluations: " << this->evaluations << std::endl;
            eval_flag = this->evaluations - 100000;
            this->population->print_MO();
            fflush(stdout);
        }

        if(this->evaluations == 0) break;
        if(bdd_getnodenum() >= (int) (0.8 * bdd_getallocnum())) bdd_gbc();
        first_individual = this->population->size;
    }
    std::cout << "Evaluations: " << this->evaluations << std::endl;
    std::cout << "Optimization Final Solution:" << std::endl;
    this->population->print_MO();
}