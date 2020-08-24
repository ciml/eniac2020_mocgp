#include "individual.tpp"

struct doubint{
    double value;
    int indv;
};

template <typename S, typename T>
class Population{
    public:
        std::vector<Individual<S, T>> individuals;     /**< Vector of individuals */
        int max_size;                                  /**< Maximum population size */
        int size;                                      /**< Population size */

        /**
		* @brief Constructor
		* @param population_size - The population size
		* @param num_outputs - The circuit's number of outputs
		* @return none
		*/
        Population(int population_size, int num_outputs);
        
        /**
		* @brief Destructor
		* @return none
		*/
        ~Population();

        /**
		* @brief Initialize all individuals
		* @param circuit - The circuit
		* @return none
		*/
        void initialize(Circuit *circuit);
        
        /**
		* @brief Apply the Point Mutation in the individuals starting at the first 
        * individual send as parameter
		* @param circuit - The circuit
		* @param first_individual - The first individual to be mutated
		* @return none
		*/
        void apply_PM(Circuit *circuit, int first_individual);

        /**
		* @brief Apply the Single Active Mutation in the individuals starting 
        * at the first individual send as parameter
		* @param circuit - The circuit
		* @param first_individual - The first individual to be mutated
		* @return none
		*/
        void apply_SAM(Circuit *circuit, int first_individual);

        /**
		* @brief Apply the Guided Active Mutation in the individuals starting 
        * at the first individual send as parameter
		* @param circuit - The circuit
		* @param first_individual - The first individual to be mutated
		* @return none
		*/
        void apply_GAM(Circuit *circuit, int first_individual);

        /**
		* @brief Apply the Single and Active Mutation in the individuals 
        * starting at the first individual send as parameter
		* @param circuit - The circuit
		* @param first_individual - The first individual to be mutated
		* @return none
		*/
        void apply_SG(Circuit *circuit, int first_individual);

        /**
		* @brief Evaluate the Sat Count for the individuals 
        * starting at the first individual send as parameter
		* @param circuit - The circuit
		* @param first_individual - The first individual to be evaluated
		* @return none
		*/
        void evaluate(Circuit *circuit, int first_individual);

        /**
		* @brief Evaluate the Sat Count, the delay and the power for the 
        * individuals starting at the first individual send as parameter
		* @param circuit - The circuit
		* @param first_individual - The first individual to be evaluated
		* @return none
		*/
        void evaluate_MO(Circuit *circuit, int first_individual);

        /**
		* @brief Finds the individual with the lower Sat Count
		* @return The position of the best individual in the vector
		*/
        int get_best_sat_count();

        /**
		* @brief Finds the individual with the lowe Sat Count and
        * the lower number of transistors
		* @return none
		*/
        int get_optimized();

        /**
		* @brief Clone the best individual
		* @param best_indv - Best individual
		* @return none
		*/
        void clone_best_individual(int best_indv);
        
        /**
		* @brief Clone the individuals from the first Pareto's Front
		* @return none
		*/
        void clone_MO();

        /**
		* @brief Print information of the individuals from the first Pareto's Front
		* @return none
		*/
        void print_MO();

        /**
		* @brief Fast Non Dominated Sort algorithm from NSGAII
		* @return none
		*/
        void fast_non_dominated_sort();

        /**
		* @brief Fast Non Dominated Sort algorithm from NSGAII adapted to 
        * an Constrained Multiobjective Optimization
		* @return none
		*/
        void constrained_fast_non_dominated_sort();

        /**
		* @brief Sort an restriction in increasing order
		* @param restriction - The restriction to be sorted
		* @param Fi - The i-th Pareto's Front
		* @return The order of individuals according to the ordered restriction value
		*/
        std::vector<doubint> sort_restriction(int restriction, std::vector<int> Fi);
        
        /**
		* @brief Crowding Distance algorithm from NSGAII
		* @param Fi - The i-th Pareto's Front
		* @return The individuals crowding distance
		*/
        std::vector<double> crowding_distance(std::vector<int> Fi);

        /**
		* @brief Crowded Comparison algorithm from NSGAII
		* @param rank - The Pareto's Front rank
		* @param num_individuals - The number of individuals to be choose 
        * to the next generation
		* @param Fi - The i-th Pareto's Front
		* @param individuals - The vector of individuals
		* @return none
		*/
        void crowded_comparison(int rank, int num_individuals, std::vector<int> Fi, std::vector<Individual<S, T>> individuals);
        
        /**
		* @brief Non Dominated Sorting Genetic Algorithm 2 selection method based on NSGAII
		* @return none
		*/
        void select_nsga2();
        
        /**
		* @brief Adaptative Population Size selection method
		* @return none
		*/
        void select_aps();
};

template <typename S, typename T>
Population<S, T>::Population(int population_size, int num_outputs){
    this->size = population_size;
    this->max_size = population_size;
    this->individuals = std::vector<Individual<S, T>>(this->size, Individual<S,T>(num_outputs));


}

template <typename S, typename T>
Population<S, T>::~Population(){
    for(int i = 0; i < this->size; i++) this->individuals[i].~Individual();
    std::vector<Individual<S, T>>().swap(this->individuals);
}

template <typename S, typename T>
void Population<S, T>::initialize(Circuit *circuit){
    for(int i = 0; i < this->size; i++){
        individuals[i].initialize(circuit);
    }
}

template <typename S, typename T>
void Population<S, T>::apply_PM(Circuit *circuit, int first_individual){
    for(int i = first_individual; i < this->size; i++){
        this->individuals[i].apply_PM(circuit);
    }
}

template <typename S, typename T>
void Population<S, T>::apply_SAM(Circuit *circuit, int first_individual){
    for (int i = first_individual; i < this->size; i++)
    {
        this->individuals[i].apply_SAM(circuit);
    }
}

template <typename S, typename T>
void Population<S, T>::apply_GAM(Circuit *circuit, int first_individual){
    for (int i = first_individual; i < this->size; i++)
    {
        this->individuals[i].apply_GAM(circuit);
    }
}

template <typename S, typename T>
void Population<S, T>::apply_SG(Circuit *circuit, int first_individual){
    int temp = (this->size - first_individual) / 2;
    for (int i = first_individual; i < (first_individual + temp); i++)
    {
        this->individuals[i].apply_SAM(circuit);
    }
    for (int i = (first_individual + temp); i < this->size; i++)
    {
        this->individuals[i].apply_GAM(circuit);
    }
}

template <typename S, typename T>
void Population<S, T>::evaluate(Circuit *circuit, int first_individual){
    for(int i = first_individual; i < this->size; i++){
        this->individuals[i].evaluate_sat_count(circuit);
    }
}

template <typename S, typename T>
void Population<S, T>::evaluate_MO(Circuit *circuit, int first_individual){
    for(int i = first_individual; i < this->size; i++){
        this->individuals[i].evaluate_sat_count(circuit);
        #ifdef MO
        this->individuals[i].evaluate_delay(circuit);
        this->individuals[i].evaluate_power(circuit);
        #endif
    }
}

template <typename S, typename T>
int Population<S, T>::get_best_sat_count(){
    std::vector<int> temparray;
    int best = 0;

    temparray.push_back(0);

    for(int i = 1; i < this->size; i++){
        if(this->individuals[i].error < this->individuals[best].error){
            best = i;
            temparray.clear();
            temparray.push_back(i);
        }
        else if(this->individuals[i].error == this->individuals[best].error){
            temparray.push_back(i);
        }
    }
    return temparray[rand() % (int)temparray.size()];
}

template <typename S, typename T>
int Population<S, T>::get_optimized(){
    std::vector<int> temparray;
    int best = 0;

    temparray.push_back(0);

    for(int i = 1; i < this->size; i++){
        if(this->individuals[i].error == this->individuals[best].error &&
        this->individuals[i].transistors <= this->individuals[best].transistors){
            best = i;
            temparray.clear();
            temparray.push_back(i);
        }
        else if(this->individuals[i].error == this->individuals[best].error &&
        this->individuals[i].transistors == this->individuals[best].transistors){
            temparray.push_back(i);
        }
    }
    return temparray[rand() % (int)temparray.size()];
}

template <typename S, typename T>
void Population<S, T>::clone_best_individual(int best_indv){
    this->individuals[best_indv].clear_active_nodes();
    for(int i = 0; i < this->size; i++){
        if(i != best_indv){
            this->individuals[i] = this->individuals[best_indv];
        }
    }
}

#ifdef MO
struct sort_crescent{
    inline bool operator() (const doubint& struct1, const doubint& struct2)
    {
        return (struct1.value < struct2.value);
    }
};

template <typename S, typename T>
void Population<S, T>::fast_non_dominated_sort(){
    std::vector<int> np(this->size, 0);
    std::vector<std::vector<int>> sp(this->size);
    std::vector<int> front;
    //std::cout << "Front 0: ";
    for(int p = 0; p < this->size; p++){
        for(int q = 0; q < this->size; q++){
            if(p != q){
                if(this->individuals[p].mre <= 0.10 && this->individuals[q].mre <= 0.10){
                    if(((this->individuals[p].error <= this->individuals[q].error)  && 
                        (this->individuals[p].delay <= this->individuals[q].delay)  && 
                        (this->individuals[p].power <= this->individuals[q].power)) &&
                        ((this->individuals[p].error < this->individuals[q].error)  || 
                        (this->individuals[p].delay < this->individuals[q].delay)   || 
                        (this->individuals[p].power < this->individuals[q].power))){
                            sp[p].push_back(q);
                    }
                    else if (((this->individuals[q].error <= this->individuals[p].error)  && 
                            (this->individuals[q].delay <= this->individuals[p].delay)  && 
                            (this->individuals[q].power <= this->individuals[p].power)) &&
                            ((this->individuals[q].error < this->individuals[p].error)  || 
                            (this->individuals[q].delay < this->individuals[p].delay)   || 
                            (this->individuals[q].power < this->individuals[p].power))){
                                np[p] += 1;
                    }
                }
                else if(this->individuals[p].mre <= 0.10 || this->individuals[q].mre <= 0.10){
                    if(this->individuals[p].mre <= 0.10) sp[p].push_back(q);
                    else if(this->individuals[q].mre <= 0.10) np[p] += 1;
                }
                else if(this->individuals[p].mre > 0.10 && this->individuals[q].mre > 0.10){
                    if(this->individuals[p].mre < this->individuals[q].mre) sp[p].push_back(q);
                    else if(this->individuals[q].mre < this->individuals[p].mre) np[p] += 1;
                }
            }
        }
        if(np[p] == 0){
            this->individuals[p].rank = 0;
            front.push_back(p);
            //std::cout << p << " ";
        }
    }
    //std::cout << std::endl;
    int i = 0;
    while((int)front.size() > 0)
    {
        std::vector<int> Q;
        //std::cout << "Front " << i+1 << ": ";
        for(int p = 0; p < (int)front.size(); p++){
            int sp_size = (int)sp[front[p]].size();
            for(int q = 0; q < sp_size; q++){
                int temp = sp[front[p]][q];
                np[temp] -= 1;

                if(np[temp] == 0){
                    //std::cout << temp << " ";
                    this->individuals[temp].rank = i + 1;
                    Q.push_back(temp);
                }
            }
        }
        //std::cout<< std::endl;
        i = i + 1;
        front = Q;
    }
}

template <typename S, typename T>
void Population<S, T>::constrained_fast_non_dominated_sort(){
    std::vector<int> np(this->size, 0);
    std::vector<std::vector<int>> sp(this->size);
    std::vector<int> front;
    //std::cout << "Front 0: ";
    for(int p = 0; p < this->size; p++){
        for(int q = 0; q < this->size; q++){
            if(p != q){
                if(this->individuals[p].error == 0 && this->individuals[q].error == 0){
                    if(((this->individuals[p].transistors <= this->individuals[q].transistors)  && 
                        (this->individuals[p].delay <= this->individuals[q].delay)  && 
                        (this->individuals[p].power <= this->individuals[q].power)) &&
                        ((this->individuals[p].transistors < this->individuals[q].transistors)  || 
                        (this->individuals[p].delay < this->individuals[q].delay)   || 
                        (this->individuals[p].power < this->individuals[q].power))){
                            sp[p].push_back(q);
                    }
                    else if (((this->individuals[q].transistors <= this->individuals[p].transistors)  && 
                            (this->individuals[q].delay <= this->individuals[p].delay)  && 
                            (this->individuals[q].power <= this->individuals[p].power)) &&
                            ((this->individuals[q].transistors < this->individuals[p].transistors)  || 
                            (this->individuals[q].delay < this->individuals[p].delay)   || 
                            (this->individuals[q].power < this->individuals[p].power))){
                                np[p] += 1;
                    }
                }
                else if(this->individuals[p].error == 0 || this->individuals[q].error == 0){
                    if(this->individuals[p].error == 0) sp[p].push_back(q);
                    else if(this->individuals[q].error == 0) np[p] += 1;
                }
                else if(this->individuals[p].error != 0 && this->individuals[q].error != 0){
                    if(this->individuals[p].error < this->individuals[q].error) sp[p].push_back(q);
                    else if(this->individuals[q].error < this->individuals[p].error) np[p] += 1;
                }
            }
        }
        if(np[p] == 0){
            this->individuals[p].rank = 0;
            front.push_back(p);
            //std::cout << p << " ";
        }
    }
    //std::cout << std::endl;
    int i = 0;
    while((int)front.size() > 0)
    {
        std::vector<int> Q;
        //std::cout << "Front " << i+1 << ": ";
        for(int p = 0; p < (int)front.size(); p++){
            int sp_size = (int)sp[front[p]].size();
            for(int q = 0; q < sp_size; q++){
                int temp = sp[front[p]][q];
                np[temp] -= 1;

                if(np[temp] == 0){
                    //std::cout << temp << " ";
                    this->individuals[temp].rank = i + 1;
                    Q.push_back(temp);
                }
            }
        }
        //std::cout<< std::endl;
        i = i + 1;
        front = Q;
    }
}

template<typename S, typename T>
std::vector<doubint> Population<S, T>::sort_restriction(int restriction, std::vector<int> Fi){
    std::vector<doubint> temp_array;
    
    for(int i = 0; i < (int)Fi.size(); i++){
        doubint temp;

        if(restriction == 1) temp.value = (double)this->individuals[Fi[i]].error;
        if(restriction == 2) temp.value = this->individuals[Fi[i]].delay;
        if(restriction == 3) temp.value = this->individuals[Fi[i]].power;
        temp.indv = Fi[i];

        temp_array.push_back(temp);
    }

    std::sort(temp_array.begin(), temp_array.end(), sort_crescent());

    // if(restriction == 1) std::cout << "Error:" << std::endl;
    // if(restriction == 2) std::cout << "Delay:" << std::endl;
    // if(restriction == 3) std::cout << "Power:" << std::endl;
    // for(int i = 0; i < (int)temp_array.size(); i++){
    //     std::cout << temp_array[i].indv << " " << temp_array[i].value << std::endl;
    // }
    return temp_array;
}

template <typename S, typename T>
std::vector<double> Population<S, T>::crowding_distance(std::vector<int> Fi){
    std::vector<double> CRD(Fi.size(), 0.0);
    std::vector<int>::iterator it;
    std::vector<doubint> temp;
    int Fi_size = (int)Fi.size();
    int index;

    // std::cout << "Initial CRD\n";
    // for(int i = 0; i < (int)CRD.size(); i++){
    //     std::cout << CRD[i] << " ";
    // }
    // std::cout << std::endl;

    for(int i = 1; i < 4; i++){
        temp = sort_restriction(i, Fi);
        
        it = std::find(Fi.begin(), Fi.end(), temp[0].indv);
        index = std::distance(Fi.begin(), it);
        CRD[index] += 1000.0;

        it = std::find(Fi.begin(), Fi.end(), temp[Fi_size - 1].indv);
        index = std::distance(Fi.begin(), it);
        CRD[index] += 1000.0;

        for(int j = 1; j < Fi_size - 1; j++){
            it = std::find(Fi.begin(), Fi.end(), temp[j].indv);
            index = std::distance(Fi.begin(), it);
            CRD[index] += (temp[j+1].value - temp[j-1].value)/(temp[Fi_size - 1].value - temp[0].value);
        }
    }

    // for(int i = 0; i < (int)CRD.size(); i++){
    //     std::cout << CRD[i] << " ";
    // }
    // std::cout << std::endl;
    return CRD;
}

template <typename S, typename T>
void Population<S, T>::crowded_comparison(int rank, int num_individuals, std::vector<int> Fi, std::vector<Individual<S, T>> individuals){
    std::vector<double> CRD;
    std::vector<doubint> temp_array;
    int last = 0;

    CRD = crowding_distance(Fi);

    for(int i = 0; i < (int)CRD.size(); i++){
        doubint temp;
        temp.value = CRD[i];
        temp.indv = Fi[i];
        
        temp_array.push_back(temp);
    }
    sort(temp_array.begin(), temp_array.end(), sort_crescent());
    
    last = (int)temp_array.size() - 1;
    for(int i = last; i > last - num_individuals; i--){
        individuals.push_back(this->individuals[temp_array[i].indv]);
    }

    this->individuals = individuals;
    this->size = this->individuals.size();
}

template <typename S, typename T>
void Population<S, T>::clone_MO(){
    for(int i = 0; i < this->size; i++){
        this->individuals.push_back(this->individuals[i]);
    }
    this->size = (int)this->individuals.size();
}

template <typename S, typename T>
void Population<S, T>::select_nsga2(){
    int pop_size = this->max_size;
    int new_pop_size = 0;
    int rank = 0;
    std::vector<Individual<S,T>> new_individuals;

    fast_non_dominated_sort();
    
    while(1){
        std::vector<int> Fi;
        for(int i = 0; i < this->size; i++){
            if(this->individuals[i].rank == rank){
                Fi.push_back(i);
            }
        }

        if(new_pop_size + (int)Fi.size() <= pop_size){
            for(int i = 0; i < (int)Fi.size(); i++){
                new_individuals.push_back(this->individuals[Fi[i]]);
            }
            new_pop_size += (int)Fi.size();
        }
        else{
            crowded_comparison(rank, pop_size - new_pop_size, Fi, new_individuals);
            break;
        }
        rank++;
    }
}

template <typename S, typename T>
void Population<S, T>::select_aps(){
    int pop_size = this->max_size;
    std::vector<Individual<S,T>> new_individuals;

    fast_non_dominated_sort();
    
    std::vector<int> Fi;
    for(int i = 0; i < this->size; i++){
        if(this->individuals[i].rank == 0){
            Fi.push_back(i);
        }
    }

    if((int)Fi.size() <= pop_size){
        for(int i = 0; i < (int)Fi.size(); i++){
            new_individuals.push_back(this->individuals[Fi[i]]);
        }
        this->individuals = new_individuals;
        this->size = (int)this->individuals.size();
    }
    else{
        crowded_comparison(0, pop_size, Fi, new_individuals);
    }
    // std::cout << "PopSize: " << this->size << "  ";
}


template <typename S, typename T>
void Population<S, T>::print_MO(){
    for(int i = 0; i < this->size; i++){
        if(this->individuals[i].rank == 0) this->individuals[i].print_MO();
    }
    std::cout << std::endl;
}

#endif