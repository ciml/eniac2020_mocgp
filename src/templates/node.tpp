#include<vector>
#define INPUTS 2

template<typename T>
class Node{
    public:
        std::vector<int> inputs;   /**< Vector of inputs */
        int function;              /**< Node function */
        bool active;               /**< Active node */
        T output;                  /**< Node function */

        #ifdef MO
        double power;              /**< Node power */
        double delay;              /**< Maximum total delay at this node path  */
        #endif

        /**
		* @brief Constructor
		* @return none
		*/
        Node();

        /**
		* @brief Destructor
		* @return none
		*/
        ~Node();
};

template<typename T>
Node<T>::Node(){
    this->inputs = std::vector<int>(INPUTS, -1);
    this->function = 0;
    this->active = false;

    #ifdef MO
    power = -1;
    delay = -1;
    #endif
}

template<typename T>
Node<T>::~Node(){
    std::vector<int>().swap(this->inputs);
}