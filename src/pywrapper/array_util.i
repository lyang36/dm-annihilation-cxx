%include <carrays.i>
%include <std_string.i>
%include <std_vector.i>


namespace std {
    %template(vectori) vector<int>;
    #%template(vectorl) vector<long int>;
    %template(vectord) vector<double>;
    %template(vectorf) vector<float>;
    %template(vectori64) vector<int64_t>;
    %template(vectoriu64) vector<uint64_t>;
    %template(vectoriu32) vector<uint32_t>;
};

%inline %{
    std::vector<long long int> getInt64Array(int64_t * array, long int num){
        std::vector<long long int> a;
        for(int i = 0; i < num; i++){
            a.push_back(array[i]);
        }
        return a;
    };
    
    std::vector<int> getInt32Array(int32_t * array, long int num){
        std::vector<int> a;
        for(int i = 0; i < num; i++){
            a.push_back(array[i]);
        }
        return a;
    };
    
    std::vector<int> getUInt32Array(uint32_t * array, long int num){
        std::vector<int> a;
        for(int i = 0; i < num; i++){
            a.push_back(array[i]);
        }
        return a;
    };
    
    std::vector<int> getUInt64Array(uint64_t * array, long int num){
        std::vector<int> a;
        for(int i = 0; i < num; i++){
            a.push_back(array[i]);
        }
        return a;
    };
    
    std::vector<double> getDoubleArray(double * array, long int num){
        std::vector<double> a;
        for(int i = 0; i < num; i++){
            a.push_back(array[i]);
        }
        return a;
    };
    
    std::vector<float> getDoubleArray(float * array, long int num){
        std::vector<float> a;
        for(int i = 0; i < num; i++){
            a.push_back(array[i]);
        }
        return a;
    };
    
    double * getDoublePointer(std::vector<double> &a){
        return a.data();
    };
    
    float * getFloatPointer(std::vector<float> &a){
        return a.data();
    };
    
    int * getIntPointer(std::vector<int> &a){
        return a.data();
    };
    
    long int * getLongPointer(std::vector<long int> &a){
        return a.data();
    };
    
    uint32_t * getU32Pointer(std::vector<uint32_t> &a){
        return a.data();
    };
    
    uint64_t * getU64Pointer(std::vector<uint64_t> &a){
        return a.data();
    };
    
    int64_t * getI64Pointer(std::vector<int64_t> &a){
        return a.data();
    };
    
%}


