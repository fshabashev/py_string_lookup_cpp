#include <iostream>
#include <math.h>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "arithmetic_coding/encode.h"
#include "arithmetic_coding/decode.h"



double pi_integral() {
    double sum = 0.0;
    long n = 100*1000*1000;
    for (long k = 0; k < n; k++){
        double x = ((double) k)/n;
        double val = sqrt(1.0 - x*x);
        sum += val;
    }
    return 4.0*sum/n;
}

double bbp_pi_calculation_no_multiplier(long k) {
    double k_f64 = (double) k;
    double val = 4.0/(8.0*k_f64+1.0) - 2.0/(8.0*k_f64 + 4.0) - 1.0/(8.0*k_f64+5.0) - 1.0/(8.0*k_f64+6.0);
    return val;
}

double bbp_pi_calcuation(long k){
    double val = bbp_pi_calculation_no_multiplier(k);
    double k_f64 = (double) k;
    double multiplier = 1.0/(pow(16, k_f64));
    return val*multiplier;
}

void calculate_pi_bbp() {
    double sum = 0.0;
    for (long i = 0; i < 100; i++){
        sum += bbp_pi_calcuation(i);
    }
    std::cout << sum << std::endl;
}

std::string gen_random(const int len) {
    static const char alphanum[] =
            "0123456789"
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            "abcdefghijklmnopqrstuvwxyz";
    std::string tmp_s;
    tmp_s.reserve(len);

    for (int i = 0; i < len; ++i) {
        tmp_s += alphanum[rand() % (sizeof(alphanum) - 1)];
    }

    return tmp_s;
}

std::vector<std::string> gen_rand_string_vector(long N){
    std::vector<std::string> string_vec;
    int P = 10;
    for (int i = 0; i<N; i++){
        string_vec.push_back(gen_random(P));
    }
    return string_vec;

}

// sort vector of pairs
template <typename T>
void sort_vector_of_pairs(std::vector<std::pair<T, T>> &vec) {
    std::sort(vec.begin(), vec.end(),
              [](const std::pair<T, T> &a, const std::pair<T, T> &b) {
                  return a.first < b.first;
              });
}

std::vector<std::pair<std::string, std::string> > gen_random_pairs(long N){
    std::vector<std::pair<std::string, std::string> > string_vec;
    int P = 10;
    for (int i = 0; i<N; i++){
        string_vec.push_back(std::make_pair(gen_random(P), gen_random(P)));
    }
    return string_vec;
}

class StringPairHolder {
public:
    std::vector<std::pair<std::string, std::string> > pair_storage;
    long N = 100*1000;
    StringPairHolder() {
        pair_storage = gen_random_pairs(N);
        sort_vector_of_pairs(pair_storage);
    }

    StringPairHolder(std::vector< std::pair<std::string, std::string> > vec) {
        pair_storage = vec;
        sort_vector_of_pairs(pair_storage);
    }

    long binary_search(std::string key) {
        long left = 0;
        long right = pair_storage.size() - 1;
        while (left <= right) {
            long mid = (left + right) / 2;
            if (pair_storage[mid].first == key) {
                return mid;
            } else if (pair_storage[mid].first < key) {
                left = mid + 1;
            } else {
                right = mid - 1;
            }
        }
        return -1;
    }

    std::string get_value(std::string key) {
        long index = binary_search(key);
        if (index == -1) {
            return "";
        } else {
            return pair_storage[index].second;
        }
    }

    bool contains(std::string key) {
        return binary_search(key) != -1;
    }
};

class StringAccumulator {
private:
    long size;
    long iter = 0;
public:
    std::vector<std::string> string_storage;

    StringAccumulator() {
        this->iter = 0;
        this->size = 0;
    }
    StringAccumulator(long n) {
        this->size = n;
        this->iter = 0;
        this->string_storage = std::vector<std::string>(n);
    }
    void add(std::string str_to_add) {
        if ((this->iter) < (this->size)) {
            string_storage[iter] = str_to_add;
            (this->iter)++;
        }
        else {
            throw std::runtime_error("StringAccumulator is full");
        }
    }

    void clean() {
        iter = 0;
        size = 0;
        this->string_storage.clear();
    }
};

class CharFrequencyHolder {
public:
    std::map<std::string, unsigned long> freq_array;

    CharFrequencyHolder() {
        this->freq_array = std::map<std::string, unsigned long>();
    }

    void incrementFreq(std::string str){
        for (int i = 0; i < str.length(); i++){
            if (freq_array.find(str.substr(i, 1)) == freq_array.end()){
                freq_array[str.substr(i, 1)] = 1;
            }
            else {
                freq_array[str.substr(i, 1)] += 1;
            }
        }
    }

    void printFreq(){
        for (auto it = freq_array.begin(); it != freq_array.end(); it++){
            std::cout << it->first << " " << it->second << std::endl;
        }
    }
};

class StringFrequencyHolder {
public:
    std::map<std::string, unsigned long> freq_array;
    
    StringFrequencyHolder(CharFrequencyHolder encoding){
        unsigned long sum = 0;
        for (auto it = encoding.freq_array.begin(); it != encoding.freq_array.end(); it++){
            sum += it->second;
            freq_array[it->first] = sum;
        }
    }
};



class StringHolder {
public:
    std::vector<std::string> string_storage;
    long N = 100000;
    StringHolder(){
        this->string_storage = gen_rand_string_vector(this->N);
        std::sort(this->string_storage.begin(), this->string_storage.end());
    }

    StringHolder(std::vector<std::string> string_storage_param){
        this->string_storage = string_storage_param;
        std::sort(this->string_storage.begin(), this->string_storage.end());
    }

    StringHolder(StringAccumulator accumulator){
        this->string_storage = accumulator.string_storage;
        std::sort(this->string_storage.begin(), this->string_storage.end());
    }

    long get_index(std::string str){
        auto v = this->string_storage;
        auto val = str;
        auto lower = std::lower_bound(v.begin(), v.end(), val);
        const bool found = lower != v.end() && *lower == val;
        if (!found) {
            return -1;
        }
        auto idx = std::distance(v.begin(), lower);
        return idx;
    }

    bool lookup(std::string str) {
        return std::binary_search(this->string_storage.begin(), this->string_storage.end(), str);
    }
};

int main() {
    std::cout << "printing PI_   " << pi_integral() << std::endl;
    calculate_pi_bbp();

    auto str_storage = StringHolder();

    for (int i =0; i<str_storage.N; i++){
        std::cout << str_storage.string_storage[i];
        std::cout << "\n";
    }
    std::cout << str_storage.lookup(str_storage.string_storage[20]) << std::endl;
    std::cout << str_storage.lookup(std::string("asdfasdf")) << std::endl;


    auto str_pair_storage = StringPairHolder();
    std::cout << "test get_value " << str_pair_storage.get_value("asdf") << std::endl;
    std::cout << "test contains " <<  str_pair_storage.contains(std::string("asdf")) << std::endl;
    std::cout << "test get_value " << str_pair_storage.get_value(str_pair_storage.pair_storage[1].first) << std::endl;
    return 0;
}
namespace py = pybind11;

float some_fn (float arg1, float arg2) {
    return arg1 + arg2;
}

std::vector<int> read_vec_from_file(std::string filename) {
    std::ifstream file;
    file.open(filename, ios::binary | ios::in);
    std::vector<int> vec;
    int i;
    while (true) {
        auto val = file.get();
        vec.push_back(val);
        if (val==EOF){
            break;
        }
    }
    file.close();
    return vec;
}

void write_vec_to_file(std::string filename, std::vector<int> vec) {
    std::ofstream file;
    file.open(filename, ios::binary | ios::out);
    for (int i = 0; i < vec.size(); i++){
        file.put(vec[i]);
    }
    file.close();
}

void print_int_vector(std::vector<int> vec){
    for (int i = 0; i < vec.size(); i++){
        std::cout << vec[i] << std::endl;
    }
}

void print_cum_frequency_table(int *arr){
    for (int i = 0; i < NO_OF_SYMBOLS; i++){
        std::cout << i << " " << arr[i] << std::endl;
    }
}

class EncoderModel {
public:
    Encode obj;

    EncoderModel(){}
    

    void train(std::vector<int> vector_from_text){
        obj.in = my_ifstream(vector_from_text);
        obj.out = my_ofstream();
        obj.updating = true;
        obj.encode_streams();
        obj.out.storage.shrink_to_fit();
    }

    std::vector<int> encoding_fixed_freq(std::vector<int> vector_from_text){
        obj.in = my_ifstream(vector_from_text);
        obj.out = my_ofstream();
        obj.updating = false;
        obj.encode_streams();
        obj.out.storage.shrink_to_fit();
        return obj.out.storage;
    }

    std::vector<int> get_encoded_vector(){
        return obj.out.storage;
    }
};

class DecoderModel {
public:
    Decode obj;
    std::vector<int> decode_fixed_freq(std::vector<int> vector_from_text){
        obj.in = my_ifstream(vector_from_text);
        obj.out = my_ofstream();
        obj.updating = false;
        obj.decode_streams();
        obj.out.storage.shrink_to_fit();
        return obj.out.storage;
    }
    DecoderModel(EncoderModel &encoder){
        for (int i = 0; i < NO_OF_SYMBOLS + 1; i++){
            obj.cum_freq[i] = encoder.obj.cum_freq[i];
        }
        for (int i = 0; i < NO_OF_SYMBOLS + 1; i++){
            obj.freq[i] = encoder.obj.freq[i];
        }
        for (int i = 0; i < NO_OF_SYMBOLS; i++){
            obj.index_to_char[i] = encoder.obj.index_to_char[i];
        }
        for (int i = 0; i < NO_OF_CHARS; i++){
            obj.char_to_index[i] = encoder.obj.char_to_index[i];
        }
    }
};


std::vector<int> encode_decode_vector(std::vector<int> vector_from_text){

    /*
    EncoderModel model = EncoderModel();

    //model.train(vector_from_text);

    auto encoded_message = model.encoding_fixed_freq(vector_from_text);

    std::cout << "size of storage " << model.obj.out.storage.size() << std::endl;

    // print cum frequency tables
    std::cout << "cum frequency table encoder" << std::endl;
    print_cum_frequency_table(model.obj.cum_freq);
    std::cout << "end cum frequency table encoder" << std::endl;

    DecoderModel decoder_model = DecoderModel(model);

    std::cout << "cum frequency table decoder"  << std::endl;
    print_cum_frequency_table(decoder_model.obj.cum_freq);
    std::cout << "end cum frequency table decoder" << std::endl;


    // decode data
    Decode obj2;

//    obj2.in = my_ifstream(model.obj.out.storage);
//    obj2.out = my_ofstream();
//    obj2.decode_streams();
//    obj2.out.storage.shrink_to_fit();
    auto decoded_text = decoder_model.decode_fixed_freq(encoded_message);
    return decoded_text;
    */

    return std::vector<int>();
    
}

void test_fun(void){
    // encode the data

    auto vector_from_text = read_vec_from_file("test_file_input.txt");

    auto uncompressed_vector = encode_decode_vector(vector_from_text);

    write_vec_to_file("test_file_output.txt", uncompressed_vector);
}


PYBIND11_MODULE(cpp_string_lookup, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
    m.def("some_fn", &some_fn, "A function which adds two numbers");
    m.def("test_fun", &test_fun, "A function which runs arithmetic encoding");
    py::class_<StringHolder>(m, "StringHolder")
        .def(py::init< std::vector<std::string> >())
        .def(py::init< StringAccumulator >())
        .def("get_index", &StringHolder::get_index)
        .def("lookup", &StringHolder::lookup);
    py::class_<StringPairHolder>(m, "StringPairHolder")
        .def(py::init< std::vector<std::pair<std::string, std::string> > >())
        .def("get_value", &StringPairHolder::get_value)
        .def("contains", &StringPairHolder::contains);
    py::class_<StringAccumulator>(m, "StringAccumulator")
        .def(py::init<long>())
        .def("add", &StringAccumulator::add)
        .def("clean", &StringAccumulator::clean);
    py::class_<CharFrequencyHolder>(m, "CharFrequencyHolder")
        .def(py::init<>())
        .def("incrementFreq", &CharFrequencyHolder::incrementFreq)
        .def("printFreq", &CharFrequencyHolder::printFreq);

}


