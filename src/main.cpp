#include <iostream>
#include <math.h>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "arithmetic_coding/encode.h"
#include "arithmetic_coding/decode.h"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <sstream>

#include "ar_enc/ArithmeticCoder.hpp"
#include "ar_enc/BitIoStream.hpp"
#include "ar_enc/FrequencyTable.hpp"



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

//    void printFreq(){
//        for (auto it = freq_array.begin(); it != freq_array.end(); it++){
//            std::cout << it->first << " " << it->second << std::endl;
//        }
//    }
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

int main_old() {
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
    std::cout << "vector size is " << vec.size() << std::endl;
    for (int i = 0; i < vec.size(); i++){
        std::cout << vec[i] << std::endl;
    }
}

void print_cum_frequency_table(int *arr){
    for (int i = 0; i < NO_OF_SYMBOLS; i++){
        std::cout << i << " " << arr[i] << std::endl;
    }
}

const int NO_OF_SYMBOLS_IN_FREQ_TABLE = 256;

SimpleFrequencyTable calculate_freq(std::vector<int> vec){
    // replace with the actual number of symbols in the input data
    SimpleFrequencyTable freqs(std::vector<uint32_t>(NO_OF_SYMBOLS_IN_FREQ_TABLE + 1, 0));
    freqs.increment(256);  // EOF symbol gets a frequency of 1
    for (int b: vec){
        freqs.increment(static_cast<uint32_t>(b));

        if (!(0 <= b && b <= NO_OF_SYMBOLS_IN_FREQ_TABLE))
            throw std::logic_error("Error char isnt in the 0 to 255 range");

    }
    return freqs;
}

typedef std::pair<SimpleFrequencyTable, std::vector<int>> FreqTableAndVec;


FreqTableAndVec calc_freqs_and_encode(std::vector<int> vec){
    auto freqs = calculate_freq(vec);
    // std::ofstream out(outputFile, std::ios::binary); // need to inherit and redefine put
    std::stringstream out;
    BitOutputStream bout(out);
    try {
        // Write frequency table
        /*
        for (uint32_t i = 0; i < 256; i++) {
            uint32_t freq = freqs.get(i);
            for (int j = 31; j >= 0; j--)
                bout.write(static_cast<int>((freq >> j) & 1));  // Big endian
        }*/
        ArithmeticEncoder enc(32, bout);
        std::vector<int> message_to_encode;
        for (int val: vec) {
            // Read and encode one byte
            int symbol = val;
            if (symbol == EOF)
                break;
            if (!(0 <= symbol && symbol <= 255))
                throw std::logic_error("Assertion error");
            message_to_encode.push_back(symbol);
        }
        // write encoded message to disk
        for (int symbol : message_to_encode){
            //std::cout << "symbol is " << symbol << std::endl;
            enc.write(freqs, static_cast<uint32_t>(symbol));
        }

        enc.write(freqs, 256);  // EOF
        enc.finish();  // Flush remaining code bits
        bout.finish();
        // std::cout << "Encoded message: " << std::endl;
        auto encoded_message = out.str();
        // std::cout << "printing stringstream" << encoded_message << endl;
        return FreqTableAndVec(freqs, std::vector<int>(encoded_message.begin(), encoded_message.end()));

    } catch (const char *msg) {
        std::cerr << msg << std::endl;
        return FreqTableAndVec(freqs, std::vector<int>());
    }
}

std::vector<int> decompress_encoded_message(SimpleFrequencyTable freqs, std::vector<int> encoded_message){
    std::stringstream in(std::string(encoded_message.begin(), encoded_message.end()));
    std::stringstream out;
    BitInputStream bin(in);
    ArithmeticDecoder dec(32, bin);
    while (true) {
        uint32_t symbol = dec.read(freqs);
        if (symbol == 256)  // EOF symbol
            break;
        int b = static_cast<int>(symbol);
        if (std::numeric_limits<char>::is_signed)
            b -= (b >> 7) << 8;
        out.put(static_cast<char>(b));
    }
    std::vector<int> decompressed_message;
    for (auto elem : out.str()){
        decompressed_message.push_back(static_cast<int>(elem));
    }
    return decompressed_message;
}

class SymbolMapper {
public:
    std::map<int, int> symbol_to_index;
    std::map<int, int> index_to_symbol;

    SymbolMapper(std::vector<int> vec){
        std::set<int> unique_symbols;
        for (int val: vec){
            unique_symbols.insert(val);
        }
        int i = 0;
        for (int val: unique_symbols){
            symbol_to_index[val] = i;
            index_to_symbol[i] = val;
            i++;
        }
    }

    std::vector<int> encode(std::vector<int> vec){
        std::vector<int> encoded_vec;
        for (int i = 0; i < vec.size(); i++){
            encoded_vec.push_back(symbol_to_index[vec[i]]);
        }
        return encoded_vec;
    }
    std::vector<int> decode(std::vector<int> vec){
        std::vector<int> decoded_vec;
        for (int i = 0; i < vec.size(); i++){
            decoded_vec.push_back(index_to_symbol[vec[i]]);
        }
        return decoded_vec;
    }

    int size(){
        return symbol_to_index.size();
    }
};



int submain(int argc, char* argv[]){
    // Handle command line arguments
    const char *inputFile  = argv[1];
    // const char *outputFile = argv[2];

    // Read input file once to compute symbol frequencies
    std::ifstream in(inputFile, std::ios::binary);

    std::vector<int> vec;
    while (true) {
        int b = in.get();
        if (b == EOF)
            break;
        if (b < 0 || b > 255)
            throw std::logic_error("Assertion error");
        vec.push_back(b);
    }

    auto freqs = calculate_freq(vec);

    try {
        auto freq_and_encoding = calc_freqs_and_encode(vec);
        auto decompressed_data = decompress_encoded_message(freq_and_encoding.first, freq_and_encoding.second);
        // print decompressed data

//        for (int i = 0; i < decompressed_data.size(); i++){
//            std::cout << decompressed_data[i];
//        }
        return EXIT_SUCCESS;

    } catch (const char *msg) {
        std::cerr << msg << std::endl;
        return EXIT_FAILURE;
    }
}


std::vector<int> string2vec(std::string str){
    std::vector<int> vec;
    for (auto elem : str){
        vec.push_back(static_cast<int>(elem));
    }
    return vec;
}

class EncodingModel{
public:
    SimpleFrequencyTable freqs;

    EncodingModel(std::vector<int> vec): freqs(calculate_freq(vec)){}


    void learn_freqs(std::vector<int> vec){
        auto freq_and_enc = calc_freqs_and_encode(vec);
        this->freqs = freq_and_enc.first;
    }

    void learn_freqs_str(std::string str){
        std::vector<int> vec = string2vec(str);
        this->learn_freqs(vec);
    }

    std::vector<int> string2vec(std::string str){
        std::vector<int> vec;
        for (auto elem : str){
            vec.push_back(static_cast<int>(elem));
        }
        return vec;
    }

    std::vector<int> compress_vec(std::vector<int> vec){
        auto freq_and_enc = calc_freqs_and_encode(vec);
        return freq_and_enc.second;
    }

    std::vector<int> compress(std::string str){
        return compress_vec(this->string2vec(str));
    }

    std::vector<int> decompress(std::vector<int> vec){
        return decompress_encoded_message(this->freqs, vec);
    }

};
void test_speed_vec(void) {
    for (int i = 0; i < 100; i++) {
        std::vector<int> vec;
        for (int j = 0; j < 1000000; j++) {
            vec.push_back(42);
        }
    }
}

class Compressor{
public:
    EncodingModel enc_model;
    SymbolMapper sym_mapper;
    Compressor(std::vector<int> vec): sym_mapper(vec), enc_model(std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}){
        auto encoded_vec = this->sym_mapper.encode(vec);
        this->enc_model.learn_freqs(encoded_vec);
    }
    Compressor(std::string str): Compressor(string2vec(str)){}
    std::vector<int> compress(std::string str){
        auto encoded_vec = sym_mapper.encode(this->enc_model.string2vec(str));
        auto compressed_vec = this->enc_model.compress_vec(encoded_vec);
        return compressed_vec;
    }
    std::string decompress(std::vector<int> vec){
        auto decompressed_vec = this->enc_model.decompress(vec);
        auto decoded_vec = sym_mapper.decode(decompressed_vec);
        std::string str;
        for (auto elem : decoded_vec){
            str += static_cast<char>(elem);
        }
        return str;
    }
};

void test_fun(void) {
    // submain(3, (char *[]){"./test", "test.txt", "test_out.txt"});

    auto enc = EncodingModel(std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10});
    // auto str = "In 2015, ывафыавфыав bioinformatician asdfl;askjfdaslkjwqeioprjak;jsdfla;jkoiwqjrek;lasmdf laskdjfowiqjer;lkamsdfl;kjas;dlkfjsa;ldfj . But in general UTF16 includes surrogate pairs, so a unicode code point cannot be represented with a single wide character. You need wide string instead. Your problem is also partly to do with printing UTF16 character in Windows console. If you use MessageBoxW to view a wide string it will work as expected ";
    auto str = "1234567890asdfasfфвавыфафы";
    // convert string to vector of ints
    auto str_in_vector = string2vec(str);
    // remap ints into the smaller range
    auto mapper = SymbolMapper(str_in_vector);
    auto encoded_str = mapper.encode(str_in_vector);
    // learn the frequencies of the encoded string

    //std::cout << "string size = " << std::string(str).size() << std::endl;
    // print string before encoding
    //std::cout << "string before encoding = " << std::endl;
//    for (int i = 0; i < str_in_vector.size(); i++){
//        std::cout << str_in_vector[i] << " ";
//    }
    //std::cout << "string after encoding = " << std::endl;
    // print encoded string
    /*
    for (int i = 0; i < encoded_str.size(); i++){
        std::cout << encoded_str[i] << " ";
    }
    */
    enc.learn_freqs(encoded_str);
    auto compressed = enc.compress_vec(encoded_str);
    //std::cout << "size of compressed is " << compressed.size() << std::endl;
    auto decompressed = enc.decompress(compressed);
    // remap back to original ints
    auto decoded_str = mapper.decode(decompressed);

    //std::cout << "decompressed is " << std::endl;

    for (auto elem : decoded_str){
        std::cout << static_cast<char> (elem);
    }
}

void test_fun2(void){
    auto str = "1234567890фывафыва";
    auto comp = Compressor(string2vec(str));
    auto compressed = comp.compress(str);
    auto decompressed = comp.decompress(compressed);
    std::cout << decompressed << std::endl;
}

int main(int argc, char* argv[]) {
    // run test_fun and measure time
    auto start = std::chrono::high_resolution_clock::now();
    test_fun2();
    //test_speed_vec();
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count() ;
    return 0;
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
        .def("incrementFreq", &CharFrequencyHolder::incrementFreq);
    py::class_<EncodingModel>(m, "EncodingModel")
        .def(py::init< std::vector<int> >())
        .def("learn_freqs", &EncodingModel::learn_freqs)
        .def("compress", &EncodingModel::compress)
        .def("decompress", &EncodingModel::decompress);
    py::class_<Compressor>(m, "Compressor")
        .def(py::init< std::vector<int> >())
        .def(py::init< std::string >())
        .def("compress", &Compressor::compress)
        .def("decompress", &Compressor::decompress);
}


