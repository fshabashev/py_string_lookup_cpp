#pragma once
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

class my_ofstream : public std::ofstream {
public:
    std::vector<int> storage;
    my_ofstream(const char *filename) : std::ofstream(filename, std::ios::binary) {}
    my_ofstream() : std::ofstream() {
        this->storage = std::vector<int>();
    }
    void put_new(int c) {
        std::cout << "put value = " << c << std::endl ;
        std::ofstream::put(c);
    }
    void put(int c) {
        storage.push_back(c);
    }
};


