#pragma once

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

class my_ifstream : public std::ifstream {
public:
    std::vector<int> storage;
    int index;
    my_ifstream(const char *filename) : std::ifstream(filename, std::ios::binary) {
        index = 0;
        storage = std::vector<int>();
    }
    my_ifstream() : std::ifstream() {
        this->storage = std::vector<int>();
        index = 0;
    }

    my_ifstream(std::vector<int> storage) : std::ifstream() {
        this->storage = storage;
        index = 0;
    }

    int get_old() {
        auto c = std::ifstream::get();
        return c;
    }

    int get() {
        auto c = storage[index];
        std::cout << " c =  " << c << " ";
        if (index < storage.size()) {
            return storage[index++];
        }
        return -1;
    }

    bool eof() {
        return index >= storage.size();
    }
};

