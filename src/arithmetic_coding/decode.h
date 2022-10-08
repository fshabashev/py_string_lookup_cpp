#pragma once
#include "compress.h"
#include "my_ifstream.h"
#include "my_ofstream.h"

class Decode : public Compress
{
	int low, high;
	int value;

	int buffer;
	int	bits_in_buf;
	bool end_decoding;

public:
    my_ifstream in;
    my_ofstream out;

    Decode(void);
	~Decode(void);

	void load_first_value(void);
	void decode(char *infile, char *outfile);
    void decode_streams(void);
	int decode_symbol(void);
	int get_bit(void);
};

