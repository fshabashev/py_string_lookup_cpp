#pragma once
#include "compress.h"
#include "my_ifstream.h"
#include "my_ofstream.h"


class Encode : public Compress
{
	int low, high;
	int opposite_bits;
	int buffer;
	int	bits_in_buf;

public:
    my_ifstream in;
    my_ofstream out;

    Encode(void);
	~Encode(void);
	
	void write_bit( int bit);
	void output_bits(int bit);
	void end_encoding(void);
	void encode_symbol(int symbol);
	void encode(char *infile, char *outfile);
    void encode_streams(void);

};

