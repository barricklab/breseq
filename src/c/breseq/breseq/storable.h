/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the
  terms the GNU General Public License as published by the Free Software
  Foundation; either version 1, or (at your option) any later version.

 *****************************************************************************/


#ifndef _BRESEQ_STORABLE_H_
#define _BRESEQ_STORABLE_H_

#include "common.h"

using namespace std;

namespace breseq
{
	class Storable
	{
	  protected:

		Storable() {};
		~Storable() {};

		static void write_to_file(ofstream& outfile, string& input)
		{
			uint32_t size = input.size();
			outfile.write(reinterpret_cast<char*>(&size), sizeof(uint32_t));
			const char* output = input.c_str();
			outfile.write(output, size * sizeof(char));
		}
		static void read_from_file(ifstream& infile, string& output)
		{
			uint32_t size;
			infile.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
			char input[size];
			infile.read(input, size * sizeof(char));
			output = string(input);
		}
		template<typename T> static void write_to_file(ofstream& outfile, T& input)
		{
			outfile.write(reinterpret_cast<char*>(&input), sizeof(T));
		}
		template<typename T> static void read_from_file(ifstream& infile, T& output)
		{
			infile.read(reinterpret_cast<char*>(&output), sizeof(T));
		}
		template<typename T, typename U> static void write_to_file(ofstream& outfile, map<T,U>& input)
		{
			uint32_t size = input.size();
			write_to_file(outfile, size);
			for (typename map<T,U>::iterator it = input.begin(); it != input.end(); it++)
			{
				T t = it->first;
				U u = it->second;
				write_to_file(outfile, t);
				write_to_file(outfile, u);
			}
		}
		template<typename T, typename U> static void read_from_file(ifstream& infile, map<T,U>& output)
		{
			uint32_t size;
			read_from_file(infile, size);
			for (uint32_t i = 0; i < size; i++)
			{
				T t; read_from_file(infile, t);
				U u; read_from_file(infile, u);
				output[t] = u;
			}
		}

	  public:

		virtual void store(string filename) = 0;
		template<typename T, typename U> static void store(map<T,U>& input, string filename)
		{
			ofstream outfile(filename.c_str());
			write_to_file(outfile, input);
			outfile.close();
		}
		virtual void retrieve(string filename) = 0;
		template<typename T, typename U> static void retrieve(map<T,U>& output, string filename)
		{
			ifstream infile(filename.c_str());
			read_from_file(infile, output);
			infile.close();
		}

	}; // class Storable

} // breseq namespace

#endif