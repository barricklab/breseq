/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2017 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the
  terms the GNU General Public License as published by the Free Software
  Foundation; either version 1, or (at your option) any later version.

 *****************************************************************************/


#ifndef _BRESEQ_STORABLE_H_
#define _BRESEQ_STORABLE_H_

#include "common.h"

#include "json.hpp"

using namespace std;
using namespace nlohmann;

namespace breseq
{
  
  // This class stores things in
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
			char input[size+1];
			infile.read(input, size * sizeof(char));
      input[size] = '\0'; // need to null terminate
			output = string(input);
		}
    
    static void write_to_file(ofstream& outfile, Storable& s)
		{
      s.serialize(outfile);
		}
		static void read_from_file(ifstream& infile, Storable& s)
		{
      s.deserialize(infile);
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

    // classes that inherit from this must define how to write and read themselves
    virtual void serialize(ofstream& f) = 0;
    virtual void deserialize(ifstream& f) = 0;
    
    void store(string filename)
    {
      ofstream outfile(filename.c_str());
      assert(!outfile.fail());
      //ASSERT(!outfile.fail(), "Error storing in file: " + filename);
      serialize(outfile);
      outfile.close();    
    }
    
    void retrieve(string filename)
    {
      ifstream infile(filename.c_str());
      assert(!infile.fail());
      //ASSERT(!infile.fail(), "Error retrieving from file: " + filename);
      deserialize(infile);
      infile.close();
    }

    // static convenience functions to store/retrieve maps
		template<typename T, typename U> static void store(map<T,U>& input, string filename)
		{
			ofstream outfile(filename.c_str());
      outfile << setprecision(10) << scientific;
			write_to_file(outfile, input);
			outfile.close();
		}
		template<typename T, typename U> static void retrieve(map<T,U>& output, string filename)
		{
			ifstream infile(filename.c_str());
			read_from_file(infile, output);
			infile.close();
		}
    


	}; // class Storable
  
  
  template <typename T,typename U> class storable_map : public Storable, public map<T,U>
  {
  public: 
    void serialize(ofstream& outfile)
		{
			uint32_t elements = map<T,U>::size();
			write_to_file(outfile, elements);
			for (typename map<T,U>::iterator it = map<T,U>::begin(); it != map<T,U>::end(); it++)
			{
				T t = it->first;
				U& u = it->second;
				write_to_file(outfile, t);
				u.serialize(outfile);
			}
		}
		void deserialize(ifstream& infile)
		{
			uint32_t size;
			read_from_file(infile, size);
			for (uint32_t i = 0; i < size; i++)
			{
				T t; read_from_file(infile, t);
				U u;
        u.deserialize(infile);
				map<T,U>::insert(pair<T,U>(t,u));
			}
		}
    
  };
  
  template <typename T> class storable_vector : public Storable, public vector<T> 
  {
  public: 
    
    void serialize(ofstream& outfile)
		{
			uint32_t elements = vector<T>::size();
			write_to_file<uint32_t>(outfile, elements);
			for (typename vector<T>::iterator it = vector<T>::begin(); it != vector<T>::end(); it++)
			{
				write_to_file(outfile, *it);
			}
		}
		void deserialize(ifstream& infile)
		{
			uint32_t size;
			read_from_file(infile, size);
      vector<T>::resize(size);
			for (uint32_t i = 0; i < size; i++)
			{
				T t; read_from_file(infile, t);
				(*this)[i] = t;
			}
		}
    
  };
  
  template <class C> class JSONStorable
  {
  public:
    
    virtual void store(const string& filename)
    {
      ofstream outfile(filename.c_str());
      assert(!outfile.fail());
      //ASSERT(!outfile.fail(), "Error storing in file: " + filename);
      json j;
      to_json(j, static_cast<C&>(*this));
      outfile << j.dump(2);
      outfile.close();
    }
    
    virtual void retrieve(const string& filename)
    {
      ifstream infile(filename.c_str());
      assert(!infile.fail());
      //ASSERT(!infile.fail(), "Error retrieving from file: " + filename);
      json j;
      infile >> j;
      from_json(j, dynamic_cast<C&>(*this));
      infile.close();
    }
  };
  
  // Summary
 // inline void to_json(json& j, const JSONStorable& s) {(void)j; (void)s;};
 // inline void from_json(const json& j, JSONStorable& s) {(void)s; (void)j;};

  
} // breseq namespace

#endif
