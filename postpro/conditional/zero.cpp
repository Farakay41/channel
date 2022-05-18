

#include <iostream>
#include <fstream> 
#include <vector> 
#include <numeric> 
#include <array>
#include <algorithm>
#include <string>
#include <stdio.h>

// This program writes out the positions of zero-crossings in a large-scale velocity signal

constexpr auto DIMZ = 2048; // TODO: READ THIS VALUE FROM A FILE! If too difficult: leave it hard coded
// TODO: num is too generic as a variable name -> call it DIMT

using Matrix_dimz = std::array<std::vector<int>, DIMZ>; // we are going to increase the size of this array at runtime, but this is good!
using Matrix_dimz_double = std::array<std::vector<double>, DIMZ>;

using Matrix_dimz_2 = std::array<Matrix_dimz, 2>;

using Matrix_2 = std::array < std::array<int, DIMZ>, 2 >;


using Pair = std::pair<Matrix_dimz_2, Matrix_2 >;
// TODO: create a new main which calls the new function



std::vector<int> read_int(std::string filename)
{

	std::ifstream infile_a(filename, std::ifstream::binary);

	infile_a.seekg(0, infile_a.end);
	int N = infile_a.tellg();
	infile_a.seekg(0, infile_a.beg);

	std::vector<int> array(N / sizeof(int));
	infile_a.read(reinterpret_cast<char*>(array.data()), static_cast<std::streamsize>(array.size()) * sizeof(int));

	return array;

}

std::vector<double> read_double(std::string filename)
{

	std::ifstream infile_a(filename, std::ifstream::binary);

	infile_a.seekg(0, infile_a.end);
	int N = infile_a.tellg();
	infile_a.seekg(0, infile_a.beg);

	std::vector<double> array(N / sizeof(double));
	infile_a.read(reinterpret_cast<char*>(array.data()), static_cast<std::streamsize>(array.size()) * sizeof(double));

	return array;

}


Pair find_zero_crossings(std::string ul_filename)

{
	// TODO: turn this into a function find_zero_crossings(Ul, time) -> dumps everything to file + returns (array_pos and array_neg, size_pos, size_neg)
	// when you dump: take the name of the input file (maybe you need pass name of the file as argument as well)
	// ul_bottom_101.out -> ul_bottom_101.npzc *I'm making this up*


		//This reads in the large scale velocity signals from a binary file and write it into an array

	std::vector<double> Ul;
	//std::vector<double> t;

	Ul = read_double(ul_filename);
	//t = read_double(time_filename);

	std::string file;

	for (int i = 0; i < ul_filename.size() - 4; ++i)
	{
		file = file + ul_filename[i];

	}

	const int DIMT = Ul.size() / DIMZ;
	//const int DIMT = t.size();

	Matrix_dimz_double ul; // ul[j][i]

	int total_size = DIMT * DIMZ;

	double* ul_flat = new double[total_size];



	for (int j = 0; j < total_size; ++j)
	{

		ul_flat[j] = Ul[j];

	}

	std::cout << "velocity" << ul_flat[5655];

	std::cout << "velocity" << Ul[5655]; // check if they are equal

#define ul(y,x) *(ul_flat + x*DIMZ + y)

	//std::cout << "HGCVL" << ul(1779, 567);


	Matrix_dimz array_pos;
	Matrix_dimz array_neg;

	//std::vector<int> size_pos; // TODO: maybe use std::array with size DIMZ
	std::array < int, DIMZ > size_pos;

	std::array < int, DIMZ > size_neg; // TODO: maybe use std::array with size DIMZ



	for (int z = 0; z < DIMZ; z++)
	{


		int i = 100; // TODO: this can be 0, I think! // changed this to 100 because of the cutting for the velocity
		//forgot to discuss this in the last meeting
		// to find a zero-crossing, you need two points!
		//positions for positive to negative is written to this array

		while (i < (DIMT - 100)) // TODO: ul[0].size() is DIMT or num
		{ 

			if (((ul(z, i)) > 0) && ((ul(z, (i + 1))) < 0))
			{
				//for every zero-crossing found, the positions is written to the array "array_pos"
				array_pos[z].push_back(i);


			}
			//std::cout << array_pos[z].size() << std::endl;
			i++;

		}

		//After this while loop, we have all the zero-crossings for positive to negative written in "array_pos"

		int j = 100; // TODO: this can be zero
		// to avoid the first "zero-crossing at time = 0"


		while (j < (DIMT - 100)) // TODO: check corrections from previous loop
		{
			if (((ul(z, j)) < 0) && ((ul(z, (j + 1))) > 0))
			{
				array_neg[z].push_back(j);

			}
			//std::cout << array_neg[z].size() << std::endl;
			j++;
		}

		//After this while loop, we have all the zero-crossings for negative to positive written in "array_pos"

		// TODO: these two arrays should not be vectros anymore
		// do not use push back
		// use size_pos[...] instead

		size_pos[z] = (array_pos[z].size());
		size_neg[z] = (array_neg[z].size());


		std::cout << "uzt" << size_pos[z] << std::endl; // debug line

	}



	// TODO: fix this whole section - find a method to write like this: array_pos_f[i][j]


	std::fstream pos_file_P;
	pos_file_P.open(file + ".pnzc", std::ios::out | std::ios::binary);
	for (int z = 0; z < DIMZ; z++)
	{
		for (int j = 0; j < array_pos[z].size(); j++)
		{

			pos_file_P.write((char*)(&array_pos[z][j]), sizeof(int));

		}
	}



	std::fstream neg_file_P;
	neg_file_P.open(file + ".npzc", std::ios::out | std::ios::binary);
	for (int z = 0; z < DIMZ; z++)
	{
		for (int j = 0; j < array_neg[z].size(); j++)
		{

			neg_file_P.write((char*)(&array_neg[z][j]), sizeof(int));

		}
	}


	std::fstream  out_neg_size;
	out_neg_size.open(file + "_size.pnzc", std::ios::out | std::ios::binary);
	for (int i = 0; i < size_pos.size(); i++)
	{
		out_neg_size.write((char*)(&size_pos[i]), sizeof(int));
	}


	std::fstream  out_pos_size;
	out_pos_size.open(file + "_size.npzc", std::ios::out | std::ios::binary);
	for (int i = 0; i < size_neg.size(); i++)
	{
		out_pos_size.write((char*)(&size_neg[i]), sizeof(int));
	}



	return { {array_pos, array_neg} , {size_pos, size_neg} }; // TODO: return array_pos and array_neg, size_pos, size_neg

}

int main()

{
	Pair ULE;

	ULE = find_zero_crossings("ul.out");


	std::cout << " check" << ULE.first[0][16].size() << std::endl;

	std::cout << "another check" << ULE.second[1].size();


	return 0;


}