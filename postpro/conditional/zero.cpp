
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


// TODO: create a new main which calls the new function


int main() // TODO: turn this into a function find_zero_crossings(Ul, time) -> dumps everything to file + returns (array_pos and array_neg, size_pos, size_neg)
// when you dump: take the name of the input file (maybe you need pass name of the file as argument as well)
// ul_bottom_101.out -> ul_bottom_101.npzc *I'm making this up*
{
	std::cout << " pleaee work";
	//This reads in the large scale velocity signals from a binary file and write it into an array
	std::ifstream infile_a("ul.out", std::ifstream::binary);

	infile_a.seekg(0, infile_a.end);
	int P = infile_a.tellg();
	infile_a.seekg(0, infile_a.beg);

	std::vector<double> Ul(P / sizeof(double));
	infile_a.read(reinterpret_cast<char*>(Ul.data()), static_cast<std::streamsize>(Ul.size()) * sizeof(double));

	std::cout << Ul.size();

	const int num = Ul.size() / DIMZ;

	//std::cout << num;
	//using Matrix = std::array<std::vector<double>, 2>;

	//td::vector < std::vector<double> > ul;
  //ul.resize((std::vector<double>(DIMZ), num));

  //ul.resize( DIMZ, (std::vector<double>(num)));

	Matrix_dimz_double ul; // ul[j][i]



  for (int j = 0; j < DIMZ; ++j)
  {
	  // debug
	  std::cout << "hvdkjkjdk" << ul[j].size() << std::endl;

  }

  // TODO: suggestion - this is pseudocode
  // define a type like timevec_type = type(vector<double num>)
  // define the array as array<timevec_type DIMZ>
  // this hopefully avoids using snippet on github

  //ul.resize((std::vector<double>(DIMZ), num));


	// transposing the array ul = Ul^T
	for (int j = 0; j < DIMZ; ++j)
	{
		for (int i = 0; i < num; ++i)
		{

			ul[j].push_back(Ul.at((i * DIMZ) + j));

			//std::cout << ul[j].at(i) << std::endl;
			//std::cout << "fhb" << ul[j].at(i) << std::endl;
		}

	}
	
	for (int i = 0; i < ul[0].size(); i++)
	{
		// debug
		std::cout << ul[100].at(i) << std::endl;
	}

	// int num_steps;
   //std::ifstream infilt_txt("num_timesteps.txt");
   //infilt_txt >> num_steps;





	//This reads in the time signals from a binary file and write it into an array
	std::ifstream infile_t("time.out", std::ifstream::binary);

	
	infile_t.seekg(0, infile_a.end);
	int T = infile_t.tellg();
	infile_t.seekg(0, infile_t.beg);

	std::vector<double> t(T / sizeof(double));
	infile_t.read(reinterpret_cast<char*>(t.data()), static_cast<std::streamsize>(t.size()) * sizeof(double));

	std::cout << t.size();
	
	
	/*std::vector<double> t;
	double dt;
	std::ifstream infilt_txt("time_step.txt");
	infilt_txt >> dt;
	//for (int i = 0; i < ul[0].size(); i++)
	for (int i = 0; i < num; i++)
	{

		t.push_back(dt * i);
	}
	*/



	int method; // the method you wish to use


	std::ifstream infile_txt("selection.txt"); // read in the method assigned

	infile_txt >> method;

	//If "2" is chosen, the threshold is also read in from a file 
	double threshold;
	std::ifstream infil_txt("criterion.txt");
	infil_txt >> threshold;

	// TODO: check if there is some __easy__ c++ library for JSON or TOML -> all inputs in one file!


	Matrix_dimz array_pos;
	Matrix_dimz array_neg;

	std::vector<int> size_pos; // TODO: maybe use std::array with size DIMZ

	std::vector<int> size_neg; // TODO: maybe use std::array with size DIMZ

	std::cout << ul[0].size(); // debug line
	for (int z = 0; z < DIMZ; z++)
	{
		

		int i = 1; // TODO: this can be 0, I think!
		// to find a zero-crossing, you need two points!
		//positions for positive to negative is written to this array

		while (i < (ul[0].size() - 1)) // TODO: ul[0].size() is DIMT or num
		{ // 'at' is for assessing an array element
			// TODO: array.at() controls that you're not running out of bounds
			// SO: array.at() is slower than array[], but safer
			// BUT: we are crazy, we do not want safe code
			// try to use [], if not working, .at() is also ok
			if (((ul[z].at(i)) > 0) && ((ul[z].at(i + 1)) < 0))
			{
				//for every zero-crossing found, the positions is written to the array "array_pos"
				array_pos[z].push_back(i);


			}
			//std::cout << array_pos[z].size() << std::endl;
			i++;

		}

		//After this while loop, we have all the zero-crossings for positive to negative written in "array_pos"

		int j = 1; // TODO: this can be zero
		// to avoid the first "zero-crossing at time = 0"



		while (j < (ul[0].size() - 1 )) // TODO: check corrections from previous loop
		{
			if (((ul[z].at(j)) < 0) && ((ul[z].at(j + 1)) > 0))
			{
				array_neg[z].push_back(j);

			}
			//std::cout << array_neg[z].size() << std::endl;
			j++;
		}
	
		// To define the method for the zero-crossing
		//1 -> Use all zero-crossings
		//2 -> implement a threshold system ( use some zero-crossings)


		// TODO: I would "move" this whole switch to the function/program that calculates the average
		// WHY: I want this program/function to calculate and save to disk all zero crossings
		// ALSO: it's maybe easier to write a code that selects only zero-crossings that are good for 2nd method inside the loop that iterates over zero-crossings

		// pseudo-code: in the "cut1" program, you have some loop over all zero crossings in z
		// loop for i in "list of zero crossings at a given z-position" {
		//		if (this_is_last_zerocrossing_at_this_z OR this_is_last_zerocrossing_at_this_z) continue;
		//		if ( delta_t_left < threshold   OR   delta_t_right < threshold ) {
		//  		continue
		//  	}
		// }
		switch (method)
		{

		case 1:

			break;

		case 2:

			
			if( array_pos[z].size() < 3 )

			{ 
				//std::cout << z << array_pos[z].size() << std::endl;
				
				array_pos[z].clear();

			
			}
			

			else if (array_neg[z].size() < 3)

			{
				array_neg[z].clear();

			}


			else
			 {

               std::vector<int> array_pos_;
			   std::vector<int> array_neg_;
				array_pos_ = array_pos[z];
				array_neg_ = array_neg[z];

				array_pos[z].clear();
				array_neg[z].clear();

				for (int i = 1; i < (array_pos_.size() - 1); i++)
				{

					double precedent_gap = t[array_pos_.at(i)] - t[array_pos_.at(i - 1)];
					//std::cout << "gu" <<  precedent_gap << std::endl;
					double successive_gap = (t[array_pos_.at(i + 1)] - t[array_pos_.at(i)]);
					//std::cout << "d" << successive_gap << std::endl;
					if ((precedent_gap > (threshold * 1.0)) && (successive_gap > (threshold * 1.0)))

					{
						
						array_pos[z].push_back(array_pos_.at(i));

					}


				}

				for (int j = 1; j < (array_neg_.size() - 1); j++)
				{
					double precedent_gap = t[array_neg_.at(j)] - t[array_neg_.at(j - 1)];
					double successive_gap = (t[array_neg_.at(j + 1)] - t[array_neg_.at(j)]);

					if ((precedent_gap > (threshold * 1.0)) && (successive_gap > (threshold * 1.0)))
					{
						array_neg[z].push_back(array_neg_.at(j));

					}

				}


			}


		}

		
		//std::cout << z << array_pos[z].size() << std::endl;


		for (int i = 0; i < array_neg.size(); i++)
		{
			// this was debug
			//std::cout << array_neg[z].at(i) << std::endl;
		}

		// TODO: these two arrays should not be vectros anymore
		// do not use push back
		// use size_pos[...] instead
		size_pos.push_back(array_pos[z].size());
		size_neg.push_back(array_neg[z].size());


	std::cout << "uzt" << size_pos.at(z) << std::endl; // debug line

	}

	for (int i = 0; i < array_pos[100].size(); i++) // debug line
	{

		//std::cout << ul[100].at(array_pos[100].at(i));
		//std::cout << ul[100].at(array_pos[100].at(i) + 1);
		//std::cout << std::endl;


	}

	//int sum = std::accumulate(size_pos.begin(), size_pos.end(), 0);


		//std::cout <<" f" << array_pos.size() << std::endl;


	// TODO: this should not be needed, if you change the code that saves the data
	std::vector<int> array_neg_f;
	std::vector<int> array_pos_f;

	for (int z = 0; z < DIMZ; z++)
	{
		for (int i = 0; i < array_neg[z].size(); i++)
		{
			array_neg_f.push_back(array_neg[z].at(i)); // TODO: do not use .at()

		}


	}


	for (int z = 0; z < DIMZ; z++)
	{
		for (int i = 0; i < array_pos[z].size(); i++)
		{
			array_pos_f.push_back(array_pos[z].at(i));
			//std::cout << array_pos[z].at(i) << std::endl;
		}

	}


	//std::cout << "sum" << sum;
	//std::cout << "sum 1" << array_pos_f.size();

	for (int z = 0; z < array_neg_f.size(); z++)

	{

		std::cout << "fdh" << array_neg_f.at(z) << std::endl;


	}

	

	/*
	// This writes array_neg to a binary file
	std::ofstream out_neg;
	out_neg.open("neg_pos.out");

	if (out_neg)
	{
		out_neg.write(reinterpret_cast<char*>(array_neg_f.data()), static_cast<std::streamsize>(array_neg_f.size()) * sizeof(int));
		out_neg.close();

	}

	*/

	// TODO: fix this whole section - find a method to write like this: array_pos_f[i][j]

	std::fstream time_file_P;
	time_file_P.open("pos_neg.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < array_pos_f.size(); i++)
	{
		time_file_P.write((char*)(&array_pos_f[i]), sizeof(int));
	}




	std::fstream time_file_N;
	time_file_N.open("neg_pos.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < array_neg_f.size(); i++)
	{
		time_file_N.write((char*)(&array_neg_f[i]), sizeof(int));
	}



	std::fstream  out_pos_size;
	out_pos_size.open("size_neg_pos.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < size_neg.size(); i++)
	{
		out_pos_size.write((char*)(&size_neg[i]), sizeof(int));
	}

	std::fstream  out_neg_size;
	out_neg_size.open("size_pos_neg.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < size_pos.size(); i++)
	{
		out_neg_size.write((char*)(&size_pos[i]), sizeof(int));
	}

	





	return 0; // TODO: return array_pos and array_neg, size_pos, size_neg

}