
#include <iostream>
#include <fstream> 
#include <vector> 
#include <array>
#include <algorithm>


using Matrix = std::vector < std::vector<double> >; // This is just an alias for 2D vector

//This is the main function for calculating the averaged velocity signals across the zero-crossings
std::vector<double> cut_velocity(std::vector <double> large, std::vector <double> small, std::vector <int> positions, int n)

{



	//std::cout << n;  //number of velocity signals each right and left) of the zero crossing

   // std::vector<double> cut_vel;// the array takes in the cut small scale velocity for each zero crossing positions

	int size = positions.size(); // This calculates the size of the positions (of zero-crossing) array, Also needed to calculate the average velocity signal later

	//vector w contains the calculated weights

	std::vector<double> w;

	for (int i = 0; i < size; ++i)

	{
		int pos = positions.at(i);

		w.push_back(large.at(pos) / (large.at(pos) - large.at(pos + 1)));
		//std::cout << w.at(i) << std::endl;
	}

	std::vector<double> cut_vel;

	//cut_vel.resize((std::vector<double>(size), ((2 * n) - 1)));

	// the matrix takes in the cut small scale velocity for each zero crossing position
	for (size_t i = 0; i < size; i++) // for each of the zero-crossing positions
	{

		int  temp = positions.at(i); // allocates a variable for each zero crossing positions


		for (size_t j = 0; j < ((2 * n) - 1); j++) // n is the number of small scale velocity signals and time cut for each zero-crossing position

		{
			int loc = (temp - (n - 1) + j);


			double vel = ((1 - w.at(i)) * small.at(loc)) + (w.at(i) * small.at(loc + 1));


			cut_vel.push_back(vel);

		}

	}

	return cut_vel;
}


std::vector <double> average_cut_velocity(std::vector<double> velocity, std::vector <int> positions, int n)
{

	std::vector <double> averaged;
	std::vector <double> vel_sum;

	int size = positions.size();
	//std::cout << velocity.size();

	for (int i = 0; i < ((2 * n) - 1); i++)
	{
		double vel = 0.0;

		for (int j = 0; j < size; j++)
		{

			vel = vel + velocity.at((((2 * n) - 1) * j) + i);


		}
		averaged.push_back(vel / size);
		// std::cout << averaged.size();

	}

	return averaged;

}



int main()
{
	//Read in the positions for positive to negative zero crossings

	std::ifstream infile("pos_neg.out", std::ifstream::binary);

	infile.seekg(0, infile.end);
	int P = infile.tellg();
	infile.seekg(0, infile.beg);

	std::vector<int> array_pos(P / sizeof(int));
	infile.read(reinterpret_cast<char*>(array_pos.data()), static_cast<std::streamsize>(array_pos.size()) * sizeof(double));


	//Read in the positions for negative to positive zero crossings

	std::ifstream infile_a("neg_pos.out", std::ifstream::binary);

	infile_a.seekg(0, infile_a.end);
	int N = infile_a.tellg();
	infile_a.seekg(0, infile_a.beg);

	std::vector<int> array_neg(N / sizeof(int));
	infile_a.read(reinterpret_cast<char*>(array_neg.data()), static_cast<std::streamsize>(array_neg.size()) * sizeof(int));


	//Read in the the small-scale velocity signals

	std::ifstream infile_s("us.out", std::ifstream::binary);

	infile_s.seekg(0, infile_s.end);
	int M = infile_s.tellg();
	infile_s.seekg(0, infile_s.beg);

	std::vector<double> us(M / sizeof(double));
	infile_s.read(reinterpret_cast<char*>(us.data()), static_cast<std::streamsize>(us.size()) * sizeof(double));


	//This reads in the large scale velocity signals from a binary file and write it into an array
	std::ifstream infile_l("ul.out", std::ifstream::binary);

	infile_l.seekg(0, infile_a.end);
	int V = infile_l.tellg();
	infile_l.seekg(0, infile_a.beg);

	std::vector<double> ul(V / sizeof(double));
	infile_l.read(reinterpret_cast<char*>(ul.data()), static_cast<std::streamsize>(ul.size()) * sizeof(double));


	//Matrix cut_us_P; // cut small velocity and time for positive to negative
	//Matrix cut_us_N; // cut small velocity and time for negative to positive
	std::vector<double> cut_us_P;
	std::vector<double> cut_us_N;
	int num_of_cuts; // number of velocity signals in the right(or left) of the zero crossing


	// the num_of_cuts variable is read in from a text file
	std::ifstream infile_txt("input.txt");

	infile_txt >> num_of_cuts;

	//call the cut_velocity function to calculate the averaged velocity signals (positive to negative)
	cut_us_P = cut_velocity(ul, us, array_pos, num_of_cuts);

	//call the cut_velocity function to calculate the averaged velocity signals (negative to positive)
	cut_us_N = cut_velocity(ul, us, array_neg, num_of_cuts);


	std::vector<double> cut_us_Pos;

	std::vector<double> cut_us_Neg;

	cut_us_Pos = average_cut_velocity(cut_us_P, array_pos, num_of_cuts);
	cut_us_Neg = average_cut_velocity(cut_us_N, array_neg, num_of_cuts);


	//write the result to a textfile
	std::ofstream fileOut_P;
	fileOut_P.open("Velocity_pos_neg.txt");

	for (int i = 0; i < cut_us_Pos.size(); i++)
	{
		fileOut_P << cut_us_Pos.at(i) << std::endl;
	}
	fileOut_P.close();


	std::ofstream fileOut_N;
	fileOut_N.open("Velocity_neg_pos.txt");

	for (int i = 0; i < cut_us_Neg.size(); i++)
	{
		fileOut_N << cut_us_Neg.at(i) << std::endl;
	}
	fileOut_N.close();

	
	//Write the unaveraged velocities into .out files

	std::ofstream out_pos;
	out_pos.open("data_pos_neg.out");

	if (out_pos)
	{
		out_pos.write(reinterpret_cast<char*>(cut_us_P.data()), static_cast<std::streamsize>(cut_us_P.size()) * sizeof(double));
		out_pos.close();

	}


	std::ofstream out_neg;
	out_pos.open("data_neg_pos.out");

	if (out_neg)
	{
		out_neg.write(reinterpret_cast<char*>(cut_us_N.data()), static_cast<std::streamsize>(cut_us_N.size()) * sizeof(double));
		out_neg.close();

	}

	/*
	std::vector<double> arr_neg;

	for (int y = 0; y < array_neg.size(); ++y) {
		for (int x = 0; x < ((2 * num_of_cuts) - 1); ++x) {
			arr_neg.push_back(cut_us_P[y][x]);
		}
	}

	//std::cout << arr[123];
	//std::cout << cut_us_P[2][3];

	std::fstream file_;
	file_.open("data_neg.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < arr_neg.size(); i++)
	{
		file_.write((char*)(&arr_neg[i]), sizeof(double));
	}

	*/

	return 0;

}