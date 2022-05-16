#include <iostream>
#include <fstream> 
#include <vector> 
#include <array>
#include <algorithm>

constexpr auto DIMZ = 2048;

using Matrix = std::vector < std::vector<double> >; // This is just an alias for 2D vector

using Matrix2 = std::array <std::vector<double>, 2>;


using Matrix_dimz_double = std::array<std::vector<double>, DIMZ>;

//This is the main function for calculating the averaged velocity signals across the zero-crossings




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




std::vector<double> cut_velocity(std::vector <double> large, std::vector <double> small, std::vector <int> positions, int n)

{
	//std::cout << positions.at(0);


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
	std::vector <double> time;
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



Matrix2 average_ave(std::vector<double> velocity, int extent, int n, double dt)
{
	std::vector <double> averaged;
	std::vector <double> vel_sum;
	
	//std::cout << velocity.size();

	for (int i = 0; i < ((2 * n) - 1); i++)
	{
		double vel = 0.0;

		for (int j = 0; j < extent; j++)
		{

			vel = vel + velocity.at((((2 * n) - 1) * j) + i);


		}
		averaged.push_back(vel / extent);
		// std::cout << averaged.size();

	}


	std::vector <double> time;

	for (int j = 0; j < n - 1; j++)
	{

		time.push_back(dt * (j - (n - 1)));
		std::cout << time.at(j) << std::endl;

	}

	for (int j = n - 1; j < (2 * n) - 1; j++)
	{

		time.push_back(dt * (j - (n - 1)));
		std::cout << time.at(j) << std::endl;

	}

	std::cout << time.size();



	return { averaged, time };

}




int main()
{
	

	std::vector<int> array_neg;
	array_neg = read_int("neg_pos.out");


	std::vector<int> array_pos;
	array_pos = read_int("pos_neg.out");


	std::vector<int> size_array_neg;
	size_array_neg = read_int("size_neg_pos.out");


	std::vector<int> size_array_pos;
	size_array_pos = read_int("size_pos_neg.out");

	std::vector<double> Us;
	Us = read_double("us.out");


	std::vector<double> Ul;
	Ul = read_double("ul.out");

	

    //std::cout << Ul.size();

	const int num = Ul.size() / DIMZ;


	//using Matrix = std::array<std::vector<double>, 2>;

	//std::vector < std::vector<double> > ul;

	//ul.resize( DIMZ, std::vector<double> (num));
	//ul.resize((std::vector<double>(DIMZ), num));

	Matrix_dimz_double ul;
		
	for (int j = 0; j < DIMZ; ++j)
	{
		for (int i = 0; i < num; ++i)
		{

			ul[j].push_back(Ul.at((i * DIMZ) + j));

		}

	}


	//std::vector < std::vector<double> > us;

	//us.resize(DIMZ, std::vector<double>(num));
	//us.resize((std::vector<double>(DIMZ), num));

	Matrix_dimz_double us;

	for (int j = 0; j < DIMZ; ++j)
	{
		for (int i = 0; i < num; ++i)
		{

			us[j].push_back(Us.at((i * DIMZ) + j));

		}

	}

	std::cout << "sgs" << us[45].at(46);


	//Matrix cut_us_P; // cut small velocity and time for positive to negative
	//Matrix cut_us_N; // cut small velocity and time for negative to positive

	int num_of_cuts; // number of velocity signals in the right(or left) of the zero crossing


	// the num_of_cuts variable is read in from a text file
	std::ifstream infile_txt("input.txt");

	infile_txt >> num_of_cuts;

	double dt_;
	std::ifstream infilt_txt("time_step.txt");
	infilt_txt >> dt_;



	int starting_size_Pos = 0;

	std::vector <double> vel_matrix_pos;

	int counter_pos = 0;

	for (int z = 0; z < DIMZ; z++)
	{
		if (size_array_pos.at(z) < 4)
		{

			starting_size_Pos = starting_size_Pos + size_array_pos[z];
		}

		else
		{

			std::vector<double> cut_us_P; // cut small velocity and time for positive to negative

			std::vector<double> cut_us_Pos;


			std::vector <int> local_Pos;

			for (int i = 0; i < size_array_pos.at(z); i++)
			{

				local_Pos.push_back(array_pos.at(starting_size_Pos + i));

			}


			starting_size_Pos = starting_size_Pos + size_array_pos[z];


			cut_us_P = cut_velocity(ul[z], us[z], local_Pos, num_of_cuts);



			//std::cout << ul[z].size() << std::endl;

			cut_us_Pos = average_cut_velocity(cut_us_P, local_Pos, num_of_cuts);
			local_Pos.clear();



			for (int i = 0; i < cut_us_Pos.size(); i++)
			{

				vel_matrix_pos.push_back(cut_us_Pos.at(i));


			}

			counter_pos++;
		}

	}



		int starting_size_Neg = 0;

		std::vector <double> vel_matrix_neg;

		int counter_neg = 0;

		for (int z = 0; z < DIMZ; z++)
		{
			if (size_array_neg.at(z) < 4)
			{

				starting_size_Neg = starting_size_Neg + size_array_neg[z];
			}

			else
			{

				std::vector<double> cut_us_N; // cut small velocity and time for positive to negative

				std::vector<double> cut_us_Neg;


				std::vector <int> local_Neg;

				for (int i = 0; i < size_array_neg.at(z); i++)
				{

					local_Neg.push_back(array_neg.at(starting_size_Neg + i));

				}


				starting_size_Neg = starting_size_Neg + size_array_neg[z];


				cut_us_N = cut_velocity(ul[z], us[z], local_Neg, num_of_cuts);



				//std::cout << ul[z].size() << std::endl;

				cut_us_Neg = average_cut_velocity(cut_us_N, local_Neg, num_of_cuts);
				local_Neg.clear();



				for (int i = 0; i < cut_us_Neg.size(); i++)
				{

					vel_matrix_neg.push_back(cut_us_Neg.at(i));


				}

				counter_neg++;
			}

		}


			Matrix2 averaged_pos_vel;
			Matrix2 averaged_neg_vel;

			averaged_pos_vel = average_ave(vel_matrix_pos, counter_pos, num_of_cuts, dt_);
			averaged_neg_vel = average_ave(vel_matrix_neg, counter_neg, num_of_cuts, dt_);





	std::fstream file_P;
	file_P.open("Velocity_pos_neg.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < averaged_pos_vel[0].size(); i++)
	{
		file_P.write((char*)(&averaged_pos_vel[0][i]), sizeof(double));
	}


	std::fstream file_N;
	file_N.open("Velocity_neg_pos.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < averaged_neg_vel[0].size(); i++)
	{
		file_N.write((char*)(&averaged_neg_vel[0][i]), sizeof(double));
	}


	std::fstream time_file_P;
	time_file_P.open("time_pos_neg.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < averaged_pos_vel[1].size(); i++)
	{
		time_file_P.write((char*)(&averaged_pos_vel[1][i]), sizeof(double));
	}

	std::fstream time_file_N;
	time_file_N.open("time_neg_pos.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < averaged_neg_vel[1].size(); i++)
	{
		time_file_N.write((char*)(&averaged_neg_vel[1][i]), sizeof(double));
	}




	/*
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
	*/


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