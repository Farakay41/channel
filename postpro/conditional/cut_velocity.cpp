#include <iostream>
#include <fstream> 
#include <vector> 
#include <array>
#include <algorithm>

// This program writes the averaged cut-velocity signals


using Matrix = std::vector < std::vector<double> >; // This is just an alias for 2D vector

using Matrix_int = std::vector < std::vector<int> >;

using Pair = std::pair<std::vector<double>, Matrix_int>;

using Matrix2 = std::array<std::vector<double>, 2>; // This is just an alias for 2D vector
//This is the main function for calculating the averaged velocity signals across the zero-crossings
Pair cut_velocity(std::vector <double> large, std::vector <double> small, std::vector <int> positions, int n)

{



	//std::cout << n;  //number of velocity signals each right and left) of the zero crossing

   // std::vector<double> cut_vel;// the array takes in the cut small scale velocity for each zero crossing positions

	int size = positions.size(); // This calculates the size of the positions (of zero-crossing) array, Also needed to calculate the average velocity signal later

	//vector w contains the calculated weights

	std::vector<double> w;

	for (int i = 1; i < (size - 1); ++i)

	{
		int pos = positions.at(i);

		w.push_back(large.at(pos) / (large.at(pos) - large.at(pos + 1)));
		//std::cout << w.at(i-1) << std::endl;
	}

	std::vector<int> diff_pre;
	std::vector<int> diff_suc;

	for (int i = 1; i < (size - 1); ++i)
	{
		diff_pre.push_back(positions.at(i) - positions.at(i - 1));

		diff_suc.push_back(positions.at(i + 1) - positions.at(i));
	}


	std::vector<int> num_pre;
	std::vector<int> num_suc;

	for (int i = 0; i < (diff_pre.size()); ++i)
	{

		num_pre.push_back(floor(diff_pre.at(i) / 2));
		num_suc.push_back(floor(diff_suc.at(i) / 2));

	}


	for (int i = 0; i < num_pre.size(); i++)
	{
		std::cout << num_pre.at(i) << std::endl;
	}

	for (int i = 0; i < num_suc.size(); i++)
	{
		std::cout << num_suc.at(i) << std::endl;
	}

	Matrix_int diff_count;

	diff_count.resize((std::vector<double>(2), num_pre.size()));

	for (int i = 0; i < (num_pre.size()); ++i)
	{

		diff_count[0].push_back(num_pre.at(i));

	}

	for (int i = 0; i < (num_suc.size()); ++i)
	{

		diff_count[1].push_back(num_suc.at(i));

	}

	std::vector<double> cut_vel;



	// the matrix takes in the cut small scale velocity for each zero crossing position
	for (int i = 1; i < (size - 1); i++) // for each of the zero-crossing positions
	{

		int  temp = positions.at(i); // allocates a variable for each zero crossing positions


		for (int j = 0; j < (num_pre.at(i - 1) + num_suc.at(i - 1)); j++) // n is the number of small scale velocity signals and time cut for each zero-crossing position

		{
			int loc = (temp - (num_pre.at(i - 1)) + j);


			double vel = ((1 - w.at(i - 1)) * small.at(loc)) + (w.at(i - 1) * small.at(loc + 1));


			cut_vel.push_back(vel);


		}



	}
	//std::cout << cut_vel.at(304066) << std::endl;
	//std::cout << cut_vel[0][0]<< std::endl;
	//std::cout << cut_vel[1][0] << std::endl;
	//std::cout << cut_vel[2][0] << std::endl;

	return { cut_vel, diff_count };
}


std::vector <double> average_cut_velocity(Pair velocity, std::vector <int> positions)
{

	//std::cout << velocity.first.at(304067) << std::endl;

	const int size = positions.size();

	std::vector<int> length_vel;

	for (int i = 1; i < (size - 1); i++)

	{
		length_vel.push_back(velocity.second[0].at(i - 1) + velocity.second[1].at(i - 1) + 2);


	}

	int max_2 = *std::max_element(length_vel.begin(), length_vel.end());

	//had issues using te vector library, had to use a dynamic array
	double* cut_velo = new double[1000000];

	std::vector <double> averaged;
	std::vector <double> vel_sum;


	std::vector<int> length;

	/*I try to make the arrays same size by padding with zero where necessary.So i use the zero - crossing position with the highest preceding velocities when cut to calculate
	how much zeros other zero crossing position needs befor the zero crossing*/

	for (int i = 1; i < (size - 1); i++)

	{
		//std::cout << velocity.second[0].at(i - 1);
		//(num_pre.at(i - 1) + num_suc.at(i - 1))
		length.push_back(velocity.second[0].at(i - 1));
		//std::cout << length.at(i-1);
		//std::cout << length.at(i - 1) << std::endl;
	}

	int max = *std::max_element(length.begin(), length.end());

	//std::cout << max;
	//std::cout << velocity[0].size();

	//std::cout << velocity.size();



	//difference


	//cut_velo.resize((std::vector<double>(size - 2), 1000));
	int total_size = 0;
	for (int i = 1; i < (size - 1); i++)

	{

		//int fill = local_size - min;


		for (int j = 0; j < (max - velocity.second[0].at(i - 1)) + 1; j++)
		{


			cut_velo[(max_2 * (i - 1)) + j] = 0.0;

		}
		//std::cout << cut_velo[i - 1].size() << std::endl;

		int local_size = velocity.second[0].at(i - 1) + velocity.second[1].at(i - 1);
		total_size += local_size;

		for (int j = 0; j < (local_size); j++)
		{

			//cut_velo[i - 1].push_back(velocity[0][i - 1].at(j));

			//cut_velo.at(((size - 2) * i) + (max - velocity.second[0].at(i - 1)) + j);
			//double a =

			cut_velo[(max_2 * (i - 1)) + (max - velocity.second[0].at(i - 1)) + 1 + j] = velocity.first.at(total_size - local_size + j);
			// std::cout << cut_velo[(max_2 * (i - 1)) + (max - velocity.second[0].at(i - 1)) + 1 + j] << std::endl;
		}


		//std::cout << cut_velo[i - 1].size() << std::endl;


	}



	//Fill the other elements with to zero to have the same size as the largest "cut-velocity array"

	for (int i = 1; i < (size - 1); i++)

	{
		int local_size = velocity.second[0].at(i - 1) + velocity.second[1].at(i - 1);

		for (int j = 0; j < max_2 - length_vel.at(i - 1) + 1; j++)
		{

			//cut_velo[(((size - 2) * i) + (max - velocity.second[0].at(i - 1)) + local_size + j)] = (0.0);

			cut_velo[(max_2 * (i - 1)) + (max - velocity.second[0].at(i - 1)) + 1 + local_size + j] = 0.0;

		}

		//std::cout << cut_velo[i-1].size();

	}



	std::vector <int> divisor;

	for (int i = 1; i < max_2; i++)

	{
		int num = 0;

		for (int j = 1; j < size - 1; j++)
		{
			double temp = cut_velo[((max_2 * (j - 1)) + i)];

			if (temp == 0.0)
			{
				num += 0;
			}

			else
			{
				num += 1;
			}

		}

		divisor.push_back(num);

	}
	std::cout << std::endl;
	


	for (int i = 1; i < max_2; i++)
	{
		double vel = 0.0;

		for (int j = 1; j < size - 1; j++)
		{

			vel = vel + cut_velo[((max_2 * (j - 1)) + i)];


		}


		averaged.push_back(vel / (divisor.at(i - 1) * 1.0));
	

	}

	delete []cut_velo;

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

	std::ifstream infile_s("us_file.out", std::ifstream::binary);

	infile_s.seekg(0, infile_s.end);
	int M = infile_s.tellg();
	infile_s.seekg(0, infile_s.beg);

	std::vector<double> us(M / sizeof(double));
	infile_s.read(reinterpret_cast<char*>(us.data()), static_cast<std::streamsize>(us.size()) * sizeof(double));


	//This reads in the large scale velocity signals from a binary file and write it into an array
	std::ifstream infile_l("ul_file.out", std::ifstream::binary);

	infile_l.seekg(0, infile_a.end);
	int V = infile_l.tellg();
	infile_l.seekg(0, infile_a.beg);

	std::vector<double> ul(V / sizeof(double));
	infile_l.read(reinterpret_cast<char*>(ul.data()), static_cast<std::streamsize>(ul.size()) * sizeof(double));

	Pair cut_us_P; // cut small velocity and time for positive to negative
	Pair cut_us_N; // cut small velocity and time for negative to positive

	int num_of_cuts; // number of velocity signals in the right(or left) of the zero crossing


	// the num_of_cuts variable is read in from a text file
	std::ifstream infile_txt("input.txt");

	infile_txt >> num_of_cuts;

	//call the cut_velocity function to calculate the averaged velocity signals (positive to negative)
	cut_us_P = cut_velocity(ul, us, array_pos, num_of_cuts);

	//call the cut_velocity function to calculate the averaged velocity signals (negative to positive)
	cut_us_N = cut_velocity(ul, us, array_neg, num_of_cuts);

	//std::cout << cut_us_N[0][0][0];

	//std::cout << cut_us_P[0].size();

	std::vector<double> cut_us_Pos;

	std::vector<double> cut_us_Neg;

	cut_us_Pos = average_cut_velocity(cut_us_P, array_pos);
	cut_us_Neg = average_cut_velocity(cut_us_N, array_neg);

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


}

/*

*/