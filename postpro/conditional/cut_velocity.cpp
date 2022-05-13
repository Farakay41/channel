
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
using Matrix3 = std::array<std::vector<double>, 3>;
//This is the main function for calculating the averaged velocity signals across the zero-crossings
Pair cut_velocity(std::vector <double> large, std::vector <double> small, std::vector <int> positions, int n, std::string type)

{



	//std::cout << n;  //number of velocity signals each right and left) of the zero crossing

   // std::vector<double> cut_vel;// the array takes in the cut small scale velocity for each zero crossing positions

	int size = positions.size(); // This calculates the size of the positions (of zero-crossing) array, Also needed to calculate the average velocity signal later

	//vector w contains the calculated weights

	std::vector<double> w;

	for (int i = 1; i < (size - 1); ++i)

	{
		int pos = positions.at(i);
		std::cout << pos << std::endl;
		w.push_back(large.at(pos) / (large.at(pos) - large.at(pos + 1)));
		//std::cout << w.at(i-1) << std::endl;
	}

	std::vector<int> diff_pre;
	std::vector<int> diff_suc;


	if (type == "pos_to_neg")

	{
		for (int i = 1; i < (size - 1); ++i)
		{
			int p = positions.at(i);
			int p_minus = positions.at(i - 1);
			int p_plus = positions.at(i + 1);
			int front_diff = p_plus - p;
			int back_diff = p - p_minus;

			int n = 0;
			while (large.at(p) > 0.0)
			{
				n++;
				p--;

			}

			diff_pre.push_back(n);

			int p_ = positions.at(i);
			int m = 0;
			while (large.at(p_) < 0.0)
			{
				m++;
				p_++;

			}

			diff_suc.push_back(m);

		}

	}

	else
	{

		for (int i = 1; i < (size - 1); ++i)
		{
			int p = positions.at(i);
			int p_minus = positions.at(i - 1);
			int p_plus = positions.at(i + 1);
			int front_diff = p_plus - p;
			int back_diff = p - p_minus;


			int n = 0;
			while (large.at(p) < 0.0)
			{
				n++;
			    
				p--;

			}
			std::cout << n << std::endl;

			diff_pre.push_back(n);

			int p_ = positions.at(i);
			int m = 0;
			while (large.at(p_) > 0.0)
			{
				m++;
				p_++;

			}
			std::cout << n << std::endl;
			diff_suc.push_back(m);

		}


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
	std::cout << "size" << cut_vel.size();
	//std::cout << cut_vel.at(304066) << std::endl;
	//std::cout << cut_vel[0][0]<< std::endl;
	//std::cout << cut_vel[1][0] << std::endl;
	//std::cout << cut_vel[2][0] << std::endl;

	return { cut_vel, diff_count };
}


Matrix2 average_cut_velocity(Pair velocity, std::vector <int> positions, double dt, std::string filename)
{

	//std::cout << velocity.first.at(304067) << std::endl;

	const int size = positions.size();

	std::vector<int> length_vel;

	for (int i = 1; i < (size - 1); i++)

	{
		length_vel.push_back(velocity.second[0].at(i - 1) + velocity.second[1].at(i - 1) + 2);


	}

	int max_2 = *std::max_element(length_vel.begin(), length_vel.end());

	std::cout << max_2;
	//had issues using te vector library, had to use a dynamic array
	double* cut_velo = new double[10000000];

	std::vector <double> averaged;
	std::vector <double> vel_sum;


	std::vector<int> length_1;
	std::vector<int> length_2;

	/*I try to make the arrays same size by padding with zero where necessary.So i use the zero - crossing position with the highest preceding velocities when cut to calculate
	how much zeros other zero crossing position needs befor the zero crossing*/

	for (int i = 1; i < (size - 1); i++)

	{
		//std::cout << velocity.second[0].at(i - 1);
		//(num_pre.at(i - 1) + num_suc.at(i - 1))
		length_1.push_back(velocity.second[0].at(i - 1));
		//std::cout << length.at(i-1);
		//std::cout << length.at(i - 1) << std::endl;
	}

	int max = *std::max_element(length_1.begin(), length_1.end());

	for (int i = 1; i < (size - 1); i++)

	{
		//std::cout << velocity.second[0].at(i - 1);
		//(num_pre.at(i - 1) + num_suc.at(i - 1))
		length_2.push_back(velocity.second[1].at(i - 1));
		//std::cout << length.at(i-1);
		//std::cout << length.at(i - 1) << std::endl;
	}

	int max_3 = *std::max_element(length_2.begin(), length_2.end());

	int array_size = max + max_3 + 2;
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


			cut_velo[(array_size * (i - 1)) + j] = -1e-12;

		}
		//std::cout << cut_velo[i - 1].size() << std::endl;

		int local_size = velocity.second[0].at(i - 1) + velocity.second[1].at(i - 1);
		total_size += local_size;

		for (int j = 0; j < (local_size); j++)
		{

			//cut_velo[i - 1].push_back(velocity[0][i - 1].at(j));

			//cut_velo.at(((size - 2) * i) + (max - velocity.second[0].at(i - 1)) + j);
			//double a =

			cut_velo[(array_size * (i - 1)) + (max - velocity.second[0].at(i - 1)) + 1 + j] = velocity.first.at(total_size - local_size + j);

			// std::cout << cut_velo[(max_2 * (i - 1)) + (max - velocity.second[0].at(i - 1)) + 1 + j] << std::endl;
		}


		//std::cout << cut_velo[i - 1].size() << std::endl;


	}



	//Fill the other elements with to zero to have the same size as the largest "cut-velocity array"

	for (int i = 1; i < (size - 1); i++)

	{
		int local_size = velocity.second[0].at(i - 1) + velocity.second[1].at(i - 1);

		for (int j = 0; j < (max_3 - velocity.second[1].at(i - 1)) + 1; j++)
		{
			//cut_velo[(((size - 2) * i) + (max - velocity.second[0].at(i - 1)) + local_size + j)] = (0.0);
			cut_velo[(array_size * (i - 1)) + (max - velocity.second[0].at(i - 1)) + 1 + local_size + j] = -1e-12;

		}

	}



	std::vector <int> divisor;

	for (int i = 1; i < array_size - 1; i++)

	{
		int num = 0;

		for (int j = 1; j < size - 1; j++)
		{
			double temp = cut_velo[((array_size * (j - 1)) + i)];

			if (temp == -1e-12)
			{
				num += 0;
			}

			else
			{
				num += 1;
			}

		}
		//std::cout << num << std::endl;
		divisor.push_back(num);

	}

	std::fstream time_file_N;
	time_file_N.open(filename, std::ios::out | std::ios::binary);
	for (int i = 0; i < divisor.size(); i++)
	{
		time_file_N.write((char*)(&divisor[i]), sizeof(int));
	}

	std::cout << "number " << divisor.size();

	std::cout << std::endl;



	for (int i = 1; i < array_size - 1; i++)
	{
		double vel = 0.0;

		for (int j = 1; j < size - 1; j++)
		{

			vel = vel + cut_velo[((array_size * (j - 1)) + i)];


		}


		averaged.push_back(vel / (divisor.at(i - 1) * 1.0));


	}

	delete[] cut_velo;

	std::vector<double> time;


	for (int j = 1; j < max; j++)
	{

		time.push_back(dt * (j - max));


	}

	//time.at(max - 1) = 0.0;



	for (int j = max; j < array_size - 1; j++)
	{

		time.push_back(dt * (j - max));


	}

	for (int j = 0; j < time.size(); j++)
	{

		//std::cout << time.at(j) << std::endl;


	}

	std::cout << averaged.size() << std::endl;
	std::cout << time.size() << std::endl;


	std::cout << " ave" << averaged.size();
	return { averaged, time };

}



int main()
{
	//Read in the positions for positive to negative zero crossings

	std::ifstream infile("pos_neg.out", std::ifstream::binary);

	infile.seekg(0, infile.end);
	int P = infile.tellg();
	infile.seekg(0, infile.beg);

	std::vector<int> array_pos(P / sizeof(int));
	infile.read(reinterpret_cast<char*>(array_pos.data()), static_cast<std::streamsize>(array_pos.size()) * sizeof(int));


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

	Pair cut_us_P; // cut small velocity and time for positive to negative
	Pair cut_us_N; // cut small velocity and time for negative to positive

	int num_of_cuts; // number of velocity signals in the right(or left) of the zero crossing


	// the num_of_cuts variable is read in from a text file
	std::ifstream infile_txt("input.txt");

	infile_txt >> num_of_cuts;

	//call the cut_velocity function to calculate the averaged velocity signals (positive to negative)
	cut_us_P = cut_velocity(ul, us, array_pos, num_of_cuts, "pos_to_neg");

	//call the cut_velocity function to calculate the averaged velocity signals (negative to positive)
	cut_us_N = cut_velocity(ul, us, array_neg, num_of_cuts, " neg_to_pos");

	//std::cout << cut_us_N[0][0][0];

	//std::cout << cut_us_P[0].size();
	Matrix2 cut_us_Pos;

	Matrix2 cut_us_Neg;

	double dt = 0.000127175055596395;
	cut_us_Pos = average_cut_velocity(cut_us_P, array_pos, dt, "ensemble_pos.out");
	cut_us_Neg = average_cut_velocity(cut_us_N, array_neg, dt, "ensemble_neg.out");

	//write the result to a textfile


	std::fstream file_P;
	file_P.open("Velocity_pos_neg.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < cut_us_Pos[0].size(); i++)
	{
		file_P.write((char*)(&cut_us_Pos[0][i]), sizeof(double));
	}


	std::fstream file_N;
	file_N.open("Velocity_neg_pos.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < cut_us_Neg[0].size(); i++)
	{
		file_N.write((char*)(&cut_us_Neg[0][i]), sizeof(double));
	}


	std::fstream time_file_P;
	time_file_P.open("time_pos_neg.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < cut_us_Pos[1].size(); i++)
	{
		time_file_P.write((char*)(&cut_us_Pos[1][i]), sizeof(double));
	}

	std::fstream time_file_N;
	time_file_N.open("time_neg_pos.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < cut_us_Neg[1].size(); i++)
	{
		time_file_N.write((char*)(&cut_us_Neg[1][i]), sizeof(double));
	}




	/*
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

	*/
}

/*
*/