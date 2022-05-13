
#include <iostream>
#include <fstream> 
#include <vector> 
#include <array>
#include <algorithm>


constexpr auto DIMZ = 2048;

// This program writes the averaged cut-velocity signals


using Matrix = std::vector < std::vector<double> >; // This is just an alias for 2D vector

using Matrix_int = std::vector < std::vector<int> >;

using Matrix2_int = std::array<std::vector<int>, 2>;

using Pair = std::pair<std::vector<double>, Matrix2_int>;

using Matrix_dimz_double = std::array<std::vector<double>, DIMZ>;


using Matrix2 = std::array<std::vector<double>, 2>; // This is just an alias for 2D vector
using Matrix3 = std::array<std::vector<double>, 3>;

 // This is just an alias for 2D vector

using Pair_T = std::pair<Matrix2, std::vector<int> >;

//This is the main function for calculating the averaged velocity signals across the zero-crossings
Pair cut_velocity(std::vector <double> large, std::vector <double> small, std::vector <int> positions, int n, std::string type, std::vector <int> limit)

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
		//std::cout << "weights at " << w.at(i - 1) << " = " << std::endl;
	}

	//the arrays below take in the distances between a particular zero crossing and preceding and successive zero crossings
	std::vector<int> diff_pre;
	std::vector<int> diff_suc;

	//to solve for positive to negative
	if (type == "pos_to_neg")

	{
		for (int i = 1; i < (size - 1); ++i)
		{
			
			//int p_minus = positions.at(i - 1);
			//int p_plus = positions.at(i + 1);

			//int front_diff = p_plus - p;
			//int back_diff = p - p_minus;
			
			int p = positions.at(i);
			int n = 0;
			while ( (large.at(p) > 0.0) && (p >= limit[0]) )  
			{
				n++;
				p--;

			}
			
			diff_pre.push_back(n);

			int p_ = positions.at(i) + 1;
			int m = 0;
			while ( (large.at(p_) < 0.0)  && (p_ <= limit[1]) )
			{
				m++;
				p_++;

			}
			
			diff_suc.push_back(m);

		}

	}

	//to solve for negative to positive
	else
	{

		for (int i = 1; i < (size - 1); ++i)
		{

			int p = positions.at(i);

			int n = 0;
			while ((large.at(p) < 0.0) && (p >= limit[2]) )
			{
				n++;

				p--;

			}
			//std::cout << n << std::endl;
			
			diff_pre.push_back(n);

			int p_ = positions.at(i) + 1 ;
			int m = 0;
			while ( (large.at(p_) > 0.0) && (p_ <= limit[3]) )
			{
				m++;
				p_++;

			}
			//std::cout << n << std::endl;
			
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
		//std::cout << num_pre.at(i) << std::endl;
	}

	for (int i = 0; i < num_suc.size(); i++)
	{
		//std::cout << num_suc.at(i) << std::endl;
	}

	Matrix2_int diff_count;

	//diff_count.resize((std::vector<double>(2), num_pre.size()));


	for (int i = 0; i < (num_pre.size()); ++i)
	{

		diff_count[0].push_back(num_pre.at(i));

	}

	for (int i = 0; i < (num_suc.size()); ++i)
	{

		diff_count[1].push_back(num_suc.at(i));

	}


	// the matrix takes in the cut small scale velocity for each zero crossing position
	std::vector<double> cut_vel;

	for (int i = 1; i < (size - 1); i++) // for each of the zero-crossing positions
	{

		int  temp = positions.at(i); 


		for (int j = 0; j < (num_pre.at(i - 1) + num_suc.at(i - 1)); j++) 

		{
			int loc = (temp - (num_pre.at(i - 1)) + j);


			double vel = ((1 - w.at(i - 1)) * small.at(loc)) + (w.at(i - 1) * small.at(loc + 1));


			cut_vel.push_back(vel);

			//std::cout<<vel << std::endl;
		}

		
	}
	
	
 return { cut_vel, diff_count };

}


Pair_T average_cut_velocity(Pair velocity, std::vector <int> positions)
{


	const int size = positions.size();

	std::vector<int> length_vel;

	for (int i = 1; i < (size - 1); i++)
	{
       //length_vel takes in the size of the cut velocities for each zero crossing + 2. 
		//we add 2 to the size so as to allow for padding with -1e-6 later on, since the cut velocities are initially of unequal sizes
		length_vel.push_back(velocity.second[0].at(i - 1) + velocity.second[1].at(i - 1) + 2);

	}


	//int max_2 = *std::max_element(length_vel.begin(), length_vel.end());

	//std::cout << max_2;
	//had issues using the vector library, had to use a dynamic array

	//This is an array which takes all the cut velocities for all zero-crossings with the -1e-6 padding

	double* cut_velo = new double[10000000];


	std::vector <double> averaged;
	std::vector <double> vel_sum;


	std::vector<int> length_1;
	std::vector<int> length_2;

	//I try to make the arrays same size by padding with -1e-6 where necessary so its easy for averaging.

	//Note: 1) velocity.second[0] = preceding size
	 //     2) velocity.second[0] = succesive size

	for (int i = 1; i < (size - 1); i++)
	{
		
		length_1.push_back(velocity.second[0].at(i - 1));

	}

	int max = *std::max_element(length_1.begin(), length_1.end());

	for (int i = 1; i < (size - 1); i++)

	{
		length_2.push_back(velocity.second[1].at(i - 1));

	}

	int max_3 = *std::max_element(length_2.begin(), length_2.end());

	int array_size = max + max_3 + 2; //



	int total_size = 0;
	for (int i = 1; i < (size - 1); i++)

	{


		for (int j = 0; j < (max - velocity.second[0].at(i - 1)) + 1; j++)

		{

			cut_velo[(array_size * (i - 1)) + j] = -1e-6;

		}
		//std::cout << cut_velo[i - 1].size() << std::endl;

		int local_size = velocity.second[0].at(i - 1) + velocity.second[1].at(i - 1);
		total_size += local_size;

		for (int j = 0; j < (local_size); j++)
		{


			cut_velo[(array_size * (i - 1)) + (max - velocity.second[0].at(i - 1)) + 1 + j] = velocity.first.at(total_size - local_size + j);


			// std::cout << cut_velo[(max_2 * (i - 1)) + (max - velocity.second[0].at(i - 1)) + 1 + j] << std::endl;
		}


		//std::cout << cut_velo[i - 1].size() << std::endl;


	}



	for (int i = 1; i < (size - 1); i++)

	{
		int local_size = velocity.second[0].at(i - 1) + velocity.second[1].at(i - 1);

		for (int j = 0; j < (max_3 - velocity.second[1].at(i - 1)) + 1; j++)
		{
			//cut_velo[(((size - 2) * i) + (max - velocity.second[0].at(i - 1)) + local_size + j)] = (0.0);
			cut_velo[(array_size * (i - 1)) + (max - velocity.second[0].at(i - 1)) + 1 + local_size + j] = -1e-6;

		}

	}



	std::vector <int> divisor;

	for (int i = 1; i < array_size - 1; i++)

	{
		int num = 0;

		for (int j = 1; j < size - 1; j++)
		{
			double temp = cut_velo[((array_size * (j - 1)) + i)];

			if (temp <= -1e-7)
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

	/*
	std::fstream time_file_N;
	time_file_N.open(filename, std::ios::out | std::ios::binary);
	for (int i = 0; i < divisor.size(); i++)
	{
		time_file_N.write((char*)(&divisor[i]), sizeof(int));
	}

	std::cout << "number " << divisor.size();

	std::cout << std::endl;

	*/

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




	std::vector<double> placement;

	placement.push_back((max - 1));

	placement.push_back(array_size - max - 1);


	//std::cout << averaged.size() << std::endl;



	//std::cout << " ave" << averaged.size();
	return { { averaged, placement }, divisor };

}




std::vector<double> average_ave(std::vector<double> velocity, std::vector<int> pred, std::vector<int> succ, double dt, std::string filename, int extent)

{


	//had issues using the vector library, had to use a dynamic array
	double* cut_velo = new double[10000000];

	std::vector <double> averaged;
	std::vector <double> vel_sum;


	std::vector<int> length_1;
	std::vector<int> length_2;


	for (int i = 0; i < extent; i++)

	{

		length_1.push_back(pred.at(i));

	}

	int max = *std::max_element(length_1.begin(), length_1.end());

	for (int i = 0; i < extent; i++)

	{

		length_2.push_back(succ.at(i));

	}

	int max_3 = *std::max_element(length_2.begin(), length_2.end());

	int array_size = max + max_3 + 2;


	int total_size = 0;
	for (int i = 0; i < extent; i++)

	{

		//int fill = local_size - min;


		for (int j = 0; j < max - pred.at(i) + 1; j++)
		{


			cut_velo[(array_size * i) + j] = -1e-6;

		}


		int local_size = pred.at(i) + succ.at(i);
		total_size += local_size;


		for (int j = 0; j < (local_size); j++)
		{

			cut_velo[(array_size * (i) ) + (max - pred.at(i)) + 1 + j] = velocity.at(total_size - local_size + j);

		}


	}


	for (int i = 0; i < extent; i++)

	{
		int local_size = pred.at(i) + succ.at(i);

		for (int j = 0; j < (max_3 - succ.at(i)) + 1; j++)
		{

			cut_velo[(array_size * (i)) + (max - pred.at(i)) + 1 + local_size + j] = -1e-6;

		}

	}

	std::vector <int> divisor;

	for (int i = 1; i < array_size - 1; i++)

	{
		int num = 0;

		for (int j = 0; j < extent; j++)
		{
			double temp = cut_velo[((array_size * (j)) + i)];

			if (temp <= -1e-7)
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



	for (int i = 1; i < array_size - 1; i++)
	{
		double vel = 0.0;

		for (int j = 0; j < extent; j++)
		{

			vel = vel + cut_velo[((array_size * (j)) + i)];


		}


		averaged.push_back(vel / (divisor.at(i - 1) * 1.0));

	}


	//std::cout << "sg" << max;
	//std::cout << "gdba " << max_3;

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


	std::fstream time_file_P;
	time_file_P.open(filename, std::ios::out | std::ios::binary);
	for (int i = 0; i < time.size(); i++)
	{
		time_file_P.write((char*)(&time[i]), sizeof(double));
	}


	return averaged;

}



int main()

{

	//Read in the positions for negative to positive zero crossings

	std::ifstream infile_a("neg_pos.out", std::ifstream::binary);

	infile_a.seekg(0, infile_a.end);
	int N = infile_a.tellg();
	infile_a.seekg(0, infile_a.beg);

	std::vector<int> array_neg(N / sizeof(int));
	infile_a.read(reinterpret_cast<char*>(array_neg.data()), static_cast<std::streamsize>(array_neg.size()) * sizeof(int));

	for (int i = 0; i < array_neg.size(); ++i)
	{
		//std::cout << array_neg.at(i) << std::endl;
	}


	std::ifstream infile_ap("pos_neg.out", std::ifstream::binary);

	infile_ap.seekg(0, infile_ap.end);
	int P = infile_ap.tellg();
	infile_ap.seekg(0, infile_ap.beg);

	std::vector<int> array_pos(P / sizeof(int));
	infile_ap.read(reinterpret_cast<char*>(array_pos.data()), static_cast<std::streamsize>(array_pos.size()) * sizeof(int));


	std::ifstream infile_s("size_pos_neg.out", std::ifstream::binary);

	infile_s.seekg(0, infile_s.end);
	int P_ = infile_s.tellg();
	infile_s.seekg(0, infile_s.beg);

	std::vector<int> size_array_pos(P_ / sizeof(int));
	infile_s.read(reinterpret_cast<char*>(size_array_pos.data()), static_cast<std::streamsize>(size_array_pos.size()) * sizeof(int));


	std::ifstream infile_k("size_neg_pos.out", std::ifstream::binary);

	infile_k.seekg(0, infile_k.end);
	int N_ = infile_k.tellg();
	infile_k.seekg(0, infile_k.beg);

	std::vector<int> size_array_neg(N_ / sizeof(int));
	infile_k.read(reinterpret_cast<char*>(size_array_neg.data()), static_cast<std::streamsize>(size_array_neg.size()) * sizeof(int));



	std::vector<int> lim;

	lim.push_back(*std::min_element(array_neg.begin(), array_neg.end()));

	lim.push_back(*std::max_element(array_neg.begin(), array_neg.end()));

	lim.push_back(*std::min_element(array_pos.begin(), array_pos.end()));

	lim.push_back(*std::max_element(array_pos.begin(), array_pos.end()));


	std::ifstream infile_st("us.out", std::ifstream::binary);

	infile_st.seekg(0, infile_st.end);
	int M = infile_st.tellg();
	infile_st.seekg(0, infile_st.beg);

	std::vector<double> Us(M / sizeof(double));
	infile_st.read(reinterpret_cast<char*>(Us.data()), static_cast<std::streamsize>(Us.size()) * sizeof(double));


	std::ifstream infile_ar("ul.out", std::ifstream::binary);

	infile_ar.seekg(0, infile_ar.end);
	int Pj = infile_ar.tellg();
	infile_ar.seekg(0, infile_ar.beg);

	std::vector<double> Ul(Pj / sizeof(double));
	infile_ar.read(reinterpret_cast<char*>(Ul.data()), static_cast<std::streamsize>(Ul.size()) * sizeof(double));




	//std::cout << Ul.size();

	const int num = Ul.size() / DIMZ;


	//using Matrix = std::array<std::vector<double>, 2>;

	//std::vector < std::vector<double> > ul;

	
   //ul.resize(DIMZ, std::vector<double>(num));

	//ul.resize(DIMZ, std::vector<double>(num));

	//ul.resize((std::vector<double>(DIMZ), num));

	Matrix_dimz_double ul;
	//ul = Ul;
	//std::cout << Ul.size();

	for (int j = 0; j < DIMZ; ++j)
	{
		for (int i = 0; i < num; ++i)
		{

			ul[j].push_back(Ul.at((i * DIMZ) + j));

		}

	}

		//std::vector < std::vector<double> > us;

		Matrix_dimz_double us;

	//us.resize((std::vector<double>(DIMZ), num));
		

		for (int j = 0; j < DIMZ; ++j)
		{
			for (int i = 0; i < num; ++i)
			{

				us[j].push_back(Us.at((i * DIMZ) + j));

			}

		}



	int starting_size_Pos = 0;

	std::vector <double> vel_matrix_pos;

	std::vector <double> ensemble_pos;

	std::vector < int> pre_pos;
	std::vector < int> suc_pos;

	//std::cout << size_array_pos.size();

	int counter_pos = 0;

	for (int z = 0; z < DIMZ; z++)
	{
		if (size_array_pos.at(z) < 4)
		{
		
			starting_size_Pos = starting_size_Pos + size_array_pos[z];
		}

		else
		{

			Pair cut_us_P; // cut small velocity and time for positive to negative

			Pair_T cut_us_Pos;

			int num_of_cuts;
			// the num_of_cuts variable is read in from a text file
			std::ifstream infile_txt("input.txt");

			infile_txt >> num_of_cuts;

			std::vector <int> local_Pos;



			for (int i = 0; i < size_array_pos.at(z); i++)
			{

				local_Pos.push_back(array_pos.at(starting_size_Pos + i));

			}


			starting_size_Pos = starting_size_Pos + size_array_pos[z];


			cut_us_P = cut_velocity(ul[z], us[z], local_Pos, num_of_cuts, "pos_to_neg", lim);

			

			//std::cout << ul[z].size() << std::endl;

			cut_us_Pos = average_cut_velocity(cut_us_P, local_Pos);
			local_Pos.clear();

			int local_pre = round(cut_us_Pos.first[1].at(0));
			int local_suc = round(cut_us_Pos.first[1].at(1));

			pre_pos.push_back(local_pre);
			suc_pos.push_back(local_suc);


			for (int i = 0; i < cut_us_Pos.first[0].size(); i++)
			{

				vel_matrix_pos.push_back(cut_us_Pos.first[0].at(i));

				ensemble_pos.push_back(cut_us_Pos.second.at(i));

			}

			counter_pos++;
		}

		//starting_size_Pos = starting_size_Pos + size_array_pos[z];
		

	}

	
	int starting_size_Neg = 0;


	std::vector <double> vel_matrix_neg;

	std::vector <double> ensemble_neg;

	std::vector < int> pre_neg;
	std::vector < int> suc_neg;

	int counter_neg = 0;

	for (int z = 0; z < DIMZ; z++)
	{
		if (size_array_neg.at(z) < 4)
		{
			starting_size_Neg = starting_size_Neg + size_array_neg[z];
		}

		else
		{

			Pair cut_us_N; // cut small velocity and time for positive to negative

			Pair_T cut_us_Neg;

			int num_of_cuts;
			// the num_of_cuts variable is read in from a text file
			std::ifstream infile_txt("input.txt");

			infile_txt >> num_of_cuts;


			std::vector <int> local_Neg;

			for (int i = 0; i < size_array_neg.at(z); i++)
			{

				local_Neg.push_back(array_neg.at(starting_size_Neg + i));

			}

			starting_size_Neg = starting_size_Neg + size_array_neg[z];

			cut_us_N = cut_velocity(ul[z], us[z], local_Neg, num_of_cuts, "neg_to_pos", lim);

			cut_us_Neg = average_cut_velocity(cut_us_N, local_Neg);


			int local_pre = round(cut_us_Neg.first[1].at(0));
			int local_suc = round(cut_us_Neg.first[1].at(1));

			pre_neg.push_back(local_pre);
			suc_neg.push_back(local_suc);

		//std::cout << "hzu" << cut_us_Neg.first[0].size();

		//std::cout << "tu" << local_pre + local_suc;

			for (int i = 0; i < cut_us_Neg.first[0].size(); i++)
			{

				vel_matrix_neg.push_back(cut_us_Neg.first[0].at(i));

				ensemble_neg.push_back(cut_us_Neg.second.at(i));
				


			}
			//std::cout << "yd" << vel_matrix_neg.size() << std::endl;
			//std::cout << "grg" << ensemble_neg.size() << std::endl;

			
			counter_neg++;
		}
	}


	


	double dt = 0.000127175055596395;



	std::vector<double> averaged_pos_vel;
	std::vector<double> averaged_neg_vel;


	std::vector<double> averaged_pos_ens;
	std::vector<double> averaged_neg_ens;


	averaged_pos_vel = average_ave(vel_matrix_pos, pre_pos, suc_pos, dt, "time_pos_neg.out", counter_pos);

	averaged_neg_vel = average_ave(vel_matrix_neg, pre_neg, suc_neg, dt, "time_neg_pos.out", counter_neg);

	averaged_pos_ens = average_ave(ensemble_pos, pre_pos, suc_pos, dt, "time_pos_neg.out", counter_pos);

	averaged_neg_ens = average_ave(ensemble_neg, pre_neg, suc_neg, dt, "time_neg_pos.out", counter_neg);



	std::fstream file_P;
	file_P.open("Velocity_pos_neg.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < averaged_pos_vel.size(); i++)
	{
		file_P.write((char*)(&averaged_pos_vel[i]), sizeof(double));
	}



	std::fstream file_N;
	file_N.open("Velocity_neg_pos.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < averaged_neg_vel.size(); i++)
	{
		file_N.write((char*)(&averaged_neg_vel[i]), sizeof(double));
	}



	std::fstream time_file_P;
	time_file_P.open("ensemble_pos.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < averaged_pos_ens.size(); i++)
	{
		time_file_P.write((char*)(&averaged_pos_ens[i]), sizeof(double));
	}


	std::fstream time_file_N;
	time_file_N.open("ensemble_neg.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < averaged_neg_ens.size(); i++)
	{
		time_file_N.write((char*)(&averaged_neg_ens[i]), sizeof(double));
	}




	

}
