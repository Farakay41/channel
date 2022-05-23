

#include "Zero_Cross.h"

#include <jsoncpp/json/json.h>


constexpr auto num_of_cuts = 70;



using Matrix = std::vector < std::vector<double> >; // This is just an alias for 2D vector

using Matrix2 = std::array <std::vector<double>, 2>;

using Matrix_dimz_double = std::array<std::vector<double>, DIMZ>;


//This is the main function for calculating the averaged velocity signals across the zero-crossings




void cut_velocity(std::vector <double> large, std::vector <double> small, std::vector <int> positions, std::string time_filename, int n, double* vel_buf, int* num_samples)

// TODO: pass cut_vel and number_of_samples_array (choose a name for that) as POINTERS
// these are now the return values of the function
// instead, we just pass a reference to these arrays and we change them
// MAYBE: you do not need to pass a pointer
// maybe arrays that are defined in the main can be directly accessed by the function

{
	//std::cout << positions.at(0);


	//std::cout << n;  //number of velocity signals each right and left) of the zero crossing

   // std::vector<double> cut_vel;// the array takes in the cut small scale velocity for each zero crossing positions

	int size = positions.size(); // This calculates the size of the positions (of zero-crossing) array, Also needed to calculate the average velocity signal later

	//vector w contains the calculated weights

	std::vector<double> w;

	for (int i = 0; i < size; ++i) // looping over zero crossings

	{
		int pos = positions.at(i);

		w.push_back(large.at(pos) / (large.at(pos) - large.at(pos + 1)));
		//std::cout << w.at(i) << std::endl;
	}


	// TODO: add the following command: cut_vel = 0
	// TODO: you need an array with the same size as vel, called number of samples, and you also set it to zero initially
	

	//cut_vel.resize((std::vector<double>(size), ((2 * n) - 1)));

	// the matrix takes in the cut small scale velocity for each zero crossing position


	std::vector<double> t;

	
	t = read_double(time_filename);
	
	

	std::ifstream ifs("inputed.json");

	Json::Reader reader;

	Json::Value obj;

	reader.parse(ifs, obj);


	int method;
	double dt, threshold;

	threshold = obj["threshold"].asDouble();

	method = obj["method"].asUInt();

	for (size_t i = 0; i < size; i++) // for each of the zero-crossing positions
	{

		// TODO: we need to add some conditions here: we cannot use all zero crossings
		// for first and second method: postions.at(i) must be > (or >=, check it please) window size n
		// for the second method: position(i+1) - position(i) > something  AND  position(i) - position(i-1) > something

		// SUGGESTION: do this by checking the negative of these statements
		// example for the first case:
		// if (positions.at(i) < something) {
		// 	 continue
		// }


		int  temp = positions.at(i); // allocates a variable for each zero crossing positions


		if ((temp <= n) || (large.size() - temp) <= n)

		{

			continue;

		}

		switch (method)
		{

		case 1:

			for (size_t j = 0; j < ((2 * n) - 1); j++) // n is the number of small scale velocity signals and time cut for each zero-crossing position

			{
				int loc = (temp - (n - 1) + j);



				double vel = ((1 - w.at(i)) * small.at(loc)) + (w.at(i) * small.at(loc + 1));


				//cut_vel.push_back(vel); // TODO: don't do this - you know the size of vel!

				// TODO: this is how you write to cut_vel


				vel_buf[j] += vel;

				num_samples[j] += 1;


			}

		case 2:

			for (size_t i = 0; i < size; i++) // for each of the zero-crossing positions
			{


				int  temp = positions[i]; // allocates a variable for each zero crossing positions


				if ((temp == positions[0]) || (temp == positions[size - 1]))

				{
					continue;
				}

				else

				{
					int temp_plus = positions[i + 1];
					int temp_neg = positions[i - 1];

					double precedent_gap = t[temp] - t[temp_neg];

					double successive_gap = (t[temp_plus] - t[temp]);

					if ((precedent_gap < (threshold * 1.0)) && (successive_gap < (threshold * 1.0)))
					{

						continue;

					}

					else
					{

						for (size_t j = 0; j < ((2 * n) - 1); j++) // n is the number of small scale velocity signals and time cut for each zero-crossing position
						{

							int loc = (temp - (n - 1) + j);

							double vel = ((1 - w.at(i)) * small.at(loc)) + (w.at(i) * small.at(loc + 1));

							vel_buf[j] += vel;

							num_samples[j] += 1;


						}


					}

				}

			}


		}


	}




}




int main()
{

	
	Matrix_dimz array_neg;
	//array_neg = read_int("neg_pos.out");


	Matrix_dimz  array_pos;
	//array_pos = read_int("pos_neg.out");


	std::vector<int> size_array_neg;
	//size_array_neg = read_int("size_neg_pos.out");


	std::vector<int> size_array_pos;
	//size_array_pos = read_int("size_pos_neg.out");

	std::vector<double> Ul;
	//Ul = read_double("ul.out");


	std::vector<double> Us;
	Us = read_double("us.out");


	Pair zero_result;

	zero_result = find_zero_crossings("ul.out");


	array_pos = zero_result.first.first[0];
	array_neg = zero_result.first.first[1];
	

	Ul = zero_result.second;


	//std::cout << Ul.size();

	const int num = Ul.size() / DIMZ;


	//using Matrix = std::array<std::vector<double>, 2>;

	//std::vector < std::vector<double> > ul;

	//ul.resize( DIMZ, std::vector<double> (num));
	//ul.resize((std::vector<double>(DIMZ), num));

	Matrix_dimz_double ul;

	// we are going to call the function from zero_crossings.cpp
	// no need to read files in both places:
	// pass the data as argument, for instance, or as result

	// read ul here
	// then pass a pointer to ul to find_zero_crossings
	// (you just need to copy paste some code around ;) )



	for (int j = 0; j < DIMZ; ++j)
	{
		ul[j].resize(num);

		//std::cout << "f" << ul[j].size() << std::endl;
		for (int i = 0; i < num; ++i)
		{
			
			ul[j][i] = (Ul.at((i * DIMZ) + j));

		}

	}


	//std::vector < std::vector<double> > us;

	//us.resize(DIMZ, std::vector<double>(num));
	//us.resize((std::vector<double>(DIMZ), num));

	Matrix_dimz_double us;

	for (int j = 0; j < DIMZ; ++j)
	{
		us[j].resize(num);

		//std::cout << "f" << ul[j].size() << std::endl;
		for (int i = 0; i < num; ++i)
		{

			us[j][i] = (Us.at((i * DIMZ) + j));

		}

	}

	
	
	std::ifstream ifs("inputed.json");

	Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);

	double dt;
    dt = obj["timestep"].asDouble();



	//int starting_size_Pos = 0;

	//std::vector <double> vel_matrix_pos;


	double average_vel_pos[(2 * num_of_cuts) - 1] = { 0 };

	double average_vel_pos_ave[(2 * num_of_cuts) - 1] = { 0 };

	int counter_pos = 0;

	//int max_pos = *std::max_element(size_array_pos.begin(), size_array_pos.end());

	for (int z = 0; z < DIMZ; z++)
	{


		//std::vector<double> cut_us_P; // cut small velocity and time for positive to negative

		//std::vector<double> cut_us_Pos;


	/*
		std::vector <int> local_Pos;

		local_Pos.resize(size_array_pos.at(z));

		for (int i = 0; i < size_array_pos.at(z); i++) // loop over zero crossings at the current z
		{

			local_Pos[i] = (array_pos.at(starting_size_Pos + i)); // TODO: CHECK ME - I think this is not needed

		}


		starting_size_Pos = starting_size_Pos + size_array_pos[z];
*/



		double cut_us_Pos[(2 * num_of_cuts) - 1] = { 0 };

		int divider_pos[(2 * num_of_cuts) - 1] = { 0 };

		//std::cout << ul[z].size() << std::endl;

		//cut_us_Pos = average_cut_velocity (cut_us_P, local_Pos, num_of_cuts);

		

		cut_velocity(ul[z], us[z], array_pos[z], "time.out", num_of_cuts, cut_us_Pos, divider_pos);

		
		if (divider_pos[0] != 0)

		{

			counter_pos++;

			//std::cout << counter_pos << std::endl;



		//cut_us_N = cut_velocity(ul[z], us[z], local_Neg, num_of_cuts);

		//local_Pos.clear();



		/*	for (int i = 0; i < cut_us_Pos.size(); i++)
			{

				vel_matrix_pos.push_back(cut_us_Pos.at(i));


			}
			*/


			for (int i = 0; i < ((2 * num_of_cuts) - 1); i++)
			{

				average_vel_pos[i] = cut_us_Pos[i] / divider_pos[i];

			}


			for (int i = 0; i < ((2 * num_of_cuts) - 1); i++)
			{

				average_vel_pos_ave[i] += average_vel_pos[i];

			}

		}
		


	}

	// TODO: fix code below as for the positive case


	//int starting_size_Neg = 0;

	//std::vector <double> vel_matrix_neg;

	

	double average_vel_neg[(2 * num_of_cuts) - 1] = { 0 };

	double average_vel_neg_ave[(2 * num_of_cuts) - 1] = { 0 };

	int counter_neg = 0;

	
	//int max_neg = *std::max_element(size_array_neg.begin(), size_array_neg.end());

	for (int z = 0; z < DIMZ; z++)
	{


		double cut_us_Neg[(2 * num_of_cuts) - 1] = { 0 };

		int divider_neg[(2 * num_of_cuts) - 1] = { 0 };


		//std::vector<double> cut_us_N; // cut small velocity and time for positive to negative

		//std::vector<double> cut_us_Neg;

		
/*
		std::vector <int> local_Neg;

		local_Neg.resize(size_array_neg[z]);

		for (int i = 0; i < size_array_neg[z]; i++)
		{

			local_Neg[i] = array_neg[(starting_size_Neg + i)];

		}

*/
		//starting_size_Neg = starting_size_Neg + size_array_neg[z];


		//cut_us_N = cut_velocity(ul[z], us[z], local_Neg, num_of_cuts);

		

		cut_velocity(ul[z], us[z], array_neg[z], "time.out", num_of_cuts, cut_us_Neg, divider_neg);

		//std::cout << ul[z].size() << std::endl;

		//cut_us_Neg = average_cut_velocity(cut_us_N, local_Neg, num_of_cuts);

		//local_Neg.clear();

		

		std::cout << counter_neg;
		if (divider_neg[0] != 0)

		{

			counter_neg++;
			std::cout << counter_neg << std::endl;




			for (int i = 0; i < ((2 * num_of_cuts) - 1); i++)
			{
				average_vel_neg[i] = cut_us_Neg[i] / divider_neg[i];
				
				//std::cout << cut_us_Neg[i] << std::endl;
			}



			for (int i = 0; i < ((2 * num_of_cuts) - 1); i++)
			{

				average_vel_neg_ave[i] += average_vel_neg[i];

			}

		}
		

	}


	for (int i = 0; i < ((2 * num_of_cuts) - 1); i++)
	{

		average_vel_neg_ave[i] = average_vel_neg_ave[i] / counter_neg;
		

	}

	for (int i = 0; i < ((2 * num_of_cuts) - 1); i++)
	{

		average_vel_pos_ave[i] = average_vel_pos_ave[i] / counter_pos;

	}



	std::vector<double> time;

	time.resize((2 * num_of_cuts) - 1);

	for (int j = 0; j < num_of_cuts - 1; j++)
	{

		time[j] = (dt * (j - (num_of_cuts - 1)));
		//std::cout << time.at(j) << std::endl;

	}


	for (int j = num_of_cuts - 1; j < (2 * num_of_cuts) - 1; j++)
	{

		time[j] = (dt * (j - (num_of_cuts - 1)));
		//std::cout << time.at(j) << std::endl;

	}

	std::fstream file_PT;
	file_PT.open("time_result.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < ((2 * num_of_cuts) - 1); i++)
	{
		file_PT.write((char*)(&time[i]), sizeof(double));
	}




	std::fstream file_P;
	file_P.open("Velocity_pos_neg.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < ((2 * num_of_cuts) - 1); i++)
	{
		file_P.write((char*)(&average_vel_pos_ave[i]), sizeof(double));
	}

	std::cout << "z";

	std::fstream file_N;
	file_N.open("Velocity_neg_pos.out", std::ios::out | std::ios::binary);
	for (int i = 0; i < ((2 * num_of_cuts) - 1); i++)
	{
		file_N.write((char*)(&average_vel_neg_ave[i]), sizeof(double));
	}


	return 0;

}
