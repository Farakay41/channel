

#include "Zero_Cross.h"
#include <jsoncpp/json/json.h>


constexpr auto num_of_cuts = 20;



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

   //double dt = 0.000127175055596395;
	
	/*
	for (int i = 0; i < large.size(); ++i)
	{
		t.push_back(i * dt);
	}
	*/

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

		//std::cout << "tw" << large.size();
		//std::cout << "tw2" << n;

		if ((temp <= n) || ( (large.size() - temp) <= n) )
		{


			//std::cout << "3e" << std::endl;

			continue;

			std::cout << "ue";

		}

		else
		{

			switch (method)
			{

			case 1:

				for (size_t j = 0; j < ((2 * n) - 1); j++) // n is the number of small scale velocity signals and time cut for each zero-crossing position

				{
					//std::cout << "��f" << temp << std::endl;
					int loc = (temp - (n - 1) + j);



					double vel = ((1 - w.at(i)) * small.at(loc)) + (w.at(i) * small.at(loc + 1));


					//cut_vel.push_back(vel); // TODO: don't do this - you know the size of vel!

					// TODO: this is how you write to cut_vel


					vel_buf[j] += vel;

					num_samples[j] += 1;


				}

			case 2:

				//for (size_t i = 0; i < size; i++) // for each of the zero-crossing positions
				//{



					//int  temp = positions[i]; // allocates a variable for each zero crossing positions


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
								//std::cout << "wer" << temp;
								int loc = (temp - (n - 1) + j);


								double vel = ((1 - w.at(i)) * small.at(loc)) + (w.at(i) * small.at(loc + 1));

								vel_buf[j] += vel;

								num_samples[j] += 1;


							}


						}
					}

				//}

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


	std::cout << "oa" << Ul.size();
	std::cout << "oe" << Us.size();

	const int num = Ul.size() / DIMZ;
	std::cout << "oe" << num;

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
			
			ul[j][i] = Ul[(i * DIMZ) + j];

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

			us[j][i] = Us[(i * DIMZ) + j];

		}

	}

	
	
	std::ifstream ifs("inputed.json");

	Json::Reader reader;
    Json::Value obj;
    reader.parse(ifs, obj);

	double dt;
    dt = obj["timestep"].asDouble();

	int p;
	p = obj["z-position"].asInt();



	double average_vel_pos[(2 * num_of_cuts) - 1] = { 0 };

	double average_vel_pos_ave[(2 * num_of_cuts) - 1] = { 0 };

	int counter_pos = 0;
	double cut_us_Pos[(2 * num_of_cuts) - 1] = { 0 };

	int divider_pos[(2 * num_of_cuts) - 1] = { 0 };



	double average_vel_neg[(2 * num_of_cuts) - 1] = { 0 };

	double average_vel_neg_ave[(2 * num_of_cuts) - 1] = { 0 };

	double cut_us_Neg[(2 * num_of_cuts) - 1] = { 0 };

	int divider_neg[(2 * num_of_cuts) - 1] = { 0 };

	

	int counter_neg = 0;

		cut_velocity(ul[p], us[p], array_pos[p], "time.out", num_of_cuts, cut_us_Pos, divider_pos);

		for (int i = 0; i < ((2 * num_of_cuts) - 1); i++)
		{

			average_vel_pos[i] = cut_us_Pos[i] / divider_pos[i];

		}


		std::fstream file_P;
		file_P.open("Velocity_pos_neg.out", std::ios::out | std::ios::binary);
		for (int i = 0; i < ((2 * num_of_cuts) - 1); i++)
		{
			file_P.write((char*)(&average_vel_pos[i]), sizeof(double));

		}


		//double average_vel_neg[(2 * num_of_cuts) - 1] = { 0 };

		//double average_vel_pos_ave[(2 * num_of_cuts) - 1] = { 0 };

		//double cut_us_Neg[(2 * num_of_cuts) - 1] = { 0 };

		//int divider_neg[(2 * num_of_cuts) - 1] = { 0 };

		cut_velocity(ul[p], us[p], array_neg[p], "time.out", num_of_cuts, cut_us_Neg, divider_neg);

		for (int i = 0; i < ((2 * num_of_cuts) - 1); i++)
		{

			average_vel_neg[i] = cut_us_Neg[i] / divider_neg[i];

		}


		std::fstream file_N;
		file_N.open("Velocity_neg_pos.out", std::ios::out | std::ios::binary);
		for (int i = 0; i < ((2 * num_of_cuts) - 1); i++)
		{
			file_N.write((char*)(&average_vel_neg[i]), sizeof(double));
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

	return 0;


}
