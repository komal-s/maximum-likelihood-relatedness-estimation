#ifndef RELATEDNESS_H
#define RELATEDNESS_H

#include <iostream>
#include <fstream>

#include <string>
#include <sstream>
#include <vector>

#include <algorithm>
#include <math.h>
#include <limits>
#include <Eigen/Dense>

#include <time.h>
#include <sys/resource.h>
#include <sys/time.h>

#define IBD_COUNT 3
#define GENOTYPE_COUNT 9
#define FILENAME "data/sim2_5.vcf"

class relatedness {

private:
		
	std::vector<std::string> header;

	std::vector<std::vector<std::string>> snp_data;	
	
	int snp_count;

	std::vector<std::string> unrelated_individuals {"1_0", "1_1", "1_2", "1_3", "1_4", "1_5", "1_6", "1_7",
								 	    		  	"2_0", "2_1", "2_2", "2_3", "2_4", "2_5", "2_6", "2_7",
										    	    "3_0", "3_1", "3_2", "3_3", "3_4", "3_5", "3_6", "3_7",
										        	"4_0", "4_1", "4_2", "4_3", "4_4", "4_5", "4_6", "4_7",
										    	    "5_0", "5_1", "5_2", "5_3", "5_4", "5_5", "5_6", "5_7",
											        "6_0", "6_1", "6_2", "6_3", "6_4", "6_5", "6_6", "6_7",
											        "7_0", "7_1", "7_2", "7_3", "7_4", "7_5", "7_6", "7_7",
										    	    "8_0", "8_1", "8_2", "8_3", "8_4", "8_5", "8_6", "8_7",
										    	    "9_0", "9_1", "9_2", "9_3", "9_4", "9_5", "9_6", "9_7",
										        	"10_0", "10_1", "10_2", "10_3", "10_4", "10_5", "10_6", "10_7"};
	std::vector<int> unrelated_individual_index;

	Eigen::VectorXd allele_frequency;

	Eigen::Array<Eigen::MatrixXd, GENOTYPE_COUNT, 1> ibs_all;	//Eigen::Array of GENOTYPE_COUNT Eigen::Matrix

	std::vector<std::pair<int,int>> pairs;

	Eigen::MatrixXd ibs_pairwise;

	Eigen::VectorXd mask_snp;

	Eigen::MatrixXd ibs_best;

public:

	void populate_data(std::string);

	void calculate_allele_frequencies();

	void calculate_ibs();

	void calculate_pairwise_likelihood(std::pair<int,int>);

	void calculate_pairwise_ibd();

	Eigen::Vector3d optimize_parameters();

	double kin(std::pair<double,double>);

	double gl_kin(std::pair<double,double>);

	Eigen::Vector3d em_optimization(Eigen::Vector3d k_values);

};

#endif