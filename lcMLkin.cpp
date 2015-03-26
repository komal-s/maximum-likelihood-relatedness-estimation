
#include <iostream>
#include <fstream>

#include <string>
#include <sstream>
#include <vector>

#include <algorithm>
#include <math.h>
#include <limits>

#include <time.h>
#include <sys/resource.h>
#include <sys/time.h>

#define FILENAME "sim2_5.vcf"
#define IBD_COUNT 3
#define GENOTYPE_COUNT 9

//Split a string with the given delimiter
std::vector<std::string> split(std::string &s, char delim) {
    
    std::vector<std::string> elements;

    std::stringstream stream(s);
    std::string element;
    while (std::getline(stream, element, delim)) {
        elements.push_back(element);
    }

    return elements;
}

//Print time elapsed in seconds
void print_time_elapsed(std::string desc, struct timeval* start, struct timeval* end) {
    
    /*
    struct timeval {
        time_t      tv_sec;
        suseconds_t tv_usec;    
    }*/
    struct timeval elapsed;
    
    if(start->tv_usec > end->tv_usec) {
        end->tv_usec += 1000000;
        end->tv_sec--;
    }
    elapsed.tv_usec = end->tv_usec - start->tv_usec;
    elapsed.tv_sec  = end->tv_sec  - start->tv_sec;
    float time_elapsed = (elapsed.tv_sec*1000000 + elapsed.tv_usec)/1000000.f;
    std::cout << desc << " Total Time Elapsed = " << time_elapsed << std::endl;

    return;
}

void transpose(double** m, double** m_transpose, int r_size, int c_size){
	
	for(int i=0; i<r_size; i++){
    	for(int j=0; j<c_size; j++){
			m_transpose[j][i]=m[i][j];
    	}
    }
}

double kin(std::pair<double,double> k12, double** IBS, int* mask){

	double k0 = 1.0-(k12.first+k12.second);
	double k1 = k12.first;
	double k2 = k12.second;

	std::vector<double> k {k0,k1,k2};
	int snp_count = sizeof(IBS)/sizeof(IBS[0]);

	double **IBS_K = new double*[snp_count];
	for(int i=0; i<snp_count; i++){
		IBS_K[i] = new double[IBD_COUNT]; 
		for(int j=0; j<IBD_COUNT; j++){
			IBS_K[i][j] = 0;
		}
	}
	double *IBS_K_SUM = new double[snp_count];
	for(int i=0; i<snp_count; i++){
		IBS_K_SUM[i]=0;
	}
		
	double IBS_SUM = 0;
	for(int i=0; i<snp_count; i++){
		for(int j=0; j<k.size();j++){
			IBS_K[i][j] = IBS[i][j]*k[j];
			IBS_K_SUM[i] +=IBS_K[i][j];
		}
		IBS_K_SUM[i] = std::log(IBS_K_SUM[i]);
		if(mask[i]==1){
			IBS_K_SUM[i]=0;
		}
		IBS_SUM+=IBS_K_SUM[i];
	}
	IBS_SUM*=-1;

	int pen=0;
	if(k0<0 || k0>1) {pen+=1;}
	if(k1<0 || k1>1) {pen+=1;}
	if(k2<0 || k2>1) {pen+=1;}    
    if(4*k2*k0>=pow(k1,2)) {pen+=1;}
    if (IBS_SUM==std::numeric_limits<double>::infinity()) {pen+=1;}
    
    if(pen>0){
        IBS_SUM=100000;
    }

    /*Clean Up*/
    for(int j=0; j<snp_count; j++){
		delete[] IBS_K[j];
	}
	delete[] IBS_K;
	delete[] IBS_K_SUM;

    return IBS_SUM;
}

//to-do: fix infinity issue
double gl_kin(std::pair<double,double> k12, double** GL, double*** IBS, int* mask){
	
	std::cout<<"gl_kin"<<std::endl;

	double k0 = 1.0-(k12.first+k12.second);
	double k1 = k12.first;
	double k2 = k12.second;

	std::vector<double> kx {k0,k1,k2};

	int snp_count = sizeof(IBS[0])/sizeof(IBS[0][0]);

	double ***IBS_K = new double**[GENOTYPE_COUNT];
	for(int i=0; i<GENOTYPE_COUNT; i++){
		IBS_K[i] = new double*[snp_count];
		for(int j=0; j<snp_count; j++){
			IBS_K[i][j] = new double[IBD_COUNT]; 
			for(int k=0; k<IBD_COUNT; k++){
				IBS_K[i][j][k] = 0;
			}
		}
	}

	double **IBS_K_SUM = new double*[GENOTYPE_COUNT];
	for(int j=0; j<GENOTYPE_COUNT; j++){
		IBS_K_SUM[j] = new double[snp_count]; 
		for(int k=0; k<snp_count; k++){
			IBS_K_SUM[j][k] = 0;
		}
	}
	double *IBS_K_PROD = new double[GENOTYPE_COUNT];
	for(int i=0; i<GENOTYPE_COUNT; i++){
		IBS_K_PROD[i]=0;
	}

	double IBS_SUM = 0;
	
	for(int i=0; i<GENOTYPE_COUNT; i++){

		for(int j=0; j<snp_count;j++){
			for(int k=0; k<IBD_COUNT;k++){
				IBS_K[i][j][k] = IBS[i][j][k]*kx[k];
				IBS_K_SUM[i][j] +=IBS_K[i][j][k];
			}	
		}
		for(int j=0; j<snp_count;j++){
			IBS_K_SUM[i][j]*=GL[i][j];
			IBS_K_PROD[i]+=IBS_K_SUM[i][j];
		}

		IBS_K_PROD[i] = std::log(IBS_K_PROD[i]);

		for(int j=0; j<snp_count;j++){
			if(mask[j]==1){
				IBS_K_PROD[i]=0;
			}
		}
			
		IBS_SUM+=IBS_K_PROD[i];
	}
	IBS_SUM*=-1;
	std::cout<<IBS_SUM<<"\t";
	int pen=0;
	if(k0<0 || k0>1) {pen+=1;}
	if(k1<0 || k1>1) {pen+=1;}
	if(k2<0 || k2>1) {pen+=1;}    
    if(4*k2*k0>=pow(k1,2)) {pen+=1;}
    if (IBS_SUM==std::numeric_limits<double>::infinity()) {pen+=1;}
    
    if(pen>0){
        IBS_SUM=100000;
    }

    //Clean Up
    for(int j=0; j<snp_count; j++){
		delete[] IBS_K_SUM[j];
	}
	delete[] IBS_K_SUM;
	delete[] IBS_K_PROD;

	for(int i=0; i<GENOTYPE_COUNT; i++){
		for(int j=0; j<snp_count; j++){
			delete[] IBS_K[i][j]; 
		}
		delete[] IBS_K[i];
	}
	delete[] IBS_K;	

    return IBS_SUM;
}

void calculate_allele_frequencies(double* allele_frequency, int snp_count, std::vector<int> unrelated_individual_index, std::vector<std::vector<std::string>> snp_data){

	for(int i=0; i<snp_count; i++){
		
		std::vector<int> counts;

		for(int id: unrelated_individual_index){
			std::string geno = split(snp_data[i][id],':')[0];
			if(geno != "./."){ 
				counts.push_back(std::stoi(std::string(1,geno.front())));
				counts.push_back(std::stoi(std::string(1,geno.back())));
			}			
		}
		int zeroes=0;
		for(auto iter=counts.begin(); iter!=counts.end(); iter++){
			if(*iter==0) zeroes++;
		}
		allele_frequency[i++] = zeroes/(double)counts.size();
	}

}

void calculate_ibs(double*** IBS_ALL, int snp_count, double* allele_frequency){	

	/*Populate matrices with actual probablities*/
	for(int i=0; i<snp_count; i++){

		double p = allele_frequency[i];
		double q = 1.0-p;

		IBS_ALL[0][i][0] = pow(p,2)*pow(q,2); 	IBS_ALL[0][i][1] = 0;			 		IBS_ALL[0][i][2] = 0; //PPQQ
    	IBS_ALL[1][i][0] = pow(q,2)*pow(p,2);	IBS_ALL[1][i][1] = 0; 					IBS_ALL[1][i][2] = 0; //QQPP
		
		IBS_ALL[2][i][0] = pow(p,2)*(2*p*q);	IBS_ALL[2][i][1] = (2*p*q)*pow(p,2); 	IBS_ALL[2][i][2] = 0; //PPPQ
		IBS_ALL[3][i][0] = (2*p*q)*pow(p,2);	IBS_ALL[3][i][1] = (p*q)*p; 			IBS_ALL[3][i][2] = 0; //PQPP
		IBS_ALL[4][i][0] = (2*p*q)*pow(q,2);	IBS_ALL[4][i][1] = (p*q)*q; 			IBS_ALL[4][i][2] = 0; //PQQQ
		IBS_ALL[5][i][0] = pow(q,2)*(2*p*q);	IBS_ALL[5][i][1] = pow(q,2)*p;			IBS_ALL[5][i][2] = 0; //QQPQ

		IBS_ALL[6][i][0] = pow(p,4); 			IBS_ALL[6][i][1] = pow(p,3);			IBS_ALL[6][i][2] = pow(p,3); //PPPP
		IBS_ALL[7][i][0] = (2*p*q)*(2*p*q); 	IBS_ALL[7][i][1] = p*q*(p+q);			IBS_ALL[7][i][2] = 2*(p*q);  //PQPQ
		IBS_ALL[8][i][0] = pow(q,4);		 	IBS_ALL[8][i][1] = pow(q,3);			IBS_ALL[8][i][2] = pow(q,2); //QQQQ

	}

}

void calculate_pairwise_likelihood(double** P_IBS, int* mask, int snp_count, std::pair<int,int> pair, double* allele_frequency, std::vector<std::vector<std::string>> snp_data){

	//Parallallize
	/*Iterate through all SNPs*/
	for(int j=0; j<snp_count; j++){

	    double p=allele_frequency[j];
	    double q=1.0-p;

		/*Mask SNPs where allele frequency is fixed*/
		if (allele_frequency[j]==1.0||allele_frequency[j]==0.0){
		    mask[j]=1;
		}
	    
	    /*Pulls out info field for first individual from VCF, pulls out the three precomputed genotype likelihoods*/
		std::vector<std::string> l1 = split(split(snp_data[j][pair.first],':')[3],',');
		std::vector<std::string> l2 = split(split(snp_data[j][pair.second],':')[3],',');

		//Assert: l1==l2?
		double* likelihood_1 = new double[l1.size()];
		double* likelihood_2 = new double[l2.size()];

		/*Convert likelihoods from strings to floating point numbers*/
		for(int k=0; k<l1.size(); k++){
			likelihood_1[k]=std::stod(l1[k]);
			likelihood_2[k]=std::stod(l2[k]);
			/*If one of those likelihoods comes out as negative (anomaly), we mask those SNPs*/
			if(likelihood_1[k]==-9.0 || likelihood_2[k]==-9.0){
				mask[j]=1;
			}
		}

		/*Calculate the probability of observing all possible two genotype combinations by multiplying their likelihoods*/
		P_IBS[j][0]=(likelihood_1[0]*likelihood_2[2]);
	    P_IBS[j][1]=(likelihood_1[2]*likelihood_2[0]);

	    P_IBS[j][2]=(likelihood_1[0]*likelihood_2[1]);
	    P_IBS[j][3]=(likelihood_1[1]*likelihood_2[0]);
	    P_IBS[j][4]=(likelihood_1[1]*likelihood_2[2]);
	    P_IBS[j][5]=(likelihood_1[2]*likelihood_2[1]);

	    P_IBS[j][6]=(likelihood_1[0]*likelihood_2[0]);
	    P_IBS[j][7]=(likelihood_1[1]*likelihood_2[1]);
	    P_IBS[j][8]=(likelihood_1[2]*likelihood_2[2]);

		/*Clean up*/
		delete[] likelihood_1;
		delete[] likelihood_2;
	}
}


int main(){

	//Calculating time elapsed
    struct timeval start, end;
    struct timezone tzp;
    gettimeofday(&start, &tzp);

	//Read the VCF file line by line
	std::ifstream vcfFile (FILENAME);  	
  	std::vector<std::string> data;
  	std::string line; 
  	while (std::getline(vcfFile, line)){
  		if(!line.empty())
 	 		data.push_back(line);
	}
  	vcfFile.close();
 	
  	std::vector<std::string> head= split(data[0],'\t');
  	std::vector<std::string> header (head.begin()+9, head.end());

  	//Fabricated unrelated individuals (removing individuals known to be related)
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
	for(int i=0; i<header.size();i++){
		if(std::find(unrelated_individuals.begin(), unrelated_individuals.end(), header[i]) != unrelated_individuals.end())
			unrelated_individual_index.push_back(i);	
	}

	std::vector<std::vector<std::string>> snp_data;	
	for(auto iter=data.begin()+1; iter!=data.end(); iter++){
		std::vector<std::string> elements = split(*iter,'\t');
  		std::copy(elements.begin()+9, elements.end(),elements.begin());
  		if(elements.size()>1)
	  		snp_data.push_back(elements);
	}

	int snp_count = snp_data.size();
	

	std::cout<<"Calculating Allelle Frequencies"<<std::endl;
	double* allele_frequency = new double[snp_count];
	for(int i=0; i<snp_count; i++){
		allele_frequency[i]=0;
	}
	calculate_allele_frequencies(allele_frequency, snp_count, unrelated_individual_index, snp_data);


	std::cout<<"Calculating Prestored IBS|IBD"<<std::endl;
	/*	For every SNP (as each has a different allele frequency) calculate the P(IBS=x|IBD=z) for all combinations. 
	We only do it for the 3 IBD values, i.e. assuming no inbreeding in the two individuals */
	/*	Make a matrix to store these values - one dimension for each of the 9 possible genotype combinations we might observe 
	One matrix for each possible genotype combination we might observe. 
	If we consider the reference as always p, there are 9 possible combinations.
	In a standard table like Mulligan or Thompson, some of these are collapsed 
	(i.e. ppqq is the same as qqpp, as p has no meaning  has no meaning in that context) */
	/* AATT,TTAA,AAAT,ATAA,ATTT,TTAT,AAAA,ATAT,TTTT */
	double ***IBS_ALL = new double**[GENOTYPE_COUNT];
	for(int i=0; i<GENOTYPE_COUNT; i++){
		IBS_ALL[i] = new double*[snp_count];
		for(int j=0; j<snp_count; j++){
			IBS_ALL[i][j] = new double[IBD_COUNT]; 
			for(int k=0; k<IBD_COUNT; k++){
				IBS_ALL[i][j][k] = 0;
			}
		}
	}
	calculate_ibs(IBS_ALL, snp_count, allele_frequency);	

	/*Pairwise Analysis (of all possible combinations)*/
	//Parallalizable
	std::vector<std::pair<int,int>> pairs;
	for(int i=0; i<header.size(); i++){
		for(int j=i+1; j<header.size(); j++){
			if(split(header[i],'_')[0]==split(header[j],'_')[0])
				pairs.push_back(std::make_pair(i,j));
		}
	}


	std::string output_file (FILENAME +std::string(".relateF_optim_cpp"));
	std::ofstream outfile (output_file);
	outfile << "Ind1\tInd2\tZ0bg\tZ1bg\tZ2bg\tPI_HATbg\tZ0ag\tZ1ag\tZ2ag\tPI_HATag\tnbSNP\n";
	outfile.close();

	std::cout<<"Starting Pairwise IBD Computations"<<std::endl;

	//Parallalizable
	/*Iterate through all pairwise computations */
	for(int i=0; i<pairs.size();i++) {
		
		/* Matrices for all possible pairs of genotypes for every SNP. 
		This will eventually store genotype likelihoods based on the calling likelihood (for example based on read depth)*/
		/*Store these pairwise likelihoods in an array
		AATT,TTAA,AAAT,ATAA,ATTT,TTAT,AAAA,ATAT,TTTT */
		double **P_IBS = new double*[snp_count];
		for(int j=0; j<snp_count; j++){
			P_IBS[j] = new double[GENOTYPE_COUNT]; 
			for(int k=0; k<GENOTYPE_COUNT; k++){
				P_IBS[j][k] = 0;
			}
		}
		
	    /*A matrix to denote SNPs we may want to mask for two reasons(see below)*/
	    int* mask = new int[snp_count];
	    for(int j=0; j<snp_count; j++){
	    	mask[j]=0;
	    }

	    calculate_pairwise_likelihood(P_IBS, mask, snp_count, pairs[i], allele_frequency, snp_data);
	
 		/*Identify the most likely genotype combination: Index of genotype for that SNP*/
 		int* genotype_best = new int[snp_count];
		
 		/*For each SNP, given the best genotype combination, pull out the appropriate P(IBS|IBD) for all three IBS possibilities*/
		double **IBS_best = new double*[snp_count];
		for(int j=0; j<snp_count; j++){
			IBS_best[j] = new double[IBD_COUNT]; 
			for(int k=0; k<IBD_COUNT; k++){
				IBS_best[j][k] = 0;
			}
		}

		for(int j=0; j<snp_count; j++){
			genotype_best[j] = std::distance(P_IBS[j], std::max_element(P_IBS[j], P_IBS[j]+GENOTYPE_COUNT)); 
			for(int k=0; k<IBD_COUNT; k++){
				IBS_best[j][k] = IBS_ALL[genotype_best[j]][j][k];		
			}
			
		}

		/*ML optimization*/
    	srand (time(0));

		for(int j=0; j<snp_count; j++){
    		bool flag = false;
    		while(flag==false){
				double k1 = rand()/(double)(RAND_MAX);
				double k2 = rand()/(double)(RAND_MAX);
				if(k1+k2<=1.0){
					double k0=1.0-(k1+k2);
					if(4*k2*k0<pow(k1,2)){
						std::pair<double,double> k12 = std::make_pair(k1,k2);
						if(kin(k12,IBS_best,mask) != 100000){
							flag = true;
						}
					}
				}
    		}
    		//to-do: minimize function
    		//to-do: opt_parameter <- parameter that minimizes the function
    	}

		double **P_IBS_T = new double*[GENOTYPE_COUNT];
		for(int j=0; j<GENOTYPE_COUNT; j++){
			P_IBS_T[j] = new double[snp_count]; 
			for(int k=0; k<snp_count; k++){
				P_IBS_T[j][k] = 0;
			}
		}
		transpose(P_IBS, P_IBS_T,snp_count,GENOTYPE_COUNT);
  	
		for(int j=0; j<snp_count; j++){
    		bool flag = false;
    		while(flag==false){
				double k1 = rand()/(double)(RAND_MAX);
				double k2 = rand()/(double)(RAND_MAX);
				if(k1+k2<=1.0){
					double k0=1.0-(k1+k2);
					if(4*k2*k0<pow(k1,2)){
						std::pair<double,double> k12 = std::make_pair(k1,k2);
						if(gl_kin(k12,P_IBS_T,IBS_ALL,mask) != 100000){
							flag = true;
						}
					}
				}
    		}
    		//to-do: minimize function
    		//to-do: opt_parameter <- parameter that minimizes the function
    	}

   		std::ofstream outfile (output_file);
		outfile << pairs[i].first << "\t" << pairs[i].second;
		outfile << "opt_parameter \t results \n";
		outfile.close();

 		/*Clean up*/
 		for(int j=0; j<snp_count; j++){
			delete[] P_IBS[j];
			delete[] IBS_best[j];
		}
		delete[] P_IBS;
 		delete[] IBS_best;
 		delete[] mask;		
 		delete[] genotype_best;
	}

	/*Clean up*/
	for(int i=0; i<GENOTYPE_COUNT; i++){
		for(int j=0; j<snp_count; j++){
			delete[] IBS_ALL[i][j]; 
		}
		delete[] IBS_ALL[i];
	}
	delete[] IBS_ALL;	
	delete[] allele_frequency;

	gettimeofday(&end, &tzp);
    print_time_elapsed("", &start, &end);
  	return 0;
}