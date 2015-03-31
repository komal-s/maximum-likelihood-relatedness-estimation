#include <iostream>

#include "relatedness.h"
#include "utils.h"

void relatedness::populate_data(std::string filename){

	//Read the VCF file line by line
	std::ifstream vcfFile (filename);  	
  	std::vector<std::string> data;
  	std::string line; 
  	while (std::getline(vcfFile, line)){
  		if(!line.empty())
 	 		data.push_back(line);
	}
  	vcfFile.close();
 	
  	std::vector<std::string> head= split(data[0],'\t');
  	header = {head.begin()+9, head.end()};

	for(int i=0; i<header.size();i++){
		if(std::find(unrelated_individuals.begin(), unrelated_individuals.end(), header[i]) != unrelated_individuals.end())
			unrelated_individual_index.push_back(i);	
	}

	for(auto iter=data.begin()+1; iter!=data.end(); iter++){
		std::vector<std::string> elements = split(*iter,'\t');
  		std::copy(elements.begin()+9, elements.end(),elements.begin());
  		if(elements.size()>1)
	  		snp_data.push_back(elements);
	}

	snp_count = snp_data.size();

}

void relatedness::calculate_allele_frequencies(){

	allele_frequency = Eigen::VectorXd::Zero(snp_count);

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
		allele_frequency(i++) = zeroes/(double)counts.size();
	}

}

void relatedness::calculate_ibs(){	

	//ibs_all = Eigen::MatrixXd::Zero(snp_count, IBD_COUNT);	
	for(int i=0; i<GENOTYPE_COUNT; i++){
		ibs_all(i) = Eigen::MatrixXd::Zero(snp_count,IBD_COUNT);
	}

	/*Populate matrices with actual probablities*/
	for(int i=0; i<snp_count; i++){

		double p = allele_frequency(i);
		double q = 1.0-p;

		ibs_all(0)(i,0) = pow(p,2)*pow(q,2); 	ibs_all(0)(i,1) = 0;			 		ibs_all(0)(i,2) = 0; //PPQQ
    	ibs_all(1)(i,0) = pow(q,2)*pow(p,2);	ibs_all(1)(i,1) = 0; 					ibs_all(1)(i,2) = 0; //QQPP
		
		ibs_all(2)(i,0) = pow(p,2)*(2*p*q);		ibs_all(2)(i,1) = (2*p*q)*pow(p,2); 	ibs_all(2)(i,2) = 0; //PPPQ
		ibs_all(3)(i,0) = (2*p*q)*pow(p,2);		ibs_all(3)(i,1) = (p*q)*p; 				ibs_all(3)(i,2) = 0; //PQPP
		ibs_all(4)(i,0) = (2*p*q)*pow(q,2);		ibs_all(4)(i,1) = (p*q)*q; 				ibs_all(4)(i,2) = 0; //PQQQ
		ibs_all(5)(i,0) = pow(q,2)*(2*p*q);		ibs_all(5)(i,1) = pow(q,2)*p;			ibs_all(5)(i,2) = 0; //QQPQ

		ibs_all(6)(i,0) = pow(p,4); 			ibs_all(6)(i,1) = pow(p,3);				ibs_all(6)(i,2) = pow(p,3); //PPPP
		ibs_all(7)(i,0) = (2*p*q)*(2*p*q); 		ibs_all(7)(i,1) = p*q*(p+q);			ibs_all(7)(i,2) = 2*(p*q);  //PQPQ
		ibs_all(8)(i,0) = pow(q,4);		 		ibs_all(8)(i,1) = pow(q,3);				ibs_all(8)(i,2) = pow(q,2); //QQQQ

	}

}

void relatedness::calculate_pairwise_ibd(){

	//Generate Pairs
	for(int i=0; i<header.size(); i++){
		for(int j=i+1; j<header.size(); j++){
			if(split(header[i],'_')[0]==split(header[j],'_')[0])
				pairs.push_back(std::make_pair(i,j));
		}
	}

	//Parallalizable
	/*Iterate through all pairwise computations */
	for(int i=0; i<pairs.size();i++) {
					
		/* Matrices for all possible pairs of genotypes for every SNP.
		This will eventually store genotype likelihoods based on the calling likelihood (for example based on read depth)*/
		/*Store these pairwise likelihoods in an array AATT,TTAA,AAAT,ATAA,ATTT,TTAT,AAAA,ATAT,TTTT */
		ibs_pairwise = Eigen::MatrixXd::Zero(snp_count,GENOTYPE_COUNT);

		/*A matrix to denote SNPs we may want to mask for two reasons(see below)*/
		mask_snp = Eigen::VectorXd::Zero(snp_count);

		calculate_pairwise_likelihood(pairs[i]);

		/*Identify the most likely genotype combination: Index of genotype for that SNP*/
		/*For each SNP, given the best genotype combination, pull out the appropriate P(IBS|IBD) for all three IBS possibilities*/
		ibs_best = Eigen::MatrixXd::Zero(snp_count,GENOTYPE_COUNT);

		for(int j=0; j<snp_count; j++){
			Eigen::MatrixXf::Index bestIndex;
			double max = ibs_pairwise.row(j).maxCoeff(&bestIndex);

			for(int k=0; k<IBD_COUNT; k++){
				ibs_best(j,k) = ibs_all(bestIndex)(j,k);		
			}	
		}

		optimize_parameters();
	}
}

void relatedness::calculate_pairwise_likelihood(std::pair<int,int> pair){

	//Parallallize
	/*Iterate through all SNPs*/
	for(int j=0; j<snp_count; j++){

	    double p=allele_frequency(j);
	    double q=1.0-p;

		/*Mask SNPs where allele frequency is fixed*/
		if (allele_frequency(j)==1.0||allele_frequency(j)==0.0){
		    mask_snp(j)=1;
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
				mask_snp[j]=1;
			}
		}

		/*Calculate the probability of observing all possible two genotype combinations by multiplying their likelihoods*/
		ibs_pairwise(j,0)=(likelihood_1[0]*likelihood_2[2]);
	    ibs_pairwise(j,1)=(likelihood_1[2]*likelihood_2[0]);

	    ibs_pairwise(j,2)=(likelihood_1[0]*likelihood_2[1]);
	    ibs_pairwise(j,3)=(likelihood_1[1]*likelihood_2[0]);
	    ibs_pairwise(j,4)=(likelihood_1[1]*likelihood_2[2]);
	    ibs_pairwise(j,5)=(likelihood_1[2]*likelihood_2[1]);

	    ibs_pairwise(j,6)=(likelihood_1[0]*likelihood_2[0]);
	    ibs_pairwise(j,7)=(likelihood_1[1]*likelihood_2[1]);
	    ibs_pairwise(j,8)=(likelihood_1[2]*likelihood_2[2]);

		/*Clean up*/
		delete[] likelihood_1;
		delete[] likelihood_2;
	}
}


void relatedness::optimize_parameters(){

	srand (time(0));

	for(int j=0; j<3; j++){
		bool flag = false;
		while(flag==false){
			double k1 = rand()/(double)(RAND_MAX);
			double k2 = rand()/(double)(RAND_MAX);
			if(k1+k2<=1.0){
				double k0=1.0-(k1+k2);
				if(4*k2*k0<pow(k1,2)){
					std::pair<double,double> k12 = std::make_pair(k1,k2);
					if(kin(k12) != 100000){
						flag = true;
					}
				}
			}
		}
		//to-do: minimize function
		//to-do: opt_parameter <- parameter that minimizes the function
	}

	for(int j=0; j<3; j++){
		bool flag = false;
		while(flag==false){
			double k1 = rand()/(double)(RAND_MAX);
			double k2 = rand()/(double)(RAND_MAX);
			if(k1+k2<=1.0){
				double k0=1.0-(k1+k2);
				if(4*k2*k0<pow(k1,2)){
					std::pair<double,double> k12 = std::make_pair(k1,k2);
					if(gl_kin(k12) != 100000){
						flag = true;
					}
				}
			}
		}
		//to-do: minimize function
		//to-do: opt_parameter <- parameter that minimizes the function
	}

}

double relatedness::kin(std::pair<double,double> k12){

	double k0 = 1.0-(k12.first+k12.second);
	double k1 = k12.first;
	double k2 = k12.second;

	Eigen::Vector3d k;
	k << k0,k1,k2;

	Eigen::VectorXd ibs_k_sum = Eigen::VectorXd::Zero(snp_count);
	double ibs_sum;
	for(int i=0; i<snp_count; i++){
		for(int j=0; j<IBD_COUNT;j++){
			ibs_k_sum(i)+=ibs_best(i,j)*k(j);
		}
		if(mask_snp(i)==1){
			ibs_k_sum(i)=0;
		}else{
			ibs_k_sum(i) = std::log(ibs_k_sum(i));
		}
		ibs_sum+=ibs_k_sum(i);
	}
	ibs_sum*=-1;

	bool flag=false;
	if(k0<0 || k0>1) {flag=true;}
	if(k1<0 || k1>1) {flag=true;}
	if(k2<0 || k2>1) {flag=true;}    
    if(4*k2*k0>=pow(k1,2)) {flag=true;}
    if (ibs_sum==std::numeric_limits<double>::infinity()) {flag=true;}
    
    if(flag==true){
        ibs_sum=100000;
    }

    return ibs_sum;
}

double relatedness::gl_kin(std::pair<double,double> k12){

	double k0 = 1.0-(k12.first+k12.second);
	double k1 = k12.first;
	double k2 = k12.second;

	Eigen::Vector3d k_values;
	k_values << k0,k1,k2;

	Eigen::MatrixXd ibs_pairwise_t = Eigen::MatrixXd::Zero(GENOTYPE_COUNT,snp_count);
	ibs_pairwise_t = ibs_pairwise.transpose();

	Eigen::MatrixXd ibs_k = Eigen::MatrixXd::Zero(GENOTYPE_COUNT,snp_count);
	Eigen::VectorXd ibs_k_sum = Eigen::VectorXd::Zero(snp_count);
	double ibs_sum;
	
	for(int i=0; i<GENOTYPE_COUNT; i++){
		for(int j=0; j<snp_count; j++){
			for(int k=0; k<IBD_COUNT;k++){
				ibs_k(i,j)+=ibs_all(i)(j,k)*k_values(k);
			}
		}
		for(int j=0; j<snp_count; j++){
			ibs_k(i,j)*=ibs_pairwise_t(i,j);
		}
	}
	for(int i=0; i<snp_count; i++){
		for(int j=0; j<GENOTYPE_COUNT; j++){
			ibs_k_sum(i)+=ibs_k(j,i);
		}
		if(mask_snp(i)==1){
			ibs_k_sum(i)=0;
		}else{
			ibs_k_sum(i) = std::log(ibs_k_sum(i));
		}
		ibs_sum+=ibs_k_sum(i);
	}
	ibs_sum*=-1;

	bool flag=false;
	if(k0<0 || k0>1) {flag=true;}
	if(k1<0 || k1>1) {flag=true;}
	if(k2<0 || k2>1) {flag=true;}    
    if(4*k2*k0>=pow(k1,2)) {flag=true;}
    if (ibs_sum==std::numeric_limits<double>::infinity()) {flag=true;}
    
    if(flag==true){
        ibs_sum=100000;
    }

    return ibs_sum;
}


int main(){

	struct timeval start, end;
    struct timezone tzp;
    gettimeofday(&start, &tzp);

    relatedness r;
    std::string filename ("data/sim2_5.vcf");

    std::cout<<"Populating Data"<<std::endl;
    r.populate_data(filename);

    std::cout<<"Calculating Allelle Frequencies"<<std::endl;
    r.calculate_allele_frequencies();

	std::cout<<"Calculating Prestored IBS|IBD"<<std::endl;
	r.calculate_ibs();    

	std::cout<<"Starting Pairwise IBD Computations"<<std::endl;
	r.calculate_pairwise_ibd();

    gettimeofday(&end, &tzp);
    print_time_elapsed("", &start, &end);

	return 0;
}