import subprocess
from subprocess import Popen
import os
import string
import math
import random
import numpy as np
import datetime
import numpy.ma as ma
import numpy.linalg
import multiprocessing
from multiprocessing import Pool, Process, Queue

# Python *sucks* at multi-threading (hence the multi-*processing*)
# module.  However, this is a "poor-man's" attempt at parallelizing
# the lcMLKin code.
class RelatednessCalculator(multiprocessing.Process):

    def __init__(self, args, nbSNPs, datum, IBS_all, AF, taskQueue, resultQueue):
        multiprocessing.Process.__init__(self)
        self.args = args
        self.nbSNPs = nbSNPs
        self.datum = datum
        self.IBS_all = IBS_all
        self.AF = AF
        self.taskQueue = taskQueue
        self.resultQueue = resultQueue

    def run(self):
        while True:
            nextTask = self.taskQueue.get()
            if nextTask is None:
                print("{}: Exiting".format(self.name))
                self.taskQueue.task_done()
                break

            id1 = nextTask[0]
            id2 = nextTask[1]
            print("processing {} vs {}".format(id1, id2))
            # matrix for all possible pairs of genotypes for every SNP.
            # This will eventually store genotype likelihoods based on
            # the calling likelihood (for example based on read depth)
            PIBS=np.zeros((self.nbSNPs,9),float)
            # matrix all possible pairs of genotypes for every SNP.
            # This will eventually store genotype likelihoods based on
            # the calling likelihood (for example based on read depth)
            PIBSg=np.zeros((self.nbSNPs,9),float)

            # A matrix to denote SNPs we may want to mask for two reasons(see below)
            mask_mat=np.zeros(self.nbSNPs,int)

            # iterate through each SNP
            for gg in range(self.nbSNPs):
                # mask SNPs where allele frequency is fixed
                if self.AF[gg]==1.0:
                    mask_mat[gg]=1
                if self.AF[gg]==0.0:
                    mask_mat[gg]=1

                # pulls out info field for first individual from VCF,
                # pulls out the three precomputed genotype likelihoods
                l1=string.split(string.split(self.datum[gg][id1],':')[3],',')
                l2=string.split(string.split(self.datum[gg][id2],':')[3],',')

                # converts likelihoods from strings to floating point numbers
                for ggg in range(len(l1)):
                    l1[ggg]=float(l1[ggg])
                    l2[ggg]=float(l2[ggg])
                    # if one of those likelihoods comes out as
                    # negative (some anomaly) we mask that SNPs
                    if l1[ggg]==-9.0:
                        mask_mat[gg]=1
                    if l2[ggg]==-9.0:
                        mask_mat[gg]=1

                # We now calculate the probability of observing all possible two
                # genotype combination by multiplying their likelihoods together
                PPQQ=(l1[0]*l2[2])
                QQPP=(l1[2]*l2[0])

                PPPQ=(l1[0]*l2[1])
                PQPP=(l1[1]*l2[0])
                PQQQ=(l1[1]*l2[2])
                QQPQ=(l1[2]*l2[1])

                PPPP=(l1[0]*l2[0])
                PQPQ=(l1[1]*l2[1])
                QQQQ=(l1[2]*l2[2])

                # Store these pairwise likelihoods in an array
                PIBS[gg]=[PPQQ,QQPP,PPPQ,PQPP,PQQQ,QQPQ,PPPP,PQPQ,QQQQ]  #AATT,TTAA,AAAT,ATAA,ATTT,TTAT,AAAA,ATAT,TTTT

            ###We transpose the matrices
            PIBSt=PIBS.transpose()
            #PIBStg=PIBSg.transpose()

            ###identify the most likely genotype combination
            BestGT=PIBS.argmax(axis=1)

            BestIBS=np.zeros((self.nbSNPs,3),float)

            ###For each SNP, given the best genotype combination, pulls out the appropriate P(IBS|IBD) for all three IBS possibilities
            for gg in range(self.nbSNPs):
                BestIBS[gg]=self.IBS_all[BestGT[gg]][gg]


            def EMOpt(k, IBS, mask_mat, GL = None, verbose = False):
                """Determine the k coefficients via the EM algorithm

                    k --- An initial guess of the k coefficient
                    IBS --- The matrix that gives the identity by site probs
                    mask_mat --- The matrix that masks SNPs that should not be
                                 considered
                    Keyword arguments:
                    GL --- The matrix of genotype likelihoods if we are
                           considering all (rather than just the best)
                           genotypes.
                    verbose --- If True, print out periodic progress of
                                the EM algorithm, if False, be silent.
                """

                useAllGenotypes = GL is not None
                # number of rows and columns depends on
                # whether IBS is a 3D or 2D matrix
                if useAllGenotypes:
                    nr = IBS[0].shape[0]
                    nc = IBS[0].shape[1]
                else:
                    nr = IBS.shape[0]
                    nc = IBS.shape[1]

                # This matrix will store the probabilities of
                # identity-by-descent for each SNP
                PIBD = np.zeros((nr, nc))
                k3 = k
                # The current norm of the difference between
                # subsequent parameter estimates
                d = 100
                # The current iteration number
                it = 0

                # Arbitrary threshold (should be a parameter?)
                while d > 1e-4:
                    # Num SNPs x 3
                    X = np.zeros((nr, nc))
                    X = IBS * k3

                    # If we use *all* genotypes
                    if useAllGenotypes:
                        for i in xrange(9):
                            # Pr(X_m=j | S_{i,m}) * k_{j}
                            # Probability of IDB = k given IBS
                            X[i, :, 0] *= GL[i,:]
                            X[i, :, 1] *= GL[i,:]
                            X[i, :, 2] *= GL[i,:]
                        X = X.sum(axis=0)

                    # Normalize the X matrix to contain probabilities
                    XS = ma.masked_array(X.sum(axis=1), mask=mask_mat)
                    X[:,0] /= XS
                    X[:,1] /= XS
                    X[:,2] /= XS
                    # Copy over into the PIBD matrix
                    PIBD = X

                    # The new parameter estimate
                    kp = np.zeros(3)
                    # The sum of probabilities at each site
                    kp[0] += np.sum(ma.masked_array(PIBD[:,0], mask=mask_mat))
                    kp[1] += np.sum(ma.masked_array(PIBD[:,1], mask=mask_mat))
                    kp[2] += np.sum(ma.masked_array(PIBD[:,2], mask=mask_mat))
                    # Normalized
                    kp /= kp.sum()
                    # Compute the difference between successive
                    # estimates to assess convergence
                    d = np.linalg.norm(np.abs(k3-kp))
                    k3 = kp
                    it += 1
                    if verbose and it % 100 == 0:
                        print(k3)
                        print("diff = {}".format(d))
                # Return what we did all that work to calculate
                return k3

            kest = np.zeros(3)
            if self.args.genotype == "best":
                k3 = np.random.rand(3)
                k3 /= k3.sum()
                kest = EMOpt(k3, BestIBS, mask_mat)
            elif self.args.genotype == "all":
                k3 = np.random.rand(3)
                k3 /= k3.sum()
                kest = EMOpt(k3, self.IBS_all, mask_mat, PIBSt)


            self.taskQueue.task_done()
            numUsed = len(mask_mat)-np.sum(mask_mat)
            self.resultQueue.put([id1, id2, kest[0], kest[1], kest[2], numUsed])
        return

def runMLKin(args):
    filenamein = args.vcf
    file = open(filenamein)
    data=file.read()
    data=string.split(data,'\n')
    file.close()

    if data[-1]=='':
        del(data[-1])

    unrels=['1_0', '1_1', '1_2', '1_3', '1_4', '1_5', '1_6', '1_7',
            '2_0', '2_1', '2_2', '2_3', '2_4', '2_5', '2_6', '2_7',
            '3_0', '3_1', '3_2', '3_3', '3_4', '3_5', '3_6', '3_7',
            '4_0', '4_1', '4_2', '4_3', '4_4', '4_5', '4_6', '4_7',
            '5_0', '5_1', '5_2', '5_3', '5_4', '5_5', '5_6', '5_7',
            '6_0', '6_1', '6_2', '6_3', '6_4', '6_5', '6_6', '6_7',
            '7_0', '7_1', '7_2', '7_3', '7_4', '7_5', '7_6', '7_7',
            '8_0', '8_1', '8_2', '8_3', '8_4', '8_5', '8_6', '8_7',
            '9_0', '9_1', '9_2', '9_3', '9_4', '9_5', '9_6', '9_7',
            '10_0', '10_1', '10_2', '10_3', '10_4', '10_5', '10_6', '10_7']


    head=string.split(data[0],'\t')[9:]
    unrel_ind=[]

    for g in range(len(head)):
        if head[g] in unrels:
            unrel_ind.append(g)

    datum=[]
    for g in range(1,len(data)):
        k=string.split(data[g],'\t')
        datum.append(k[9:])

    data=[] #free memory of initial file
    del data

    nbSNPs=len(datum)

    AF=np.zeros(nbSNPs,float)

    print('calculating allele frequencies')
    for g in range(len(datum)):
        counts=[]
        for gg in range(len(unrel_ind)):
            geno=string.split(datum[g][unrel_ind[gg]],':')[0]
            if geno != './.':
                counts.append(geno[0])
                counts.append(geno[-1])
        AF[g]=counts.count('0')/float(len(counts))


    print 'calculating prestored IBS|IBD'
    ## For every SNP (as each has a different allele frequency) calculate the
    ## P(IBS=x|IBD=z) for all combinations. We only do it for the 3 IBD values,
    ## i.e. assuming no inbreeding in the two individuals

    ## Makes a matrix to store these values. One matrix for each possible
    ## genotype combination we might observe. If we consider the reference as
    ## always p, there are 9 possible combinations.

    ## In a standard table like Mulligan or Thompson, some of these are
    ## collapsed (i.e. ppqq is the same as qqpp, as p has no meaning  has no
    ## meaning in that context)


    #IBS0=np.zeros((nbSNPs,3),float)   #IBD0,IBD1,IBD2
    ppqq=np.zeros((nbSNPs,3),float)
    qqpp=np.zeros((nbSNPs,3),float)

    #IBS1=np.zeros((nbSNPs,3),float)   #IBD0,IBD1,IBD2
    pppq=np.zeros((nbSNPs,3),float)
    pqpp=np.zeros((nbSNPs,3),float)
    pqqq=np.zeros((nbSNPs,3),float)
    qqpq=np.zeros((nbSNPs,3),float)

    #IBS2=np.zeros((nbSNPs,3),float)   #IBD0,IBD1,IBD2
    pppp=np.zeros((nbSNPs,3),float)
    pqpq=np.zeros((nbSNPs,3),float)
    qqqq=np.zeros((nbSNPs,3),float)


    #populates matrix with actual probabilities
    for g in range(len(AF)):
        p=AF[g]
        q=1.0-p
        #IBS0[g]=[2.0*(p**2)*(q**2),0,0]
        ppqq[g]=[(p**2)*(q**2),0,0]
        qqpp[g]=[(q**2)*(p**2),0,0]

        #IBS1[g]=[4*(p**3)*q+4*p*(q**3),2*(p**2)*q+2*p*(q**2),0]
        pppq[g]=[(p**2)*(2*p*q),(p**2)*q,0]
        pqpp[g]=[(2*p*q)*(p**2),(p*q)*p,0]
        pqqq[g]=[(2*p*q)*(q**2),(p*q)*q,0]
        qqpq[g]=[(q**2)*(2*p*q),(q**2)*p,0]


        #IBS2[g]=[p**4+q**4+4*(p**2)*(q**2),p**3+q**3+(p**2)*q+p*(q**2),1]
        pppp[g]=[(p**4),p**3,p**2]
        pqpq[g]=[(2*p*q)*(2*p*q),p*q*(p+q),2*(p*q)]
        qqqq[g]=[(q**4),q**3,q**2]


    #Store in a single matrix, one dimension for each of the 9 possible genotype combinations we might observe
    IBS_all=np.array([
            ppqq, #AATT
            qqpp, #TTAA
            pppq, #AAAT
            pqpp, #ATAA
            pqqq, #ATTT
            qqpq, #TTAT
            pppp, #AAAA
            pqpq, #ATAT
            qqqq])#TTTT

    #pairwise analysis

    pw=[[0,9],[0,10],[0,11],[0,12],[0,13],[0,14],[0,15]]  #list of pairwise comparisons to make
    pw=[]   ##here we do every combination
    for g in range(len(head)):
        for gg in range(g+1,len(head)):
            if string.split(head[g],'_')[0]==string.split(head[gg],'_')[0]:
                pw.append([g,gg])


    filenameout = args.out
    file = open(filenameout,'w')
    out = 'Ind1\tInd2\tk0_hat\tk1_hat\tk2_hat\tpi_HAT\tnbSNP\n'
    file.write(out)

    print 'starting pairwise IBD computations\n'
    print out[:-1]

    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    numConsumers = multiprocessing.cpu_count() * 2
    consumers = [ RelatednessCalculator(args, nbSNPs, datum, IBS_all, AF, tasks, results)
                  for i in xrange(numConsumers) ]

    for w in consumers:
        w.start()

    numJobs = 0
    #iterate through each pairwise comparison
    for g in range(len(pw)):
        ## if we're testing
        if args.testing and (pw[g][0] > 15 or pw[g][1] > 16):
            continue

        tasks.put((pw[g][0], pw[g][1]))
        numJobs += 1
        pass
    for i in xrange(numConsumers):
        tasks.put(None)

    tasks.join()
    estRes = []
    while numJobs:
        estRes.append(results.get())
        print(estRes[-1])
        numJobs -= 1

    estRes = sorted(estRes)
    for res in estRes:
        outL = [ str(res[0]+1), str(res[1]+1) ]
        outL += [ str(round(res[2],2)),
                  str(round(res[3],2)),
                  str(round(res[4],2)),
                  str(round(0.5 * res[3] + res[4], 3)) ]
        outL += [ str(res[5]) ]
        out = '\t'.join(outL)
        print(out)
        file.write(out+'\n')
    file.close()


#        # matrix for all possible pairs of genotypes for every SNP.
#        # This will eventually store genotype likelihoods based on
#        # the calling likelihood (for example based on read depth)
#        PIBS=np.zeros((nbSNPs,9),float)
#        # matrix all possible pairs of genotypes for every SNP.
#        # This will eventually store genotype likelihoods based on
#        # the calling likelihood (for example based on read depth)
#        PIBSg=np.zeros((nbSNPs,9),float)
#
#        # A matrix to denote SNPs we may want to mask for two reasons(see below)
#        mask_mat=np.zeros(nbSNPs,int)
#
#        # iterate through each SNP
#        for gg in range(nbSNPs):
#            # mask SNPs where allele frequency is fixed
#            if AF[gg]==1.0:
#                mask_mat[gg]=1
#            if AF[gg]==0.0:
#                mask_mat[gg]=1
#
#            # pulls out info field for first individual from VCF,
#            # pulls out the three precomputed genotype likelihoods
#            l1=string.split(string.split(datum[gg][pw[g][0]],':')[3],',')
#            l2=string.split(string.split(datum[gg][pw[g][1]],':')[3],',')
#
#            # converts likelihoods from strings to floating point numbers
#            for ggg in range(len(l1)):
#                l1[ggg]=float(l1[ggg])
#                l2[ggg]=float(l2[ggg])
#                # if one of those likelihoods comes out as
#                # negative (some anomaly) we mask that SNPs
#                if l1[ggg]==-9.0:
#                    mask_mat[gg]=1
#                if l2[ggg]==-9.0:
#                    mask_mat[gg]=1
#
#            # We now calculate the probability of observing all possible two
#            # genotype combination by multiplying their likelihoods together
#            PPQQ=(l1[0]*l2[2])
#            QQPP=(l1[2]*l2[0])
#
#            PPPQ=(l1[0]*l2[1])
#            PQPP=(l1[1]*l2[0])
#            PQQQ=(l1[1]*l2[2])
#            QQPQ=(l1[2]*l2[1])
#
#            PPPP=(l1[0]*l2[0])
#            PQPQ=(l1[1]*l2[1])
#            QQQQ=(l1[2]*l2[2])
#
#            # Store these pairwise likelihoods in an array
#            PIBS[gg]=[PPQQ,QQPP,PPPQ,PQPP,PQQQ,QQPQ,PPPP,PQPQ,QQQQ]  #AATT,TTAA,AAAT,ATAA,ATTT,TTAT,AAAA,ATAT,TTTT
#
#
#            # The same as above, but this time we are also weighting the
#            # two genptype likelihoods by the probability of observing
#            # that genotype combination in the population
#            # Include probability of genotype
#            PPQQg=l1[0]*l2[2]*(p**2)*(q**2)
#            QQPPg=l1[2]*l2[0]*(q**2)*(p**2)
#
#            PPPQg=l1[0]*l2[1]*(p**2)*(2*p*q)
#            PQPPg=l1[1]*l2[0]*(2*p*q)*(p**2)
#            PQQQg=l1[1]*l2[2]*(2*p*q)*(q**2)
#            QQPQg=l1[2]*l2[1]*(q**2)*(2*p*q)
#
#            PPPPg=l1[0]*l2[0]*(p**4)
#            PQPQg=l1[1]*l2[1]*(2*p*q)*(2*p*q)
#            QQQQg=l1[2]*l2[2]*(q**2)*(q**2)
#
#            PIBSg[gg]=[PPQQg,QQPPg,PPPQg,PQPPg,PQQQg,QQPQg,PPPPg,PQPQg,QQQQg]  #AATT,TTAA,AAAT,ATAA,ATTT,TTAT,AAAA,ATAT,TTTT
#
#        ###We transpose the matrices
#        PIBSt=PIBS.transpose()
#        PIBStg=PIBSg.transpose()
#
#        ###identify the most likely genotype combination
#        BestGT=PIBS.argmax(axis=1)
#
#        BestIBS=np.zeros((nbSNPs,3),float)
#
#        ###For each SNP, given the best genotype combination, pulls out the appropriate P(IBS|IBD) for all three IBS possibilities
#        for gg in range(nbSNPs):
#            BestIBS[gg]=IBS_all[BestGT[gg]][gg]
#
#
#        def EMOpt(k, IBS, mask_mat, GL = None, verbose = False):
#            """Determine the k coefficients via the EM algorithm
#
#                k --- An initial guess of the k coefficient
#                IBS --- The matrix that gives the identity by site probs
#                mask_mat --- The matrix that masks SNPs that should not be
#                             considered
#                Keyword arguments:
#                GL --- The matrix of genotype likelihoods if we are
#                       considering all (rather than just the best)
#                       genotypes.
#                verbose --- If True, print out periodic progress of
#                            the EM algorithm, if False, be silent.
#            """
#
#            useAllGenotypes = GL is not None
#            # number of rows and columns depends on
#            # whether IBS is a 3D or 2D matrix
#            if useAllGenotypes:
#                nr = IBS[0].shape[0]
#                nc = IBS[0].shape[1]
#            else:
#                nr = IBS.shape[0]
#                nc = IBS.shape[1]
#
#            # This matrix will store the probabilities of
#            # identity-by-descent for each SNP
#            PIBD = np.zeros((nr, nc))
#            k3 = k
#            # The current norm of the difference between
#            # subsequent parameter estimates
#            d = 100
#            # The current iteration number
#            it = 0
#
#            # Arbitrary threshold (should be a parameter?)
#            while d > 1e-4:
#                # Num SNPs x 3
#                X = np.zeros((nr, nc))
#                X = IBS * k3
#
#                # If we use *all* genotypes
#                if useAllGenotypes:
#                    for i in xrange(9):
#                        # Pr(X_m=j | S_{i,m}) * k_{j}
#                        # Probability of IDB = k given IBS
#                        X[i, :, 0] *= GL[i,:]
#                        X[i, :, 1] *= GL[i,:]
#                        X[i, :, 2] *= GL[i,:]
#                    X = X.sum(axis=0)
#
#                # Normalize the X matrix to contain probabilities
#                XS = ma.masked_array(X.sum(axis=1), mask=mask_mat)
#                X[:,0] /= XS
#                X[:,1] /= XS
#                X[:,2] /= XS
#                # Copy over into the PIBD matrix
#                PIBD = X
#
#                # The new parameter estimate
#                kp = np.zeros(3)
#                # The sum of probabilities at each site
#                kp[0] += np.sum(ma.masked_array(PIBD[:,0], mask=mask_mat))
#                kp[1] += np.sum(ma.masked_array(PIBD[:,1], mask=mask_mat))
#                kp[2] += np.sum(ma.masked_array(PIBD[:,2], mask=mask_mat))
#                # Normalized
#                kp /= kp.sum()
#                # Compute the difference between successive
#                # estimates to assess convergence
#                d = np.linalg.norm(np.abs(k3-kp))
#                k3 = kp
#                it += 1
#                if verbose and it % 100 == 0:
#                    print(k3)
#                    print("diff = {}".format(d))
#            # Return what we did all that work to calculate
#            return k3
#
#        kest = np.zeros(3)
#        if args.genotype == "best":
#            k3 = np.random.rand(3)
#            k3 /= k3.sum()
#            kest = EMOpt(k3, BestIBS, mask_mat)
#        elif args.genotype == "all":
#            k3 = np.random.rand(3)
#            k3 /= k3.sum()
#            kest = EMOpt(k3, IBS_all, mask_mat, PIBSt)
#
#        tasks.join()
#        while numJobs:
#            result = results.get()
#            numJobs -= 1
#        outL = [str(pw[g][0]+1), str(pw[g][1]+1)]
#        outL += [ str(round(kest[0],2)),
#                  str(round(kest[1],2)),
#                  str(round(kest[2],2)),
#                  str(round(0.5 * kest[1] + kest[2], 3)) ]
#        outL += [str(len(mask_mat)-np.sum(mask_mat))]
#        out = '\t'.join(outL)
#
#        print(out)
#        file.write(out+'\n')
#
#    file.close()



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", help="VCF input file", required=True)
    parser.add_argument("--out", help="File where parameter estimates are written", required=True)
    parser.add_argument("-t", "--testing", help="Only compute relationship for first 16 individuals",
                                  action="store_true")
    parser.add_argument("--genotype", nargs="?",
                                  help="Type of inference to perform",
                                  choices = ["all", "best"],
                                  const="all", default="all")
    args = parser.parse_args()

    runMLKin(args)
