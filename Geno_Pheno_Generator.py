IAF = .5
var = 0
da = 0


import math
import random
import numpy as np
import scipy.stats as ss
from scipy.stats import norm

def numInd(var, da, IAF):
    """Generates genotype and phenotype data based on QTL Variance, dominance additivity, and increaser allele frequency"""
    
    #beginning calculations
    NAF = 1 - IAF
    pq2 = (2*IAF*NAF) 
    HWMin = min(IAF**2, pq2, NAF**2)
    numInd = max((10/HWMin), 1000)

    #Additive Term
    at = math.sqrt(var/(2*IAF*(1-IAF)*(1+da*(1-2*IAF))**2+4*IAF**2*(1-IAF)**2*da**2))
    #Dominance Term
    dt = at*da
    #Mean
    mean = -(IAF**2*at+2*IAF*(1-IAF)*dt-((1-IAF)**2)*at)

    #Mean + Additive Term - genotype 2
    mat = mean + at
    #Mean + Dominance Term - genotype 1
    mdt = mean + dt
    #Mean - Additive Term - genotype 0
    negmat = mean - at

    #Residual Variance
    rVar = 1-var
    #Residual Standard Deviation
    rsd = math.sqrt(rVar)

    #### Assignment of Genotype and Phenotype Data ####
    #a number is randomly generated and is sorted into a genotype based on IAF
    #a phenotype is then randomly generated based on its genotypic mean and the residual standard deviation
    randGeno = random.random()
    if (randGeno <= (NAF**2)):
        geno = 0
        pheno = ss.norm(negmat,rsd).rvs(1)
    #end of geno0 assignment

    if  (randGeno >= (1-IAF**2)) :
        geno = 2
        pheno = ss.norm(negmat,rsd).rvs(1)
    #end of geno2 assignment

    if (randGeno > NAF**2 and randGeno < (1-IAF**2)): 
        geno = 1
        pheno = ss.norm(mdt,rsd).rvs(1)
    #end of geno1 assignment

    #returns a tuple of the genotype and phenotype
    return_mat = np.array([geno, pheno])
    return(return_mat)

    
(geno, pheno) = numInd(0,0,.5)
