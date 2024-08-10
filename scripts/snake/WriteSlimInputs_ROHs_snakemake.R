make.slim.input <- function(filename.start,
                            rand.seed="1234567890", pop.size=10000,
                            burnin.gens=50000, samp.size=100, demog="onepop",
                            mating.sys="outc",
                            mate.perc, rep,
                            neutral.only = FALSE,
                            mu){
    
    pop.size <- as.numeric(pop.size)
    rep <- as.numeric(rep)
    samp.size <- as.numeric(samp.size)
    samp.size <- min(samp.size,pop.size)
    mate.perc <- as.numeric(mate.perc)
    mate.perc <- mate.perc/100

    neutralString <- paste(neutral.only)

    file.slim <- sprintf('infile.txt')
    
    
    cat('', file=file.slim, append = FALSE)
                                        # mycat = function(string) {
                                        #   cat(string, file=file.slim, append=TRUE)
                                        # }
                                        # 
    cat('initialize() { \n', file=file.slim, append=TRUE)
    
    cat(sprintf('   setSeed(%s);\n', rand.seed), file=file.slim, append=TRUE)
    cat(sprintf('   initializeSLiMOptions(preventIncidentalSelfing=T);\n'), file=file.slim, append=TRUE)
    cat(sprintf('   initializeMutationRate(%s);\n',paste(mu)), file=file.slim, append=TRUE)
    
    cat(sprintf('   initializeMutationType("m1", 0.5, "f", 0.0);			// neutral\n'), file=file.slim, append=TRUE)
    cat(sprintf('   m1.convertToSubstitution = T;\n'), file=file.slim, append=TRUE)
    cat(sprintf('   m1.mutationStackGroup = -1;\n'), file=file.slim, append=TRUE)
    cat(sprintf('   m1.mutationStackPolicy = "l"; \n\n'), file=file.slim, append=TRUE)
    
    cat(sprintf('   initializeMutationType("m2", 0.3,"e", -0.001);		// delet, exponential effect\n'), file=file.slim, append=TRUE)
    cat(sprintf('   m2.convertToSubstitution = T;\n'), file=file.slim, append=TRUE)
    cat(sprintf('   m2.mutationStackGroup = -1;\n'), file=file.slim, append=TRUE)
    cat(sprintf('   m2.mutationStackPolicy = "l";\n\n'), file=file.slim, append=TRUE)

    cat(sprintf('   initializeMutationType("m3", 0.02,"f", -1);		// delet, lethal effect\n'), file=file.slim, append=TRUE)
    cat(sprintf('   m3.convertToSubstitution = T;\n'), file=file.slim, append=TRUE)
    cat(sprintf('   m3.mutationStackGroup = -1;\n'), file=file.slim, append=TRUE)
    cat(sprintf('   m3.mutationStackPolicy = "l";\n\n'), file=file.slim, append=TRUE)

    cat(sprintf('   initializeMutationType("m4", 0.3,"e", 0.01);		// ben, exponential effect\n'), file=file.slim, append=TRUE)
    cat(sprintf('   m4.convertToSubstitution = T;\n'), file=file.slim, append=TRUE)
    cat(sprintf('   m4.mutationStackGroup = -1;\n'), file=file.slim, append=TRUE)
    cat(sprintf('   m4.mutationStackPolicy = "l";\n\n'), file=file.slim, append=TRUE)
    
    cat(sprintf('   // set mutation proportions\n'), file=file.slim, append=TRUE)
    if (!neutral.only) { # default case with selection
        cat(sprintf('   initializeGenomicElementType("g1", c(m1,m2,m3,m4), c(0.25, 0.649, 0.1, 0.001));\n\n'),
            file=file.slim, append=TRUE)
    } else { # neutral only 
        cat(sprintf('   initializeGenomicElementType("g1", c(m1,m2,m3,m4), c(1, 0, 0, 0));\n\n'),
            file=file.slim, append=TRUE)
    }
    
    cat(sprintf('   // make total genome 10Mbp, with one recomb rate\n'), file=file.slim, append=TRUE)
    cat(sprintf('   initializeGenomicElement(g1, 0, 9999999);       // one chrom only\n'), file=file.slim, append=TRUE)
    cat(sprintf('   initializeRecombinationRate(1e-8, 9999999);\n'), file=file.slim, append=TRUE)
    cat(sprintf('\n'), file=file.slim, append=TRUE)

    cat(sprintf('}\n\n'), file=file.slim, append=TRUE)
    
    
                                        # NOW THE DEMOGRAPHIC SCENARIO OPTIONS
    cat(sprintf('   // set DEMOGRAPHY\n'), file=file.slim, append=TRUE)
    end.sample.point <- burnin.gens + 5000  

    if(demog=="onepop"){
        cat(sprintf('1 {\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sim.addSubpop("p1", %i);\n', pop.size), file=file.slim, append=TRUE)
        if(mating.sys=="outc"){
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        if(mating.sys=="self"){
            cat(sprintf('       p1.setSelfingRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        if(mating.sys=="asex"){
            cat(sprintf('       p1.setCloningRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        cat(sprintf('%i late() { \n', end.sample.point), file=file.slim, append=TRUE)
        cat(sprintf('       subsampDiploids = sim.subpopulations.individuals;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sampledIndividuals = sample(subsampDiploids, %i);\n', samp.size), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf = sampledIndividuals.genomes;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf.outputVCF(filePath="out.vcf");\n'), file=file.slim, append = TRUE)   
        cat(sprintf('}\n'), file=file.slim, append=TRUE)
    }
    
    if(demog=="bneckShort"){
        bneck.size <- pop.size/100
        samp.size <- min(bneck.size,samp.size)
        
        cat(sprintf('1 {\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sim.addSubpop("p1", %i);\n', pop.size), file=file.slim, append=TRUE)
        if(mating.sys=="outc"){
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        if(mating.sys=="self"){
            cat(sprintf('       p1.setSelfingRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        if(mating.sys=="asex"){
            cat(sprintf('       p1.setCloningRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        ## make the bottleneck happen from gens 50,000 to 50,500  
        cat(sprintf('%i {\n', end.sample.point-5000), file=file.slim, append=TRUE)
        cat(sprintf('       p1.setSubpopulationSize(%i);\n', bneck.size), file=file.slim, append=TRUE)
        cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        ## sample for the contraction scenario, before going back up to big pop size
        cat(sprintf('%i late() { \n', end.sample.point-4501), file=file.slim, append=TRUE)
        cat(sprintf('       subsampDiploids = sim.subpopulations.individuals;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sampledIndividuals = sample(subsampDiploids, %i);\n', samp.size), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf = sampledIndividuals.genomes;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf.outputVCF(filePath="out_contr.vcf");\n'), file=file.slim, append = TRUE)   
        cat(sprintf('}\n'), file=file.slim, append=TRUE)
        cat(sprintf('%i {\n', end.sample.point-4500), file=file.slim, append=TRUE)
        cat(sprintf('       p1.setSubpopulationSize(%i);\n', pop.size), file=file.slim, append=TRUE)
        cat(sprintf('}\n'), file=file.slim, append=TRUE)	  
        cat(sprintf('%i late() { \n', end.sample.point), file=file.slim, append=TRUE)
        cat(sprintf('       subsampDiploids = sim.subpopulations.individuals;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sampledIndividuals = sample(subsampDiploids, %i);\n', samp.size), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf = sampledIndividuals.genomes;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf.outputVCF(filePath="out.vcf");\n'), file=file.slim, append = TRUE)   
        cat(sprintf('}\n'), file=file.slim, append=TRUE)
    }
    
    if(demog=="bneckLong"){
  	bneck.size <- pop.size/100
        samp.size <- min(bneck.size,samp.size)

        cat(sprintf('1 {\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sim.addSubpop("p1", %i);\n', pop.size), file=file.slim, append=TRUE)
        if(mating.sys=="outc"){
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        if(mating.sys=="self"){
            cat(sprintf('       p1.setSelfingRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        if(mating.sys=="asex"){
            cat(sprintf('       p1.setCloningRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        ## make the bottleneck happen from gens 50,000 to 52,500  
        cat(sprintf('%i {\n', end.sample.point-5000), file=file.slim, append=TRUE)
        cat(sprintf('       p1.setSubpopulationSize(%i);\n', bneck.size), file=file.slim, append=TRUE)
        cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        ## sample for the contraction scenario, before going back up to big pop size
        cat(sprintf('%i late() { \n', end.sample.point-2501), file=file.slim, append=TRUE)
        cat(sprintf('       subsampDiploids = sim.subpopulations.individuals;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sampledIndividuals = sample(subsampDiploids, %i);\n', samp.size), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf = sampledIndividuals.genomes;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf.outputVCF(filePath="out_contr.vcf");\n'), file=file.slim, append = TRUE)   
        cat(sprintf('}\n'), file=file.slim, append=TRUE)
        cat(sprintf('%i {\n', end.sample.point-2500), file=file.slim, append=TRUE)
        cat(sprintf('       p1.setSubpopulationSize(%i);\n', pop.size), file=file.slim, append=TRUE)
        cat(sprintf('}\n'), file=file.slim, append=TRUE)	  
        cat(sprintf('%i late() { \n', end.sample.point), file=file.slim, append=TRUE)
        cat(sprintf('       subsampDiploids = sim.subpopulations.individuals;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sampledIndividuals = sample(subsampDiploids, %i);\n', samp.size), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf = sampledIndividuals.genomes;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf.outputVCF(filePath="out.vcf");\n'), file=file.slim, append = TRUE)   
        cat(sprintf('}\n'), file=file.slim, append=TRUE)
    }  
    
    if(demog=="bneckGrad"){
        bneck.size <- pop.size/100
        samp.size <- min(bneck.size,samp.size)
        
        cat(sprintf('   function (integer)linReg(i N1, i N2, i t1, i t2, i tnow) {\n'), file=file.slim, append=TRUE)
        cat(sprintf('       newSize = asInteger(round((N2-N1)/(t2-t1)*tnow+((t2*N1-t1*N2)/(t2-t1))));\n'), file=file.slim, append=TRUE)
        cat(sprintf('       return(newSize);\n\n'), file=file.slim, append=TRUE)
        cat(sprintf('}\n'), file=file.slim, append=TRUE)
        cat(sprintf('1 {\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sim.addSubpop("p1", %i);\n', pop.size), file=file.slim, append=TRUE)
        if(mating.sys=="outc"){
            cat(sprintf('}\n'), file=file.slim, append=TRUE)
        }
        if(mating.sys=="self"){
            cat(sprintf('       p1.setSelfingRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('}\n'), file=file.slim, append=TRUE)
        }
        if(mating.sys=="asex"){
            cat(sprintf('       p1.setCloningRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('}\n'), file=file.slim, append=TRUE)
        }
        cat(sprintf('%i:%i { // tds:tdf timepoints \n', end.sample.point-5000,end.sample.point-2500), file=file.slim, append=TRUE)
        cat(sprintf('    newSize = linReg(%i, %i, %i, %i, sim.generation);\n', pop.size, bneck.size, end.sample.point-5000,end.sample.point-2500), file=file.slim, append=TRUE)
        cat(sprintf('    p1.setSubpopulationSize(newSize);\n'), file=file.slim, append=TRUE)
        cat(sprintf('}\n'), file=file.slim, append=TRUE)
        ## sample for the gradual contraction scenario, before going back up to big pop size
        cat(sprintf('%i late() { \n', end.sample.point-4501), file=file.slim, append=TRUE)
        cat(sprintf('       subsampDiploids = sim.subpopulations.individuals;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sampledIndividuals = sample(subsampDiploids, %i);\n', samp.size), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf = sampledIndividuals.genomes;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf.outputVCF(filePath="out_contr.vcf");\n'), file=file.slim, append = TRUE)
        cat(sprintf('}\n'), file=file.slim, append=TRUE)
        cat(sprintf('%i:%i { // tgs:tgf timepoints \n', end.sample.point-2499,end.sample.point), file=file.slim, append=TRUE)
        cat(sprintf('    newSize = linReg(%i, %i, %i, %i, sim.generation);\n', bneck.size, pop.size, end.sample.point-2499,end.sample.point), file=file.slim, append=TRUE)
        cat(sprintf('    p1.setSubpopulationSize(newSize);\n'), file=file.slim, append=TRUE)
        cat(sprintf('}\n'), file=file.slim, append=TRUE)
        cat(sprintf('%i late() { \n', end.sample.point), file=file.slim, append=TRUE)
        cat(sprintf('       subsampDiploids = sim.subpopulations.individuals;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sampledIndividuals = sample(subsampDiploids, %i);\n', samp.size), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf = sampledIndividuals.genomes;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf.outputVCF(filePath="out.vcf");\n'), file=file.slim, append = TRUE)
        cat(sprintf('}\n'), file=file.slim, append=TRUE)
    }

    if(demog=="twoPop"){
        cat(sprintf('1 {\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sim.addSubpop("p1", %i);\n', pop.size), file=file.slim, append=TRUE)
        cat(sprintf('       sim.addSubpop("p2", %i);\n', pop.size), file=file.slim, append=TRUE)	  
        cat(sprintf('       p1.setMigrationRates(p2, 0.0);\n'), file=file.slim, append=TRUE)
        cat(sprintf('       p2.setMigrationRates(p1, 0.0);\n'), file=file.slim, append=TRUE)
                                        # no initial migration between the two pops
        if(mating.sys=="outc"){
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        if(mating.sys=="self"){
            cat(sprintf('       p1.setSelfingRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('       p2.setSelfingRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        if(mating.sys=="asex"){
            cat(sprintf('       p1.setCloningRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('       p2.setCloningRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        ## make migration begin at gen 50,000  
        cat(sprintf('%i {\n', end.sample.point-5000), file=file.slim, append=TRUE)
        cat(sprintf('       p1.setMigrationRates(p2, 0.1);\n'), file=file.slim, append=TRUE)
        cat(sprintf('       p2.setMigrationRates(p1, 0.1);\n'), file=file.slim, append=TRUE)
        cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        ## sample for the early admixture scenario, before migration goes on for too long
        cat(sprintf('%i late() { \n', end.sample.point-4501), file=file.slim, append=TRUE)
        cat(sprintf('       subsampDiploids = sim.subpopulations.individuals;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sampledIndividuals = sample(subsampDiploids, %i);\n', samp.size), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf = sampledIndividuals.genomes;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf.outputVCF(filePath="out_EarlyAdmix.vcf");\n'), file=file.slim, append = TRUE)   
        ## same timepoint, sample only p1
        cat(sprintf('       subsampDiploids = p1.individuals;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sampledIndividuals = sample(subsampDiploids, %i);\n', samp.size), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf = sampledIndividuals.genomes;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf.outputVCF(filePath="out_EarlyAdmixP1.vcf");\n'), file=file.slim, append = TRUE)   
        cat(sprintf('}\n'), file=file.slim, append=TRUE)  
        cat(sprintf('%i late() { \n', end.sample.point), file=file.slim, append=TRUE)
        cat(sprintf('       subsampDiploids = sim.subpopulations.individuals;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sampledIndividuals = sample(subsampDiploids, %i);\n', samp.size), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf = sampledIndividuals.genomes;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf.outputVCF(filePath="out.vcf");\n'), file=file.slim, append = TRUE)   
        ## same timepoint, sample only p1
        cat(sprintf('       subsampDiploids = p1.individuals;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sampledIndividuals = sample(subsampDiploids, %i);\n', samp.size), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf = sampledIndividuals.genomes;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf.outputVCF(filePath="out_P1.vcf");\n'), file=file.slim, append = TRUE)   
        cat(sprintf('}\n'), file=file.slim, append=TRUE)
    } 

    if(demog=="twoPopIsol"){
        cat(sprintf('1 {\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sim.addSubpop("p1", %i);\n', pop.size), file=file.slim, append=TRUE)
        cat(sprintf('       sim.addSubpop("p2", %i);\n', pop.size), file=file.slim, append=TRUE)	  
        cat(sprintf('       p1.setMigrationRates(p2, 0.0);\n'), file=file.slim, append=TRUE)
        cat(sprintf('       p2.setMigrationRates(p1, 0.0);\n'), file=file.slim, append=TRUE)
                                        # no initial migration between the two pops
        if(mating.sys=="outc"){
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        if(mating.sys=="self"){
            cat(sprintf('       p1.setSelfingRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('       p2.setSelfingRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        if(mating.sys=="asex"){
            cat(sprintf('       p1.setCloningRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('       p2.setCloningRate(%f);\n', mate.perc), file=file.slim, append=TRUE)
            cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        }
        ## make migration begin at gen 50,000  
        cat(sprintf('%i {\n', end.sample.point-5000), file=file.slim, append=TRUE)
        cat(sprintf('       p1.setMigrationRates(p2, 0.1);\n'), file=file.slim, append=TRUE)
        cat(sprintf('       p2.setMigrationRates(p1, 0.1);\n'), file=file.slim, append=TRUE)
        cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        ## make migration stop at gen 50,500  
        cat(sprintf('%i {\n', end.sample.point-4500), file=file.slim, append=TRUE)
        cat(sprintf('       p1.setMigrationRates(p2, 0.0);\n'), file=file.slim, append=TRUE)
        cat(sprintf('       p2.setMigrationRates(p1, 0.0);\n'), file=file.slim, append=TRUE)
        cat(sprintf('}\n'), file=file.slim, append=TRUE) 
        cat(sprintf('%i late() { \n', end.sample.point), file=file.slim, append=TRUE)
        cat(sprintf('       subsampDiploids = sim.subpopulations.individuals;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sampledIndividuals = sample(subsampDiploids, %i);\n', samp.size), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf = sampledIndividuals.genomes;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf.outputVCF(filePath="out.vcf");\n'), file=file.slim, append = TRUE)   
        ## same timepoint, sample only p1
        cat(sprintf('       subsampDiploids = p1.individuals;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       sampledIndividuals = sample(subsampDiploids, %i);\n', samp.size), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf = sampledIndividuals.genomes;\n'), file=file.slim, append=TRUE)
        cat(sprintf('       outvcf.outputVCF(filePath="out_P1.vcf");\n'), file=file.slim, append = TRUE)   
        cat(sprintf('}\n'), file=file.slim, append=TRUE)
    }   

    return(invisible())
}



#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

args<-commandArgs(trailingOnly=TRUE)


# args[1]: -- input directory
# args[2]: -- replicate
# args[3]: -- population size
# args[4]: -- demography
# args[5]: -- mating system
# args[6]: -- % mating
# args[7]: -- neutral.only


# MAKE INPUT FILE FOR THIS SIM



setwd(args[1])

# always do 100 replicates # could speed it up by doing one rand num in R then generate these rand nums
set.seed(as.numeric(args[9]))
rand.seeds <- sample.int(1e9,as.numeric(args[10]))

replicate <- as.numeric(args[2])
pop.size <- as.numeric(args[3])
demog <- args[4]
mating.sys <-args[5]
mate.perc <-args[6]
isNeutral <- as.logical(args[7]) # FALSE: default with selection, TRUE: only neutral mutations
mu <- args[8]                    # mutation rate


inds.to.sample <- 100


make.slim.input(filename.start=args[1], rand.seed=rand.seeds[replicate],
                pop.size=pop.size, burnin.gens=50000, samp.size=inds.to.sample,
                demog=demog, mating.sys=mating.sys,mate.perc=mate.perc,
                rep=replicate, neutral.only=isNeutral, mu=mu)

