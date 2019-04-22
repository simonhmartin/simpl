#load simulation functions
source("simpl.R")

################################################################################
#######################  set simulation parameters  ############################
################################################################################

# number of loci to simulate
n_loci <- 101

# each locus needs a map position (in centimorgan).
# These can either be specified directly:

#map_pos <- 1:n_loci

# or computed based on inter-locus recombination distances (should total n_loci - 1):
rec = c(rep(0.8,20),
        rep(1.6,20),
        rep(0.2,20),
        rep(1.6,20),
        rep(0.8,20))

map_pos = c(0,cumsum(rec))

# number of individuals in the population
n_ind <- 4000

# initialise the population
# population is always an array with dimensions n_loci, 2 (diploid) and n_individuals
# all loci are given the ancestral (local) allele (1)
pop = array(1, dim = c(n_loci, 2, n_ind))


# now designate some individuals has immigrants. 10% of the population in this case
immigrants = 1:(n_ind*0.1)

# give the immigrants a derived (foreign) allele (2)
pop[,,immigrants] <- 2

# now give the immigrants a scattering of deleterious alleles (3)
pop[seq(1,n_loci,5),,immigrants] <- 3

# add a benefficial allele to one locus in the immigrants
pop[75,,immigrants] <- 4


# we have designated four different allele types
# give the alleles selection coefficients relative to the most beneficial allele (which gets 0)
# local allele (1) has the same as the most common introgressed allele (2)
# The deleterious alleles get a larger cost and the beneficial allele gets no cost
s <- c(0.6, 0.6, 0.8, 0)

################################################################################
#############################  run simulation  #################################
################################################################################

# before running the simulation, we decide which generations to store
store_gens = c(10,100,300)

#we can also define the colours we want in the progress plots
cols = c("gray90","#FFC5C5","#FF9999", "#FF0000")

#run the simulation and store population snapshots
sims <- run.simulation(300, pop, map_pos, s, store = store_gens, cores=5,
                       progress_plot_frequency=10, progress_plot_N=30, col=cols)


################################################################################
#############################  summary plot  ###################################
################################################################################

#make a summary plot showing a few individuals along with the frequency of the introgressed allels (2, 3, and 4)
par(mfrow=c(length(sims)*2,1))

for (i in 1:length(sims)){
    par(mar = c(0,4,3,1), xpd=NA)
    plot.individuals(sims[[i]][,,1:10], col = cols, main = paste("Generation", store_gens[i]))

    par(mar = c(4,4,1,1), xpd=NA)
    plot.frequency(sims[[i]], alleles=c(2,3,4), col=cols[2])
    }
    
