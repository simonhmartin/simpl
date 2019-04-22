library(parallel)

# function that takes a diploid chromosome and makes a haploid recombinant version
# genotypes has two dimensions first has length n_loci, second is has length 2 (for each chromatid)
# positions is position of each locus in centimorgans
new.gamete <- function(genotypes, positions, units="centimorgan"){
    if (units == "centimorgan") scale <- 0.01
    else scale <- 1
    n_loci = length(positions)
    crossProbs <- (positions[-1] - positions[-n_loci]) * scale #probability of a crossover between each pair of loci
    crossOvers <- rpois(length(crossProbs), crossProbs) # sample a number of crossovers for each interval
    switches = which(crossOvers %% 2 == 1) # switchs occur where there ar odd numbers of crossovers
    chroms = sample(c(1,2),2) # randomise chromosomes so we can start with either
    if (length(switches) > 0){    
        donors <- rep_len(chroms, length(switches)+1) #the donor chromosome alternates at each switch
        seg_lengths <- c(switches,n_loci) - c(0,switches) # number of loci between switches
        donors <- rep(donors, seg_lengths) #vector giving the donoor for each locus
        gam <- sapply(1:n_loci, function(x) genotypes[x,donors[x]]) #extract new sequence
        }
    else gam <- genotypes[,chroms[1]]
    gam
    }

# function to make a new individual.
# Parents are picked from the population with probability proportional to relative fitness
new.ind <- function(population, positions, fitnesses = NULL){
    N <- dim(population)[3] #number of parents we're choosing from
    parents <- sample(1:N, 2, replace=TRUE, prob=fitnesses)
    cbind(new.gamete(population[,,parents[1]],positions),new.gamete(population[,,parents[2]],positions))
    }

# function to make a new generation of the same size as population.
# relative fitnesses are defined as the sum of 1-s for each allele
# this version of the fiunction uses a single core, without the parallel package
new.generation <- function(pop, positions, s = NULL){
    N <- dim(pop)[3]
    fitnesses <- sapply(1:N, function(x) ifelse(is.null(s)==TRUE, 1, sum(1-s[pop[,,x]])))
    sapply(1:N, function(x) new.ind(pop,positions,fitnesses), simplify="array")
    }

#copy of function above, but using the package parallel to use multiple cores so that new individuals can be made simultaneously
new.generation.parallel <- function(pop, positions, s = NULL, cores){
    N <- dim(pop)[3]
    fitnesses <- sapply(1:N, function(x) ifelse(is.null(s)==TRUE, 1, sum(1-s[pop[,,x]])))
    list.of.matrices.to.array(mclapply(1:N, function(x) new.ind(pop,positions,fitnesses), mc.cores=cores))
    }

#just a function to convert the list of matrices from the mcapply command into a 3D array
list.of.matrices.to.array <- function(lst){
    a <- array(NA, dim = c(nrow(lst[[1]]),ncol(lst[[1]]),length(lst)),
               dimnames = list(rownames(lst[[1]]), colnames(lst[[1]]), names(lst)))
    for(k in 1:length(lst))  a[,,k] <- lst[[k]]
    a
    }

# this is the wrapper function for new.generation.
# it runs it multiple times and can store the output whenever desired
# it can also make progress plots where it shows a subset of individuals' chromosomes
run.simulation <- function(n_gen, population, map_pos, s = NULL, store = NULL, cores=1,
                           progress_plot_frequency=0, progress_plot_N=20, col=NULL){
    if (is.null(store)==TRUE) store <- n_gen
    stored <- list(length=length(store))
    
    store_idx <- 1
    
    for (n in 1:n_gen){
        print(n)
        
        if (cores >= 1) population <- new.generation.parallel(population, map_pos, s=s, cores=cores)
        else population <- new.generation(population, map_pos, s=s)
        
        if (n %in% store){
            stored[[store_idx]] <- population
            store_idx = store_idx+1
            }
        
        if (progress_plot_frequency >= 1 & n %% progress_plot_frequency == 0) {
            plot.individuals(population[,,1:progress_plot_N], col = cols, main = paste("Generation", n))
            }
        }
    stored
    }

#function for plotting indivial chreonmosomes for a group of diploid individuals with loci coloured
plot.individuals <- function(individuals, col = c("red","blue"), chrom_width=0.2, main=NULL){
    N = dim(individuals)[3]
    plot(0, cex = 0, bty = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n",
            ylim = c(1,N+chrom_width*3), xlim = c(1,dim(individuals)[1]), main = main)
    for (i in 1:N) {
        draw.chrom(individuals[,1,i], origin = c(0,i), chrom_width = chrom_width, loc_len=1, col=col)
        draw.chrom(individuals[,2,i], origin = c(0,i+chrom_width*2), chrom_width = chrom_width, loc_len=1, col=col)
        }
    }

#function for drawing individual chromosome using rectangles
draw.chrom <- function(chrom, col = c("red","blue"), origin = c(0,0), chrom_width = 0.2, loc_len = 1){
    nloci = length(chrom)
    loc_ends = cumsum(rep(loc_len,nloci))
    loc_starts = loc_ends - loc_len
    rect(origin[1] + loc_starts, rep(origin[2],nloci), origin[1] + loc_ends, rep(origin[2]+chrom_width,nloci), border = NA, col = col[chrom])
    }

#function to plot the frequency of one or more alleles as a filled polygon
plot.frequency <- function(population, alleles, col, border=NA){
    freq <- apply(population, 1, function(x) mean(x %in% alleles))
    pos = 1:length(freq)
    plot(0, cex=0, bty = "l", ylim = c(0,1), xlim = c(1,length(freq)), ylab = "Frequency", xlab = "Position")
    polygon(c(0, interleave(pos-1,pos), tail(pos,1)), c(0,rep(freq,each=2),0), col = col, border = border)
    }

#function to plot the recombination rate between loci as a line
plot.rec <- function(rec,chrom_pos=NULL, lty=par("lty"), add=FALSE){
    if(is.null(chrom_pos) == TRUE) chrom_pos = 1:(length(rec)+1)
    if (add==FALSE) plot(interleave(chrom_pos-1,chrom_pos), c(rec[1], rep(rec,each=2), tail(rec,1)), lty = lty, type="l")
    else lines(interleave(chrom_pos-1,chrom_pos), c(rec[1], rep(rec,each=2), tail(rec,1)), lty = lty)
    }

#just a generic data manipulation function
interleave <- function(x1,x2){
    output <- vector(length= length(x1) + length(x2))
    output[seq(1,length(output),2)] <- x1
    output[seq(2,length(output),2)] <- x2
    output
    }

