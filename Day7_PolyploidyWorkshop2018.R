#title: "Day7-PolyploidyWorkshop2018"
#author: "Patrick Monnahan"
#date: "11/23/2018"

require(assertthat)
require(ggplot2)
require(dplyr)
require(tidyr)
require(PopGenome)
require(stringr)
require(data.table)


## Useful Functions

#Fitness and allele frequency trajectories

# Generate HW genotype frequencies for dips and tets
genoFreqs = function(freqs){
  AA = freqs ^ 2
  Aa = 2 * freqs * (1 - freqs)
  aa = (1 - freqs) ^ 2
  AAAA = freqs ^ 4
  AAAa = 4 * (freqs ^ 3) * (1 - freqs)
  AAaa = 6 * (freqs ^ 2) * (1 - freqs) ^ 2
  Aaaa = 4 * freqs * (1 - freqs) ^ 3
  aaaa = (1 - freqs) ^ 4
  dfreq = data.frame("p" = freqs, "AA" = AA, "Aa" = Aa, "aa" = aa, "Ploidy" = rep(2, length(freqs)))
  tfreq = data.frame("p" = freqs, "AAAA" = AAAA, "AAAa" = AAAa, "AAaa" = AAaa, "Aaaa" = Aaaa, "aaaa" = aaaa, "Ploidy" = rep(4, length(freqs)))
  DF = rbind(melt(dfreq, id.var = c("p", "Ploidy")), melt(tfreq, id.var = c("p", "Ploidy")))
  return(DF)
}

# Generate allele frequency trajectory based on deterministic recursion equations   
dipTraj = function(s, h, start_freq, end_freq, max_gens){
  df = data.frame("s" = s, "h" = h, "gen" = 0, "freq" = start_freq, "dp1" = 0, "dp2" = 0, "w.bar" = 0, "var.w" = 0, "w.bar.p" = 0, "h1" = NA, "h3" = NA, "ploidy" = as.factor(2))
  p = start_freq
  gen = 0
  fits = c(1 + s, 1 + (s * h), 1) / (1 + s)
  while (p < end_freq & gen < max_gens){
    q = 1 - p
    Gfreqs = c(p ^ 2, 2 * p * q, q ^ 2)
    num = (Gfreqs[1] * fits[1]) + (p * q * fits[2])
    w.bar.p = (p * fits[1]) + (q * fits[2])
    w.bar = sum(Gfreqs * fits)
    p_prime = num / w.bar
    dp1 = p_prime - p
    dp2 = (p * (w.bar.p - w.bar)) / w.bar
    p = p_prime
    var.w = sum((Gfreqs * (fits - w.bar) ^ 2) / length(Gfreqs))
    gen = gen + 1
    
    df = rbind(df, c(s, h, gen, p, dp1, dp2, w.bar, var.w, w.bar.p, NA, NA, 2))
  }
  return(df[-1,])
}
tetTraj = function(s, h1, h2, h3, start_freq, end_freq, max_gens){
  assert_that(h3 >= h1)
  df = data.frame("s" = s, "h" = h2, "gen" = 0, "freq" = start_freq, "dp1" = 0, "dp2" = 0, "w.bar" = 0, "var.w" = 0, "w.bar.p" = 0, "h1" = h1, "h3" = h3, "ploidy" = as.factor(4))
  p = start_freq
  gen = 0
  freq = c(p)
  fits = c(1 + s, 1 + (s * h3), 1 + (s * h2), 1 + (s * h1), 1) / (1 + s)
  while (p < end_freq & gen < max_gens){
    q = 1 - p
    Gfreqs = c(p ^ 4, 4 * p^3 * q, 6 * p^2 * q^2, 4 * p * q^3, q ^ 4)
    num = (Gfreqs[1] * fits[1]) + (3 * p ^ 3 * q * fits[2]) + (3 * p ^ 2 * q ^ 2 * fits[3]) + (p * q ^ 3 * fits[4])
    w.bar.p = ((p ^ 3) * fits[1]) + (3 * (p ^ 2) * q * fits[2]) + (3 * p * (q ^ 2) * fits[3]) + ((q ^ 3) * fits[4])
    w.bar = sum(Gfreqs * fits)
    p_prime = num / w.bar
    dp1 = p_prime - p
    dp2 = (p * (w.bar.p - w.bar)) / w.bar
    p = p_prime
    var.w = sum((Gfreqs * (fits - w.bar) ^ 2) / length(Gfreqs))
    gen = gen + 1
    df = rbind(df, c(s, h2, gen, p, dp1, dp2, w.bar, var.w, w.bar.p, h1, h3, 4))
  }
  return(df[-1,])
}

# Function to loop over all possible combinations of s and h.
simTraj = function(s = c(0.1, 0.01, 0.001), h1 = c(0.25, 0, 1), h2 = c(0.5, 0, 1), h3 = c(0.75, 0, 1), start_freq =  0.05, end_freq = 1.0, maxGen = 9999){
  traj = data.frame()
  for (i in 1:length(s)){
    for (j in 1:length(h1)){
      print(paste("Starting; s =", s[i], ", h =", h1[j]))
      jj = dipTraj(s[i], h2[j], start_freq, end_freq, maxGen)
      kk = tetTraj(s[i], h1[j], h2[j], h3[j], start_freq, end_freq, maxGen)
      traj = rbind(traj, jj)
      traj = rbind(traj, kk)
    }
  }
  return(traj)
}

# Generate fitness stats as a function of allele frequency
getFits = function(freqs, s, h1, h, h3){
  df = data.frame("s" = 0, "h" = 0, "freq" = 0, "w.bar" = 0, "var.w" = 0, "w.bar.p" = 0, "h1" = NA, "h3" = NA, "ploidy" = 2)
  assert_that(length(h1) == length(h))
  assert_that(length(h) == length(h3))
  for (j in 1:length(s)){
    for(k in 1:length(h)){
      for (i in 1:length(freqs)){
        p = freqs[i]
        q = 1 - p
        Tfits = c(1 + s[j], 1 + (s[j] * h3[k]), 1 + (s[j] * h[k]), 1 + (s[j] * h1[k]), 1) / (1 + s[j])
        TGfreqs = c(p ^ 4, 4 * p^3 * q, 6 * p^2 * q^2, 4 * p * q^3, q ^ 4)
        Tw.bar = sum(TGfreqs * Tfits)
        Dfits = c(1 + s[j], 1 + (s[j] * h[k]), 1) / (1 + s[j])
        DGfreqs = c(p ^ 2, 2 * p * q, q ^ 2)
        Dw.bar = sum(DGfreqs * Dfits)
        Tvar.w = sum((TGfreqs * (Tfits - Tw.bar) ^ 2) / length(TGfreqs))
        Dvar.w = sum((DGfreqs * (Dfits - Dw.bar) ^ 2) / length(DGfreqs))
        Dw.bar.p = (p * Dfits[1]) + (q * Dfits[2])
        Tw.bar.p = ((p ^ 3) * Tfits[1]) + (3 * (p ^ 2) * q * Tfits[2]) + (3 * p * (q ^ 2) * Tfits[3]) + ((q ^ 3) * Tfits[4])
        df = rbind(df, c(s[j], h[k], p, Dw.bar, Dvar.w, Dw.bar.p, NA, NA, 2))
        df = rbind(df, c(s[j], h[k], p, Tw.bar, Tvar.w, Tw.bar.p, h1[k], h3[k], 4))
      }
    }
  }
  return(df[-1,])
}

# Generate fixation probabilities for diploids for a given initial frequency, selection coefficient and population size
getFix = function(p, s, N){
  num = 1 - exp(-4*N*s*p)
  den = 1 - exp(-4*N*s)
  return(num / den)
}

#Functions for simulating selective sweep trajectories and running coalescent simulations

# Stochastic simulation of selection for a beneficial allele
PloidyForSim = function(ploidy,  N, s, h, start_freq, end_freq, maxGen, maxTries){
  
  assert_that(ploidy - length(h) == 1)
  attempts = 0
  
  #Sample individuals to reproduce
  # taken from https://stats.stackexchange.com/questions/67911/how-to-sample-from-a-discrete-distribution
  sampleDist = function(PLOIDY, n, genoFreqs) { 
    sample(x = seq(0, PLOIDY), n, replace = T, prob = genoFreqs)
  }
  
  #Initalize population
  fits = sort( c(1, 1 + h * s, 1 + s) / (1 + s))
  Pop = rbinom(n = N, size = ploidy, prob = start_freq)
  p = sum(Pop) / (ploidy * N)
  freqs = c(p)
  
  while(p < end_freq & length(freqs) < maxGen & attempts < maxTries){
    #Retry if stochastic loss of beneficial allele in previous generation
    if ( p == 0 ){
      Pop = rbinom(n = N, size = ploidy, prob = start_freq)
      p = sum(Pop) / (ploidy * N)
      freqs = c(p)
      attempts = attempts + 1
    }
    
    geno_freq = c( table(factor(Pop, levels = as.character(seq(0, ploidy)))) / N )
    g_prime = ( fits * geno_freq ) / sum(fits * geno_freq) 
    # Make individuals for next generation
    p1 = sampleDist(length(g_prime) - 1, N, g_prime)
    p2 = sampleDist(length(g_prime) - 1, N, g_prime)
    Pop = c()
    for (i in 1:N){
      g1 = sample(c(rep(0, ploidy - p1[i]), rep(1, p1[i])), ploidy / 2, replace = F)
      g2 = sample(c(rep(0, ploidy - p2[i]), rep(1, p2[i])), ploidy / 2, replace = F)
      Pop = c(Pop, sum(g1) + sum(g2))
    }
    p = sum(Pop) / (ploidy * N)
    print(p)
    freqs = c(freqs, p)
  }
  return(freqs)
}

writeTraj = function(file, traj, numPops, selPops, timeScale, trajName = "rep1", startGen = 1){
  gen = seq(startGen, startGen + length(traj)) / timeScale
  fileConn = file(file)
  writeLines(c(paste("ntraj:", "1"), paste("npop:", numPops), paste("n:", length(traj), trajName)), con = fileConn)
  dat = data.frame("gen" = gen, "traj" = c(rev(traj), 0), "anc" = rep(0, length(traj) + 1))
  write.table(dat, file, quote = F, col.names = F, row.names = F, sep = "\t", append = T)
  close(fileConn)
}

msselRun = function(N, n, trajectory, name, L = 1000000, mu = 1e-8, r = 1e-8, ploidy = 2, numWindows = 200, slideRate = 0.5, selPos = 0.5, npop = 2, selPop = 1, fuseGen = 1, sampleGen = 1, ms = "/Users/pmonnahan/Documents/Research/code/dmc/mssel_modified/mssel"){ 
  options(scipen=999)
  p = ploidy
  time_scale = 2 * N * p
  theta = 2 * N * p * mu * L 
  rho = 2 * N * p * r * L
  fuse_time =  (length(trajectory) + sampleGen + fuseGen) / time_scale
  mig_str = paste(npop, "0", n * p, n * p, "0 -ej", fuse_time, "2 1")
  
  # Prepare input/output
  writeTraj(paste(name, ".traj.txt", sep = ""), trajectory, npop, selPop, time_scale, name, sampleGen)
  out1 = paste(name, ".ms.out", sep = "")
  
  # Format argument strings
  args1 = paste((2 * p * n), 1, n * p, n * p, traj_file, L * selPos, "-r", rho, L, "-t", theta, "-I", mig_str)
  
  # Run mssel and analysis script
  print(paste(ms, args1, ">", out1))
  cmd1 = system2(ms, args1, stdout = out1)
  return(cmd1)
}

msselCalc <- function(in_file, numWindows, outgroup, slideRate=0.5, selPop=1, linkage_stats = c("Kelly.Z_nS"), neutrality_stats = c("Tajima.D", "Fay.Wu.H", "Zeng.E")){
  # Define some necessary functions
  read.ms.output2 <- function(txt=NA, file.ms.output=NA, MSMS=FALSE) {
    
    if( !is.na(file.ms.output) ) txt <- scan(file=file.ms.output,
                                             what=character(0), sep="\n", quiet=TRUE)
    if( is.na(txt[1]) ){
      print("Usage: read.ms.output(txt), or read.ms.output(file=filename)")
      return()
    }
    
    
    if(MSMS[1]==FALSE){
      nsam   <- as.integer(strsplit(txt[1], split=" ")[[1]][2] )
      ndraws <- as.integer(strsplit(txt[1], split=" ")[[1]][3] )
    }
    
    #print(strsplit(txt[1], split=" "))
    
    h         <- numeric()
    result    <- list()
    gamlist   <- list()
    positions <- list()
    
    #marker <- grep("prob",txt)
    #probs <- sapply(strsplit(txt[marker], split=":"), function(vec) as.numeric(vec[2]))
    #marker <- grep("time",txt)
    #times <- sapply(strsplit(txt[marker], split="\t"), function(vec){ as.numeric(vec[2:3])} )
    times <- NaN
    probs <- NaN
    
    ## THE OUTPUT TEXT FOR EACH DRAW SHOULD CONTAIN THE WORD "segsites"
    
    
    marker <- grep("segsites", txt)
    
    if(MSMS[1]!=FALSE){ndraws <- length(marker);nsam <- MSMS$nsam} # MSMS
    
    
    stopifnot(length(marker) == ndraws)
    
    
    
    
    ## GET NUMBERS OF SEGREGATING SITES IN EACH DRAW
    segsites <- sapply(strsplit(txt[marker], split=" "), function(vec) as.integer(vec[2]) )
    
    
    for(draw in seq(along=marker)) {
      # if(!(draw %% 100)) cat(draw, " ")
      if(segsites[draw] > 0) {
        tpos <- strsplit(txt[marker[draw]+1], split=" ")
        positions[[draw]] <- as.numeric( tpos[[1]][ 2:(segsites[draw]+1) ] )
        
        haplotypes <- txt[(marker[draw] + 2):(marker[draw] + 2 + nsam - 1)]
        
        haplotypes <- strsplit(haplotypes, split="")
        
        h <- sapply(haplotypes, function(el) c(as.integer(el)))
        
        ## IF THERE'S 1 SEGREGATING SITE, THIS WON'T BE A MATRIX 
        
        if(segsites[draw] == 1) h <- as.matrix(h)
        ## OTHERWISE, IT NEEDS TO BE TRANSPOSED
        else h <- t(h)
        
        
      }
      else {
        h <- matrix(nrow=nsam, ncol=0)
        positions[[draw]] <- NA	
      }
      
      gamlist[[draw]] <- h
      stopifnot(all(dim(h) == c(nsam, segsites[draw]))) 
    }
    
    list(segsites=segsites, gametes=gamlist, probs=probs, times=t(times), positions=positions, nsam=nsam, nreps=ndraws ) 
  }
  makePops <- function(sampSize_list){
    popList = list()
    start = 1
    running_tot = 0
    for (i in 1:length(sampSize_list)){
      running_tot = running_tot + sampSize_list[i]
      popList[[i]] = seq(start,running_tot)
      start = start + sampSize_list[i]
    }
    return(popList)
  }
  munge <- function(inp, Nsites, positions, selPop, numWindows){
    coords = c("rep","snp.start", "snp.end")
    DF = data.frame()
    parseInfo = function(info){
      df = as.data.frame(info)
      df$snp.start=as.numeric(str_split_fixed(rownames(df)," ",4)[,2])
      df$snp.end=as.numeric(str_split_fixed(str_split_fixed(rownames(df)," ",4)[,4]," ",2)[,1])
      df$rep=str_split_fixed(str_split_fixed(rownames(df),"_",3)[,3], "[.]",2)[,1]
      return(df)
    }
    for (i in 1:length(inp@populations)){
      div = parseInfo(get.diversity(inp)[[i]])
      ld = parseInfo(get.linkage(inp)[[i]])
      neutrality_tests = parseInfo(get.neutrality(inp)[[i]])
      df = merge(div, ld, by = coords)
      df = merge(df, neutrality_tests[, c(coords, neutrality_stats)], by = coords)
      df$pop = as.character(i)
      DF=rbind(DF,df)
    }
    vars = c(coords, "nuc.diversity.within", linkage_stats, neutrality_stats, "pop")
    DF = melt(DF[,vars], id.vars = c(coords, "pop"))
    DF = dcast(DF, rep + snp.start + snp.end ~ variable + pop, value.var = "value")
    pos = melt(positions)
    colnames(pos)[2] = "rep"
    pos$rep=as.character(pos$rep)
    ends = pos %>% group_by(rep) %>% mutate(snp.end=row_number()) %>% slice(unique(DF$snp.end)) %>% as.data.frame()
    starts = pos %>% group_by(rep) %>% mutate(snp.start=row_number()) %>% slice(unique(DF$snp.start)) %>% as.data.frame()
    colnames(ends)[1] = "bp.end"
    colnames(starts)[1] = "bp.start"
    ends$bp.end = ends$bp.end * Nsites
    starts$bp.start = starts$bp.start * Nsites
    DF = merge(DF,starts,by=c("snp.start","rep"))
    DF = merge(DF,ends,by=c("snp.end","rep"))
    
    DF = DF[order(DF$rep, DF$snp.start),]
    
    dxy = dcast(melt(inp@nuc.diversity.between), Var2 ~ Var1, value.var="value")
    colnames(dxy)[2] = 'dxy'
    fst = dcast(melt(inp@nuc.F_ST.pairwise), Var2 ~ Var1, value.var="value")
    colnames(fst)[2] = 'fst'
    DF = cbind(DF, dxy)[,-2]
    DF = cbind(DF, fst)[,-2]
    DF$selPop = selPop
    return(DF)
  }
  
  # Read output of mssel simulation
  print("Reading mssel input...")
  sim = readMS(in_file)
  positions = read.ms.output2(file.ms.output = in_file)
  positions = positions$positions
  winSize = round(length(positions[[1]]) / numWindows, digits = -1)
  print("Done")
  
  # Define populations...number of anc and der specified in command line?
  pops = makePops(samp_sizes)
  sim=set.populations(sim,pops) #Only works for two populations currently
  sim=diversity.stats(sim)
  sim=diversity.stats.between(sim)
  
  ## EXPERIMENTAL ##
  sim = set.outgroup(sim, new.outgroup = outgroup)
  ## ##########
  
  # Divide simulated data into windows (type=1: based on SNP counts; type=2: based on nucleotide counts)
  cat(paste("","Creating windows...", "", sep = "\n"))
  sim.slide = sliding.window.transform(sim, width = winSize, jump = winSize * slideRate, type = 1, whole.data = FALSE)
  cat(paste("Done", "Calculating metrics...", "", sep = "\n"))
  sim.slide = diversity.stats(sim.slide)
  sim.slide = diversity.stats.between(sim.slide)
  sim.slide = linkage.stats(sim.slide)
  sim.slide = F_ST.stats(sim.slide)
  sim.slide = neutrality.stats(sim.slide, detail=TRUE)
  cat(paste("Done", "Munging data...", "", sep = "\n"))
  df = munge(sim.slide, Nsites, positions, selPop)
  df$mid = (df$bp.end + df$bp.start)/2
  df$bp.len = (df$bp.end - df$bp.start)
  df$Pi.1 = df$nuc.diversity.within_1 / df$bp.len
  df$Pi.2 = df$nuc.diversity.within_2 / df$bp.len
  df$dxy.1.2 = df[['dxy']] / df$bp.len
  print("Done")
  
  return(df)
}


## Class example data

# Calculate diploid/tetraploid genotype frequencies across entire range of allele frequency values
freqs = genoFreqs(seq(0,1,0.01))

# Plot HW genotype frequencies
ggplot(freqs, aes(x = p, y = value, color = as.factor(Ploidy), group = variable)) + geom_line() + scale_color_manual(name="Ploidy", values = c('red','blue')) + annotate("text", x = 0.25, y = 0.625, label="aa", color = "red") + annotate("text", x = 0.5, y = 0.55, label="Aa", color = "red") + annotate("text", x = 0.75, y = 0.625, label="AA", color = "red") + annotate("text", x = 0.085, y = 0.6, label="aaaa", color = "blue") + annotate("text", x = 0.085, y = 0.365, label="Aaaa", color = "blue") + annotate("text", x = 0.5, y = 0.395, label="AAaa", color = "blue") + annotate("text", x = 0.915, y = 0.365, label="AAAa", color = "blue") + annotate("text", x = 0.915, y = 0.6, label="AAAA", color = "blue") + ylab("Genotype Frequency") + xlab("Allele Frequency") + theme_bw()

# Calculate diploid/tetraploid fitness metrics across entire range of allele frequency values
fits = getFits(seq(0, 1, 0.01), c(0.1, 0.01, 0.001), c(0.25, 1, 0), c(0.5, 1, 0), c(0.75, 1, 0))

# Plot mean fitness 
ggplot(fits[fits$s==0.1,], aes(x=freq, y=w.bar, color = as.factor(ploidy), linetype=as.factor(h))) + geom_line() + scale_color_manual(name="Ploidy",values=c("red","blue")) + scale_linetype_discrete(name="Dominance", labels=c("Recessive","Additive","Dominant")) + xlab("Allele Frequency") + ylab("Mean Fitness") + theme_bw()

# Plot variance in fitness
ggplot(fits, aes(x=freq, y=var.w, color = as.factor(ploidy), linetype=as.factor(h))) + geom_line() + scale_color_manual(name="Ploidy",values=c("red","blue")) + scale_linetype_discrete(name="Dominance", labels=c("Recessive","Additive","Dominant")) + facet_wrap(~s, scale = "free_y")+ xlab("Allele Frequency") + ylab("Variance in Fitness") + theme_bw()

#Generate allele frequency trajectory for diploids and tetraploids given selection strength and dominance
s = 0.1
h = 0.5
h1 = 0.25
h2 = 0.5
h3 = 0.75

traj1 = rbind(dipTraj(s, h, 0.05, 0.99, 9999), tetTraj(s, h1, h2, h3, 0.05, 0.99, 9999))

# Plot allele frequency change over time for additive beneficial allele
ggplot(traj1, aes(x=gen, y=freq, color=ploidy)) + geom_line() + scale_color_manual(name="Ploidy",values=c('red','blue')) + theme_bw() + xlab("Generation") + ylab("Allele Frequency")

#Generate allele frequency trajectory for a range of basic values of s and h
traj = simTraj()

# Plot allele frequency change over time for additive, dominant, and recessive beneficial allele
ggplot(traj[traj$s==0.1,], aes(x=gen, y=freq, linetype=as.factor(h), color=ploidy)) + geom_line() + scale_color_manual(name="Ploidy",values=c("red","blue")) + scale_linetype_discrete(name="Dominance", labels=c("Recessive","Additive","Dominant")) + xlab("Generation") + ylab("Allele Frequency") + theme_bw() + xlim(0,300)

# Perform stochastic simulations for selection on beneficial allele
ploidy = 2
pop_size = 10000
selection_coeff = 0.1
dominance = 0.5 # Must be vector of 3 numbers if ploidy = 4

new_traj = PloidyForSim(ploidy, pop_size, selection_coeff, dominance, 0.05, 0.95, 999, 10)

# Write allele frequency trajectory to file
time_scale = 20000 # this should be 2*N for diploids or 4*N for tetraploids
traj_file = "/Users/pmonnahan/Documents/Research/PloidySim/testing.txt"
writeTraj(traj_file, new_traj, 2, 1, time_scale)

# Run mssel and calculate population genetic metrics in sliding windows across simulated region
samp_num = 10
dat = msselRun(pop_size, samp_num, traj_file, ploidy = ploidy)
