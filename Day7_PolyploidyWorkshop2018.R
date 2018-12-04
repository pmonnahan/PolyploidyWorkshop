#title: "Day7-PolyploidyWorkshop2018"
#author: "Patrick Monnahan"
#date: "11/23/2018"

####### Begin: Define Useful Functions ##########

# Copy and Paste this entire section (scroll down until you see "End: Define Useful Functions") into the R console

require(assertthat)
require(ggplot2)
require(dplyr)
require(tidyr)
require(PopGenome)
require(stringr)
require(data.table)

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
dipTraj = function(s, h, start_freq, N = -9, end_freq = 0.99, maxGens = 1){
  df = data.frame("s" = s, "h" = h, "gen" = 0, "freq" = start_freq, "dp1" = 0, "dp2" = 0, "w.bar" = 0, "var.w" = 0, "w.bar.p" = 0, "h1" = NA, "h3" = NA, "ploidy" = as.factor(2))
  p = start_freq
  gen = 0
  fits = c(1 + s, 1 + (s * h), 1) / (1 + s)
  
  while (p < end_freq & p > 0 & gen < maxGens){
    q = 1 - p
    Gfreqs = c(p ^ 2, 2 * p * q, q ^ 2)
    num = (Gfreqs[1] * fits[1]) + (p * q * fits[2])
    w.bar.p = (p * fits[1]) + (q * fits[2])
    w.bar = sum(Gfreqs * fits)
    p_prime = num / w.bar
    if (N != -9){
      p_prime = sum(rbinom(N, 2, p_prime)) / (2 * N) 
    }
    dp1 = p_prime - p
    dp2 = (p * (w.bar.p - w.bar)) / w.bar
    var.w = sum((Gfreqs * (fits - w.bar) ^ 2) / length(Gfreqs))
    df = rbind(df, c(s, h, gen, p, dp1, dp2, w.bar, var.w, w.bar.p, NA, NA, 2))
    p = p_prime
    gen = gen + 1
  }
  df = rbind(df, c(s, h, gen, p, dp1, dp2, w.bar, var.w, w.bar.p, NA, NA, 2))
  return(df[-1,])
}

tetTraj = function(s, h1, h2, h3, start_freq, N = -9, end_freq = 0.99, maxGens = 1){
  assert_that(h3 >= h1)
  df = data.frame("s" = s, "h" = h2, "gen" = 0, "freq" = start_freq, "dp1" = 0, "dp2" = 0, "w.bar" = 0, "var.w" = 0, "w.bar.p" = 0, "h1" = h1, "h3" = h3, "ploidy" = as.factor(4))
  p = start_freq
  gen = 0
  freq = c(p)
  fits = c(1 + s, 1 + (s * h3), 1 + (s * h2), 1 + (s * h1), 1) / (1 + s)
  while (p < end_freq & p > 0 & gen < maxGens){
    q = 1 - p
    Gfreqs = c(p ^ 4, 4 * p^3 * q, 6 * p^2 * q^2, 4 * p * q^3, q ^ 4)
    num = (Gfreqs[1] * fits[1]) + (3 * p ^ 3 * q * fits[2]) + (3 * p ^ 2 * q ^ 2 * fits[3]) + (p * q ^ 3 * fits[4])
    w.bar.p = ((p ^ 3) * fits[1]) + (3 * (p ^ 2) * q * fits[2]) + (3 * p * (q ^ 2) * fits[3]) + ((q ^ 3) * fits[4])
    w.bar = sum(Gfreqs * fits)
    p_prime = num / w.bar
    if (N != -9){
      p_prime = sum(rbinom(N, 4, p_prime)) / (4 * N) 
    }
    dp1 = p_prime - p
    dp2 = (p * (w.bar.p - w.bar)) / w.bar
    var.w = sum((Gfreqs * (fits - w.bar) ^ 2) / length(Gfreqs))
    df = rbind(df, c(s, h2, gen, p, dp1, dp2, w.bar, var.w, w.bar.p, h1, h3, 4))
    p = p_prime
    gen = gen + 1
  }
  df = rbind(df, c(s, h2, gen, p, dp1, dp2, w.bar, var.w, w.bar.p, h1, h3, 4))
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
getFix = function(p, s, N, ploidy){
  num = 1 - exp(-2 * ploidy * N * s * p)
  den = 1 - exp(-2 * ploidy * N * s)
  fp = num / den
  t = (2 * log(ploidy * N - 1)) / s
  return(data.frame("fix.prob" = fp, "fix.time" = t))
}


#Functions for simulating selective sweep trajectories and running coalescent simulations
# Stochastic simulation of selection for a beneficial allele
PloidyForSim = function(ploidy,  N, s, h, start_freq = 0.05, end_freq = 0.95, maxGen = 9999, maxTries = 10){
  
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
  
  while(p < end_freq & length(freqs) <= maxGen & attempts <= maxTries){
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
    print(paste("Frequency = ", p))
    freqs = c(freqs, p)
  }
  if (attempts > maxTries){
    print(paste("Beneficial allele was lost due to drift for", maxTries, "consecutive attempts"))
  }
  if (length(freqs) > maxGen){
    print(paste("Beneficial allele did not fix before the maximum number of generations (", maxGen, ")."))
  }
  return(freqs)
}

PloidyForSim2 = function(ploidy,  N, s, h, start_freq = 0.05, end_freq = 0.95, maxGen = 9999, maxTries = 10){
  
  assert_that(ploidy - length(h) == 1)
  attempts = 0
  freqs = 0
  while (attempts < maxTries & freqs[length(freqs)] == 0){
    if (ploidy == 2){
      freqs = dipTraj(s, h, start_freq, N, 0.99, maxGen)$freq
    }
    if (ploidy == 4){
      freqs = tetTraj(s, h[1], h[2], h[3], start_freq, N, 0.99, maxGen)$freq
    }
    if (freqs[length(freqs)] == 0){
      attempts = attempts + 1
    }
  }

  if (attempts > maxTries){
    print(paste("Beneficial allele was lost due to drift for", maxTries, "consecutive attempts"))
  }
  if (length(freqs) > maxGen){
    print(paste("Beneficial allele did not fix before the maximum number of generations (", maxGen, ")."))
  }
  return(freqs)
}

writeTraj = function(file, traj, numPops, selPops, timeScale, trajName = "rep1", startGen = 1){
  assert_that(length(traj) > 10, msg = "Error with trajectory")
  gen = seq(startGen, startGen + length(traj)) / timeScale
  fileConn = file(file)
  writeLines(c(paste("ntraj:", "1"), paste("npop:", numPops), paste("n:", length(traj) + 1, trajName)), con = fileConn)
  dat = data.frame("gen" = gen, "traj" = c(rev(traj), 0), "anc" = rep(0, length(traj) + 1))
  # dat = data.frame("gen" = gen, "traj" = c(rev(traj), 0))
  write.table(dat, file, quote = F, col.names = F, row.names = F, sep = "\t", append = T)
  close(fileConn)
}

msselRun = function(N, n, trajectory, outfile = -9, L = 1000000, mu = 1e-8, r = 1e-8, ploidy = 2, selPos = 0.5, npop = 2, selPop = 1, fuseGen = 1, sampleGen = 1, ms = "/Users/pmonnahan/Documents/Research/code/dmc/mssel_modified/mssel"){ 
  options(scipen=999)
  if (outfile == -9){
    name = paste(getwd(), "/mssel.", sample(1:99999999, 1), sep = "")
    outfile = paste(name, ".ms.out", sep = "")
    traj_file = paste(name, ".traj.txt", sep = "")
  } else {
    name = tools::file_path_sans_ext(outfile)
    traj_file = paste(name, ".traj.txt", sep = "")
  }
  p = ploidy
  time_scale = 2 * N * p
  theta = 2 * N * p * mu * L 
  rho = 2 * N * p * r * L
  fuse_time =  (length(trajectory) + sampleGen + fuseGen) / time_scale
  mig_str = paste(npop, "0", n * p, n * p, "0 -ej", fuse_time, "2 1")
  
  # Prepare input/output
  writeTraj(traj_file, trajectory, npop, selPop, time_scale, name, sampleGen)
  
  # Format argument strings
  args1 = paste((2 * p * n), 1, n * p, n * p, traj_file, L * selPos, "-r", rho, L, "-t", theta, "-I", mig_str, ">", outfile)
  
  # Run mssel
  print(paste(ms, args1))
  cmd1 = system2(ms, args1)
  return(outfile)
}

msselCalc <- function(in_file, numWindows, samp_sizes, Nsites, outgroup = 21, slideRate=0.5, selPop=1, linkage_stats = c("Kelly.Z_nS"), neutrality_stats = c("Tajima.D", "Fay.Wu.H", "Zeng.E")){
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

####### End: Define Useful Functions ##########



####### Begin: Class examples ###########

# Q84
dipTraj(s = 0.5, h = 0, start_freq = 0.5)
dipTraj(s = 0.5, h = 0.5, start_freq = 0.5)
dipTraj(s = 0.5, h = 1, start_freq = 0.5)
dipTraj(s = 0.5, h = 0, start_freq = 0.05)
dipTraj(s = 0.5, h = 0.5, start_freq = 0.05)
dipTraj(s = 0.5, h = 1, start_freq = 0.05)

# Q86
getFits(freqs = 0.5, s = 0.1, h1 = 0.25, h = 0.5, h3 = 0.75) %>% select(ploidy, w.bar)

# Q87
# Calculate diploid/tetraploid fitness metrics across entire range of allele frequency values
fits = getFits(seq(0, 1, 0.01), c(0.1, 0.01, 0.001), c(0.25, 1, 0), c(0.5, 1, 0), c(0.75, 1, 0))
# Plot mean fitness 
mean_fit = ggplot(fits[fits$s==0.1,], aes(x=freq, y=w.bar, color = as.factor(ploidy), linetype=as.factor(h))) + geom_line() + scale_color_manual(name="Ploidy",values=c("red","blue")) + scale_linetype_discrete(name="Dominance", labels=c("Recessive","Additive","Dominant")) + xlab("Allele Frequency") + ylab("Mean Fitness") + theme_bw()
mean_fit

# Q88
dipTraj(s = 0.1, h = 0.5, start_freq = 0.5) %>% select(dp1)
tetTraj(s = 0.1, h1 = 0.25, h2 = 0.5, h3 = 0.75, start_freq = 0.5) %>% select(dp1)
dipTraj(s = 0.1, h = 1, start_freq = 0.5) %>% select(dp1)
tetTraj(s = 0.1, h1 = 1, h2 = 1, h3 = 1, start_freq = 0.5) %>% select(dp1)
dipTraj(s = 0.1, h = 0, start_freq = 0.5) %>% select(dp1)
tetTraj(s = 0.1, h1 = 0, h2 = 0, h3 = 0, start_freq = 0.5) %>% select(dp1)

# Q89
#Generate allele frequency trajectory for diploids and tetraploids given selection strength and dominance
dip_add_traj = dipTraj(s = 0.1, h = 0.5, start_freq = 0.05, end_freq = 0.99, maxGens = 9999)
tet_add_traj = tetTraj(s = 0.1, h1 = 0.25, h2 = 0.5, h3 = 0.75, start_freq = 0.05, end_freq = 0.99, maxGens = 9999)

traj_add = rbind(dip_add_traj, tet_add_traj)

# Plot allele frequency change over time for additive beneficial allele
traj_plot_add = ggplot(traj_add, aes(x=gen, y=freq, color=ploidy)) + geom_line() + scale_color_manual(name="Ploidy",values=c('red','blue')) + theme_bw() + xlab("Generation") + ylab("Allele Frequency")

# Q92 
getFix(p = 1/200, s = 0.01, N = 100, ploidy = 2)
getFix(p = 1/2000, s = 0.01, N = 1000, ploidy = 2)
getFix(p = 1/20000, s = 0.01, N = 10000, ploidy = 2)

# Q93
getFix(p = 1/200, s = 0.01, N = 100, ploidy = 2)
getFix(p = 1/400, s = 0.01, N = 100, ploidy = 4)
getFix(p = 1/2000, s = 0.01, N = 1000, ploidy = 2)
getFix(p = 1/4000, s = 0.01, N = 1000, ploidy = 4) 

# Q94
getFix(p = 1/2000, s = 0.5 * 0.01, N = 1000, ploidy = 2)
getFix(p = 1/4000, s = 0.25 * 0.01, N = 1000, ploidy = 4)

# Q96
# Plot the mean fitness over time for an additive beneficial allele. 
fit_plot_add = ggplot(traj_add, aes(x=gen, y=w.bar, color=ploidy)) + geom_line() + scale_color_manual(name="Ploidy",values=c('red','blue')) + theme_bw() + xlab("Generation") + ylab("Mean Fitness")
fit_plot_add

# Generate allele frequency trajectory for diploids and tetraploids given selection strength and dominance Q97
dip_traj_dom = dipTraj(s = 0.1, h = 1, start_freq = 0.05, end_freq = 0.99, maxGens = 300)
tet_traj_dom = tetTraj(s = 0.1, h1 = 1, h2 = 1, h3 = 1, start_freq = 0.05, end_freq = 0.99, maxGens = 300)
traj_dom = rbind(dip_traj_dom, tet_traj_dom)

# Plot allele frequency change over time for additive beneficial allele
traj_plot_dom = ggplot(traj_dom, aes(x=gen, y=freq, color=ploidy)) + geom_line() + scale_color_manual(name="Ploidy",values=c('red','blue')) + theme_bw() + xlab("Generation") + ylab("Allele Frequency")
# Plot the mean fitness over time for a dominant beneficial allele. 
fit_plot_dom = ggplot(traj_dom, aes(x=gen, y=w.bar, color=ploidy)) + geom_line() + scale_color_manual(name="Ploidy",values=c('red','blue')) + theme_bw() + xlab("Generation") + ylab("Mean Fitness")
traj_plot_dom
fit_plot_dom

# Q100
dip_traj1 = dipTraj(s = 0.01, h = 0, start_freq = 0.00316, end_freq = 0.99, maxGens = 9999)
tet_traj1 = tetTraj(s = 0.01, h1 = 0, h2 = 0.5, h3 = 0.75, start_freq = 0.0562, end_freq = 0.99, maxGens = 9999)
traj1 = rbind(dip_traj1, tet_traj1)
traj1_plot = ggplot(traj1, aes(x=gen, y=freq, color=ploidy)) + geom_line() + scale_color_manual(name="Ploidy",values=c('red','blue')) + theme_bw() + xlab("Generation") + ylab("Allele Frequency")
traj1_plot

dip_traj2 = dipTraj(s = 0.01, h = 0.5, start_freq = 0.00002, end_freq = 0.99999999, maxGens = 9999)
tet_traj2 = tetTraj(s = 0.01, h1 = 0.25, h2 = 0.5, h3 = 0.75, start_freq = 0.00004, end_freq = 0.999999999, maxGens = 9999)
traj2 = rbind(dip_traj2, tet_traj2)
traj2_plot = ggplot(traj2, aes(x=gen, y=freq, color=ploidy)) + geom_line() + scale_color_manual(name="Ploidy",values=c('red','blue')) + theme_bw() + xlab("Generation") + ylab("Allele Frequency")
traj2_plot

####### Part 4: Linked Selection ########
# Set your parameters here
ploidy = 2
pop_size = 10000 # Keep this value above 100 and below 1000000 (computation time will increase with increasing pop_size)
selection_coeff = 0.1 # Keep this between 0 and 1
dominance = 0.5 # Must be vector of 3 numbers if ploidy = 4. For example, c(0.25, 0.5, 0.75) for an additive allele.
seq_len = 1000000 # Length of sequence that we will simulate with mssel.  Increasing this value will increase computation time.
mutation_rate = 1e-8 # per-base mutation rate; mu
recomb_rate = 1e-8 # per-base recombination rate; r
samp_num = 10 # number of individuals to sample

### HERE YOU MUST ENTER THE PATH TO MSSEL
path_to_mssel = "XXXXXXX"
path_to_mssel = "/Users/pmonnahan/Documents/Research/code/dmc/mssel_modified/mssel"
 
# Perform stochastic simulations for selection on beneficial allele
new_traj = PloidyForSim(ploidy, pop_size, selection_coeff, dominance)

# Run mssel
infile = msselRun(N = pop_size, n = samp_num, new_traj, L = seq_len, mu = mutation_rate, r = recomb_rate, ploidy = ploidy, ms = path_to_mssel)

# calculate population genetic metrics in sliding windows across simulated region
dat = msselCalc(infile, numWindows = 200, rep(ploidy * samp_num, 2), Nsites = seq_len)

# Plotting
ggplot() + geom_line(data = dat, aes(x = bp.end, y = Pi.1), color = "red") + geom_line(data = dat, aes(x = bp.end, y = Pi.2), color = "blue") + xlab("Position (bp)") + ylab("Diversity")

ggplot() + geom_line(data = dat, aes(x = bp.end, y = Kelly.Z_nS_1), color = "red") + geom_line(data = dat, aes(x = bp.end, y = Kelly.Z_nS_2), color = "blue") + xlab("Position (bp)")

ggplot() + geom_line(data = dat, aes(x = bp.end, y = Pi.1), color = "red") + geom_line(data = dat, aes(x = bp.end, y = Pi.2), color = "blue") + xlab("Position (bp)")

ggplot() + geom_line(data = dat, aes(x = bp.end, y = fst)) + xlab("Position (bp)")

####### END: Class examples #######

####### EXTRA #######
# Calculate diploid/tetraploid genotype frequencies across entire range of allele frequency values
freqs = genoFreqs(seq(0,1,0.01))

# Plot HW genotype frequencies
ggplot(freqs, aes(x = p, y = value, color = as.factor(Ploidy), group = variable)) + geom_line() + scale_color_manual(name="Ploidy", values = c('red','blue')) + annotate("text", x = 0.25, y = 0.625, label="aa", color = "red") + annotate("text", x = 0.5, y = 0.55, label="Aa", color = "red") + annotate("text", x = 0.75, y = 0.625, label="AA", color = "red") + annotate("text", x = 0.085, y = 0.6, label="aaaa", color = "blue") + annotate("text", x = 0.085, y = 0.365, label="Aaaa", color = "blue") + annotate("text", x = 0.5, y = 0.395, label="AAaa", color = "blue") + annotate("text", x = 0.915, y = 0.365, label="AAAa", color = "blue") + annotate("text", x = 0.915, y = 0.6, label="AAAA", color = "blue") + ylab("Genotype Frequency") + xlab("Allele Frequency") + theme_bw()

# Plot variance in fitness
var_fit = ggplot(fits, aes(x=freq, y=var.w, color = as.factor(ploidy), linetype=as.factor(h))) + geom_line() + scale_color_manual(name="Ploidy",values=c("red","blue")) + scale_linetype_discrete(name="Dominance", labels=c("Recessive","Additive","Dominant")) + facet_wrap(~s, scale = "free_y")+ xlab("Allele Frequency") + ylab("Variance in Fitness") + theme_bw()

#Generate allele frequency trajectory for a range of basic values of s and h
traj = simTraj()

# Plot allele frequency change over time for additive, dominant, and recessive beneficial allele
traj_all_plot = ggplot(traj[traj$s==0.1,], aes(x=gen, y=freq, linetype=as.factor(h), color=ploidy)) + geom_line() + scale_color_manual(name="Ploidy",values=c("red","blue")) + scale_linetype_discrete(name="Dominance", labels=c("Recessive","Additive","Dominant")) + xlab("Generation") + ylab("Allele Frequency") + theme_bw() + xlim(0,300)

ploidy = 2
pop_size = 10000 # Keep this value above 100 and below 1000000 (computation time will increase with increasing pop_size)
selection_coeff = 0.1 # Keep this between 0 and 1
dominance = 0.5 # Must be vector of 3 numbers if ploidy = 4. For example, c(0.25, 0.5, 0.75) for an additive allele.
seq_len = 1000000 # Length of sequence that we will simulate with mssel.  Increasing this value will increase computation time.
mutation_rate = 1e-8 # per-base mutation rate; mu
recomb_rate = 1e-8 # per-base recombination rate; r
samp_num = 10 # number of individuals to sample
path_to_mssel = "/Users/pmonnahan/Documents/Research/code/dmc/mssel_modified/mssel"

dat = data.frame()
num_reps = 10
for (i in 1:num_reps){
  new_traj = PloidyForSim2(ploidy, pop_size, selection_coeff, dominance)
  
  # Run mssel
  infile = msselRun(N = pop_size, n = samp_num, new_traj, L = seq_len, mu = mutation_rate, r = recomb_rate, ploidy = ploidy, ms = path_to_mssel)
  
  # calculate population genetic metrics in sliding windows across simulated region
  ndat = msselCalc(infile, numWindows = 200, rep(ploidy * samp_num, 2))
  ndat['rep'] = i
  dat = rbind(dat, ndat)
}

