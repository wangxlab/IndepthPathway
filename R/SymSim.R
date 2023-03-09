#' Make a gct file format from expression data matrix
#'
#' @param ExpMatrix Gene expression data matrix. Row names are gene symbols and column names are sample names.
#'
#' @return GCT file format of gene expression matrix
#' @export
#'
#' @examples MakeGCT(ExpMatrix)
MakeGCT<-function(ExpMatrix) {
  #tpmlog2.matrix<-trueData.TPMlog2
  newRowName<-c("#1.2",ncol(ExpMatrix), rownames(ExpMatrix))
  header<-matrix(0,ncol=ncol(ExpMatrix), nrow=2)
  colnames(header)<-colnames(ExpMatrix)
  gctheader.tpmlog.matrix<-rbind(header,ExpMatrix)
  rownames(gctheader.tpmlog.matrix)<-newRowName

  ExpMatrix.Descript<-tibble::rownames_to_column(as.data.frame(gctheader.tpmlog.matrix),"Description")
  rownames(ExpMatrix.Descript)<-ExpMatrix.Descript$Description
  ExpMatrix.Name.Descript<- tibble::rownames_to_column(as.data.frame(ExpMatrix.Descript),"Name")  #
  ExpMatrix.Name.Descript[1:2,2:ncol(ExpMatrix.Name.Descript)]<-""
  ExpMatrix.Name.Descript[2,1:2]<-c(nrow(ExpMatrix), ncol(ExpMatrix))

  ExpData.gct.cellname<-rbind(colnames(ExpMatrix.Name.Descript), ExpMatrix.Name.Descript)
  ExpData.gct.rowoder<-ExpData.gct.cellname[c(2,3,1,4:nrow(ExpData.gct.cellname)),]
  return(ExpData.gct.rowoder)
}


#' Title Simulate observed count matrix given technical biases and the true counts
#'
#' @param SE Default is NULL
#' @param true_counts Gene expression count data.frame
#' @param protocol a string, can be "nonUMI" or "UMI"
#' @param alpha_mean the mean of rate of subsampling of transcripts during capture step, default at 10 percent efficiency
#' @param alpha_sd the std of rate of subsampling of transcripts
#' @param lenslope amount of length bias
#' @param nbins number of bins for gene length
#' @param gene_len a vector with lengths of all genes
#' @param amp_bias_limit range of amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high, default is 0.8
#' @param nPCR1 the number of PCR cycles in "pre-amplification" step, default is 16
#' @param nPCR2 the number of PCR cycles used after fragmentation.
#' @param LinearAmp if linear amplification is used for pre-amplification step, default is FALSE
#' @param LinearAmp_coef the coeficient of linear amplification, that is, how many times each molecule is amplified by
#' @param depth_mean mean of sequencing depth
#' @param depth_sd std of sequencing depth
#' @param nbatch
#'
#' @return simulated observed count matrix
#' @export
#'
#' @examples True2ObservedCounts(true_counts=ExpDataProc_RmvNA, protocol=protocol, alpha_mean=MyAlpha, alpha_sd=0.02, gene_len=GenelengthDataMatchProc, nPCR1= MyP, depth_mean=5e4, depth_sd=3e3)
True2ObservedCounts <- function(SE=NULL,true_counts,protocol,alpha_mean=0.1,alpha_sd=0.002,  # meta_cell
                                lenslope=0.02,nbins=20,gene_len,amp_bias_limit=c(-0.2, 0.2),
                                rate_2PCR=0.8,nPCR1=16, nPCR2=10, LinearAmp=F, LinearAmp_coef=2000,
                                depth_mean, depth_sd, nbatch=1){
  print(paste0("protocol: ", protocol))
  print(paste0("alpha_mean: ", alpha_mean))
  print(paste0("nPCR1: ", nPCR1))
  if(!is.null(SE)){
    #meta_cell <- colData(SE)
    true_counts <- assays(SE)$count
  }
  ngenes <- dim(true_counts)[1]; ncells <- dim(true_counts)[2]
  amp_bias <- cal_amp_bias(lenslope, nbins, gene_len, amp_bias_limit)
  rate_2cap_lb <- 0.0005; depth_lb <- 200 # lower bound for capture efficiency and sequencing depth
  rate_2cap_vec <- rnorm_trunc(n=ncells, mean = alpha_mean, sd=alpha_sd, a=rate_2cap_lb, b=Inf)
  depth_vec <- rnorm_trunc(n=ncells, mean = depth_mean, sd=depth_sd,a=depth_lb,b=Inf)

  observed_counts <- lapply(c(1:ncells),function(icell){
    amplify_1cell(true_counts_1cell =  true_counts[, icell], protocol=protocol,
                  rate_2cap=rate_2cap_vec[icell], gene_len=gene_len, amp_bias = amp_bias,
                  rate_2PCR=rate_2PCR, nPCR1=nPCR1, nPCR2=nPCR2, LinearAmp = LinearAmp,
                  LinearAmp_coef = LinearAmp_coef, N_molecules_SEQ = depth_vec[icell])
  })
  ## assign random batch ID to cells
  batchIDs <- sample(1:nbatch, ncells, replace = TRUE)
  #meta_cell2 <- data.frame(alpha=rate_2cap_vec,depth=depth_vec, batch=batchIDs)
  #meta_cell <- cbind(meta_cell, meta_cell2)

  if (protocol=="UMI"){
    UMI_counts <- do.call(cbind, lapply(observed_counts, "[[", 1))
    nreads_perUMI <- lapply(observed_counts, "[[", 2)
    nUMI2seq <- sapply(observed_counts, "[[", 3)
    observed_counts <- UMI_counts
  } else
    observed_counts <- do.call(cbind,observed_counts)

  ## add batch effects to observed counts
  # use different mean and same sd to generate the multiplicative factor for different gene in different batch
  if (nbatch>1){
    mean_matrix <- matrix(0, ngenes, nbatch)
    batch_effect_size <- 2
    gene_mean <- rnorm(ngenes, 0, 1)
    temp <- lapply(1:ngenes, function(igene) {
      return(runif(nbatch, min = gene_mean[igene]-batch_effect_size, max = gene_mean[igene]+batch_effect_size))
    })
    mean_matrix <- do.call(rbind, temp)

    batch_factor <- matrix(0, ngenes, ncells)
    for (igene in 1:ngenes){
      for (icell in 1:ncells){
        batch_factor[igene, icell] <- rnorm(n=1, mean=mean_matrix[igene, batchIDs[icell]], sd=0.01)
      }
    }
    observed_counts <- 2^(log2(observed_counts)+batch_factor)
  }

  if(is.null(SE)){
    if (protocol=="UMI"){return(list(counts=observed_counts, nreads_perUMI=nreads_perUMI,   # cell_meta=meta_cell,
                                     nUMI2seq=nUMI2seq))} else
                                       return(list(counts=observed_counts))  #cell_meta=meta_cell
  } else{
    assays(SE)$observed_counts <- observed_counts
    #colData(SE)<-meta_cell
    return(SE)
  }
}


#' Simulate technical biases
#' @param lenslope amount of length bias. This value sould be less than 2*amp_bias_limit[2]/(nbins-1)
#' @param nbins number of bins for gene length
#' @param gene_len transcript length of each gene
#' @param amp_bias_limit range of amplification bias for each gene, a vector of length ngenes
#'
#' @return simulated technicall biases
#' @export
#'
#' @examples cal_amp_bias(lenslope, nbins, gene_len, amp_bias_limit)
cal_amp_bias <- function(lenslope, nbins, gene_len, amp_bias_limit){

  ngenes <- length(gene_len)
  len_bias_bin <- (-c(1:nbins))*lenslope
  len_bias_bin <- len_bias_bin-median(len_bias_bin)
  if (max(len_bias_bin) > amp_bias_limit[2]) {
    stop("The lenslope parameter is too large.")
  }
  max_rand_bias <- amp_bias_limit[2] - max(len_bias_bin)

  rand_bias <- rnorm(ngenes, mean=0, sd=max_rand_bias)
  rand_bias[rand_bias > max_rand_bias] <- max_rand_bias
  rand_bias[rand_bias < -max_rand_bias] <- -max_rand_bias
  #rand_bias <- runif(ngenes, -max_rand_bias,  max_rand_bias)

  binsize <- floor(ngenes/nbins)
  genes_in_bins <- vector("list", nbins)
  bin4genes <- numeric(ngenes)
  for (ibin in 1:(nbins-1)){
    genes_in_bins[[ibin]] <- order(gene_len)[((ibin-1)*binsize+1) : (ibin*binsize)]
    bin4genes[genes_in_bins[[ibin]]] <- ibin
  }
  genes_in_bins[[nbins]] <- order(gene_len)[((nbins-1)*binsize+1) : ngenes]
  bin4genes[genes_in_bins[[nbins]]] <- nbins

  len_bias <- numeric(ngenes); len_bias <- len_bias_bin[bin4genes]
  amp_bias <- rand_bias+len_bias
  return(amp_bias)
}


#' sample from truncated normal distribution
#' @param a the minimum value allowed
#' @param b the maximum value allowed
#'
#' @return vector of turncated normal distribution
#' @export
#'
#' @examples rnorm_trunc(n, mean, sd, a, b)
rnorm_trunc <- function(n, mean, sd, a, b){
  vec1 <- rnorm(n, mean = mean, sd=sd)
  beyond_idx <- which(vec1 < a | vec1 > b)
  if (length(beyond_idx) > 0) { # for each value < rate_2cap_lb
    substi_vec <- sapply(1:length(beyond_idx), function(i){
      while (TRUE){
        temp <- rnorm(1, mean = mean, sd=sd)
        if (temp > a & temp < b) {break}}
      return(temp)} )
    vec1[beyond_idx] <- substi_vec
  }
  return(vec1)
}


#' This function simulates the amplification, library prep, and the sequencing processes.
#' @param true_counts_1cell the true transcript counts for one cell (one vector)
#' @param protocol a string, can be "nonUMI" or "UMI"
#' @param rate_2cap the capture efficiency for this cell
#' @param gene_len gene lengths for the genes/transcripts, sampled from real human transcript length
#' @param amp_bias amplification bias for each gene, a vector of length ngenes
#' @param rate_2PCR PCR efficiency, usually very high
#' @param nPCR1 the number of PCR cycles
#' @param LinearAmp if linear amplification is used for pre-amplification step, default is FALSE
#' @param LinearAmp_coef the coeficient of linear amplification, that is, how many times each molecule is amplified by
#' @param N_molecules_SEQ number of molecules sent for sequencing; sequencing depth
#'
#' @return read counts (if protocol="nonUMI") or UMI counts (if protocol="UMI)
#' @export
#'
#' @examples amplify_1cell(true_counts_1cell, protocol, rate_2cap, gene_len, amp_bias, rate_2PCR, nPCR1, nPCR2, LinearAmp, LinearAmp_coef, N_molecules_SEQ)
amplify_1cell <- function(true_counts_1cell, protocol, rate_2cap, gene_len, amp_bias,
                          rate_2PCR, nPCR1, nPCR2, LinearAmp, LinearAmp_coef, N_molecules_SEQ){
  ngenes <- length(gene_len)
  if (protocol=="nonUMI"){data(len2nfrag)} else
    if(protocol=="UMI"){ } else
    {stop("protocol input should be nonUMI or UMI")}
  inds <- vector("list",2)
  # expand the original vector and apply capture efficiency
  # maintain a transcript index vector: which transcript the molecule belongs to
  expanded_res <- expand2binary(c(true_counts_1cell,1))
  expanded_vec <- expanded_res[[1]]; trans_idx <- expanded_res[[2]]

  inds[[1]] <- which(expanded_vec > 0); expanded_vec <- expanded_vec[inds[[1]]]
  trans_idx <- trans_idx[inds[[1]]]

  captured_vec <- expanded_vec; captured_vec[runif(length(captured_vec)) > rate_2cap] <- 0
  if (sum(captured_vec[1:(length(captured_vec)-1)]) < 1) {return(rep(0, ngenes))}
  captured_vec[length(captured_vec)] <- 1
  inds[[2]] <- which(captured_vec > 0); captured_vec <- captured_vec[inds[[2]]]
  trans_idx <- trans_idx[inds[[2]]]

  amp_rate <- c((rate_2PCR+amp_bias[trans_idx[1:(length(trans_idx)-1)]]),1)

  # pre-amplification:
  if (LinearAmp){
    PCRed_vec <- captured_vec*LinearAmp_coef
  } else {
    temp <- runif(length(captured_vec)) < amp_rate
    temp <- temp*2+captured_vec-temp
    for (iPCR in 2:nPCR1){
      eff <- runif(length(temp))*amp_rate
      v1 <- temp*(1-eff)
      round_down <- (v1-floor(v1)) < runif(length(v1))
      v1[round_down] <- floor(v1[round_down]); v1[!round_down] <- ceiling(v1[!round_down])
      temp <- v1 + 2*(temp-v1)
    }
    PCRed_vec <- temp
  }

  if (protocol=="nonUMI"){ # add fragmentation step here
    temp_vec <- PCRed_vec
    for (i in seq(2,1,-1)){
      temp_vec1 <- numeric(); temp_vec1[inds[[i]]] <- temp_vec;
      temp_vec <- temp_vec1; temp_vec[is.na(temp_vec)] <- 0
    }
    recovered_vec <- temp_vec[1:(length(temp_vec)-1)]
    amp_mol_count=numeric(ngenes);
    GI=c(0, cumsum(true_counts_1cell));
    for (i in which(true_counts_1cell>0)){
      x=recovered_vec[(GI[i]+1):GI[i+1]]
      amp_mol_count[i]=sum(x)
    }

    # for every copy of each transcript, convert it into number of fragments
    frag_vec <- numeric(ngenes)
    for (igene in which(amp_mol_count>0)){
      #print(paste0("current igene: ", igene))
      if(! (as.character(gene_len[igene])) %in% rownames(len2nfrag)) {  # I fixed here.
        #print(paste0("skipped igene: ", igene))
        next;
      }
      dim(len2nfrag)# 8118 500; head(rownames(len2nfrag)) #  "8"  "9"  "11" "13" "15" "16".
      ## gene_len[igene]: 3474    len2nfrag rownames doesn't have '3474'
      frag_vec[igene] <- sum(sample(len2nfrag[as.character(gene_len[igene]),], amp_mol_count[igene], replace = TRUE))
    }
    # another 8 rounds of amplification to the fragments (fragmentation bias gets amplified)
    for (iPCR in 1:2){
      frag_vec <- frag_vec + sapply(frag_vec, function(x) rbinom(n=1, x, prob = rate_2PCR))
    }
    for (iPCR in 3:nPCR2){
      frag_vec <- frag_vec + round(frag_vec*rate_2PCR)
    }
    SEQ_efficiency=N_molecules_SEQ/sum(frag_vec)
    if (SEQ_efficiency >= 1) {read_count <- frag_vec} else{
      read_count <- sapply(frag_vec,function(Y){rbinom(n=1,size=Y,prob=SEQ_efficiency)}) }
    return(read_count)
  } else if (protocol=="UMI"){

    prob_vec <- sapply(gene_len[trans_idx[1:(length(trans_idx)-1)]], get_prob)
    # fragmentation:
    frag_vec <- sapply(1:(length(PCRed_vec)-1), function(igene)
    {return(rbinom(n=1, size = PCRed_vec[igene], prob = prob_vec[igene] ))})

    # another 10 rounds of amplification to the fragments (fragmentation bias gets amplified)
    for (iPCR in 1:2){
      frag_vec <- frag_vec + sapply(frag_vec, function(x) rbinom(n=1, x, prob = rate_2PCR))
    }

    frag_vec <- round(frag_vec * (1+rate_2PCR)^(nPCR2-1))

    SEQ_efficiency <- N_molecules_SEQ/sum(frag_vec)
    if (SEQ_efficiency >= 1){sequenced_vec <- frag_vec} else {
      sequenced_vec <- sapply(frag_vec,function(Y){rbinom(n=1,size=Y,prob=SEQ_efficiency)})}

    temp_vec <- c(sequenced_vec,1)
    for (i in seq(2,1,-1)){
      temp_vec1 <- numeric(); temp_vec1[inds[[i]]] <- temp_vec;
      temp_vec <- temp_vec1; temp_vec[is.na(temp_vec)] <- 0
    }
    recovered_vec <- temp_vec[1:(length(temp_vec)-1)]

    UMI_counts=numeric(ngenes);
    GI=c(0, cumsum(true_counts_1cell));
    for (i in which(true_counts_1cell>0)){
      x=recovered_vec[(GI[i]+1):GI[i+1]];
      UMI_counts[i]=sum(x>0);
    }

    return(list(UMI_counts, sequenced_vec, sum(frag_vec>0)))
  }
}




#' Title Expand transcript counts to a vector of binaries of the same length of as the number of transcripts
#'
#' @param true_counts_1cell number of transcript in one cell
#'
#' @return list of expanded binary data.
#' @export
#'
#' @examples expand2binary(true_counts_1cell)
expand2binary <- function(true_counts_1cell){
  expanded_vec <- rep(1, sum(true_counts_1cell))
  trans_idx <- sapply(which(true_counts_1cell>0),
                      function(igene){return(rep(igene, true_counts_1cell[igene]))})
  trans_idx <- unlist(trans_idx)
  return(list(expanded_vec, trans_idx))
}


