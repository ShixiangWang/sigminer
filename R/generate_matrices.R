
# SBS ---------------------------------------------------------------------

generate_matrix_SBS <- function(query, ref_genome, genome_build = "hg19", add_trans_bias = FALSE) {
  ## TODO make 4 transcriptional bias categories all right
  query <- query[query$Variant_Type == "SNP"]
  if (nrow(query) == 0) {
    stop("Zero SNPs to analyze!")
  }

  extract.tbl <- data.table::data.table(
    Chromosome = query$Chromosome, Start = query$Start_Position, End = query$End_Position,
    Reference_Allele = query$Reference_Allele, Tumor_Seq_Allele2 = query$Tumor_Seq_Allele2,
    Tumor_Sample_Barcode = query$Tumor_Sample_Barcode, upstream = query$Start_Position - 20,
    downstream = query$End_Position + 20
  )

  message("=> Extracting 5' and 3' adjacent bases")
  ss <- BSgenome::getSeq(
    x = ref_genome,
    names = extract.tbl$Chromosome,
    start = extract.tbl$Start - 2,
    end = extract.tbl$End + 2,
    as.character = TRUE
  )

  message("=> Extracting +/- 20bp around mutated bases for background C>T estimation")
  updwn <- BSgenome::getSeq(
    x = ref_genome, names = extract.tbl$Chromosome, start = extract.tbl$upstream,
    end = extract.tbl$downstream, as.character = FALSE
  )
  updwn.alphFreq <- data.table::as.data.table(BSgenome::alphabetFrequency(x = updwn))[, c("A", "T", "G", "C")] # Nucleotide frequency
  updwn.tnmFreq <- data.table::as.data.table(Biostrings::trinucleotideFrequency(x = updwn, step = 1))

  extract.tbl$pentanucleotide <- as.character(ss)
  extract.tbl$trinucleotide <- substr(extract.tbl$pentanucleotide, 2, 4)
  extract.tbl$updown <- as.character(updwn)

  extract.tbl <- cbind(extract.tbl, updwn.alphFreq)
  extract.tbl <- cbind(extract.tbl, updwn.tnmFreq[, c("TCA", "TCT", "AGA", "TGA")])
  extract.tbl[, tcw := rowSums(extract.tbl[, c("TCA", "TCT")])]
  extract.tbl[, wga := rowSums(extract.tbl[, c("TGA", "AGA")])]

  ## All combinations
  extract.tbl[, Substitution := paste(extract.tbl$Reference_Allele, extract.tbl$Tumor_Seq_Allele2, sep = ">")]
  extract.tbl$SubstitutionMotif <- paste0(
    substr(x = as.character(extract.tbl$pentanucleotide), 1, 2),
    "[", extract.tbl$Substitution, "]",
    substr(as.character(extract.tbl$pentanucleotide), 4, 5)
  )
  extract.tbl$TriSubstitutionMotif <- substr(extract.tbl$SubstitutionMotif, 2, 8)

  # substitutions are referred to by the pyrimidine of the mutated Watson-Crick base pair
  conv <- c("T>C", "T>C", "C>T", "C>T", "T>A", "T>A", "T>G", "T>G", "C>A", "C>A", "C>G", "C>G")
  names(conv) <- c("A>G", "T>C", "C>T", "G>A", "A>T", "T>A", "A>C", "T>G", "C>A", "G>T", "C>G", "G>C")
  complement <- c("A", "C", "G", "T")
  names(complement) <- c("T", "G", "C", "A")


  extract.tbl$SubstitutionType <- conv[extract.tbl$Substitution]
  # need to reverse-complement triplet for mutated purines (not just the middle base)
  extract.tbl$should_reverse <- extract.tbl$Substitution != extract.tbl$SubstitutionType
  extract.tbl$SubstitutionTypeMotif <- ifelse(extract.tbl$should_reverse,
    paste0(
      complement[substr(x = extract.tbl$pentanucleotide, 5, 5)],
      complement[substr(x = extract.tbl$pentanucleotide, 4, 4)],
      "[", extract.tbl$SubstitutionType, "]",
      complement[substr(x = extract.tbl$pentanucleotide, 2, 2)],
      complement[substr(x = extract.tbl$pentanucleotide, 1, 1)]
    ),
    paste0(
      substr(x = as.character(extract.tbl$pentanucleotide), 1, 2),
      "[", extract.tbl$SubstitutionType, "]",
      substr(as.character(extract.tbl$pentanucleotide), 4, 5)
    )
  )
  extract.tbl$TriSubstitutionTypeMotif <- substr(extract.tbl$SubstitutionTypeMotif, 2, 8)


  # Possible substitution types after being referred to by the pyrimidine of the mutated Watson-Crick base pair

  # penta_comb <- expand.grid(
  #   complement,
  #   complement,
  #   "[",
  #   unique(as.character(conv)),
  #   "]",
  #   complement,
  #   complement,
  #   stringsAsFactors = FALSE
  # ) %>%
  #   apply(1, paste0, collapse = "") %>%
  #   unique()

  penta_comb <- vector_to_combination(
    complement,
    complement,
    "[",
    unique(as.character(conv)),
    "]",
    complement,
    complement
  )

  tri_comb <- substr(penta_comb, 2, 8) %>%
    unique()

  tri_comb2 <- vector_to_combination(
    complement,
    "[",
    unique(c(as.character(conv), names(conv))),
    "]",
    complement
  )

  # Set levels for type (mainly component)
  extract.tbl$SubstitutionType <- factor(extract.tbl$SubstitutionType, levels = unique(as.character(conv)))
  extract.tbl$TriSubstitutionTypeMotif <- factor(extract.tbl$TriSubstitutionTypeMotif, levels = tri_comb)
  extract.tbl$SubstitutionTypeMotif <- factor(extract.tbl$SubstitutionTypeMotif, levels = penta_comb)

  extract.tbl$TriSubstitutionMotif <- factor(extract.tbl$TriSubstitutionMotif, levels = tri_comb2)
  extract.tbl$Substitution <- factor(extract.tbl$Substitution, levels = unique(c(as.character(conv), names(conv))))

  # Compile data
  ## This is nucleotide frequcny and motif frequency across 41 bp in C>T and C>G context.
  apobecSummary <- extract.tbl[
    as.character(extract.tbl$SubstitutionType) %in% c("C>T", "C>G"),
    .(
      A = sum(A), T = sum(T), G = sum(G), C = sum(C), tcw = sum(tcw),
      wga = sum(wga), bases = sum(A, T, G, C)
    ), Tumor_Sample_Barcode
  ]
  ## The by operation may remove some samples here
  ## should it be added??

  ## This is per sample conversion events
  sub.tbl <- extract.tbl[, .N, list(Tumor_Sample_Barcode, Substitution)]
  sub.tbl <- data.table::dcast(data = sub.tbl, formula = Tumor_Sample_Barcode ~ Substitution, fill = 0, value.var = "N", drop = FALSE)
  sub.tbl[, n_A := rowSums(sub.tbl[, c("A>C", "A>G", "A>T")], na.rm = TRUE)][, n_T := rowSums(sub.tbl[, c("T>A", "T>C", "T>G")], na.rm = TRUE)][, n_G := rowSums(sub.tbl[, c("G>A", "G>C", "G>T")], na.rm = TRUE)][, n_C := rowSums(sub.tbl[, c("C>A", "C>G", "C>T")], na.rm = TRUE)]
  sub.tbl[, n_mutations := rowSums(sub.tbl[, c("n_A", "n_T", "n_G", "n_C")], na.rm = TRUE)]
  sub.tbl[, "n_C>G_and_C>T" := rowSums(sub.tbl[, c("C>G", "G>C", "C>T", "G>A")], na.rm = TRUE)] # number of APOBEC type mutations (C>G and C>T type)

  ## This is per substitution type events
  subType.tbl <- extract.tbl[, .N, .(Tumor_Sample_Barcode, TriSubstitutionMotif)]
  subType.tbl <- data.table::dcast(data = subType.tbl, formula = Tumor_Sample_Barcode ~ TriSubstitutionMotif, fill = 0, value.var = "N", drop = FALSE)

  ### tCw events
  subType.tbl[, tCw_to_A := rowSums(subType.tbl[, .(`T[C>A]A`, `T[C>A]T`)], na.rm = TRUE)]
  subType.tbl[, tCw_to_G := rowSums(subType.tbl[, .(`T[C>G]A`, `T[C>G]T`)], na.rm = TRUE)]
  subType.tbl[, tCw_to_T := rowSums(subType.tbl[, .(`T[C>T]A`, `T[C>T]T`)], na.rm = TRUE)]
  subType.tbl[, tCw := rowSums(subType.tbl[, .(tCw_to_A, tCw_to_G, tCw_to_T)], na.rm = TRUE)]

  ### wGa events
  subType.tbl[, wGa_to_C := rowSums(subType.tbl[, .(`A[G>C]A`, `T[G>C]A`)], na.rm = TRUE)]
  subType.tbl[, wGa_to_T := rowSums(subType.tbl[, .(`A[G>T]A`, `T[G>T]A`)], na.rm = TRUE)]
  subType.tbl[, wGa_to_A := rowSums(subType.tbl[, .(`A[G>A]A`, `T[G>A]A`)], na.rm = TRUE)]
  subType.tbl[, wGa := rowSums(subType.tbl[, .(wGa_to_C, wGa_to_T, wGa_to_A)], na.rm = TRUE)]

  ## tCw_to_G+tCw_to_T
  subType.tbl[, "tCw_to_G+tCw_to_T" := rowSums(subType.tbl[, .(`T[C>G]T`, `T[C>G]A`, `T[C>T]T`, `T[C>T]A`, `T[G>C]A`, `A[G>C]A`, `T[G>A]A`, `A[G>A]A`)], na.rm = TRUE)]

  ### Merge data
  sub.tbl <- merge(sub.tbl, subType.tbl[, .(
    tCw_to_A, tCw_to_T, tCw_to_G, tCw, wGa_to_C, wGa_to_T, wGa_to_A,
    wGa, `tCw_to_G+tCw_to_T`, Tumor_Sample_Barcode
  )],
  by = "Tumor_Sample_Barcode"
  )
  sub.tbl <- merge(sub.tbl, apobecSummary, by = "Tumor_Sample_Barcode")

  ### Estimate APOBEC enrichment
  sub.tbl[, APOBEC_Enrichment := (`tCw_to_G+tCw_to_T` / `n_C>G_and_C>T`) / ((tcw + wga) / (C + G))]
  sub.tbl[, non_APOBEC_mutations := n_mutations - `tCw_to_G+tCw_to_T`]
  sub.tbl[, fraction_APOBEC_mutations := round((n_mutations - non_APOBEC_mutations) / n_mutations, digits = 3)]
  data.table::setDF(sub.tbl)

  message("=> Estimating APOBEC enrichment scores")
  apobec.fisher.dat <- sub.tbl[, c(19, 28, 32, 33, 34)]
  if (nrow(apobec.fisher.dat) == 1) {
    apobec.fisher.dat <- t(as.matrix(apply(X = apobec.fisher.dat, 2, as.numeric)))
  } else {
    apobec.fisher.dat <- apply(X = apobec.fisher.dat, 2, as.numeric)
  }

  ### One way Fisher test to estimate over representation og APOBEC associated tcw mutations
  message("=> Performing one-way Fisher's test for APOBEC enrichment")
  sub.tbl <- cbind(sub.tbl, data.table::rbindlist(apply(X = apobec.fisher.dat, 1, function(x) {
    xf <- fisher.test(matrix(c(x[2], sum(x[3], x[4]), x[1] - x[2], x[3] - x[4]), nrow = 2), alternative = "g")
    data.table::data.table(fisher_pvalue = xf$p.value, or = xf$estimate, ci.up = xf$conf.int[1], ci.low = xf$conf.int[2])
  })))

  data.table::setDT(sub.tbl)
  colnames(sub.tbl)[29:35] <- paste0("n_bg_", colnames(sub.tbl)[29:35])
  sub.tbl <- sub.tbl[order(sub.tbl$fisher_pvalue)]

  ## Choosing APOBEC Enrichment scores > 2 as cutoff
  sub.tbl$APOBEC_Enriched <- ifelse(test = sub.tbl$APOBEC_Enrichment > 2, yes = "yes", no = "no")
  sub.tbl[, fdr := p.adjust(sub.tbl$fisher_pvalue, method = "fdr")] # Adjusted p-values

  message(
    paste0("=> APOBEC related mutations are enriched in "),
    round(nrow(sub.tbl[APOBEC_Enriched %in% "yes"]) / nrow(sub.tbl) * 100, digits = 3),
    "% of samples (APOBEC enrichment score > 2; ",
    nrow(sub.tbl[APOBEC_Enriched %in% "yes"]), " of ", nrow(sub.tbl), " samples)"
  )

  message("=> Creating SBS sample-by-component matrice")

  message("==> Creating SBS-6 matrix")
  SBS_6 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "SubstitutionType")
  SBS_6 <- SBS_6[, c("T>C", "C>T", "T>A", "T>G", "C>A", "C>G")] %>% as.matrix()

  message("==> Creating SBS-96 matrix")
  SBS_96 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "TriSubstitutionTypeMotif")
  SBS_96 <- SBS_96[, tri_comb] %>% as.matrix()

  message("==> Creating SBS-1536 matrix")
  SBS_1536 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "SubstitutionTypeMotif")
  SBS_1536 <- SBS_1536[, penta_comb] %>% as.matrix()

  if (add_trans_bias) {
    t_labels <- c("T:", "U:", "B:", "N:")

    message("==> Creating SBS-24 (6x4) matrix")
    SBS_24 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "SubstitutionType",
      add_trans_bias = TRUE, build = genome_build
    )
    SBS_24 <- SBS_24[, vector_to_combination(
      t_labels,
      c("T>C", "C>T", "T>A", "T>G", "C>A", "C>G")
    )] %>% as.matrix()

    message("==> Creating SBS-384 (96x4) matrix")
    SBS_384 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "TriSubstitutionTypeMotif",
      add_trans_bias = TRUE, build = genome_build
    )
    SBS_384 <- SBS_384[, vector_to_combination(t_labels, tri_comb)] %>% as.matrix()

    message("==> Creating SBS-6144 (1536x4) matrix")
    SBS_6144 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "SubstitutionTypeMotif",
      add_trans_bias = TRUE, build = genome_build
    )
    SBS_6144 <- SBS_6144[, vector_to_combination(t_labels, penta_comb)] %>% as.matrix()
  }

  if (add_trans_bias) {
    message("=> Return SBS-384 as major matrix")
    res <- list(
      nmf_matrix = SBS_384,
      all_matrice = list(
        SBS_6 = SBS_6,
        SBS_24 = SBS_24,
        SBS_96 = SBS_96,
        SBS_384 = SBS_384,
        SBS_1536 = SBS_1536,
        SBS_6144 = SBS_6144
      ),
      APOBEC_scores = sub.tbl
    )
  } else {
    message("=> Return SBS-96 as major matrix")
    res <- list(
      nmf_matrix = SBS_96,
      all_matrice = list(
        SBS_6 = SBS_6,
        SBS_96 = SBS_96,
        SBS_1536 = SBS_1536
      ),
      APOBEC_scores = sub.tbl
    )
  }
  res
}


records_to_matrix <- function(dt, samp_col, component_col, add_trans_bias = FALSE, build = "hg19") {
  if (add_trans_bias) {
    transcript_dt <- get(paste0("transcript.", build), envir = as.environment("package:sigminer"))
    data.table::setkey(transcript_dt, chrom, start, end)

    loc_dt <- dt[, .(Chromosome, Start, End)]
    loc_dt$MutIndex <- 1:nrow(loc_dt)
    colnames(loc_dt)[1:3] <- c("chrom", "start", "end")
    m_dt <- data.table::foverlaps(loc_dt, transcript_dt, type = "any")[
      , .(MatchCount = .N, strand = paste0(unique(strand), collapse = "/")),
      by = MutIndex
    ]
    ## Actually, the MatchCount should only be 0, 1, 2
    if (any(m_dt$MatchCount > 2)) {
      stop("More than 2 regions counted, please report your data and code to developer!")
    }

    dt$transcript_bias_label <- ifelse(
      m_dt$MatchCount == 2, "B:",
      ifelse(m_dt$MatchCount == 0, "N:",
        ifelse(xor(dt$should_reverse, m_dt$strand == "-"), ## Need expand the logical to DBS and ID
          # (dt$should_reverse & m_dt$strand == "+") | (!dt$should_reverse & m_dt$strand == "-"),
          "T:", "U:"
        )
      ) ## If should reverse, the base switch to template strand from coding strand
    )
    new_levels <- vector_to_combination(
      c("T:", "U:", "B:", "N:"),
      levels(dt[[component_col]])
    )
    dt[[component_col]] <- paste0(dt$transcript_bias_label, dt[[component_col]])
    dt[[component_col]] <- factor(dt[[component_col]], levels = new_levels)
  }

  dt.summary <- dt[, .N, by = c(samp_col, component_col)]
  mat <- as.data.frame(data.table::dcast(dt.summary,
    formula = as.formula(paste(samp_col, "~", component_col)),
    fill = 0, value.var = "N", drop = FALSE
  ))

  rownames(mat) <- mat[, 1]
  mat <- mat[, -1]
  mat
}



# DBS ---------------------------------------------------------------------

generate_matrix_DBS <- function(query, ref_genome, genome_build = "hg19", add_trans_bias = FALSE) {
  query <- query[query$Variant_Type == "SNP"]
  if (nrow(query) == 0) {
    stop("Zero DBSs to analyze!")
  }
  sample <- as.data.frame(table(query$Tumor_Sample_Barcode))
  sample$Var1 <- as.character(sample$Var1)
  matrix <- data.frame(mutation_type = c(
    "AC>CA", "AC>CG", "AC>CT", "AC>GA", "AC>GG", "AC>GT",
    "AC>TA", "AC>TG", "AC>TT", "AT>CA", "AT>CC", "AT>CG",
    "AT>GA", "AT>GC", "AT>TA", "CC>AA", "CC>AG", "CC>AT",
    "CC>GA", "CC>GG", "CC>GT", "CC>TA", "CC>TG", "CC>TT",
    "CG>AT", "CG>GC", "CG>GT", "CG>TA", "CG>TC", "CG>TT",
    "CT>AA", "CT>AC", "CT>AG", "CT>GA", "CT>GC", "CT>GG",
    "CT>TA", "CT>TC", "CT>TG", "GC>AA", "GC>AG", "GC>AT",
    "GC>CA", "GC>CG", "GC>TA", "TA>AT", "TA>CG", "TA>CT",
    "TA>GC", "TA>GG", "TA>GT", "TC>AA", "TC>AG", "TC>AT",
    "TC>CA", "TC>CG", "TC>CT", "TC>GA", "TC>GG", "TC>GT",
    "TG>AA", "TG>AC", "TG>AT", "TG>CA", "TG>CC", "TG>CT",
    "TG>GA", "TG>GC", "TG>GT", "TT>AA", "TT>AC", "TT>AG",
    "TT>CA", "TT>CC", "TT>CG", "TT>GA", "TT>GC", "TT>GG"
  ), reverse = c(
    "GT>TG", "GT>CG", "GT>AG", "GT>TC", "GT>CC", "GT>AC",
    "GT>TA", "GT>CA", "GT>AA", "AT>TG", "AT>GG", "AT>CG",
    "AT>TC", "AT>GC", "AT>TA", "GG>TT", "GG>CT", "GG>AT",
    "GG>TC", "GG>CC", "GG>AC", "GG>TA", "GG>CA", "GG>AA",
    "CG>AT", "CG>GC", "CG>AC", "CG>TA", "CG>GA", "CG>AA",
    "AG>TT", "AG>GT", "AG>CT", "AG>TC", "AG>GC", "AG>CC",
    "AG>TA", "AG>GA", "AG>CA", "GC>TT", "GC>CT", "GC>AT",
    "GC>TG", "GC>CG", "GC>TA", "TA>AT", "TA>CG", "TA>AG",
    "TA>GC", "TA>CC", "TA>AC", "GA>TT", "GA>CT", "GA>AT",
    "GA>TG", "GA>CG", "GA>AG", "GA>TC", "GA>CC", "GA>AC",
    "CA>TT", "CA>GT", "CA>AT", "CA>TG", "CA>GG", "CA>AG",
    "CA>TC", "CA>GC", "CA>AC", "AA>TT", "AA>GT", "AA>CT",
    "AA>TG", "AA>GG", "AA>CG", "AA>TC", "AA>GC", "AA>CC"
  ), stringsAsFactors = F)
  matrix_sb <- matrix %>%
    mutate(u = "U", N = "N", T = "T", B = "B", Q = "Q") %>%
    pivot_longer(c("u", "N", "T", "B", "Q"), names_to = "type", values_to = "value") %>%
    select(-type) %>%
    mutate(mutation_type = paste(value, mutation_type, sep = ":"), reverse = paste(value, reverse, sep = ":")) %>%
    filter((!substr(mutation_type, 3, 4) %in% c("TC", "TT", "CT", "CC") & substr(mutation_type, 1, 1) == "Q") |
      substr(mutation_type, 3, 4) %in% c("TC", "TT", "CT", "CC") & substr(mutation_type, 1, 1) != "Q") %>%
    select(-value)

  query <- query %>% mutate(
    index = paste(Chromosome, Start_Position, Tumor_Sample_Barcode, sep = ":"),
    index_up = paste(Chromosome, (Start_Position - 1), Tumor_Sample_Barcode, sep = ":"),
    index_down = paste(Chromosome, (Start_Position + 1), Tumor_Sample_Barcode, sep = ":")
  )



  ## 找上游的双碱基突变
  query <- query %>%
    mutate(
      source_up_p = chmatch(index_up, index),
      source_down_p = chmatch(index_down, index),
      source_up = ifelse(is.na(source_up_p), NA, paste(query$Reference_Allele[source_up_p], Reference_Allele, sep = "")),
      mutation_up = ifelse(is.na(source_up_p), NA, paste(query$Tumor_Seq_Allele2[source_up_p], Tumor_Seq_Allele2, sep = "")),
      source_down = ifelse(is.na(source_down_p), NA, paste(Reference_Allele, query$Reference_Allele[source_down_p], sep = "")),
      mutation_down = ifelse(is.na(source_down_p), NA, paste(Tumor_Seq_Allele2, query$Tumor_Seq_Allele2[source_down_p], sep = ""))
    )

  ## 将匹配的上游的突变位点打个标签，后面找下游突变的时候就不找有这些标签的突变
  tag_p <- query$source_up_p
  tag_p <- na.omit(tag_p)
  query$index[tag_p] <- "link"

  ## 去掉下游突变中和上游重复的
  query$source_down <- ifelse(query$index == "link", NA, query$source_down)
  query$mutation_down <- ifelse(query$index == "link", NA, query$mutation_down)
  query <- query %>% filter(!(is.na(source_up) & is.na(source_down) & is.na(mutation_down) & is.na(mutation_up)))
  query <- query %>%
    mutate(
      mutation_type_up = ifelse(is.na(source_up), NA, paste(source_up, mutation_up, sep = ">")),
      mutation_type_down = ifelse(is.na(source_down), NA, paste(source_down, mutation_down, sep = ">"))
    )

  ## 将互补的碱基转化

  query$mutation_type_up <- sapply(
    query$mutation_type_up,
    function(x) {
      if (x %in% matrix$reverse) {
        a <- which(matrix$reverse == x)
        return(matrix$mutation_type[a])
      } else {
        return(x)
      }
    }
  )
  query$mutation_type_down <- sapply(
    query$mutation_type_down,
    function(x) {
      if (x %in% matrix$reverse) {
        a <- which(matrix$reverse == x)
        return(matrix$mutation_type[a])
      } else {
        return(x)
      }
    }
  )
  ## 合并
  query <- setDT(query)
  a <- query[, .N, by = .(mutation_type_up, Tumor_Sample_Barcode)] %>% as.data.frame()
  colnames(a) <- c("mutation_type", "sample", "n")
  b <- query[, .N, by = .(mutation_type_down, Tumor_Sample_Barcode)] %>% as.data.frame()
  colnames(b) <- c("mutation_type", "sample", "n")
  all <- rbind(a, b)
  all <- na.omit(all)

  matrix[, 3:(nrow(sample) + 2)] <- NA
  colnames(matrix)[3:(nrow(sample) + 2)] <- sample$Var1

  ## 生成矩阵
  for (i in 3:(nrow(sample) + 2)) {
    a <- which(all$sample == colnames(matrix)[i])
    b <- all[a, ]
    if (length(b) == 0) {
      matrix[, i] <- 0
    }
    else {
      matrix[, i] <- sapply(
        matrix$mutation_type,
        function(x) {
          c <- which(b$mutation_type == x)
          if (length(c) == 0) {
            return(0)
          }
          else {
            return(sum(b[c, 3]))
          }
        }
      )
    }
  }

  ## 加入strand bias
  transcript_dt <- get(paste0("transcript.", genome_build), envir = as.environment("package:sigminer"))
  dt <- query[, .(x = Chromosome, pos = Start_Position)] %>% distinct(x, pos)

  plus_strand <- transcript_dt[strand == "+", .(x = chrom, start, end, strand)]
  minus_strand <- transcript_dt[strand == "-", .(x = chrom, start, end, strand)]
  overlap_plus <- dt[plus_strand, .(x, pos = x.pos, start, end, strand), on = .(x, pos >= start, pos <= end), nomatch = F, allow.cartesian = TRUE][, index := paste(x, pos, sep = ":")] %>%
    distinct(index, .keep_all = T)
  overlap_minus <- dt[minus_strand, .(x, pos = x.pos, start, end, strand), on = .(x, pos >= start, pos <= end), nomatch = F, allow.cartesian = TRUE][, index := paste(x, pos, sep = ":")] %>%
    distinct(index, .keep_all = T)
  intersect_position <- intersect(overlap_minus$index, overlap_plus$index)
  query <- query[, index_p := paste(Chromosome, Start_Position, sep = ":")]

  query <- query %>%
    mutate(
      bias_type_up = ifelse(index_p %in% overlap_plus$index & substr(mutation_type_up, 1, 2) %in% c("TC", "TT", "CT", "CC"), "U",
        ifelse(index_p %in% overlap_plus$index & !(substr(mutation_type_up, 1, 2) %in% c("TC", "TT", "CT", "CC")), "Q",
          ifelse(index_p %in% overlap_minus$index & substr(mutation_type_up, 1, 2) %in% c("TC", "TT", "CT", "CC"), "T",
            ifelse(index_p %in% overlap_minus$index & !(substr(mutation_type_up, 1, 2) %in% c("TC", "TT", "CT", "CC")), "Q", "N")
          )
        )
      ),
      bias_type_down = ifelse(index_p %in% overlap_plus$index & substr(mutation_type_down, 1, 2) %in% c("TC", "TT", "CT", "CC"), "U",
        ifelse(index_p %in% overlap_plus$index & !(substr(mutation_type_down, 1, 2) %in% c("TC", "TT", "CT", "CC")), "Q",
          ifelse(index_p %in% overlap_minus$index & substr(mutation_type_down, 1, 2) %in% c("TC", "TT", "CT", "CC"), "T",
            ifelse(index_p %in% overlap_minus$index & !(substr(mutation_type_down, 1, 2) %in% c("TC", "TT", "CT", "CC")), "Q", "N")
          )
        )
      )
    )

  query <- query %>%
    mutate(
      bias_type_up = ifelse(index_p %in% intersect_position & substr(mutation_type_up, 1, 2) %in% c("TC", "TT", "CT", "CC"), "B", bias_type_up),
      bias_type_down = ifelse(index_p %in% intersect_position & substr(mutation_type_down, 1, 2) %in% c("TC", "TT", "CT", "CC"), "B", mutation_type_down)
    )
  query <- query %>%
    mutate(
      bias_type_up = ifelse(is.na(bias_type_up), NA, paste(bias_type_up, mutation_type_up, sep = ":")),
      bias_type_down = ifelse(is.na(bias_type_down), NA, paste(bias_type_down, mutation_type_down, sep = ":"))
    )

  ## 合并
  query <- setDT(query)
  c <- query[, .N, by = .(bias_type_up, Tumor_Sample_Barcode)] %>% as.data.frame()
  colnames(c) <- c("mutation_type", "sample", "n")
  d <- query[, .N, by = .(bias_type_down, Tumor_Sample_Barcode)] %>% as.data.frame()
  colnames(d) <- c("mutation_type", "sample", "n")
  all2 <- rbind(c, d)
  all2 <- na.omit(all2)

  matrix_sb[, 3:(nrow(sample) + 2)] <- NA
  colnames(matrix_sb)[3:(nrow(sample) + 2)] <- sample$Var1

  ## 生成矩阵
  for (i in 3:(nrow(sample) + 2)) {
    e <- which(all2$sample == colnames(matrix_sb)[i])
    k <- all2[e, ]
    if (length(k) == 0) {
      matrix_sb[, i] <- 0
    }
    else {
      matrix_sb[, i] <- sapply(
        matrix_sb$mutation_type,
        function(x) {
          j <- which(k$mutation_type == x)
          if (length(j) == 0) {
            return(0)
          }
          else {
            return(sum(k[j, 3]))
          }
        }
      )
    }
  }

  matrix <- matrix[-2]
  matrix_sb <- matrix_sb[-2]

  if (add_trans_bias) {
    all_m <- list(DBS_78 = matrix, DBS_186 = matrix_sb)
  } else {
    all_m <- matrix
  }

  all_m
}


# INDELs (ID) -------------------------------------------------------------

generate_matrix_INDEL <- function(query, ref_genome, genome_build = "hg19", add_trans_bias = FALSE) {
  query <- query[query$Variant_Type != "SNP"]
  if (nrow(query) == 0) {
    stop("Zero INDELs to analyze!")
  }
}

# Utils -------------------------------------------------------------------

vector_to_combination <- function(...) {
  expand.grid(
    ...,
    stringsAsFactors = FALSE
  ) %>%
    apply(1, paste0, collapse = "") %>%
    unique()
}


utils::globalVariables(
  c(
    "A", "APOBEC_Enriched", "APOBEC_Enrichment",
    "A[G>A]A", "A[G>C]A", "A[G>T]A", "G",
    "Substitution",
    "T[C>A]A", "T[C>A]T",
    "T[C>G]A", "T[C>G]T",
    "T[C>T]A", "T[C>T]T", "T[G>A]A",
    "T[G>C]A", "T[G>T]A",
    "TriSubstitutionMotif", "Tumor_Sample_Barcode", "fdr",
    "fraction_APOBEC_mutations", "n_A", "n_C", "n_C>G_and_C>T",
    "n_G", "n_T", "n_mutations", "non_APOBEC_mutations",
    "tCw", "tCw_to_A", "tCw_to_G", "tCw_to_G+tCw_to_T",
    "tCw_to_T", "tcw", "wGa", "wGa_to_A", "wGa_to_C", "wGa_to_T", "wga",
    "Start", "End", "MutIndex"
  )
)
