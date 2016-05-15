Transcriptome Assembly and Evaluation
================
Nate Olson & Leann Biancani
May 2, 2016

Background
----------

#### Transcriptomics

Transcriptomics is the study of the complete set of RNA transcripts in a given sample. Various aims of transcriptomic studies include cataloging all species of transcripts (mRNA, non-coding RNAs, smRNAs), determining gene structure, studying splicing patterns, identifying post-transcriptional modifications, or investigating changes in expression patterns (Wang, Gerstein, and Snyder 2009). Most RNA sequencing (RNA-seq) methods utilize short-read sequencing technologies to sequence a library of cDNA fragments generated from an isolation of RNA molecules. In order to characterize transcripts, these short reads sequences mush be assembled into contiguous sequences (contigs). In order to assemble a transcriptome, reads may either be aligned to a known reference genome or transcriptome sequence, or assembled *de novo*.

#### De novo Transcriptome assembly

*De novo* transcriptome assembly algorithms build contigs from RNA-seq data without the use of a reference genome. These assemblies, however, can be highly variable. There are many different algorithms available and a variety of parameters that can be specified for each method, resulting in the ability to generate multiple assemblies from a single set of reads (Smith-Unna et al. 2015). As such, it is necessary to be able to evaluate the quality and accuracy of a transcriptome assembly.

-   Variation in assemblies - same read data can generate very different assemblies (Smith-Unna et al. 2015)
    -   within assembly method: many parameters & heuristics (to accommodate variation in library construction, coverage depth, organisms)
    -   between assembly methods: many different approaches to assembly problem
        - Oases (Schulz et al. 2012 - not in bib yet)
        - SOAPdenovo-Trans (Xie et al. 2014 - not in bib yet)
        - Trinity (M. G. Grabherr et al. 2011, Haas et al. (2013))

#### Transcriptome assembly evaluation

      Two basic methods exist for evaluating transcriptomes. Reference based methods compare the assembly to known genome or transcriptome sequences, while read based methods map the RNA-seq reads back to the assembly.

      Two basic methods exist for evaluating transcriptomes. Reference based methods compare the assembly to known genome or transcriptome sequences, while read based methods map the RNA-seq reads back to the assembly.

-   Reference based
-   Read based
-   Transrate - contig level quality score (Smith-Unna et al. 2015)

#### C. elegans transcriptome assembly and evaluation (what we did)

Assembly
--------

### Source data

Sequence data identified by searching SRA for *C. elegans* and limiting results to RNA, identified a study for an Illumina HiSeq 2000, 94 bp paired end run. Accession SRR2969230 was used to generate and evaluate the assembly.

> Caenorhabditis elegans Transcriptome or Gene expression BioProject Abstract: Insulin/IGF-1 signaling (IIS) is a critical regulator of the most important biological decisions, from rates of growth, development, and metabolism, to reproduction and longevity. However, the neuron-specific IIS/FOXO targets that regulate these functions are still unknown. By isolating C. elegans neurons for transcriptional profiling, we identified both the wild-type and IIS/FOXO neuronal transcriptomes. Here we deposit the sequences for wt adult neurons and IIS mutants, as well as for isolated adult mechanosensory neurons.

Relevant *C. elegans* Reference Sequences

-   Transcriptome: <ftp://ftp.ensemblgenomes.org/pub/release-31/metazoa/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.31.gtf.gz>
-   Genome: <ftp://ftp.ensemblgenomes.org/pub/release-31/metazoa/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.31.dna.chromosome.*.fa.gz>
-   Proteome: <ftp://ftp.ensemblgenomes.org/pub/release-31/metazoa/fasta/caenorhabditis_elegans/pep/Caenorhabditis_elegans.WBcel235.31.pep.all.fa.gz>

### How to run trinity

-   Link to website
-   Parameter description
-   Command used to generate *C. elegans* assembly

### Assembly results

Assembly metrics calculated using Transrate, see <http://hibberdlab.com/transrate/metrics.html> for a description of the different metrics.
**Metrics** \* n\_seq - number of contigs in the assembly
\* smallest - size of the shortest contig (bp)
\* largest - size of the longest contig (bp)
\* n50 - largest contig size at which at least 50% of bases are contained in contigs at least this length. \* n\_with\_orf - number of contigs with a open reading frame

``` r
assembly_metrics_df <- list(no_mods = "data/no_mods/", 
                            corrected = "data/corrected/",
                            trimmed = "data/trimmed/", 
                            trimmedcorrected = "data/trimmedcorrected/") %>% 
      map(paste0, "assemblies.csv") %>%  
      map_df(read_csv, .id = "read_set") %>% 
      select(-assembly) %>% 
      gather("metric","value", -read_set) %>% 
      mutate(read_set = factor(read_set, levels = c("no_mods","corrected","trimmed","trimmedcorrected")))
```

``` r
assembly_metrics_df %>% 
      filter(metric %in% c("n_seqs", "smallest", "largest", 
                           "n50", "n_with_orf")) %>%
      spread(metric, value) %>% 
      select(read_set, n_seqs, smallest, n50, largest, n_with_orf) %>% 
      kable(digits = 2,
            caption = "Trinity _C. elegans_ transcriptome assembly summary metrics.")
```

| read\_set        |  n\_seqs|  smallest|  n50|  largest|  n\_with\_orf|
|:-----------------|--------:|---------:|----:|--------:|-------------:|
| no\_mods         |    63353|       224|  465|     5446|         12426|
| corrected        |    61059|       224|  467|     4364|         12075|
| trimmed          |    42438|       224|  531|     7533|         10552|
| trimmedcorrected |    42851|       224|  532|     6816|         10645|

Assembly evaluation
-------------------

### Running transrate

#### Installation

-   Precompiled binaries available from website (<http://hibberdlab.com/transrate/>), along with instructions to build from source.
-   Command to install dependencies `transrate --install-deps type` where type can be either `all`, `read`, or `ref`.
-   See Transrate website for additional information for install.

#### Read Based Evaluation

For the read based evaluation Transrate maps the raw read data to the transcriptome assembly then uses a set of four metrics to calculate, contig and assembly quality scores.

**Contig Score:** the product of the following metrics. 1. Edit distance - characterization of contig base accuracy, calculated from read-contig edit distance.
2. Coverage - characterization of based being part of a transcript, calculated from contig based read coverage.
3. Chimera - probability contig is derived from a single transcript, calculated based on whether the contigs fits a single or double Dirichlet distribution. 4. Complete - probability contig is complete and correct, calculated using read pair alignments.

**Assembly Score** is calculated from the mean contig score and proprotion of read pairs mapping to the assembly.

**Command Line**

    transrate \
      --left=SRR2969230_1.fastq \
      --right=SRR2969230_2.fastq \
      --assembly=Trinity_sequences.Trinity.fixed.fasta \
      --output=transrate/

#### Reference Base Evaluation

Transrate uses a reciprical blast approach for reference based assembly evaluation. A sequence is classified as the reference for a contig if the reference sequence and contig are the top BLAST hits when comparing the reference sequence against all contigs and comparing the contig againts all reference sequences. The reciprical BLAST method is described in Aubry et al. (2014).

**Command Line**

    transrate \
      --reference=Caenorhabditis_elegans.WBcel235.31.pep.all.fa \
      --assembly=Trinity_sequences.Trinity.fixed.fasta \
      --output=transrate/

### Transrate results

The unmodified read set was passed as input to Transrate, potentally biasing some metric results towards the unmodifed assembly.

#### Assembly Scores

Read trimming and error correction resulted in the highest weighted assembly score, while having the lowest contig score cutoff, but a lower optimal score.

Weighted assembly score - takes into consideration expression level, can bias assemblies with high proportions of low expressed poorly assembled contigs.

``` r
assembly_metrics_df %>% 
      filter(metric %in% c("score", "optimal_score", "cutoff", "weighted")) %>%
      spread(metric, value) %>% 
      kable(digits = 2,
            caption = "Trinity _C. elegans_ transcriptome Transrate score summary.")
```

| read\_set        |  cutoff|  optimal\_score|  score|  weighted|
|:-----------------|-------:|---------------:|------:|---------:|
| no\_mods         |    0.19|            0.07|   0.02|      0.54|
| corrected        |    0.24|            0.07|   0.02|      0.65|
| trimmed          |    0.22|            0.07|   0.03|      0.67|
| trimmedcorrected |    0.07|            0.06|   0.03|      0.78|

#### Assembly Score Optimization

``` r
assembly_score_opt <- results_list %>% map(paste0,"assembly_score_optimisation.csv") %>% 
      map_df(read_csv, .id = "read_set") %>% 
      mutate(read_set = factor(read_set, levels = c("no_mods","corrected","trimmed","trimmedcorrected")))
```

Relationship between the cutoff for contig score and assembly score.

``` r
assembly_score_opt %>% 
      ggplot() + geom_path(aes(x = cutoff, y = assembly_score, color = read_set)) +
            theme_bw() +
            labs(x = "Contig Score Threshold", y = "Assembly Score")
```

    ## Warning: Removed 2 rows containing missing values (geom_path).

![](README_files/figure-markdown_github/unnamed-chunk-5-1.png)

### Contigs

``` r
contig_stat <- results_list %>% map(paste0,"contigs.csv") %>% 
      map_df(read_csv, .id = "read_set")

contig_cutoff <- assembly_metrics_df %>% filter(metric == "cutoff") %>% 
      select(-metric) %>% rename(cutoff = value)
contig_stat <- contig_stat %>% left_join(contig_cutoff) %>% 
      mutate(contig_filt = cutoff < score)
```

    ## Joining by: "read_set"

    ## Warning in left_join_impl(x, y, by$x, by$y): joining factor and character
    ## vector, coercing into character vector

**Contig Score Distribution** Trimming reads resulted in a higher proportion of contigs with scores greater than 0.25 then assemblies using unmodified reads or only error corrected reads.

``` r
contig_stat %>% 
      ggplot() + geom_density(aes(x = score, color = read_set, fill = read_set), 
                              alpha = 0.25) + theme_bw()
```

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)

Relationship between contig score and length.

``` r
contig_stat %>% 
      ggplot() + geom_hex(aes(x = length, y = score)) +
            geom_hline(aes(yintercept = cutoff), 
                       color = "grey60", linetype = 2) +
            theme_bw() +
            labs(x = "Contig Length (bp)", y = "Contig Score") +
            facet_wrap(~read_set)
```

![2D histogram of contig score and legnth, with color indicating abundance. Grey dotted lines are the contig score cutoff for optimal assembly score.](README_files/figure-markdown_github/score_length-1.png)

Proportion of bases in reads map to the transcript that support the transcript contig.

``` r
contig_stat %>% 
      ggplot() + geom_histogram(aes(x = p_good)) +
            theme_bw() +
            labs(x = "Proportion of Bases Aggree with Transcript",
                 y = "Count") +
            facet_grid(contig_filt~read_set) +
            theme(legend.position = "bottom")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![Histogram of the proportion of bases in agreement with contig and contig length.](README_files/figure-markdown_github/p_good_density-1.png)

Proportion of bases in the contig covered, at least one read aligned to the contig coverage that position.

``` r
contig_stat %>% 
      ggplot() + geom_histogram(aes(x = p_bases_covered)) +
            theme_bw() +
            labs(x = "Proportion of Covered Bases", 
                 y = "Count") +
            facet_grid(contig_filt~read_set) +
            theme(legend.position = "bottom")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](README_files/figure-markdown_github/unnamed-chunk-8-1.png)

Mean per-base contig edit distance, average number of base pair differences between contig and reads.

``` r
contig_stat %>% 
      ggplot() + geom_histogram(aes(x = p_seq_true)) +
            theme_bw() +
            labs(x = "Mean Per-Base Edit Distance", 
                 y = "Count") +
            facet_grid(contig_filt~read_set) +
            theme(legend.position = "bottom")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](README_files/figure-markdown_github/unnamed-chunk-9-1.png)

Probability that a contig is not a chimeric misassembly, where sequences from two different contigs are incorrectly assembled together.

``` r
contig_stat %>% 
      ggplot() + geom_histogram(aes(x = p_not_segmented)) +
            theme_bw() +
            labs(x = "Proportion of Non-Segmented Contigs", 
                 y = "Count") +
            facet_grid(contig_filt~read_set) +
            theme(legend.position = "bottom")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](README_files/figure-markdown_github/unnamed-chunk-10-1.png)

The individual parameters used to calculate the contig score are weakly correlated. Correlation values are similar to thoes presented in the Transrate publication.

``` r
for(i in names(results_list)){
      pl <- contig_stat %>% filter(read_set == i) %>% 
            select(p_good, p_bases_covered, p_seq_true, p_not_segmented) %>% 
            ggcorr(method = c("pairwise","spearman"),label = TRUE) + 
                  ggtitle(i)
      print(pl)
}
```

![](README_files/figure-markdown_github/unnamed-chunk-11-1.png)![](README_files/figure-markdown_github/unnamed-chunk-11-2.png)![](README_files/figure-markdown_github/unnamed-chunk-11-3.png)![](README_files/figure-markdown_github/unnamed-chunk-11-4.png)

#### Reference Based Results

**TODO** Describe CRB method

``` r
ref_contig <- results_list %>% map(paste0,"contigs.csv") %>% 
      map_df(read_csv, .id = "read_set")
```

\_\_TODO\_\_Summary table with assembly, p contigs with CRBB (proportion of contigs iwth CRB-Blast), p\_ref crbb proprotion ref with CRBB, ref coverage - proportion of reference bases covered by CRB hits, collapse factor - high chimeras, p covX - proportion protein coverage, Comprison of contigs with and with out CRB hits. Assemblies using untrimmed read set had a larger number of contigs, but little impact of the proportion of CRB hits.

``` r
ggplot(ref_contig) + geom_bar(aes(x = read_set)) + facet_wrap(~has_crb)
```

![](README_files/figure-markdown_github/unnamed-chunk-13-1.png)

Comprison of contig length and ORF length for contigs with and without CRB hits. Contigs with CRB hits were long and contined larger predicted ORFs.

``` r
ggplot(ref_contig) + geom_boxplot(aes(x = has_crb, y = length)) + facet_wrap(~read_set)
```

![](README_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
ggplot(ref_contig) + geom_boxplot(aes(x = has_crb, y = orf_length)) + facet_wrap(~read_set)
```

![](README_files/figure-markdown_github/unnamed-chunk-15-1.png)

Exploring contigs with CRB hits.

``` r
crb_contig <- ref_contig %>% filter(has_crb == "true")
```

``` r
ggplot(crb_contig) + geom_hex(aes(x = length, y = reference_coverage)) + 
      facet_wrap(~read_set)
```

![](README_files/figure-markdown_github/unnamed-chunk-17-1.png)

``` r
ggplot(crb_contig) + geom_hex(aes(x = orf_length, y = reference_coverage)) +
      facet_wrap(~read_set)
```

![](README_files/figure-markdown_github/unnamed-chunk-18-1.png)

### Contig Score and Hit Coverage

``` r
contig_crb_score <- crb_contig %>% left_join(contig_stat)
```

    ## Joining by: c("read_set", "contig_name", "length", "prop_gc", "gc_skew", "at_skew", "cpg_count", "cpg_ratio", "orf_length", "linguistic_complexity_6")

**TODO** Better way to present relationship

``` r
ggplot(contig_crb_score) + 
      geom_hex(aes(x = reference_coverage, y = score)) +
      facet_wrap(~read_set)
```

![](README_files/figure-markdown_github/unnamed-chunk-20-1.png)

**TODO** Upper limit? correlation between reference coverage and score Looking at tpm

Conclusions
-----------

-   How did Trinity do?
-   How did transrate compare to other methods?
-   Recommendations for transcriptome assembly and evaluation

Session Information
-------------------

**Session Info**

``` r
s_info <- devtools::session_info()
print(s_info$platform)
```

    ##  setting  value                       
    ##  version  R version 3.2.4 (2016-03-10)
    ##  system   x86_64, darwin13.4.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/New_York            
    ##  date     2016-05-15

``` r
kable(s_info$packages)
```

| package    | \*  | version    | date       | source                            |
|:-----------|:----|:-----------|:-----------|:----------------------------------|
| assertthat |     | 0.1        | 2013-12-06 | CRAN (R 3.2.0)                    |
| colorspace |     | 1.2-6      | 2015-03-11 | CRAN (R 3.2.0)                    |
| DBI        |     | 0.4-1      | 2016-05-08 | CRAN (R 3.2.5)                    |
| devtools   |     | 1.11.1     | 2016-04-21 | CRAN (R 3.2.5)                    |
| digest     |     | 0.6.9      | 2016-01-08 | CRAN (R 3.2.3)                    |
| dplyr      | \*  | 0.4.3.9000 | 2016-02-11 | Github (<hadley/dplyr@9bae2aa>)   |
| evaluate   |     | 0.9        | 2016-04-29 | CRAN (R 3.2.5)                    |
| formatR    |     | 1.4        | 2016-05-09 | CRAN (R 3.2.4)                    |
| GGally     | \*  | 1.0.1      | 2016-01-14 | CRAN (R 3.2.3)                    |
| ggplot2    | \*  | 2.1.0      | 2016-03-01 | CRAN (R 3.2.4)                    |
| gtable     |     | 0.2.0      | 2016-02-26 | CRAN (R 3.2.3)                    |
| hexbin     | \*  | 1.27.1     | 2015-08-19 | CRAN (R 3.2.0)                    |
| highr      |     | 0.6        | 2016-05-09 | CRAN (R 3.2.4)                    |
| htmltools  |     | 0.3.5      | 2016-03-21 | CRAN (R 3.2.4)                    |
| knitr      | \*  | 1.13       | 2016-05-09 | CRAN (R 3.2.4)                    |
| labeling   |     | 0.3        | 2014-08-23 | CRAN (R 3.2.0)                    |
| lattice    |     | 0.20-33    | 2015-07-14 | CRAN (R 3.2.4)                    |
| lazyeval   |     | 0.1.10     | 2015-01-02 | CRAN (R 3.2.0)                    |
| magrittr   |     | 1.5        | 2014-11-22 | CRAN (R 3.2.0)                    |
| memoise    |     | 1.0.0      | 2016-01-29 | CRAN (R 3.2.3)                    |
| munsell    |     | 0.4.3      | 2016-02-13 | CRAN (R 3.2.3)                    |
| plyr       |     | 1.8.3      | 2015-06-12 | CRAN (R 3.2.0)                    |
| purrr      | \*  | 0.2.1      | 2016-02-13 | CRAN (R 3.2.3)                    |
| R6         |     | 2.1.2      | 2016-01-26 | CRAN (R 3.2.3)                    |
| Rcpp       |     | 0.12.4     | 2016-03-26 | CRAN (R 3.2.4)                    |
| readr      | \*  | 0.2.2      | 2015-10-22 | CRAN (R 3.2.0)                    |
| reshape    |     | 0.8.5      | 2014-04-23 | CRAN (R 3.2.0)                    |
| reshape2   |     | 1.4.1      | 2014-12-06 | CRAN (R 3.2.0)                    |
| rmarkdown  |     | 0.9.6      | 2016-05-01 | CRAN (R 3.2.5)                    |
| scales     |     | 0.4.0      | 2016-02-26 | CRAN (R 3.2.3)                    |
| stringi    |     | 1.0-1      | 2015-10-22 | CRAN (R 3.2.1)                    |
| stringr    | \*  | 1.0.0.9000 | 2015-11-23 | Github (<hadley/stringr@a67f8f0>) |
| tidyr      | \*  | 0.4.1      | 2016-02-05 | CRAN (R 3.2.3)                    |
| withr      |     | 1.0.1      | 2016-02-04 | CRAN (R 3.2.3)                    |
| yaml       |     | 2.1.13     | 2014-06-12 | CRAN (R 3.2.0)                    |

References
----------

Aubry, Sylvain, Steven Kelly, Britta MC Kümpers, Richard D Smith-Unna, and Julian M Hibberd. 2014. “Deep Evolutionary Comparison of Gene Expression Identifies Parallel Recruitment of Trans-Factors in Two Independent Origins of c 4 Photosynthesis.” *PLoS Genet* 10 (6). Public Library of Science: e1004365.

Grabherr, Manfred G, Brian J Haas, Moran Yassour, Joshua Z Levin, Dawn A Thompson, Ido Amit, Xian Adiconis, et al. 2011. “Full-length transcriptome assembly from RNA-Seq data without a reference genome.” *Nature Biotechnology* 29 (7): 644–52.

Haas, Brian J, Alexie Papanicolaou, Moran Yassour, Manfred Grabherr, Philip D Blood, Joshua Bowden, Matthew Brian Couger, et al. 2013. “De Novo Transcript Sequence Reconstruction from RNA-Seq Using the Trinity Platform for Reference Generation and Analysis.” *Nature Protocols* 8 (8). Nature Publishing Group: 1494–1512.

Smith-Unna, Richard D, Chris Boursnell, Rob Patro, Julian M Hibberd, and Steven Kelly. 2015. “TransRate: Reference Free Quality Assessment of de-Novo Transcriptome Assemblies.” *BioRxiv*. Cold Spring Harbor Labs Journals, 021626.

Wang, Zhong, Mark Gerstein, and Michael Snyder. 2009. “RNA-Seq: a revolutionary tool for transcriptomics.” *Nature Reviews. Genetics* 10 (1): 57–63.
