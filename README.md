# Trinity and Transrate Tutorial
Nate Olson & Leann Biancani  
May 2, 2016  




## Background 
* Transcriptome
      * RNAseq & transcriptomics [@Wang:2009di]
* De novo Transcriptome assembly
      * Trinity [@Grabherr:2011jb, @haas2013novo] 
      * Mention other de-novo assembly algorithms  
            * Oases (Schulz et al. 2012 - not in bib yet)  
            * SOAPdenovo-Trans (Xie et al. 2014 - not in bib yet)  
      * Variation in assemblies - same read data can generate very different assemblies [@smith2015transrate]  
            * withing assembly methods: many parameters & heuristics to accommodated  
* Transcriptome assembly evaluation
      * Reference based
      * Read based
      * Transrate - contig level quality score [@smith2015transrate]
* C. elegans transcriptome assembly and evaluation (what we did)

## Assembly
### Source data  
Sequence data identified by searching SRA for _C. elegans_ and limiting results to RNA, identified a study for an Illumina HiSeq 2000, 94 bp paired end run. Accession SRR2969230 was used to generate and evaluate the assembly.

> 
Caenorhabditis elegans Transcriptome or Gene expression BioProject
Abstract: Insulin/IGF-1 signaling (IIS) is a critical regulator of the most important biological decisions, from rates of growth, development, and metabolism, to reproduction and longevity. However, the neuron-specific IIS/FOXO targets that regulate these functions are still unknown. By isolating C. elegans neurons for transcriptional profiling, we identified both the wild-type and IIS/FOXO neuronal transcriptomes. Here we deposit the sequences for wt adult neurons and IIS mutants, as well as for isolated adult mechanosensory neurons.

Relevant _C. elegans_ Reference Sequences  

* Transcriptome: ftp://ftp.ensemblgenomes.org/pub/release-31/metazoa/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.31.gtf.gz   
* Genome: ftp://ftp.ensemblgenomes.org/pub/release-31/metazoa/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.31.dna.chromosome.*.fa.gz   
* Proteome: ftp://ftp.ensemblgenomes.org/pub/release-31/metazoa/fasta/caenorhabditis_elegans/pep/Caenorhabditis_elegans.WBcel235.31.pep.all.fa.gz    

### How to run trinity  
* Link to website  
* Parameter description
* Command used to generate _C. elegans_ assembly

### Assembly results  
Assembly metrics calculated using Transrate, see http://hibberdlab.com/transrate/metrics.html for a description of the different metrics.  


```r
assembly_metrics_df <- list(no_mods = "data/no_mods/", 
                            corrected = "data/corrected/",
                            trimmed = "data/trimmed/", 
                            trimmedcorrected = "data/trimmedcorrected/") %>% 
      map(paste0, "assemblies.csv") %>%  
      map_df(read_csv, .id = "read_set") %>% 
      select(-assembly) %>% 
      gather("metric","value", -read_set) 
```

__TODO__ Need to workout `read_set` order

```r
assembly_metrics_df %>% 
      filter(metric %in% c("n_seqs", "smallest", "largest", 
                           "n50", "n_with_orf")) %>%
      spread(metric, value) %>% 
      select(read_set, n_seqs, smallest, n50, largest, n_with_orf) %>% 
      kable(digits = 2,
            caption = "Trinity _C. elegans_ transcriptome assembly summary metrics.")
```



Table: Trinity _C. elegans_ transcriptome assembly summary metrics.

read_set            n_seqs   smallest   n50   largest   n_with_orf
-----------------  -------  ---------  ----  --------  -----------
corrected            61059        224   467      4364        12075
no_mods              63353        224   465      5446        12426
trimmed              42438        224   531      7533        10552
trimmedcorrected     42851        224   532      6816        10645

## Assembly evaluation
### Running transrate
* Installation
      * Precompiled binaries available from website (http://hibberdlab.com/transrate/), along with instructions to build from source
      * also command to install dependencies `transrate --install-deps type` where type can be either `all`, `read`, or `ref`.
      * See Transrate website for additional information for install. 


__Command Line__  

_Read Based_

```
transrate \
  --left=SRR2969230_1.fastq \
  --right=SRR2969230_2.fastq \
  --assembly=Trinity_sequences.Trinity.fixed.fasta \
  --output=transrate/
```

_Reference Based_
```
transrate \
  --reference=Caenorhabditis_elegans.WBcel235.31.pep.all.fa \
  --assembly=Trinity_sequences.Trinity.fixed.fasta \
  --output=transrate/
```

### Transrate results
The unmodified read set was passed as input to Transrate, 
potentally biasing some metric results towards the unmodifed assembly.

#### run time
__TODO__ Add runtime for other assemblies and reference based evaluations  

- real    890m58.016s
- user    6665m50.869s
- sys     11m5.890s

#### Assembly Scores
Read trimming and error correction resulted in the highest weighted assembly score, 
while having the lowest contig score cutoff, 
but a lower optimal score.  


```r
assembly_metrics_df %>% 
      filter(metric %in% c("score", "optimal_score", "cutoff", "weighted")) %>%
      spread(metric, value) %>% 
      kable(digits = 2,
            caption = "Trinity _C. elegans_ transcriptome Transrate score summary.")
```



Table: Trinity _C. elegans_ transcriptome Transrate score summary.

read_set            cutoff   optimal_score   score   weighted
-----------------  -------  --------------  ------  ---------
corrected             0.24            0.07    0.02       0.65
no_mods               0.19            0.07    0.02       0.54
trimmed               0.22            0.07    0.03       0.67
trimmedcorrected      0.07            0.06    0.03       0.78


#### Assembly Score Optimization

```r
assembly_score_opt <- results_list %>% map(paste0,"assembly_score_optimisation.csv") %>% 
      map_df(read_csv, .id = "read_set")
```

Relationship between the cutoff for contig score and assembly score.

```r
assembly_score_opt %>% 
      ggplot() + geom_path(aes(x = cutoff, y = assembly_score, color = read_set)) +
            theme_bw() +
            labs(x = "Contig Score Threshold", y = "Assembly Score")
```

```
## Warning: Removed 2 rows containing missing values (geom_path).
```

![](README_files/figure-html/unnamed-chunk-5-1.png)

### Contigs

```r
contig_stat <- results_list %>% map(paste0,"contigs.csv") %>% 
      map_df(read_csv, .id = "read_set")

contig_cutoff <- assembly_metrics_df %>% filter(metric == "cutoff") %>% 
      select(-metric) %>% rename(cutoff = value)
contig_stat <- contig_stat %>% left_join(contig_cutoff) %>% 
      mutate(contig_filt = cutoff < score)
```

```
## Joining by: "read_set"
```


__Contig Score Distribution__
Trimming reads resulted in a higher proportion of contigs with scores greater than 0.25 then assemblies using unmodified reads or only error corrected reads.  


```r
contig_stat %>% 
      ggplot() + geom_density(aes(x = score, color = read_set, fill = read_set), 
                              alpha = 0.25) + theme_bw()
```

![](README_files/figure-html/unnamed-chunk-7-1.png)


Relationship between contig score and length. 

```r
contig_stat %>% 
      ggplot() + geom_hex(aes(x = length, y = score)) +
            geom_hline(aes(yintercept = cutoff), 
                       color = "grey60", linetype = 2) +
            theme_bw() +
            labs(x = "Contig Length (bp)", y = "Contig Score") +
            facet_wrap(~read_set)
```

<div class="figure">
<img src="README_files/figure-html/score_length-1.png" alt="2D histogram of contig score and legnth, with color indicating abundance. Grey dotted lines are the contig score cutoff for optimal assembly score."  />
<p class="caption">2D histogram of contig score and legnth, with color indicating abundance. Grey dotted lines are the contig score cutoff for optimal assembly score.</p>
</div>

Proportion of bases in reads map to the transcript that agree with the transcript contig. 


```r
contig_stat %>% 
      ggplot() + geom_histogram(aes(x = p_good)) +
            theme_bw() +
            labs(x = "Proportion of Bases Aggree with Transcript",
                 y = "Count") +
            facet_grid(contig_filt~read_set) +
            theme(legend.position = "bottom")
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

<div class="figure">
<img src="README_files/figure-html/p_good_density-1.png" alt="Histogram of the proportion of bases in agreement with contig and contig length."  />
<p class="caption">Histogram of the proportion of bases in agreement with contig and contig length.</p>
</div>

Proportion of bases coverage by mapped reads.

```r
contig_stat %>% 
      ggplot() + geom_histogram(aes(x = p_bases_covered)) +
            theme_bw() +
            labs(x = "Proportion of Covered Bases", 
                 y = "Count") +
            facet_grid(contig_filt~read_set) +
            theme(legend.position = "bottom")
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](README_files/figure-html/unnamed-chunk-8-1.png)


```r
contig_stat %>% 
      ggplot() + geom_histogram(aes(x = p_seq_true)) +
            theme_bw() +
            labs(x = "Proportion Seq Correct", 
                 y = "Count") +
            facet_grid(contig_filt~read_set) +
            theme(legend.position = "bottom")
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](README_files/figure-html/unnamed-chunk-9-1.png)


```r
contig_stat %>% 
      ggplot() + geom_histogram(aes(x = p_not_segmented)) +
            theme_bw() +
            labs(x = "Proportion of Non-Segmented Contigs", 
                 y = "Count") +
            facet_grid(contig_filt~read_set) +
            theme(legend.position = "bottom")
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](README_files/figure-html/unnamed-chunk-10-1.png)

The individual parameters used to calculate the contig score are weakly correlated. 
Correlation values are similar to thoes presented in the Transrate publication. 


```r
for(i in names(results_list)){
      pl <- contig_stat %>% filter(read_set == i) %>% 
            select(p_good, p_bases_covered, p_seq_true, p_not_segmented) %>% 
            ggcorr(method = c("pairwise","spearman"),label = TRUE) + 
                  ggtitle(i)
      print(pl)
}
```

![](README_files/figure-html/unnamed-chunk-11-1.png)![](README_files/figure-html/unnamed-chunk-11-2.png)![](README_files/figure-html/unnamed-chunk-11-3.png)![](README_files/figure-html/unnamed-chunk-11-4.png)

#### Reference Based Results
__TODO__ Describe CRB method




```r
ref_contig <- results_list %>% map(paste0,"contigs.csv") %>% 
      map_df(read_csv, .id = "read_set")
```

Comprison of contigs with and with out CRB hits. Assemblies using untrimmed read set had a larger number of contigs, but little impact of the proportion of CRB hits.

```r
ggplot(ref_contig) + geom_bar(aes(x = read_set)) + facet_wrap(~has_crb)
```

![](README_files/figure-html/unnamed-chunk-13-1.png)

Comprison of contig length and ORF length for contigs with and without CRB hits. 
Contigs with CRB hits were long and contined larger predicted ORFs.  


```r
ggplot(ref_contig) + geom_boxplot(aes(x = has_crb, y = length)) + facet_wrap(~read_set)
```

![](README_files/figure-html/unnamed-chunk-14-1.png)


```r
ggplot(ref_contig) + geom_boxplot(aes(x = has_crb, y = orf_length)) + facet_wrap(~read_set)
```

![](README_files/figure-html/unnamed-chunk-15-1.png)

Exploring contigs with CRB hits.  

```r
crb_contig <- ref_contig %>% filter(has_crb == "true")
```


```r
ggplot(crb_contig) + geom_hex(aes(x = length, y = reference_coverage)) + 
      facet_wrap(~read_set)
```

![](README_files/figure-html/unnamed-chunk-17-1.png)


```r
ggplot(crb_contig) + geom_hex(aes(x = orf_length, y = reference_coverage)) +
      facet_wrap(~read_set)
```

![](README_files/figure-html/unnamed-chunk-18-1.png)


### Contig Score and Hit Coverage

```r
contig_crb_score <- crb_contig %>% left_join(contig_stat)
```

```
## Joining by: c("read_set", "contig_name", "length", "prop_gc", "gc_skew", "at_skew", "cpg_count", "cpg_ratio", "orf_length", "linguistic_complexity_6")
```

__TODO__ Better way to present relationship

```r
ggplot(contig_crb_score) + 
      geom_hex(aes(x = reference_coverage, y = score)) +
      facet_wrap(~read_set)
```

![](README_files/figure-html/unnamed-chunk-20-1.png)



## Conclusions
* How did Trinity do?
* How did transrate compare to other methods?
* Recommendations for transcriptome assembly and evaluation

## Session Information
__Session Info__

```r
s_info <- devtools::session_info()
print(s_info$platform)
```

```
##  setting  value                       
##  version  R version 3.2.3 (2015-12-10)
##  system   x86_64, darwin13.4.0        
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  tz       America/New_York            
##  date     2016-05-09
```

```r
kable(s_info$packages)
```



package      *    version   date         source                             
-----------  ---  --------  -----------  -----------------------------------
assertthat        0.1       2013-12-06   CRAN (R 3.2.0)                     
colorspace        1.2-6     2015-03-11   CRAN (R 3.2.0)                     
DBI               0.3.1     2014-09-24   CRAN (R 3.2.0)                     
devtools          1.9.1     2015-09-11   CRAN (R 3.2.0)                     
digest            0.6.9     2016-01-08   CRAN (R 3.2.3)                     
dplyr        *    0.4.3     2015-09-01   CRAN (R 3.2.0)                     
evaluate          0.8       2015-09-18   CRAN (R 3.2.0)                     
formatR           1.2.1     2015-09-18   CRAN (R 3.2.0)                     
GGally       *    1.0.1     2016-01-14   CRAN (R 3.2.3)                     
ggplot2      *    2.1.0     2016-03-01   CRAN (R 3.2.4)                     
gtable            0.1.2     2012-12-05   CRAN (R 3.2.0)                     
hexbin       *    1.27.1    2015-08-19   CRAN (R 3.2.0)                     
highr             0.5.1     2015-09-18   CRAN (R 3.2.0)                     
htmltools         0.3       2015-12-29   CRAN (R 3.2.3)                     
knitr        *    1.12.3    2016-01-22   CRAN (R 3.2.3)                     
labeling          0.3       2014-08-23   CRAN (R 3.2.0)                     
lattice           0.20-33   2015-07-14   CRAN (R 3.2.3)                     
lazyeval          0.1.10    2015-01-02   CRAN (R 3.2.0)                     
magrittr          1.5       2014-11-22   CRAN (R 3.2.0)                     
memoise           0.2.1     2014-04-22   CRAN (R 3.2.0)                     
munsell           0.4.2     2013-07-11   CRAN (R 3.2.0)                     
plyr              1.8.3     2015-06-12   CRAN (R 3.2.0)                     
purrr        *    0.2.1     2016-02-13   CRAN (R 3.2.3)                     
R6                2.1.1     2015-08-19   CRAN (R 3.2.0)                     
Rcpp              0.12.2    2015-11-15   CRAN (R 3.2.2)                     
readr        *    0.2.2     2015-10-22   CRAN (R 3.2.0)                     
reshape           0.8.5     2014-04-23   CRAN (R 3.2.0)                     
reshape2          1.4.1     2014-12-06   CRAN (R 3.2.0)                     
rmarkdown         0.9.5     2016-01-25   Github (rstudio/rmarkdown@e7a1e6b) 
scales            0.3.0     2015-08-25   CRAN (R 3.2.0)                     
stringi           1.0-1     2015-10-22   CRAN (R 3.2.0)                     
stringr      *    1.0.0     2015-04-30   CRAN (R 3.2.0)                     
tidyr        *    0.4.1     2016-02-05   CRAN (R 3.2.3)                     
yaml              2.1.13    2014-06-12   CRAN (R 3.2.0)                     
## References
