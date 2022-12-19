#' Oral and lung microbiome replicate samples data
#' 
#' Samples, taxa, and read counts for the Charlson-Bittinger study, "Assessing
#' bacterial populations in the lung by replicate analysis of samples from the
#' upper and lower respiratory tracts."
#' @references Charlson ES, Bittinger K, Chen J, Diamond JM, Li H, Collman RG,
#'   Bushman FD. Assessing bacterial populations in the lung by replicate
#'   analysis of samples from the upper and lower respiratory tracts. PLoS One.
#'   2012;7(9):e42786. doi: 10.1371/journal.pone.0042786. Epub 2012 Sep 6.
#'   PMID: 22970118; PMCID: PMC3435383.
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/22970118/}
#' @name lungreplicate

#' @rdname lungreplicate
#' @format \code{lungreplicate_samples} is a data frame with 42 rows and 5
#'   variables:
#' \describe{
#'   \item{sample_id}{unique identifiers for each sample}
#'   \item{subject_id}{unique identifiers for each study participant, matching
#'     those in the publication}
#'   \item{sample_type}{the type of sample, "Oral wash" or "Bronchoalveolar
#'     lavage"}
#'   \item{replicate_id}{unique identifiers for each replicate sample}
#'   \item{study_group}{the study group into which the participant was
#'     recruited, either "Healthy" or "Lung transplant"}
#' }
"lungreplicate_samples"

#' @rdname lungreplicate
#' @format \code{lungreplicate_samples} is a data frame with 42 rows and 5
#'   variables:
#' \describe{
#'   \item{otu_id}{unique identifiers for each operational taxonomic unit (OTU)}
#'   \item{sample_id}{unique identifiers for each sample}
#'   \item{reads}{number of sequence reads attributed to each OTU in each
#'     sample}
#'   \item{assignment}{the taxonomic assignment of the OTU}
#' }
"lungreplicate_reads"
