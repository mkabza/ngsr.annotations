#' Prepare reference genome data directory using Ensembl data 
#'
#' Prepares a directory containing reference genome files based on Ensembl FASTA and GTF files.
#' The output directory contains following files:
#' \itemize{
#'   \item genome.fasta - genome sequence
#'   \item annotations.gtf - genome annotations
#'   \item transcriptome.fasta - transcriptome sequence
#'   \item genes.rds - gene GRanges
#'   \item genes.csv - gene descriptions
#'   \item transcripts.rds - transcript GRanges
#'   \item transcripts.csv - transcript descriptions
#' }
#'
#' @param output_dir Output directory path
#' @param fasta A path or an URL to an Ensembl FASTA file contaning genome sequence
#' @param gtf A path or an URL to an Ensembl GTF file contaning genome annotations
#' @export
prepare_reference_ensembl <- function(output_dir, fasta, gtf) {

	# Check the arguments
	invisible({
		assertthat::assert_that(assertthat::is.string(output_dir))
		assertthat::assert_that(!file.exists(output_dir))
		assertthat::assert_that(assertthat::is.writeable(dirname(output_dir)))
		assertthat::assert_that(assertthat::is.string(fasta))
		assertthat::assert_that(assertthat::is.string(gtf))
	})

	# Create the output directory
	dir.create(output_dir, recursive = TRUE)

	# Prepare genome sequence
	genome_seq <- Biostrings::readDNAStringSet(fasta, format = "fasta")
	names(genome_seq) <- sapply(strsplit(names(genome_seq), "\\s+"), "[", 1)
	Biostrings::writeXStringSet(genome_seq, filepath = file.path(output_dir, "genome.fasta"), format = "fasta")

	# Prepare genome annotations
	txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")
	txdb_to_gtf(txdb, file = file.path(output_dir, "annotations.gtf"))

	# Prepare transcriptome sequence
	transcript_seq <- GenomicFeatures::extractTranscriptSeqs(genome_seq, txdb, use.names = TRUE)
	Biostrings::writeXStringSet(transcript_seq, filepath = file.path(output_dir, "transcriptome.fasta"), format = "fasta")

	# Prepare gene and transcript GRanges
	granges <- rtracklayer::import(gtf, format = "gtf")
	S4Vectors::mcols(granges)$location <- as.character(granges)

	gene_ranges <- granges[S4Vectors::mcols(granges)$type == "gene"]
	gene_columns <- c("gene_id", "gene_name", "gene_biotype", "location")
	S4Vectors::mcols(gene_ranges) <- S4Vectors::mcols(gene_ranges)[, gene_columns]
	saveRDS(gene_ranges, file = file.path(output_dir, "genes.rds"))
 	utils::write.csv(S4Vectors::mcols(gene_ranges), file = file.path(output_dir, "genes.csv"), row.names = FALSE)

	transcript_ranges <- granges[S4Vectors::mcols(granges)$type == "transcript"]
	transcript_columns <- c("transcript_id", "transcript_name", "transcript_biotype", "gene_id", "gene_name", "gene_biotype", "location")
	S4Vectors::mcols(transcript_ranges) <- S4Vectors::mcols(transcript_ranges)[, transcript_columns]
	saveRDS(transcript_ranges, file = file.path(output_dir, "transcripts.rds"))
 	utils::write.csv(S4Vectors::mcols(transcript_ranges), file = file.path(output_dir, "transcripts.csv"), row.names = FALSE)

	# Return nothing
	invisible()
}
