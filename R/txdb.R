#' Export a TxDb object to a clean GTF file
#'
#' Exports a TxDb object to a clean GTF file (i.e. only 'exon' and 'CDS' features 
#' and 'gene_id' and 'transcript_id' tags)
#'
#' @param txdb A TxDb object
#' @param file Output file path
#' @export
txdb_to_gtf <- function(txdb, file) {

	# Check the arguments
	invisible({
		assertthat::assert_that(class(txdb) == "TxDb")
		assertthat::assert_that(assertthat::is.string(file))
		assertthat::assert_that(assertthat::is.writeable(dirname(file)))
	})

	# Prepare exon data
	txdb_exons <- GenomicFeatures::exonsBy(txdb, by = "tx")
	txdb_exons <- S4Vectors::stack(txdb_exons, "tx")
	txdb_exons_mcols <- AnnotationDbi::select(txdb, keys = as.character(S4Vectors::mcols(txdb_exons)$tx), 
                                              keytype = "TXID", columns = c("TXNAME", "GENEID"))
	S4Vectors::mcols(txdb_exons)$gene_id <- txdb_exons_mcols$GENEID
	S4Vectors::mcols(txdb_exons)$tx_id <- txdb_exons_mcols$TXNAME

	# Prepare CDS data
	txdb_cds <- GenomicFeatures::cdsBy(txdb, by = "tx")
	txdb_cds <- S4Vectors::stack(txdb_cds, "tx")
	txdb_cds_mcols <- AnnotationDbi::select(txdb, keys = as.character(S4Vectors::mcols(txdb_cds)$tx), 
                                            keytype = "TXID", columns = c("TXNAME", "GENEID"))
	S4Vectors::mcols(txdb_cds)$gene_id <- txdb_cds_mcols$GENEID
	S4Vectors::mcols(txdb_cds)$tx_id <- txdb_cds_mcols$TXNAME

	# Convert TxDb object to GTF format
	txdb_features <- BiocGenerics::sort(c(txdb_exons, txdb_cds))
	gene_ids <- S4Vectors::mcols(txdb_features)$gene_id
	transcript_ids <- S4Vectors::mcols(txdb_features)$tx_id
	annotations_df <- data.frame(
		seqname = as.character(GenomeInfoDb::seqnames(txdb_features)), 
		source = "TxDb", 
		feature = ifelse(is.na(S4Vectors::mcols(txdb_features)$cds_id), "exon", "CDS"), 
		start = BiocGenerics::start(txdb_features), 
		end = BiocGenerics::end(txdb_features), 
		score = ".", 
		strand = as.character(BiocGenerics::strand(txdb_features)), 
		frame = ".", 
		attributes = glue::glue('gene_id "{gene_ids}"; transcript_id "{transcript_ids}"')
	)

	# Save the GTF file
 	utils::write.table(annotations_df, file = file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

	# Return nothing
	invisible()
}
