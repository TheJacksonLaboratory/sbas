#' Compile feature, phenotype and assay data into one Biobase::ExpressionSet object
#'
#' @description Creates an ExpressionSet object from 3 individual dataframes with assay data, phenotype and feature metadata information
#' 
#' @param assay_dataframe A dataframe with N columns = samples + 1, the unique feature id
#' @param pheno_metadata A dataframe with N rows >= samples present in assay_dataframe
#' @param features_metadata A dataframe with at least one column, same the unique feature id present in assay_dataframe
#' @param assay_data_feature_id_variable The unique feature id present in assay_dataframe and feature_dataframe
#' @param pheno_data_sample_id_variable The unique sample id. Must be a column in pheno_data, the values of this columns must be included in assay_dataframe
#'
#' @return A Biobase ExpressionSet object
#' @export
#'
#' @examples  \dontrun{
#'     createExpressionSet(assay_dataframe   = ijc.iso.counts.mem, 
#'                         features_metadata = fromGTF.SE, 
#'                         pheno_metadata    = meta.data, 
#'                         assay_data_feature_id_variable = "ID", 
#'                         pheno_data_sample_id_variable  = "Run" ))
#'     }
#'     
createExpressionSet <- function ( assay_dataframe   = assay_dataframe,
                                  pheno_metadata    = pheno_data,
                                  features_metadata =  features_metadata,
                                  assay_data_feature_id_variable = "ID", 
                                  pheno_data_sample_id_variable  = "Run") {
                                  
                                  # Name appropriately the rows of the daataframes to facilitate merging
                                  rownames(assay_dataframe)   <- assay_dataframe[[assay_data_feature_id_variable]]
                                  rownames(features_metadata) <- features_metadata[[assay_data_feature_id_variable]]
                                  rownames(pheno_metadata)    <- pheno_metadata[[pheno_data_sample_id_variable]]
                                  
                                  # Subset as required to satisfy the conditions related to dimensions of each object
                                  toKeepSamples      <- rownames(pheno_metadata)[ rownames(pheno_metadata) %in% colnames(assay_dataframe)  ] 
                                  assay_dataframe    <- assay_dataframe[, toKeepSamples]
                                  pheno_metadata     <- pheno_metadata[rownames(pheno_metadata) %in% toKeepSamples, ]
                                  features_metadata  <- features_metadata[rownames(features_metadata) %in% rownames(assay_dataframe) , ]

                                  # Create ExpressionSet object
                                  expression_set<- Biobase::ExpressionSet(assayData   = as.matrix(assay_dataframe),
                                                                          phenoData   = AnnotatedDataFrame(pheno_metadata),
                                                                          featureData = AnnotatedDataFrame(features_metadata))
                                  return(expression_set) 
}
