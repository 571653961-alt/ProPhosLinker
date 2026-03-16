#' Parse a semicolon-separated string of group comparisons into a tidy sample list data frame.
#'
#' This function takes an input string that encodes multiple group comparisons in a specific format,
#' parses each comparison, extracts group names and associated sample identifiers, and returns
#' a unified data frame listing all samples along with their group and comparison labels.
#'
#' @param input_str A character string encoding comparisons in the format:
#'   "GroupA-VS-GroupB@(SampleA1,SampleA2,...)/(SampleB1,SampleB2,...);..."
#'   Multiple comparisons are separated by semicolons (`;`).
#' @return A data.frame with columns: `group`, `sample`, and `comparison`.
#' @examples
#' input <- "PCOS-VS-Control@(PCOS_1,PCOS_2)/(Control_1,Control_2);Obese-VS-Lean@(Obese_1)/(Lean_1,Lean_2)"
#' run_samplelist(input)
#' @export
run_samplelist <- function(input_str) {
  # 1. Split the input string into individual comparison entries
  comparisons <- unlist(strsplit(input_str, ";"))
  
  # 2. Initialize an empty list to store parsed data frames for each comparison
  samplelist <- list()
  
  # 3. Process each comparison entry
  for (comp in comparisons) {
    
    # Split into group label and sample specification using "@" as delimiter
    parts <- unlist(strsplit(comp, "@"))
    
    # Extract the comparison label (e.g., "PCOS-VS-Control")
    group_str <- parts[1]
    
    # Split the group label into two group names using "-VS-" as delimiter
    group_parts <- unlist(strsplit(group_str, "-VS-"))
    group1 <- group_parts[1]
    group2 <- group_parts[2]
    
    # Extract the sample string (e.g., "(PCOS_1,...)/(Control_1,...)")
    sample_str <- parts[2]
    
    # Remove leading and trailing parentheses from the sample string
    sample_str <- gsub("^\\(", "", sample_str)
    sample_str <- gsub("\\)$", "", sample_str)
    
    # Split into two parts: samples for group1 and group2, separated by "/"
    samples <- unlist(strsplit(sample_str, "/"))
    
    # Helper function to clean a comma-separated sample list:
    # - split by commas
    # - remove any remaining parentheses or whitespace
    # - drop empty entries
    clean_samples <- function(sample_list) {
      samples_vec <- unlist(strsplit(sample_list, ","))
      samples_clean <- gsub("[\\(\\)\\s]", "", samples_vec)
      samples_clean <- samples_clean[samples_clean != ""]
      return(samples_clean)
    }
    
    # Clean and extract samples for group1 and group2
    group1_samples <- clean_samples(samples[1])
    group2_samples <- clean_samples(samples[2])
    
    # Create data frames for each group with metadata
    df_group1 <- data.frame(
      group = group1,
      sample = group1_samples,
      comparison = group_str,
      stringsAsFactors = FALSE
    )
    
    df_group2 <- data.frame(
      group = group2,
      sample = group2_samples,
      comparison = group_str,
      stringsAsFactors = FALSE
    )
    
    # Combine both groups into one data frame for this comparison
    samplelist[[group_str]] <- rbind(df_group1, df_group2)
  }
  
  # 4. Validate that at least one comparison was parsed
  if (length(samplelist) == 0) {
    stop("No valid comparison groups found in input string.")
  }
  
  # 5. Combine all comparison-specific data frames into a single data frame
  final_samplelist <- do.call(rbind, samplelist)
  
  # 6. Reset row names for cleanliness
  rownames(final_samplelist) <- NULL
  
  return(final_samplelist)
}