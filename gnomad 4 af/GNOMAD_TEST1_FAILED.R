# ============================================================================
# DIAGNOSTIC: Test gnomAD API with a known common variant
# ============================================================================

test_gnomad_query <- function() {
  cat("Testing gnomAD API...\n\n")
  
  # Test with a well-known common variant (rs429358 - APOE e4)
  # This should definitely be in gnomAD
  test_variants <- list(
    list(chrom = "19", pos = 44908684, ref = "T", alt = "C", name = "APOE e4 variant"),
    list(chrom = "chr19", pos = 44908684, ref = "T", alt = "C", name = "APOE e4 variant (with chr)")
  )
  
  # Try different dataset names
  datasets <- c("gnomad_r4", "gnomad_r4_0", "gnomad_v4")
  
  for (ds in datasets) {
    cat(sprintf("\n=== Testing dataset: %s ===\n", ds))
    
    for (variant in test_variants) {
      variant_id <- sprintf("%s-%s-%s-%s", variant$chrom, variant$pos, variant$ref, variant$alt)
      cat(sprintf("\nTesting: %s (%s)\n", variant_id, variant$name))
      
      query <- sprintf('{
        variant(variant_id: "%s", dataset: %s) {
          variant_id
          joint {
            af
            ac
            an
          }
        }
      }', variant$chrom, variant$pos, variant$ref, variant$alt, ds)
      
      tryCatch({
        response <- POST(
          "https://gnomad.broadinstitute.org/api",
          body = list(query = query),
          encode = "json",
          timeout(10)
        )
        
        cat(sprintf("Status code: %d\n", status_code(response)))
        
        if (status_code(response) == 200) {
          data <- content(response, "parsed")
          
          # Print the full response for debugging
          cat("Raw response:\n")
          print(str(data))
          
          if (!is.null(data$errors)) {
            cat("\nERRORS FOUND:\n")
            print(data$errors)
          }
          
          if (!is.null(data$data$variant)) {
            cat("\n✓ SUCCESS! Variant found:\n")
            cat(sprintf("  AF: %s\n", data$data$variant$joint$af))
            cat(sprintf("  AC: %s\n", data$data$variant$joint$ac))
            cat(sprintf("  AN: %s\n", data$data$variant$joint$an))
            return(list(success = TRUE, dataset = ds, format = variant$chrom))
          } else {
            cat("✗ Variant not found in response\n")
          }
        } else {
          cat("✗ Bad status code\n")
          cat(content(response, "text"))
        }
        
      }, error = function(e) {
        cat(sprintf("✗ Error: %s\n", e$message))
      })
    }
  }
  
  cat("\n\n=== DIAGNOSIS COMPLETE ===\n")
  cat("If all tests failed, the API might be down or changed.\n")
  cat("Check: https://gnomad.broadinstitute.org/\n")
}

# Run the diagnostic
test_gnomad_query()


# ============================================================================
# ALTERNATIVE: Try the REST API instead of GraphQL
# ============================================================================

test_rest_api <- function() {
  cat("\n\n=== Testing REST API Alternative ===\n")
  
  # gnomAD also has a simpler REST endpoint
  variant_id <- "19-44908684-T-C"
  url <- sprintf("https://gnomad.broadinstitute.org/api/variant/%s", variant_id)
  
  tryCatch({
    response <- GET(url, timeout(10))
    cat(sprintf("Status: %d\n", status_code(response)))
    
    if (status_code(response) == 200) {
      data <- content(response, "parsed")
      cat("REST API response:\n")
      print(str(data))
    } else {
      cat("REST API not available or different format\n")
    }
  }, error = function(e) {
    cat(sprintf("REST API Error: %s\n", e$message))
  })
}

test_rest_api()