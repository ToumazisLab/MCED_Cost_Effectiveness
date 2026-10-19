create_params_list <- function(final_result_list, additional_vars) {
  # Extract values from final_icer_results
  final_icer_results <- final_result_list[["Overall"]][["final_icer_result"]]
  
  # Extract Incremental Cost, Utility, and ICER
  params_list <- list(
    Incremental_cost = formatC(
      as.numeric(final_icer_results[final_icer_results$Variable == "Incremental_Cost", "Value"]),
      format = "e", digits = 2
    ),
    Incremental_utility = formatC(
      as.numeric(final_icer_results[final_icer_results$Variable == "Incremental_Utility", "Value"]),
      format = "e", digits = 2
    ),
    ICER = formatC(
      as.numeric(final_icer_results[final_icer_results$Variable == "ICER", "Value"]),
      format = "e", digits = 2
    )
  )
  
  # Add additional variables to the params list
  params_list <- c(params_list, additional_vars)
  
  return(params_list)
}

create_comparison_note <- function(old_params, new_params) {
  # Function to compare parameters and highlight differences
  highlight_difference <- function(param_name) {
    original <- old_params[[param_name]]
    new <- new_params[[param_name]]
    if (!identical(original, new)) {
      glue("{param_name} = NEW: {new}")
    } else {
      glue("{param_name} = {new}")
    }
  }
  
  # Get all parameter names
  param_names <- union(names(old_params), names(new_params))
  
  # Prepare the note with HTML formatting
  note <- paste(
    "Note: ",
    paste(
      sapply(param_names, highlight_difference), 
      collapse = ", "
    ),
    sep = " "
  )
  
  
  return(note)
}






