run_markov_model_SoC <- function(a_P_SoC, v_names_states, num_people) {
  n_cycles <- dim(a_P_SoC)[3]  # Should be 51 for ages 50 to 100
  n_states <- length(v_names_states)
  m_M_SoC <- matrix(0, nrow = n_cycles + 1, ncol = n_states)
  m_M_SoC[1, which(v_names_states == "H")] <- num_people  # Initialize the cohort in the "H" state
  colnames(m_M_SoC) <- v_names_states
  
  for (t in 1:n_cycles) {
    # Ensure no negative values during transitions
    m_M_SoC[t + 1, ] <- pmax(0, m_M_SoC[t, ] %*% a_P_SoC[, , t])
  }
  
  # Add age labels (ages 49 to 100, where row 1 is the initial state at age 49)
  ages <- 39:(39 + n_cycles)
  m_M_SoC <- as.data.frame(m_M_SoC)
  m_M_SoC$Age <- ages
  
  return(m_M_SoC)
}
