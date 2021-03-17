testarg <- c(
  "-q", "TRUE",                      # quiet
  "-d", "1234567890",                # seed
  "-p", "test",                      # project_name
  "-b", "1",                         # batch_ID
  "-i", "tests/sample_info_test.txt", # sample_info_file
  "-g", "tests/genome_info_test.txt", # genome_info_file
  "-f", "400",                       # num_of_gen_in_forw_sim
  "-w", "8",                         # num_of_periods_forw
  "-l", "2","1.465967","26","30",    # generation_length_prior_params
  "-s", "3",                         # num_of_sims
  "-n", "10","200",                  # population_size_prior_params
  "-u", "0.00000005","0.5"           # mutation_rate_prior_params
)

