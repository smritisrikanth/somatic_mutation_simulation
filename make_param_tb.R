library(tidyverse)

param_tb = as_tibble(expand.grid(mutation_rate = c(1e-8,1e-9,1e-10),
                                 sampling = c('fixed', 'proportional'), sample_size = c(50,100,200), num_sim_2 = 1:2, num_sim_1 = 1:20))

param_tb$input_file = paste('sample_size_', param_tb$sample_size, '_', param_tb$sampling, '_', floor((as.numeric(rownames(param_tb))+2)/3), '.rda_', param_tb$mutation_rate, '_', param_tb$num_sim_2, '.rda', sep = "")

write.table(as.data.frame(param_tb),'/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/output2/param_tb_2.txt',col.names = T,row.names = F,quote = F)

# results <- data.frame(matrix(nrow=0, ncol = 8))
# colnames(results) <- c('sample_size', 'sampling', 'mutation_rate', 'input_file', 'job_id', 'num_mut', 'KC0', 'KC1')

# write.table(results,'/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/output2/results_2.txt',col.names = T,row.names = F,quote = F)