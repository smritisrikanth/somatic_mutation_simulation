param_tb = as_tibble(expand.grid(mutation_rate = c(1e-8,1e-9,1e-10),
                                 sampling = c('fixed', 'proportional'), sample_size = c(50,100,200), num_sim_2 = 1:2, num_sim_1 = 1:20))

param_tb$input_file = paste('sample_size_', param_tb$sample_size, '_', param_tb$sampling, '_', floor((as.numeric(rownames(param_tb))-1)/3), '.rda', sep = "")

write.table(as.data.frame(param_tb),'/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/output/param_tb.txt',col.names = T,row.names = F,quote = F)