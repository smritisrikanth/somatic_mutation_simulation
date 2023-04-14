sink(stdout(), type="message")

library(tidyverse)
library(r2r)
library(qfm)
library(furrr)
#library(rgl)
#library(treespace)
library(TreeDist)
library(igraph)
library(ggraph)
setwd('/home/ssrikan2/data-kreza1/smriti/qfm2')
devtools::load_all()

setwd('/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo')

#job_id = 49
job_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# construct parameter table 

param_tb <- read.table('/home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/correctness_output/param_tb_2.txt',header = T)
filename <- paste('./input/yi_output/', param_tb$input_file[job_id], sep = "")
if (file.exists(filename)) {
    load(filename)
} else {
    print('file not found')
    stop()
}

#functions

compare_lists <- function (l1, l2) {
  if ((length(intersect(l1,l2)) == length(l1)) && (length(intersect(l1,l2)) == length(l2))) {
    print(length(intersect(l1,l2)))
    print(length(l1))
    print(length(l2))
    return(TRUE)
  }
  return(FALSE)
}

get_partitions <- function(gr) {
        obj = list_dd_and_tips_mod2(gr)
        dd_list = obj$dd
        tips_list = obj$tips
        gr_tips = gr$tip.label
        names(gr_tips) = gr_tips
        tips_list = c(tips_list, gr_tips)
        part_list = map(dd_list, function(dd_vec) {
                tips_list[dd_vec]
        })
        part_list
}

compare_partitions <- function (p1, p2) {
  overall_count = 0
  for (node1 in p1) {
    for (node2 in p2) {
      if (length(node1) != length(node2)) {
        break
      }
      count = 0
      for (i in 1:length(node1)) {
        for (j in 1:length(node2)) {
          l1 <- node1[[i]]
          l2 <- node2[[j]]
          
          if (compare_lists(l1,l2)) {
            print(compare_lists(l1,l2))
            count <- count + 1
            print(l1)
            print(l2)
          }
        }
      }
      
      if (count == length(node1)) {
        overall_count <- overall_count + 1
      }
    }
  }
  return(c(overall_count, length(p1), length(p2)))
}

phy = readRDS("/home/ssrikan2/data-kreza1/smriti/qfm2/intermediate_data/gast_phylo.rds")

#calculations
p <- get_partitions(phy)
p0 <- get_partitions(res0$gr)
p1 <- get_partitions(res1$gr)

r0 <- compare_partitions(p, p0)
r1 <- compare_partitions(p, p1)

print(r0[1])
print(r1[1])

r0_correct <- r0[1]/r0[3]
r1_correct <- r1[1]/r1[3]
r0_complete <- r0[1]/r0[2]
r1_complete <- r1[1]/r1[2]



#save results
result0 <- paste(param_tb$sample_size[job_id], param_tb$sampling[job_id], filename, 'res0', job_id, r0_correct, r0_complete, sep = '  ')
result1 <- paste(param_tb$sample_size[job_id], param_tb$sampling[job_id], filename, 'res1', job_id, r1_correct, r1_complete, sep = '  ')

system(paste("echo ",result0,' >> /home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/correctness_output/results_2.txt', sep = ""))
system(paste("echo ",result1,' >> /home/ssrikan2/data-kreza1/smriti/somatic_mut_sim/git_repo/correctness_output/results_2.txt', sep = ""))


