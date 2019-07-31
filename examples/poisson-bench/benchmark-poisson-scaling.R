#! /bin/env Rscript

library(dplyr)
library(ggplot2)
library(rstan)
library(readr)
theme_set(theme_bw(10))

binary  <- "./poisson-hierarchical-scale"


if(!file.exists(binary)) {
    ## executed once
    ## possibly adapt to your compiler (note that make/local is wiped out)
    cat("CXXFLAGS+=-DSTAN_THREADS\n", file="../../make/local", append=FALSE)
    cat("CXX=clang++\n", file="../../make/local", append=TRUE)
    cat("STANCFLAGS=--allow_undefined\n", file="../../make/local", append=TRUE)
    system("cd ../..; make -j6 build")
    system("cd ../..; make examples/poisson-bench/poisson-hierarchical-scale")
}

reduce_sum_timing <- function(threads, groups, terms, method, grainsize=0) {
    pars <-  as.integer(c(threads, groups, terms, method, grainsize))
    env <- paste(c("STAN_NUM_THREADS", "GROUPS", "TERMS", "method", "grainsize"), pars, sep="=")
    cat("Problem:", paste(env, collapse=","), "\n")
    outf <- paste0("output-", paste(env, collapse="-"), ".csv")
    data_file  <- tempfile("benchmark", fileext=".R")
    bench_data  <- list(G=as.integer(groups),
                        N=as.integer(terms),
                        grainsize=as.integer(grainsize),
                        method=as.integer(method))
    stan_rdump(names(bench_data), data_file, envir=list2env(bench_data))
    runtime  <- system.time(system2(binary, c("sample num_samples=100 num_warmup=100", "data", paste0("file=", data_file), "random seed=1 ", paste0("output file=", outf)), env=env))[["elapsed"]]
    file.remove(data_file)
    runtime
}

## test things
reduce_sum_timing(2, 10L, 10L, 0)

## maybe adapt the number of maximal cores

bench <- expand.grid(threads=c(1,2,4,6), groups=c(1,5) * 100, terms=c(50, 100),
                     grainfrac=c(0))

bench <- expand.grid(threads=c(1,2,4,6,8,12), groups=c(5) * 100, terms=c(50), method=0:3,
                     grainfrac=c(0))

bench  <- bench %>% mutate(grainsize=as.integer(groups * terms * grainfrac)) %>%
    filter(!(method==3 & threads>1))

bench

bench$runtime <- apply(bench, 1, function(run) reduce_sum_timing(run[["threads"]], run[["groups"]], run[["terms"]], run[["method"]], run[["grainsize"]]))

##total_problem <- as.integer(100000000)

prefix <- format(Sys.time(), "bench-%Y-%m-%d_%H-%M")

labels <- unique(bench[c("groups", "terms")]) %>%
    mutate(problem_size=groups*terms) %>%
    arrange(problem_size) %>%
    mutate(label=paste0(paste(as.numeric(groups), as.numeric(terms), sep="-"), ": ", groups*terms))
labels$problem <- factor(labels$label, labels$label, labels$label)
labels$problem

benchPl <- bench %>% left_join(labels) %>%
    mutate(grainfrac=factor(round(grainfrac, 3))) %>%
    group_by(method, grainsize, label) %>%
    mutate(method_speedup=runtime[threads==1]/runtime) %>%
    ungroup() %>%
    mutate(serial_speedup=runtime[method==3&threads==1]/runtime,
           method_label=factor(method, 0:3, labels=c("TBB reduce", "TBB map", "map_rect", "serial")))

benchPl

write_csv(benchPl, paste0(prefix, ".csv"))

ggplot(benchPl, aes(threads, method_speedup, colour=method_label, shape=method_label)) +
    geom_point() + geom_line() +
    geom_abline(slope=1, intercept=0, linetype=2) +
    ggtitle("Speedup vs 1 core of hierarchical Poisson likelihood reduce", "Each curve shows method specific speedup relative to 1 core of respective method") +
    ylab("Speedup vs 1 core") +
    xlab("Threads") +
    scale_x_log10(breaks=c(1, 2, 4, 6, 8, 12)) +
    scale_y_log10(breaks=c(1, 2, 3, 4, 5, 6, 8, 12)) +
    coord_fixed(xlim=c(1,12), ylim=c(1,12)) +
    theme(legend.position=c(0.2,0.825))

ggsave(paste0(prefix, "-method-relative-scale-poisson-static.png"), width=6, height=7, dpi=120)

ggplot(benchPl, aes(threads, serial_speedup, colour=method_label, shape=method_label)) +
    geom_point() + geom_line() +
    geom_abline(slope=1, intercept=0, linetype=2) +
    ggtitle("Speedup vs 1 core of hierarchical Poisson likelihood reduce", "Each curve shows speedup relative to 1 core of serial method") +
    ylab("Speedup vs 1 core") +
    xlab("Threads") +
    scale_x_log10(breaks=c(1, 2, 4, 6, 8, 12)) +
    scale_y_log10(breaks=c(1, 2, 3, 4, 5, 6, 8, 12)) +
    coord_fixed(xlim=c(1,12), ylim=c(1,12)) +
    theme(legend.position=c(0.2,0.825))

ggsave(paste0(prefix, "-method-serial-scale-poisson-static.png"), width=6, height=7, dpi=120)


library(readr)

steve  <- read_csv("~/Downloads/bench-2019-07-30_16-36.csv")

cmp  <- steve %>% filter(threads%in%c(1,16), method_label=="TBB map") %>% mutate(runtime_h=runtime/60/60)

names(cmp)

library(knitr)
kable(cmp[c("label", "runtime_h")])
