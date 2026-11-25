Control logs of the operational-rows port (gmantegna/reimagining-ra #204), della job stdout.
Scripts: gmantegna/reimagining-ra branch base/ops-port-build, mwe/della/port_*.sbatch.
tdd_18376120.txt            test-first record: each commit's test file on its parent's src (pre)
                            and on its own src (post); tdd_18376120_logs.txt = the per-run logs.
tdd_18375616_cancelled.txt  first attempt, cancelled: the thermal test unit was named "T", which
                            CSV.jl reads as Bool true (renamed TH before the record above).
mutants_18377168.txt        mutants 1-13 on the head src (sed, named test file, restore).
mutants_18377435.txt        mutants 14-26. Mutants 2 and 14 break GenX's own precompile workload
                            (src/startup/genx_startup.jl), so they are killed by a load error,
                            not by an assertion; 26 replaces 2.
absent_means_off_18376174.txt  three existing test cases on benders-generalized 544cc5cc and on
                            the head: model size, objective, cmp of the LP files.
suite_*.txt                 whole test/runtests.jl: 18378046 head, 18378273 baseline 544cc5cc,
                            18378305 / 18378650 / 18378828 the intermediate commits 1/4, 2/4, 3/4.
