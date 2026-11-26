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

Fix round 1 (review of GenX_Robust PR 14; scripts mwe/della/portfix_*.sbatch and port_suite.sbatch):
fix1_tdd_18386422.txt       the new tests L10, L11 (test_ops_loaders.jl) and W1
                            (test_ops_floors_budgets.jl) on the src of 32bd4675 (pre: L10 and L11
                            fail, nothing thrown; W1 passes), then each fix commit's own tree.
fix1_mutants.txt            job 18386423: the review's surviving mutants r05, r11, r06, r16 and two
                            on the new refusals, on the tree of fix commit 3/3; each killed by W1,
                            L10 or L11.
fix1_suite_18386424.txt     whole test/runtests.jl on the tree of fix commit 3/3.
fix1_absent_means_off_18386425.txt  the three cases of absent_means_off_18376174.txt, 544cc5cc
                            against fix commit 3/3: LP files identical.
fix1_gate_co2off_18386426.txt  the port on C_fix with CO2Cap 0 at fix commit 3/3 against d58's
                            CO2-off objective (read from the predecessor's out_d58_co2off).
The commit adding these logs changes no file outside test/controls/.

Port 2 (StorSimultaneousCharge, import CO2 on lines, Benders cut slopes; scripts mwe/della/port2_*.sbatch
and port_suite.sbatch; branch base/ops-port2):
port2_tdd_18391895.txt      the new test files on the src of 4d602268 (pre); simcharge and cut_slopes
                            are the committed files; its co2_lines run is of a first version.
port2_tdd_18392423.txt      the final test_ops_co2_lines.jl on 4d602268 (pre), then 06a34d36 (B2)
                            and d4b95f69 (B3: all three new files) (post).
port2_tdd_18391905.txt      post-change of c6d38e96 (B1); its other two runs are of first versions.
port2_mutants.txt           jobs 18392424 and 18392425: 18 mutants on the tree of d4b95f69.
port2_suite_18392426.txt    whole test/runtests.jl on the tree of d4b95f69.
port2_absent_means_off_18392427.txt  the three cases, 4d602268 against d4b95f69: LP files identical.
port2_ref_f63283ac_18391530.txt  base/d58-rungb f63283ac on C_fix: StorSimultaneousCharge 1 and 0
                            (CO2Cap 1), and 1 with CO2Cap 0 (the references of the parity runs).
port2_gate_18392428.txt     the port at d4b95f69 on C_fix: CO2Cap 1; StorSimultaneousCharge 1
                            (CO2Cap 1 and 0); CO2Cap 0.
