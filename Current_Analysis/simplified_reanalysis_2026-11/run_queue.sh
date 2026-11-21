#!/bin/bash
# Runs T5 to T9 in order after T4 finishes. Run from the repository root.
export LC_ALL=en_US.UTF-8
A=Current_Analysis/simplified_reanalysis_2026-11
until grep -q "T4 COMPLETE\|Error\|Execution halted" $A/logs/04_adjusted_dream.log; do sleep 10; done
grep -q "T4 COMPLETE" $A/logs/04_adjusted_dream.log || { echo "T4 FAILED"; exit 1; }
Rscript $A/09_external_summary.R > $A/logs/09_external_summary.log 2>&1 || echo "T9 FAILED"
Rscript $A/05_matched_null.R > $A/logs/05_matched_null.log 2>&1 || echo "T5 FAILED"
Rscript $A/06_zenith_sets.R > $A/logs/06_zenith_sets.log 2>&1 || echo "T6 FAILED"
Rscript $A/08_within_infected.R > $A/logs/08_within_infected.log 2>&1 || echo "T8 FAILED"
Rscript $A/07_random_effects.R > $A/logs/07_random_effects.log 2>&1 || echo "T7a FAILED"
Rscript $A/07b_joint_permutation.R 200 > $A/logs/07b_joint_permutation.log 2>&1 || echo "T7b FAILED"
echo QUEUE DONE
