#!/bin/bash
cat ../data/symmetry_failure/files_to_cp/Sb.castep > Sb.castep
cat ../data/symmetry_failure/files_to_cp/Sb.0001.err > Sb.0001.err
sleep 1
exit 1
