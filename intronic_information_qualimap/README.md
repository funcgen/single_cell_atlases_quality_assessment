In order to properly run the script format_intronic_info.py you first have to extract the intronic content line from qualimap report. For the smart-seq2 assay qualimap report, you can do that easily using bash tools:
$ for i in `cat cells_id.txt`; do cat $i/rnaseq_qc_results.txt | grep "intronic = "; done > intronic_info.txt
In the line of code above, cells_id.txt should contain the cells id, one per line.
