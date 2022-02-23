read_lef keccak/test.lef
read_def keccak/test.def

bin_graph_extract -num_rows 5 -max_layer 5
bin_graph_labeling -drc_rpt_file test.drc.rpt
save_bin_graph_file -prefix bin_graph

analyze_congestion


puts "star model end"
exit

