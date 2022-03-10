read_lef keccak/test.lef
read_def keccak/test.def

construct_gcell_grid -num_rows 5 -max_route_layer 7
read_routing_report -in_file keccak/test.drc.rpt

save_map_images -save_dir ./img
#bin_graph_extract -num_rows 5 -max_layer 5
#bin_graph_labeling -drc_rpt_file test.drc.rpt
#save_bin_graph_file -prefix bin_graph
#analyze_congestion

exit

