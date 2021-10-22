set design aes_cipher_top
read_lef NangateOpenCellLibrary.lef
read_def $design.def

#read_liberty NangateOpenCellLibrary_typical.lib
#read_sdc $design.sdc

puts "clique model end"

#graph_extract -out_file extract_clique_a.txt -graph_model clique -edge_weight_model a 500000 500000 550000 550000 
#graph_extract -out_file extract_clique_b.txt -graph_model clique -edge_weight_model b 500000 500000 550000 550000 
#graph_extract -out_file extract_star_a.txt -graph_model star -edge_weight_model a 500000 500000 550000 550000 
#graph_extract -out_file extract_star_b.txt -graph_model star -edge_weight_model b 500000 500000 550000 550000 
bin_graph_extract -num_rows 5
bin_graph_labeling -drc_rpt_file routing.drc.rpt
save_bin_graph_file -prefix bin_graph



puts "star model end"
exit

