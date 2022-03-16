
sta::define_cmd_args "graph_extract" {
  [-graph_model star/clique]\
  [-edge_weight_model weight]\
  [-out_file fileName]
}

sta::define_cmd_args "bin_graph_extract" {
  [-out_file fileName]\
  [-num_rows numRows]
}





proc graph_extract { args } {
  sta::parse_key_args "graph_extract" args \
    keys {-graph_model -edge_weight_model -out_file} flags {}

  # default model is star
  set graph_model "star"
  if { [info exists keys(-graph_model)] } {
    set graph_model $keys(-graph_model)
    set_graph_model_cmd $graph_model
  }

  # default weight_model is xxx
  set edge_weight_model "a"
  if { [info exists keys(-edge_weight_model)] } {
    set edge_weight_model $keys(-edge_weight_model)
    set_edge_weight_model_cmd $edge_weight_model
  }

  if { ![info exists keys(-out_file)] } {
    puts "ERROR: -out_file must be used"
    return
  } else {
    set file_name $keys(-out_file)
    set_graph_extract_save_file_name_cmd $file_name
  }

  if { [ord::db_has_rows] } {
    graph_extract_init_cmd 
    set clist [split $args]
    graph_extract_cmd [lindex $clist 0] [lindex $clist 1] [lindex $clist 2] [lindex $clist 3]
    
    # graph_extract_cmd args[0] args[1] args[2] args[3]
    graph_extract_clear_cmd
  }
}


proc bin_graph_extract { args } {
    sta::parse_key_args "bin_graph_extract" args \
        keys { -num_rows -max_layer } flags {}

    
   


    if { ![info exists keys(-num_rows)]} {
        puts "ERROR: -num_rows must be used"
        return
    } else {
        set num_rows $keys(-num_rows)
    }

    if { ![info exists keys(-max_layer)] } {
        puts "ERROR: -max_layer must be used"
        return
    } else {
        set max_layer $keys(-max_layer)
    }

    graph_extract_init_cmd 
    bin_graph_extract_cmd $num_rows $max_layer

}

proc bin_graph_labeling { args } {
    sta::parse_key_args "bin_graph_labeling" args \
        keys { -drc_rpt_file } flags {}

    
    if { ![info exists keys(-drc_rpt_file)] } {
	} else {
        set drc_rpt_file $keys(-drc_rpt_file)
		if {[file exists $drc_rpt_file]} {
        	bin_graph_labeling_cmd $drc_rpt_file
		}
    }
}

proc save_bin_graph_file { args } {
    sta::parse_key_args "save_bin_graph_file" args \
        keys { -prefix } flags {}

    if { ![info exists keys(-prefix)] } {
    
    } else {
        set prefix $keys(-prefix)
        set_graph_extract_save_file_prefix_cmd $prefix
    }
    save_bin_graph_file_cmd
}

########################################################
proc read_routing_report { args } {
    sta::parse_key_args "read_routing_report" args \
        keys { -in_file }  flags {}

    if { [info exists keys(-in_file)] } {
        set in_file $keys(-in_file)
        read_routing_report_cmd $in_file
    }
}

proc construct_gcell_grid { args } {
    sta::parse_key_args "construct_gcell_grid" args \
        keys { -num_rows -max_route_layer } flags {}

    set num_rows 5 
    set max_route_layer 7

    if { [info exists keys(-num_rows)] } {
        set num_rows $keys(-num_rows)
    } 
    if { [info exists keys(-max_route_layer)]} {
        set max_route_layer $keys(-max_route_layer)
    }

    construct_gcell_grid_cmd $num_rows $max_route_layer
}

proc save_map_images { args } {
    sta::parse_key_args "save_map_imges" args \
        keys { -save_dir } flags {}

    set save_dir "./"
    if { [info exists keys(-save_dir)] } {
        set save_dir $keys(-save_dir)
    }

    save_map_images_cmd $save_dir
}

proc save_file { args } {
    sta::parse_key_args "save_file" args \
        keys { -save_dir } flags {}

    set save_dir "./"
    if { [info exists keys(-save_dir)] } {
        set save_dir $keys(-save_dir)
    }

    save_file_cmd $save_dir
}


proc analyze_congestion { } {
    analyze_congestion_cmd
}

