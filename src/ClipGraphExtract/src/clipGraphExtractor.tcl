
#sta::define_cmd_args "graph_extract" {
#  [-graph_model star/clique]\
#  [-edge_weight_model weight]\
#  [-out_file fileName]
#}
#
#sta::define_cmd_args "bin_graph_extract" {
#  [-out_file fileName]\
#  [-num_rows numRows]
#}
#
#
#
#proc graph_extract { args } {
#  sta::parse_key_args "graph_extract" args \
#    keys {-graph_model -edge_weight_model -out_file} flags {}
#
#  # default model is star
#  set graph_model "star"
#  if { [info exists keys(-graph_model)] } {
#    set graph_model $keys(-graph_model)
#    set_graph_model_cmd $graph_model
#  }
#
#  # default weight_model is xxx
#  set edge_weight_model "a"
#  if { [info exists keys(-edge_weight_model)] } {
#    set edge_weight_model $keys(-edge_weight_model)
#    set_edge_weight_model_cmd $edge_weight_model
#  }
#
#  if { ![info exists keys(-out_file)] } {
#    puts "ERROR: -out_file must be used"
#    return
#  } else {
#    set file_name $keys(-out_file)
#    set_graph_extract_save_file_name_cmd $file_name
#  }
#
#  if { [ord::db_has_rows] } {
#    graph_extract_init_cmd 
#    set clist [split $args]
#    graph_extract_cmd [lindex $clist 0] [lindex $clist 1] [lindex $clist 2] [lindex $clist 3]
#    
#    # graph_extract_cmd args[0] args[1] args[2] args[3]
#    graph_extract_clear_cmd
#  }
#}
#

########################################################


proc parse_drc_report { args } {
    sta::parse_key_args "parse_drc_report" args \
        keys { -in_file }  flags {}

    if { [info exists keys(-in_file)] } {
        set in_file $keys(-in_file)
        parse_drc_report_cmd $in_file
    }
}

proc graph_extract_init { args } {
    sta::parse_key_args "graph_extract_init" args \
        keys { -num_rows -max_route_layer -drc_report } flags {}

    
    if { [info exists keys(-num_rows)] } {
        set num_rows $keys(-num_rows)
        set_gcell_size_cmd $num_rows
    } 
    if { [info exists keys(-max_route_layer)]} {
        set max_route_layer $keys(-max_route_layer)
        set_max_route_layer_cmd $max_route_layer
    }
    if { [info exists keys(-drc_report)] } {
        set drc_report $keys(-drc_report)
        set_drc_report_cmd $drc_report
    }

    graph_extract_init_cmd
}

proc graph_extract {} {
    graph_extract_cmd
}


proc save_grid_images { args } {
    sta::parse_key_args "save_map_imges" args \
        keys { -save_dir -prefix } flags {}

    set save_dir "./"
    set prefix ""
    if { [info exists keys(-save_dir)] } {
        set save_dir $keys(-save_dir)
    }
    if { [info exists keys(-prefix)] } {
        set prefix $keys(-prefix)
    }

    save_grid_images_cmd $save_dir $prefix
}

proc save_features { args } {
    sta::parse_key_args "save_features" args \
        keys { -save_dir } flags {}
    set save_dir "./"
    if { [info exists keys(-save_dir)] } {
        set save_dir $keys(-save_dir)
        if { ![file isdirectory $save_dir] } {
            puts "ERROR: -save_dir ${save_dir} does not exist"
            return
        }       
    }
    save_features_cmd $save_dir
}

#proc save_graphs { args } {
#    sta::parse_key_args "save_graphs" args \
#        keys { -save_dir } flags {}
#    set save_dir "./"
#    if { [info exists keys(-save_dir)] } {
#        set save_dir $keys(-save_dir)
#        if { ![file isdirectory $save_dir] } {
#            puts "ERROR: -save_dir ${save_dir} does not exist"
#            return
#        }       
#    }
#    save_graphs_cmd $save_dir
#}

proc save_labels { args } {
    sta::parse_key_args "save_labels" args \
        keys { -save_dir } flags {}
    set save_dir "./"
    if { [info exists keys(-save_dir)] } {
        set save_dir $keys(-save_dir)
        if { ![file isdirectory $save_dir] } {
            puts "ERROR: -save_dir ${save_dir} does not exist"
            return
        }       
    }
    save_labels_cmd $save_dir
}

#proc save_inst_features { args } {
#    sta::parse_key_args "save_inst_features" args \
#        keys { -save_dir } flags {}
#    set save_dir "./"
#    if { [info exists keys(-save_dir)] } {
#        set save_dir $keys(-save_dir)
#        if { ![file isdirectory $save_dir] } {
#            puts "ERROR: -save_dir ${save_dir} does not exist"
#            return
#        }       
#    }
#    save_inst_features_cmd $save_dir
#}
#
#proc save_inst_labels { args } {
#    sta::parse_key_args "save_inst_labels" args \
#        keys { -save_dir } flags {}
#    set save_dir "./"
#    if { [info exists keys(-save_dir)] } {
#        set save_dir $keys(-save_dir)
#        if { ![file isdirectory $save_dir] } {
#            puts "ERROR: -save_dir ${save_dir} does not exist"
#            return
#        }       
#    }
#    save_inst_labels_cmd $save_dir
#}
