def refParser(file):
    f = open(file, "r")
    lines = f.readlines()
    f.close()

    return lines


def writeScr(lines):
    Wlines = lines
    for line in lines:
        vocas = line.split(" ")
        
        if vocas[0] == "set":
            if vocas[1] == "LIB_PATH":

            elif vocas[1] == "LEF_PATH":
            
            elif vocas[1] == "DEF_PATH":
            
            elif vocas[1] == "SPEF_PATH":
            
            elif vocas[1] == "RPT_PATH":
            
            elif vocas[1] == "SDC_PATH":
            
            elif vocas[1] == "NUM_ROWS":
            
            elif vocas[1] == "GRAPH_MODEL":
            
            elif vocas[1] == "EDGE_WEIGHT_MODEL":
            
            elif vocas[1] == "RMAX_LYR":
            
            elif vocas[1] == "IMG_HOME":
            
            elif vocas[1] == "LABEL_HOME":
            
            elif vocas[1] == "GRAPH_HOME":
            
            elif vocas[1] == "FEATURE_HOME":



lines = refParser("2_label.tcl")


