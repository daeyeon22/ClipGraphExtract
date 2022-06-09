#include "clip_graph_ext/clipGraphExtractor.h"
#include <iostream>
#include <string>
#include <fstream>
#include "grid.h"


namespace ClipGraphExtract {

using namespace std;
using namespace odb;

void ClipGraphExtractor::saveGraphs(const char* dirPath) {
    // TODO
    Grid* grid = (Grid*) grid_;

    for(Gcell* tarGcell : grid->getGcells()) {
        int col = tarGcell->getCol();
        int row = tarGcell->getRow();
        string fileName = "Col_" + to_string(col) + "_Row_" + to_string(row);
        tarGcell->saveGraph(string(dirPath), fileName);
    }
}

void ClipGraphExtractor::saveFeatures(const char* dirPath) {

    // TODO
    Grid* grid = (Grid*) grid_;
    int maxTechLayer = db_->getTech()->getRoutingLayerCount(); // Intell22nm = 0~8
    ofstream outFile;
    string attrFileName = string(dirPath) + "/GcellFeature.csv";
    outFile.open(attrFileName, std::ios_base::out);
    int numCols=grid->getNumCols();
    int numRows=grid->getNumRows();
    int numGcells = numCols*numRows;


    outFile << "Col" << ","
            << "Row" << ","
            << "RelPosX" << ","
            << "RelPosY" << ","
            << "RelArea" << ","
            << "CellDen" << ","
            << "PinDen" << ","
            << "RUDY" << ","
            << "LNetRUDY" << ","
            << "GNetRUDY" << ","
            << "SNetRUDY" << ","
            << "WireDen(RSMT)" << ","
            << "LNetDen(RSMT)" << ","
            << "GNetDen(RSMT)" << ","
            << "ChanDen(RSMT)" << ","
            << "ChanDenV(RSMT)" << ","
            << "ChanDenH(RSMT)" << ","
            << "Wire1Den(EGR)" << ","
            << "Wire2Den(EGR)" << ","
            << "Wire3Den(EGR)" << ","
            << "Wire4Den(EGR)" << ","
            << "Wire5Den(EGR)" << ","
            << "Wire6Den(EGR)" << ","
            << "Wire7Den(EGR)" << ","
            << "Wire8Den(EGR)" << ","
            << "Chan1Den(EGR)" << ","
            << "Chan2Den(EGR)" << ","
            << "Chan3Den(EGR)" << ","
            << "Chan4Den(EGR)" << ","
            << "Chan5Den(EGR)" << ","
            << "Chan6Den(EGR)" << ","
            << "Chan7Den(EGR)" << ","
            << "Chan8Den(EGR)" << ","
            << "Via12Den(EGR)" << ","
            << "Via23Den(EGR)" << ","
            << "Via34Den(EGR)" << ","
            << "Via45Den(EGR)" << ","
            << "Via56Den(EGR)" << ","
            << "Via67Den(EGR)" << ","
            << "Via78Den(EGR)" << ","
            << "LNetDen(EGR)" << ","
            << "GNetDen(EGR)" << ","
            << "ChanDenV(EGR)" << ","
            << "ChanDenH(EGR)" << ","
            << "AvgTerms" << ","
            << "NumInsts" << ","
            << "NumTerms" << ","
            << "NumNets" << ","
            << "NumGNets" << ","
            << "NumLNets" << ","
            << "ClkRatio" << "," 
            << "TNS" << endl;

    for(Gcell* tarGcell : grid->getGcells()) {
        outFile << tarGcell->getCol() << ","
                << tarGcell->getRow() << ","
                << 1.0 * tarGcell->getCol() / numCols << "," 
                << 1.0 * tarGcell->getRow() / numRows << ","
                << 1.0 / numGcells << ","
                << tarGcell->getCellUtil() << ","
                << tarGcell->getPinUtil() << ","
                << tarGcell->getRUDY() << ","
                << tarGcell->getLNetRUDY() << ","
                << tarGcell->getGNetRUDY() << ","
                << tarGcell->getSNetRUDY() << ","
                << tarGcell->getWireUtil(ModelType::TREE) << ","
                << tarGcell->getLNetUtil(ModelType::TREE) << ","
                << tarGcell->getGNetUtil(ModelType::TREE) << ","
                << tarGcell->getChanUtil(ModelType::TREE) << ","
                << tarGcell->getChanUtilV(ModelType::TREE) << ","
                << tarGcell->getChanUtilH(ModelType::TREE) << ",";
        for(int layer = 1; layer < maxTechLayer; layer++)
            outFile << tarGcell->getWireUtil(layer, ModelType::ROUTE) << ",";
        for(int layer = 1; layer < maxTechLayer; layer++)
            outFile << tarGcell->getChanUtil(layer, ModelType::ROUTE) << ",";
        for(int layer = 1; layer < maxTechLayer-1; layer++)
            outFile << tarGcell->getViaUtil(layer) << ",";
        outFile << tarGcell->getLNetUtil(ModelType::ROUTE) << ","
                << tarGcell->getGNetUtil(ModelType::ROUTE) << ","
                << tarGcell->getChanUtilV(ModelType::ROUTE) << ","
                << tarGcell->getChanUtilH(ModelType::ROUTE) << ","
                << tarGcell->getAvgTerms() << ","
                << tarGcell->getNumInsts() << ","
                << tarGcell->getNumTerms() << ","
                << tarGcell->getNumNets() << ","
                << tarGcell->getNumGNets() << ","
                << tarGcell->getNumLNets() << ","
                << tarGcell->getClkRatio() << ","
                << tarGcell->getTNS() << endl;
    }
    outFile.close();
    cout << "End writing file." << endl;
}


// USE
void ClipGraphExtractor::saveLabels(const char* dirPath) {
    // TODO
    Grid* grid = (Grid*) grid_;
    ofstream outFile;
    string fileName = "x" + to_string(numRows_) + ".csv";
    string filePath = string(dirPath) + "/" + fileName;
    outFile.open(filePath, std::ios_base::out);

    vector<string> header{
        "col", "row", "cell_den_dr", "pin_den_dr", 
        "wire_den_dr", "lnet_den_dr", "gnet_den_dr", 
        "chan_den_dr", "chan_den_v_dr", "chan_den_h_dr",
        "buf_den_dr", 
        "num_drvs", "num_lnet_drvs", "num_gnet_drvs", 
        "num_inst_drvs", "num_net_drvs" ,"clk_ratio", "wns", "tns"
    };
    
    for(int i=0; i < header.size(); i++) {
        outFile << header[i];
        if(i!=header.size()-1)
            outFile << ",";
        else
            outFile << endl;
    }
    
    int numLNetMarkers, numGNetMarkers, numInstMarkers;
    for(Gcell* tarGcell : grid->getGcells()) {
        // get #drvs
        tarGcell->getNumMarkers(numLNetMarkers, numGNetMarkers, numInstMarkers);

        outFile << tarGcell->getCol() << ","
                << tarGcell->getRow() << ","
                << tarGcell->getCellUtil() << ","
                << tarGcell->getPinUtil() << ","
                << tarGcell->getWireUtil(ModelType::ROUTE) << ","
                << tarGcell->getLNetUtil(ModelType::ROUTE) << ","
                << tarGcell->getGNetUtil(ModelType::ROUTE) << ","
                << tarGcell->getChanUtil(ModelType::ROUTE) << ","
                << tarGcell->getChanUtilV(ModelType::ROUTE) << ","
                << tarGcell->getChanUtilH(ModelType::ROUTE) << ","
                << tarGcell->getBufferUtil() << ","
                << tarGcell->getNumMarkers() << ","
                << numLNetMarkers << "," 
                << numGNetMarkers << ","
                << numInstMarkers << ","
                << numLNetMarkers + numGNetMarkers << ","
                << tarGcell->getClkRatio() << ","
                << tarGcell->getWNS() << ","
                << tarGcell->getTNS() << endl;

    }
    outFile.close();
    cout << "End writing file." << endl;
}


void writeHeader(ofstream& outFile, int numHops) {

    vector<string> fieldId { 
        "col", "row"
    };
    vector<string> fieldData {
        "rel_pos_x", "rel_pos_y", "rel_area",
        "cell_den" , "pin_den", 
        "rudy", "lnet_rudy", "gnet_rudy", "snet_rudy", 
        "wire_den_rsmt", "lnet_den_rsmt", "gnet_den_rsmt",
        "chan_den_rsmt", "chan_den_v_rsmt", "chan_den_h_rsmt",
        "wire_den_egr", "lnet_den_egr", "gnet_den_egr",
        "chan_den_egr", "chan_den_v_egr", "chan_den_h_egr",
        "avg_terms", "num_insts", "num_terms", "num_nets",
        "num_gnets", "num_lnets", "clk_ratio", 
        "wns", "tns"
    };


    
    for(int i=0; i < fieldId.size(); i++) {
        string fieldName = fieldId[i];
        outFile << fieldName << ",";
    }

    for(int x=-numHops;x<=numHops;x++) {
        for(int y=-numHops;y<=numHops;y++) {
            string prefix = "";
            if(y < 0)
                prefix+="s"+to_string(abs(y));
            else if(y>0)
                prefix+="n"+to_string(abs(y));
            if(x < 0)
                prefix+="w"+to_string(abs(x));
            else if(x > 0)
                prefix+="e"+to_string(abs(x));
            for(int i=0; i < fieldData.size(); i++) {
                string fieldName = (prefix=="")? fieldData[i] : prefix + "_" + fieldData[i];
                outFile << fieldName;
                if(!(x==numHops && y==numHops && i==fieldData.size()-1))
                        outFile << ",";
            }
        }
    }
  
    outFile << endl;
}


void writeData(ofstream& outFile, Grid* tarGrid, int tarCol, int tarRow, int numHops) {

    int numCols = tarGrid->getNumCols();
    int numRows = tarGrid->getNumRows();
    int numGcells = numCols*numRows;

    outFile << tarCol << "," << tarRow << ",";


    for(int dx=-numHops;dx<=numHops;dx++) {
        for(int dy=-numHops;dy<=numHops;dy++) {
            int col = tarCol + dx;
            int row = tarRow + dy;
            Gcell* tarGcell = tarGrid->getGcell(col, row);
            if(tarGcell==NULL) {
                outFile << 0 << "," // 1
                    << 0 << "," // 2
                    << 0 << "," // 3
                    << 0 << "," // 4
                    << 0 << "," // 5
                    << 0 << "," // 6
                    << 0 << "," // 7
                    << 0 << "," // 8
                    << 0 << "," // 9
                    << 0 << "," // 10
                    << 0 << "," // 11
                    << 0 << "," // 12
                    << 0 << "," // 13
                    << 0 << "," // 14
                    << 0 << "," // 15 
                    << 0 << "," // 16
                    << 0 << "," // 17
                    << 0 << "," // 18
                    << 0 << "," // 19
                    << 0 << "," // 20
                    << 0 << "," // 21
                    << 0 << "," // 22
                    << 0 << "," // 23
                    << 0 << "," // 24
                    << 0 << "," // 25
                    << 0 << "," // 26
                    << 0 << "," // 27
                    << 0 << "," // 28
                    << 0 << "," // 29
                    << 0; // 30
            } else {
                outFile << 1.0 * tarGcell->getCol() / numCols << "," 
                    << 1.0 * tarGcell->getRow() / numRows << ","
                    << 1.0 / numGcells << ","
                    << tarGcell->getCellUtil() << ","
                    << tarGcell->getPinUtil() << ","
                    << tarGcell->getRUDY() << ","
                    << tarGcell->getLNetRUDY() << ","
                    << tarGcell->getGNetRUDY() << ","
                    << tarGcell->getSNetRUDY() << ","
                    << tarGcell->getWireUtil(ModelType::TREE) << ","
                    << tarGcell->getLNetUtil(ModelType::TREE) << ","
                    << tarGcell->getGNetUtil(ModelType::TREE) << ","
                    << tarGcell->getChanUtil(ModelType::TREE) << ","
                    << tarGcell->getChanUtilV(ModelType::TREE) << ","
                    << tarGcell->getChanUtilH(ModelType::TREE) << ","
                    << tarGcell->getWireUtil(ModelType::ROUTE) << ","
                    << tarGcell->getLNetUtil(ModelType::ROUTE) << ","
                    << tarGcell->getGNetUtil(ModelType::ROUTE) << ","
                    << tarGcell->getChanUtil(ModelType::ROUTE) << ","
                    << tarGcell->getChanUtilV(ModelType::ROUTE) << ","
                    << tarGcell->getChanUtilH(ModelType::ROUTE) << ","
                    << tarGcell->getAvgTerms() << ","
                    << tarGcell->getNumInsts() << ","
                    << tarGcell->getNumTerms() << ","
                    << tarGcell->getNumNets() << ","
                    << tarGcell->getNumGNets() << ","
                    << tarGcell->getNumLNets() << ","
                    << tarGcell->getClkRatio() << ","
                    << tarGcell->getWNS() << ","
                    << tarGcell->getTNS(); 
            }
        
            if(dx!=numHops || dy!=numHops)
                outFile << ",";
        }
    }

    outFile << endl;
}

// USE
void ClipGraphExtractor::saveFeatures(const char* dirPath, int numHops) {

    // TODO
    Grid* grid = (Grid*) grid_;
    ofstream outFile;
    string fileName = "x" + to_string(numRows_) + ".csv";
    string filePath = string(dirPath) + "/" + fileName;
    outFile.open(filePath, std::ios_base::out);
    writeHeader(outFile, numHops);
    
    int numCols=grid->getNumCols();
    int numRows=grid->getNumRows();

    for(Gcell* tarGcell : grid->getGcells()) {
        int tarCol = tarGcell->getCol();
        int tarRow = tarGcell->getRow();
        writeData(outFile, grid, tarCol, tarRow, numHops);
//    for(int tarCol=0; tarCol<numCols;tarCol++) {
//        for(int tarRow=0; tarRow<numRows;tarRow++) {
//        }
    }

    outFile.close();
    cout << "End writing file." << endl;
}

// USE
void ClipGraphExtractor::saveInstFeatures(const char* dirPath) {


    string filePath = string(dirPath) + "/instFeature.csv";
    ofstream outFile;
    outFile.open(filePath, std::ios_base::out);


    outFile << "inst_name" << ","
            //<< "rel_pos_x" << ","
            //<< "rel_pos_y" << ","
            << "col" << ","
            << "row" << ","
            << "x_coord" << ","
            << "y_coord" << ","
            << "abs_slack" << ","
            << "rel_slack" << ","
            << "inst_acc_points" << ","
            << "inst_blk_points" << ","
            << "inst_bnd_points" << ","
            << "avg_acc_points" << ","
            << "power_via_distance" << ","
            << "white_space_l" << ","
            << "white_space_r" << ","
            << "white_space_d" << ","
            << "white_space_t" << ","
            << "white_space_h" << ","
            << "white_space_v" << ","
            << "swire_overlap" << ","
            << "stn_bbox" << ","
            << "cell_type" << ","
            << "cell_size" << ","
            << "is_clocked" << ","
            << "num_cut_edges" << ","
            << "num_in_edges" << ","
            << "num_out_edges" << ","
            << "num_edges" << endl;


    Grid* grid = (Grid*)grid_;
    dbBlock* block = db_->getChip()->getBlock();
    for(dbInst *tarInst : block->getInsts()) {
        outFile << tarInst->getName() << ","
                //<< relPosX_[tarInst] << ","
                //<< relPosY_[tarInst] << ","
                << col_[tarInst] << ","
                << row_[tarInst] << ","
                << xCoord_[tarInst] << ","
                << yCoord_[tarInst] << ","
                << absSlack_[tarInst] << ","
                << relSlack_[tarInst] << ","
                << instAccPoints_[tarInst] << ","
                << instBlkPoints_[tarInst] << ","
                << instBndPoints_[tarInst] << ","
                << 1.0*instAccPoints_[tarInst] / (numInEdges_[tarInst]+1) << ","
                << powerViaDistance_[tarInst] << ","
                << whiteSpaceL_[tarInst] << ","
                << whiteSpaceR_[tarInst] << ","
                << whiteSpaceD_[tarInst] << ","
                << whiteSpaceT_[tarInst] << ","
                << (whiteSpaceL_[tarInst] + whiteSpaceR_[tarInst] ) / 2 << ","
                << (whiteSpaceT_[tarInst] + whiteSpaceD_[tarInst] ) / 2 << ","
                << sWireOverlap_[tarInst] << ","
                << stnBBox_[tarInst] << ","
                << cellType_[tarInst] << ","
                << cellSize_[tarInst] << ","
                << isClocked_[tarInst] << ","
                << numCutEdges_[tarInst] << ","
                << numInEdges_[tarInst] << ","
                << numOutEdges_[tarInst] << ","
                << numEdges_[tarInst] << endl;
    }
}


// USE
void ClipGraphExtractor::saveInstLabels(const char* dirPath) {
    string filePath = string(dirPath) + "/instLabel.csv";
    ofstream outFile;
    outFile.open(filePath, std::ios_base::out);

    outFile << "inst_name" << ","
            << "num_drvs" << endl;

    dbBlock* block = db_->getChip()->getBlock();
    Grid* grid = (Grid*)grid_;
    for(dbInst *tarInst : block->getInsts()) {
        outFile << tarInst->getName() << ","
                << numDrvs_[tarInst] << endl;
    }

}



};



