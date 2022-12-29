#include "clip_graph_ext/clipGraphExtractor.h"
#include <iostream>
#include <string>
#include <fstream>
#include "grid.h"


namespace ClipGraphExtract {

using namespace std;
using namespace odb;

//void ClipGraphExtractor::saveGraphs(const char* dirPath) {
//    // TODO
//    Grid* grid = (Grid*) grid_;
//
//    for(Gcell* tarGcell : grid->getGcells()) {
//        int col = tarGcell->getCol();
//        int row = tarGcell->getRow();
//        string fileName = "Col_" + to_string(col) + "_Row_" + to_string(row);
//        tarGcell->saveGraph(string(dirPath), fileName);
//    }
//}

void ClipGraphExtractor::saveFeatures(const char* dirPath) {

    // TODO

    Grid* grid = (Grid*) grid_;
    ofstream outFile;
    string fileName = "GCELL_x" + to_string(numRows_) + ".csv";
    string filePath = string(dirPath) + "/" + fileName;
    outFile.open(filePath, std::ios_base::out);
    //writeHeader(outFile);
    //ofstream outFile;
    //string attrFileName = string(dirPath) + "/GcellFeature.csv";
    //outFile.open(attrFileName, std::ios_base::out);
    int numCols=grid->getNumCols();
    int numRows=grid->getNumRows();
    int numGcells = numCols*numRows;
    int maxTechLayer = db_->getTech()->getRoutingLayerCount(); // Intell22nm = 0~8
    outFile << "col" << ","
            << "row" << ","
            << "cell_den" << ","
            << "pin_den" << ","
            << "rudy" << ","
            << "lnet_rudy" << ","
            << "gnet_rudy" << ","
            << "snet_rudy" << ","
            << "wire_den_rsmt" << ","
            << "lnet_den_rsmt" << ","
            << "gnet_den_rsmt" << ","
            << "chan_den_rsmt" << ","
            << "chan_den_v_rsmt" << ","
            << "chan_den_h_rsmt" << ",";
    for(int layer = 1; layer < maxTechLayer; layer++)
        outFile << "wire_den_egr" + to_string(layer) << ",";
    for(int layer = 1; layer < maxTechLayer; layer++)
        outFile << "chan_den_egr" + to_string(layer) << ",";
    for(int layer = 1; layer < maxTechLayer-1; layer++)
        outFile << "via_den_egr" + to_string(layer) << ",";
    outFile << "lnet_den_egr" << ","
            << "gnet_den_egr" << "," 
            << "chan_den_egr" << "," 
            << "chan_den_v_egr" << "," 
            << "chan_den_h_egr" << ","
            << "avg_terms" << "," 
            << "num_insts" << "," 
            << "num_terms" << "," 
            << "num_nets" << ","
            << "num_gnets" << "," 
            << "num_lnets" << "," 
            << "clk_ratio" << ","
            << "wns" << ","
            << "tns" << "," 
            << "num_drvs" << endl;

    //vector<string> columns {
    //    "col", "row", "cell_den" , "pin_den", 
    //    "rudy", "lnet_rudy", "gnet_rudy", "snet_rudy", 
    //    "wire_den_rsmt", "lnet_den_rsmt", "gnet_den_rsmt",
    //    "chan_den_rsmt", "chan_den_v_rsmt", "chan_den_h_rsmt",
    //    "wire_den_egr1", "wire_den_egr2", "wire_den_egr3", "wire_den_egr4", 
    //    "wire_den_egr5", "wire_den_egr6", "wire_den_egr7", "wire_den_egr8", 
    //    "via_den_egr1", "via_den_egr2", "via_den_egr3", 
    //    "via_den_egr4", "via_den_egr5", "via_den_egr6", "via_den_egr7", 
    //    "lnet_den_egr", "gnet_den_egr", "chan_den_egr", "chan_den_v_egr", "chan_den_h_egr",
    //    "avg_terms", "num_insts", "num_terms", "num_nets","num_gnets", "num_lnets", "clk_ratio", 
    //    "wns", "tns"
    //};
    //for(string column : columns)
    //    outFile << column << ",";
    //outFile << endl;
   
    for(Gcell* tarGcell : grid->getGcells()) {
        outFile << tarGcell->getCol() << ","
                << tarGcell->getRow() << ","
                //<< 1.0 * tarGcell->getCol() / numCols << "," 
                //<< 1.0 * tarGcell->getRow() / numRows << ","
                //<< 1.0 / numGcells << ","
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
                << tarGcell->getTNS() << "," 
                << tarGcell->getNumMarkers() << endl;
    }
    outFile.close();
    cout << "End writing file." << endl;
}

// USE
void ClipGraphExtractor::saveLabels(const char* dirPath) {
    // TODO
    Grid* grid = (Grid*) grid_;
    ofstream outFile;
    string fileName = "GCELL_x" + to_string(numRows_) + ".csv";
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


/*
void writeHeader(ofstream& outFile) {

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

    for(int i=0; i < fieldData.size(); i++) {
        string fieldName = fieldData[i];
        outFile << fieldName;
        if(i!=fieldData.size()-1)
            outFile << ",";
    }

    outFile << endl;
}


void writeData(ofstream& outFile, Grid* tarGrid, int tarCol, int tarRow) {
    int numCols = tarGrid->getNumCols();
    int numRows = tarGrid->getNumRows();
    int numGcells = numCols*numRows;
    outFile << tarCol << "," << tarRow << ",";


    Gcell* tarGcell = tarGrid->getGcell(tarCol, tarRow);
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

    outFile << endl;
}


// USE
void ClipGraphExtractor::saveFeatures(const char* dirPath) {

    // TODO
    Grid* grid = (Grid*) grid_;
    ofstream outFile;
    string fileName = "GCELL_x" + to_string(numRows_) + ".csv";
    string filePath = string(dirPath) + "/" + fileName;
    outFile.open(filePath, std::ios_base::out);
    writeHeader(outFile);
    
    int numCols=grid->getNumCols();
    int numRows=grid->getNumRows();

    for(Gcell* tarGcell : grid->getGcells()) {
        int tarCol = tarGcell->getCol();
        int tarRow = tarGcell->getRow();
        writeData(outFile, grid, tarCol, tarRow);
    }

    outFile.close();
    cout << "End writing file." << endl;
}


// USE
void ClipGraphExtractor::saveInstFeatures(const char* dirPath) {

    string fileName = "INST_x" + to_string(numRows_) + ".csv";
    string filePath = string(dirPath) + "/" + fileName;
    ofstream outFile;
    outFile.open(filePath, std::ios_base::out);

    outFile << "inst_name" << ","
            << "col" << ","
            << "row" << ","
            << "x_coord" << ","
            << "y_coord" << ","
            // INST
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
            << "num_edges" << ",";
    int maxTechLayer = db_->getTech()->getRoutingLayerCount(); // Intell22nm = 0~8
    outFile << "cell_den" << ","
            << "pin_den" << ","
            << "rudy" << ","
            << "lnet_rudy" << ","
            << "gnet_rudy" << ","
            << "snet_rudy" << ","
            << "wire_den_rsmt" << ","
            << "lnet_den_rsmt" << ","
            << "gnet_den_rsmt" << ","
            << "chan_den_rsmt" << ","
            << "chan_den_v_rsmt" << ","
            << "chan_den_h_rsmt" << ",";
    for(int layer = 1; layer < maxTechLayer; layer++)
        outFile << "wire_den_egr" + to_string(layer) << ",";
    for(int layer = 1; layer < maxTechLayer; layer++)
        outFile << "chan_den_egr" + to_string(layer) << ",";
    for(int layer = 1; layer < maxTechLayer-1; layer++)
        outFile << "via_den_egr" + to_string(layer) << ",";
    outFile << "lnet_den_egr" << ","
            << "gnet_den_egr" << "," 
            << "chan_den_egr" << "," 
            << "chan_den_v_egr" << "," 
            << "chan_den_h_egr" << ","
            << "avg_terms" << "," 
            << "num_insts" << "," 
            << "num_terms" << "," 
            << "num_nets" << ","
            << "num_gnets" << "," 
            << "num_lnets" << "," 
            << "clk_ratio" << ","
            << "wns" << ","
            << "tns" << endl;

    Grid* grid = (Grid*)grid_;
    dbBlock* block = db_->getChip()->getBlock();
    for(dbInst *tarInst : block->getInsts()) {
        Gcell* tarGcell = gcell_[tarInst];
        
        outFile << tarInst->getName() << ","
                //<< relPosX_[tarInst] << ","
                //<< relPosY_[tarInst] << ","
                << tarGcell->getCol() << ","
                << tarGcell->getRow() << ","
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
                << numEdges_[tarInst] << ","
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
                << tarGcell->getTNS() << endl;
                //<< tarGcell->getCellUtil() << ","
                //<< tarGcell->getPinUtil() << ","
                //<< tarGcell->getRUDY() << ","
                //<< tarGcell->getLNetRUDY() << ","
                //<< tarGcell->getGNetRUDY() << ","
                //<< tarGcell->getSNetRUDY() << ","
                //<< tarGcell->getWireUtil(ModelType::TREE) << ","
                //<< tarGcell->getLNetUtil(ModelType::TREE) << ","
                //<< tarGcell->getGNetUtil(ModelType::TREE) << ","
                //<< tarGcell->getChanUtil(ModelType::TREE) << ","
                //<< tarGcell->getChanUtilV(ModelType::TREE) << ","
                //<< tarGcell->getChanUtilH(ModelType::TREE) << ","
                //<< tarGcell->getWireUtil(ModelType::ROUTE) << ","
                //<< tarGcell->getLNetUtil(ModelType::ROUTE) << ","
                //<< tarGcell->getGNetUtil(ModelType::ROUTE) << ","
                //<< tarGcell->getChanUtil(ModelType::ROUTE) << ","
                //<< tarGcell->getChanUtilV(ModelType::ROUTE) << ","
                //<< tarGcell->getChanUtilH(ModelType::ROUTE) << ","
                //<< tarGcell->getAvgTerms() << ","
                //<< tarGcell->getNumInsts() << ","
                //<< tarGcell->getNumTerms() << ","
                //<< tarGcell->getNumNets() << ","
                //<< tarGcell->getNumGNets() << ","
                //<< tarGcell->getNumLNets() << ","
                //<< tarGcell->getClkRatio() << ","
                //<< tarGcell->getWNS() << ","
                //<< tarGcell->getTNS() << endl;
    }
}


// USE
void ClipGraphExtractor::saveInstLabels(const char* dirPath) {
    string fileName = "INST_x" + to_string(numRows_) + ".csv";
    string filePath = string(dirPath) + "/" + fileName;
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
*/


};



