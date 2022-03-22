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
    ofstream outFile;
    string attrFileName = string(dirPath) + "/GcellFeature.csv";
    outFile.open(attrFileName, std::ios_base::out);
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
            << "WireDen(EGR)" << ","
            << "LNetDen(EGR)" << ","
            << "GNetDen(EGR)" << ","
            << "ChanDen(EGR)" << ","
            << "ChanDenV(EGR)" << ","
            << "ChanDenH(EGR)" << ","
            << "AvgTerms" << ","
            << "NumInsts" << ","
            << "NumTerms" << ","
            << "NumNets" << ","
            << "NumGNets" << ","
            << "NumLNets" << ","
            << "ClkRatio" << endl;

    int numCols = grid->getNumCols();
    int numRows = grid->getNumRows();
    int numGcells = numCols*numRows;

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
                << tarGcell->getClkRatio() << endl;
    }
    outFile.close();
    cout << "End writing file." << endl;
}

void ClipGraphExtractor::saveLabels(const char* dirPath) {
    // TODO
    Grid* grid = (Grid*) grid_;
    ofstream outFile;
    string attrFileName = string(dirPath) + "/GcellLabel.csv";
    outFile.open(attrFileName, std::ios_base::out);
    outFile << "Col" << ","
            << "Row" << ","
            << "CellDen(DR)" << ","
            << "PinDen(DR)" << ","
            << "WireDen(DR)" << ","
            << "LNetDen(DR)" << ","
            << "GNetDen(DR)" << ","
            << "ChanDen(DR)" << ","
            << "ChanDenV(DR)" << ","
            << "ChanDenH(DR)" << ","
            << "BufDen(DR)" << ","
            << "NumDRVs" << ","
            << "NumLNetDRVs" << ","
            << "NumGNetDRVs" << ","
            << "NumInstDRVs" << ","
            << "ClkRatio" << ","
            << "TNS" << endl;
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
                << tarGcell->getClkRatio() << ","
                << tarGcell->getTNS() << endl;

    }
    outFile.close();
    cout << "End writing file." << endl;
}




};



