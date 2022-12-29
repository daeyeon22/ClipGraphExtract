#ifndef __GRAPH__EXTRACTOR__
#define __GRAPH__EXTRACTOR__

#include <iostream>
#include <unordered_map>

namespace odb {
class dbDatabase;
class dbInst;
}

namespace sta {
class dbSta;
}


namespace ClipGraphExtract {


class Gcell;

class ClipGraphExtractor {
  public:
    void setDb(odb::dbDatabase* db);
    void setSta(sta::dbSta* sta);
    void init();
   

    void clear();
    void extract();
    //void label();
    void extract(int lx, int ly, int ux, int uy);
    void setSaveFileName (const char* fileName);
    void setSaveFilePrefix(const char* prefix);

    void setDrcReport(const char* fileName);


    // 
    void setGcellSize(int numRows);
    void setMaxRouteLayer(int maxLayer);

    //void setGraphModel(const char* graphModel);
    //void setEdgeWeightModel(const char* edgeWeightModel);

    // defined in label.cpp
    void parseDrcReport(const char* fileName);
    void parseDrcReport_(const char* fileName);
    // defined in plot.cpp
    void saveGridImages(const char* dirPath, const char* prefix);
    // defined in writer.cpp
    //void saveGraphs(const char* dirPath);
    void saveFeatures(const char* dirPath);
    void saveLabels(const char* dirPath);
    //void saveInstFeatures(const char* dirPath);
    //void saveInstLabels(const char* dirPath);

    ClipGraphExtractor();
    ~ClipGraphExtractor();
    odb::dbDatabase* getDb() { return db_; }
    sta::dbSta* getSta() { return sta_; } 
  private:
    odb::dbDatabase* db_;
    sta::dbSta* sta_;
    //GraphModel graphModel_;
    //EdgeWeightModel edgeWeightModel_;
    std::string fileName_;
    std::string prefix_;
    std::string drcRpt_;


    int numRows_;       // GCELL SIZE (= n * height of site row)
    int maxRouteLayer_; // MAX ROUTE LAYER (need to figure out routing capacity)

    // for Def
    //void* wireRtree_;
    //void* instRtree_;
    // for Grid
    //void* rsmtRtree_;
    //void* gcellRtree_;
    // for Drc
    //void* markerRtree_;

    void* grid_;

    // for initialization
    //void initGrid();
    //void initGraph();
    //void initRtree1();
    //void initRtree2();

    // For instance feature
    std::unordered_map<odb::dbInst*, bool> isCritical_;
    std::unordered_map<odb::dbInst*, bool> isClocked_;
    std::unordered_map<odb::dbInst*, double> absSlack_;
    std::unordered_map<odb::dbInst*, double> relSlack_;
    std::unordered_map<odb::dbInst*, int> instAccPoints_;
    std::unordered_map<odb::dbInst*, int> instBlkPoints_;
    std::unordered_map<odb::dbInst*, int> instBndPoints_;
    std::unordered_map<odb::dbInst*, double> powerViaDistance_;
    std::unordered_map<odb::dbInst*, double> whiteSpaceL_;
    std::unordered_map<odb::dbInst*, double> whiteSpaceR_;
    std::unordered_map<odb::dbInst*, double> whiteSpaceT_;
    std::unordered_map<odb::dbInst*, double> whiteSpaceD_;
    std::unordered_map<odb::dbInst*, double> sWireOverlap_;
    std::unordered_map<odb::dbInst*, double> stnBBox_;
    std::unordered_map<odb::dbInst*, double> cellSize_;
    std::unordered_map<odb::dbInst*, int> numCutEdges_;
    std::unordered_map<odb::dbInst*, int> numInEdges_;
    std::unordered_map<odb::dbInst*, int> numOutEdges_;
    std::unordered_map<odb::dbInst*, int> numEdges_;
    std::unordered_map<odb::dbInst*, int> cellType_;
    std::unordered_map<odb::dbInst*, double> relPosX_;
    std::unordered_map<odb::dbInst*, double> relPosY_;
    std::unordered_map<odb::dbInst*, int> numDrvs_;
    std::unordered_map<odb::dbInst*, int> row_;
    std::unordered_map<odb::dbInst*, int> col_;
    std::unordered_map<odb::dbInst*, double> xCoord_;
    std::unordered_map<odb::dbInst*, double> yCoord_;
    std::unordered_map<odb::dbInst*, Gcell*> gcell_;


};
};

//namespace GraphExtract {
//
//class GraphExtractor {
//  public:
//
//    void setNumSites(int numSites);
//    void setNumRows(int numRows);
//    void setOutFile(const char* fileName);
//
//    void init();
//    void extract();
//    void labeling(const char* fileName);
//
//  private:
//    std::string fileName_;
//
//    int numSites_;
//    int numRows_;
//};


#endif
