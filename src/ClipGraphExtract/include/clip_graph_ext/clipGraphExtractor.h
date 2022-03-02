#ifndef __GRAPH__EXTRACTOR__
#define __GRAPH__EXTRACTOR__

#include <iostream>

struct rudy_value; // This struct is already described in ../src/bingraph.h

namespace odb {
class dbDatabase;
class dbInst;
}

namespace sta {
class dbSta;
}


namespace ClipGraphExtract {

enum GraphModel {
  Star, Clique, Hybrid
};

enum EdgeWeightModel {
  A, B, C, D, E
};

class ClipGraphExtractor {
  public:
    void setDb(odb::dbDatabase* db);
    void setSta(sta::dbSta* sta);
    void init();
    void clear();
	void updateCongRUDY(rudy_value& rudyValue);
    void extract(int lx, int ly, int ux, int uy);
    void setSaveFileName (const char* fileName);
    void setSaveFilePrefix(const char* prefix);

    void setGraphModel(const char* graphModel);
    void setEdgeWeightModel(const char* edgeWeightModel);

    void extractBinGraph(int numRows);
    void extractBinGraph(int numRows, int maxLayer);
    void labelingBinGraph(const char* invRoutingReport);
    void saveBinGraph();

    void showCongestionMap();

    // defined in grid.cpp
    void initGcellGrid(int numRows, int maxLayer);
    void readRoutingReport(const char* fileName);
    void saveMapImages(const char* dirPath);

    GraphModel getGraphModel() { return graphModel_; }
    EdgeWeightModel getEdgeWeightModel() { return edgeWeightModel_; }

    ClipGraphExtractor();
    ~ClipGraphExtractor();
    odb::dbDatabase* getDb() { return db_; }
  
  private:
    odb::dbDatabase* db_;
    sta::dbSta* sta_;
    void* rTree_;
    void* inst_rTree_;
    void* wire_rTree_;
    void* via_rTree_;
    void* pin_rTree_;
    void* rudy_rTree_;
    void* drc_rTree_;
    
    GraphModel graphModel_;
    EdgeWeightModel edgeWeightModel_;
    std::string fileName_;
    std::string prefix_;
    void* binGraph_;
   
    void* grid_;




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
