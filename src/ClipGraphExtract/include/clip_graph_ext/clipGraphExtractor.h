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
    void extract();
    void extract(int lx, int ly, int ux, int uy);
    void setSaveFileName (const char* fileName);
    void setSaveFilePrefix(const char* prefix);
    
    // 
    void setGcellSize(int numRows);
    void setMaxRouteLayer(int maxLayer);

    void setGraphModel(const char* graphModel);
    void setEdgeWeightModel(const char* edgeWeightModel);

    // defined in label.cpp
    void parseDrcReport(const char* fileName);
    // defined in plot.cpp
    void saveGridImages(const char* dirPath);
    // defined in writer.cpp
    void saveGraphs(const char* dirPath);
    void saveFeatures(const char* dirPath);
    void saveLabels(const char* dirPath);


    GraphModel getGraphModel() { return graphModel_; }
    EdgeWeightModel getEdgeWeightModel() { return edgeWeightModel_; }

    ClipGraphExtractor();
    ~ClipGraphExtractor();
    odb::dbDatabase* getDb() { return db_; }
  
  private:
    odb::dbDatabase* db_;
    sta::dbSta* sta_;
    GraphModel graphModel_;
    EdgeWeightModel edgeWeightModel_;
    std::string fileName_;
    std::string prefix_;
    
    int numRows_;       // GCELL SIZE (= n * height of site row)
    int maxRouteLayer_; // MAX ROUTE LAYER (need to figure out routing capacity)


    void* grid_;

    void initGrid();
    void initGraph();



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
