#ifndef __INST_GRAPH__
#define __INST_GRAPH__

#include <clip_graph_ext/clipGraphExtractor.h>
#include <map>
#include <set>
#include <unordered_map>


namespace odb {
  class dbInst;
  class dbDatabase;
}

namespace sta {
  class dbSta;
};

// only consider two-way connections
namespace ClipGraphExtract{

class Edge;
// vertex is equal to dbInst
class Vertex {
public:
  Vertex();
  Vertex(odb::dbInst* inst, 
      int id, float weight);

  odb::dbInst* inst() const;
  const std::vector<Edge*> & inEdges() const; 
  const std::vector<Edge*> & outEdges() const; 

  float weight() const; 

  void setInst(odb::dbInst* inst);
  void setWeight(float weight);

  void addInEdge(Edge* edge);
  void addOutEdge(Edge* edge);

  void setId(int id);
  int id() const;


    // 


    void setSlack(double clockPeriod, double slack);
    void setNumAccPoints(int numAccPoints);
    void setNumBlkPoints(int numBlkPoints);
    void setNumBndPoints(int numBndPoints);
    void setViaFeature(double powerViaDistance);
    void setWhiteSpaceL(double whiteSpace);
    void setWhiteSpaceR(double whiteSpace);
    void setWhiteSpaceT(double whiteSpace);
    void setWhiteSpaceD(double whiteSpace);
    void setRelPos(double relPosX, double relPosY);
    
    void setBBoxSize(double bboxSize);
    void setCellType(int cellType);
    void setSWireOverlap(double overlap);
    void setIsCrit(bool isCrit);
    void setCutEdges(int cutEdges);
    void setNumDrvs(int numDrvs);
    


    // For node feature
    double getRelPosX();
    double getRelPosY();
    double getAbsSlack();
    double getRelSlack();

    int getNumAccPoints();
    int getNumBlkPoints();
    int getNumBndPoints();
    double getWhiteSpaceL();
    double getWhiteSpaceR();
    double getWhiteSpaceD();
    double getWhiteSpaceT();

    double getBBoxSize();
    double getSWireOverlap();
    bool getIsCrit();
    int getCutEdges();
    int getCellType();


    bool isClocked();
    int getSize();
    int getDegree();
    int getNumInEdges();
    int getNumOutEdges();
    int getNumDrvs();
private:
  odb::dbInst* inst_;
  std::vector<Edge*> inEdges_;
  std::vector<Edge*> outEdges_;
  int id_;
  float weight_;

    // Node features
    double absSlack_;
    double relSlack_;

    int numAccPoints_;
    int numBlkPoints_;
    int numBndPoints_;
    int cutEdges_;
    int cellType_;
    int numDrvs_;
    bool isCrit_;
    double bboxSize_;
    double sWireOverlap_;
    double powerViaDistance_;
    double whiteSpaceL_;
    double whiteSpaceR_;
    double whiteSpaceD_;
    double whiteSpaceT_;
    double relPosX_, relPosY_;
};


// edge is inst1-inst2 connections
class Edge {
public:
  Edge();
  Edge(Vertex* from, Vertex* to, float weight);

  Vertex* from() const;
  Vertex* to() const;
  float weight() const;

  void setFrom(Vertex* inst);
  void setTo(Vertex* inst);
  void setWeight(float weight);

private:
  Vertex* from_;
  Vertex* to_;
  float weight_;
};

inline Vertex* Edge::from() const {
  return from_;
}

inline Vertex* Edge::to() const {
  return to_;
}

inline float Edge::weight() const {
  return weight_;
}

class Graph {
  public:   
    Graph();
    ~Graph();

    void saveFile(std::string fileName);
    
    void saveNodeFeaFile(std::string fileName);
    void saveEdgeIdxFile(std::string fileName);
    void saveEdgeAttFile(std::string fileName);
    
    void setDb(odb::dbDatabase* db);
    void setSta(sta::dbSta* sta);



    // 
    void setNumDrvs(std::unordered_map<odb::dbInst*, int> &numDrvs);

    void setSlack(double clockPeriod, std::unordered_map<odb::dbInst*, double> &minSlack);
    void setNumPoints(
        std::unordered_map<odb::dbInst*, int> &numAccPoints,
        std::unordered_map<odb::dbInst*, int> &numBlkPoints,
        std::unordered_map<odb::dbInst*, int> &numBndPoints
    );

    void setViaFeature(std::unordered_map<odb::dbInst*, double> &powerViaDistance);
    void setWhiteSpace(
        std::unordered_map<odb::dbInst*, double> &whiteSpaceL,
        std::unordered_map<odb::dbInst*, double> &whiteSpaceR,
        std::unordered_map<odb::dbInst*, double> &whiteSpaceD,
        std::unordered_map<odb::dbInst*, double> &whiteSpaceT
    );


    void setRelPos(
        std::unordered_map<odb::dbInst*, double> &relPosX, 
        std::unordered_map<odb::dbInst*, double> &relPosY
    );


    void setBBoxSize(std::unordered_map<odb::dbInst*, double> &bboxSize);
    void setCellType(std::unordered_map<odb::dbInst*, int> &cellType);
    void setSWireOverlap(std::unordered_map<odb::dbInst*, double> &sWireOverlap);
    void setIsCrit(std::unordered_map<odb::dbInst*, bool> &isCrit);
    void setCutEdges(std::unordered_map<odb::dbInst*, int> &cutEdges);


    void init(std::set<odb::dbInst*> &insts);

    void init(std::set<odb::dbInst*> & insts, GraphModel gModel,
            EdgeWeightModel eModel);
    void printEdgeList();
    void setGraphModel(GraphModel graphModel);
    void setEdgeWeightModel(EdgeWeightModel edgeWeightModel);
    Vertex* dbToGraph(odb::dbInst* inst);
    void print();


  private:

    GraphModel graphModel_;
    EdgeWeightModel edgeWeightModel_;

    odb::dbDatabase* db_;
    sta::dbSta* sta_;
    std::vector<Vertex> vertices_;
    std::vector<Edge> edges_;
    std::map<odb::dbInst*, Vertex*> vertexMap_;

    
    


    
    void updateVertsFromEdges();
};

}

#endif
