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

    void setMinSlack(double slack);
    void setNumAccPoints(int numAccPoints);
    void setNumBlkPoints(int numBlkPoints);
    void setNumBndPoints(int numBndPoints);
    void setWhiteSpaceL(int whiteSpaceL);
    void setWhiteSpaceR(int whiteSpaceL);


    // For node feature
    double getMinSlack();
    int getNumAccPoints();
    int getNumBlkPoints();
    int getNumBndPoints();
    int getWhiteSpaceL();
    int getWhiteSpaceR();

    bool isClocked();
    int getSize();
    int getDegree();
    int getNumInEdges();
    int getNumOutEdges();
private:
  odb::dbInst* inst_;
  std::vector<Edge*> inEdges_;
  std::vector<Edge*> outEdges_;
  int id_;
  float weight_;

    double minSlack_;
    int numAccPoints_;
    int numBlkPoints_;
    int numBndPoints_;

    int whiteSpaceL_;
    int whiteSpaceR_;

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
    void setMinSlack(std::unordered_map<odb::dbInst*, double> &minSlack);
    void setNumAccPoints(std::unordered_map<odb::dbInst*, int> &numAccPoints);
    void setNumBlkPoints(std::unordered_map<odb::dbInst*, int> &numBlkPoints);
    void setNumBndPoints(std::unordered_map<odb::dbInst*, int> &numBndPoints);

    void setWhiteSpaceL(std::unordered_map<odb::dbInst*, int> &whiteSpaceL);
    void setWhiteSpaceR(std::unordered_map<odb::dbInst*, int> &whiteSpaceL);



    void init(std::set<odb::dbInst*> &insts);

    void init(std::set<odb::dbInst*> & insts, GraphModel gModel,
            EdgeWeightModel eModel);
    void printEdgeList();
    void setGraphModel(GraphModel graphModel);
    void setEdgeWeightModel(EdgeWeightModel edgeWeightModel);
    Vertex* dbToGraph(odb::dbInst* inst);

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
