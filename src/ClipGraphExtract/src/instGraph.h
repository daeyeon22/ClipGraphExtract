#ifndef __INST_GRAPH__
#define __INST_GRAPH__

#include <clip_graph_ext/clipGraphExtractor.h>
#include <map>
#include <set>

namespace odb {
  class dbInst;
  class dbDatabase;
}

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

private:
  odb::dbInst* inst_;
  std::vector<Edge*> inEdges_;
  std::vector<Edge*> outEdges_;
  int id_;
  float weight_;
};

inline odb::dbInst* Vertex::inst() const {
  return inst_;
}

inline const std::vector<Edge*> & Vertex::inEdges() const {
  return inEdges_;
}

inline const std::vector<Edge*> & Vertex::outEdges() const {
  return outEdges_;
}

inline float Vertex::weight() const {
  return weight_;
}

inline void Vertex::setId(int id) {
  id_ = id;
}

inline int Vertex::id() const {
  return id_;
}

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
  void setDb(odb::dbDatabase* db);

  void init(std::set<odb::dbInst*> & insts, GraphModel gModel,
      EdgeWeightModel eModel);
  void printEdgeList();
  
  Vertex* dbToGraph(odb::dbInst* inst);

private:
  odb::dbDatabase* db_;
  std::vector<Vertex> vertices_;
  std::vector<Edge> edges_;
  std::map<odb::dbInst*, Vertex*> vertexMap_;
  void updateVertsFromEdges();
};

}

#endif
