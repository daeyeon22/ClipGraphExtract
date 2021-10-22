




#ifndef __BIN_GRAPH__
#define __BIN_GRAPH__

#include <unordered_map>
#include <set>
#include <vector>



namespace odb {
  class dbInst;
  class dbDatabase;
}

// only consider two-way connections
namespace bingraph {

typedef std::pair<int,int> edge_key_type;

class Edge;
// vertex is equal to dbInst
class Vertex {
  public:
    Vertex(int id, int lx, int ly, int ux, int uy, std::vector<odb::dbInst*> insts);

    std::vector<odb::dbInst*> getInsts();
    std::vector<Edge*>  getInEdges();
    std::vector<Edge*>  getOutEdges();

    void addInst(odb::dbInst* inst);
    void addInEdge(Edge* edge);
    void addOutEdge(Edge* edge);

    void setLabel(int label);
    


    int getLabel();
    
    // for node feature (.x)
    double getUtilization() const;
    double getAvgInEdges() const;
    double getAvgOutEdges() const;
    double getSequentialRatio() const;


    int getId() const;
    int getLx() const;
    int getLy() const;
    int getUx() const;
    int getUy() const;
    

private:
    std::vector<odb::dbInst*> insts_;
    std::vector<Edge*> inEdges_;
    std::vector<Edge*> outEdges_;
    int lx_, ly_, ux_, uy_;
    int id_;
    int label_;

};

// edge is inst1-inst2 connections
class Edge {
  public:
    Edge();
    ~Edge();
    Edge(Vertex* from, Vertex* to, float weight);

    Vertex* getFrom() const;
    Vertex* getTo() const;
    float getWeight() const;

    void setFrom(Vertex* inst);
    void setTo(Vertex* inst);
    void setWeight(float weight);
    void incrWeight(float weight);
  private:
    Vertex* from_;
    Vertex* to_;
    float weight_;
};

struct pair_hash{
    template<class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& p) const{
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

class Graph {
  public:
    Graph();
    ~Graph();

    std::vector<Vertex*> getVertices();
    std::vector<Edge*> getEdges();

    void addVertex(int lx, int ly, int ux, int uy, std::vector<odb::dbInst*> insts);
    void initEdges();
    void saveFile(const char* prefix);
    
  
  private:
    odb::dbDatabase* db_;
    std::vector<Vertex> vertices_;
    std::vector<Edge> edges_;
    std::unordered_map<odb::dbInst*, Vertex*> inst2vertex_;
    std::unordered_map<std::pair<int,int>, Edge*, pair_hash> pair2edge_;

};

}

#endif
