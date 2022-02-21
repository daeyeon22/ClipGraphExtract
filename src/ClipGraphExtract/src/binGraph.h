#ifndef __BIN_GRAPH__
#define __BIN_GRAPH__

#include <unordered_map>
#include <set>
#include <vector>
#include <unordered_map>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// for easier coding with boost
typedef bg::model::point<int, 2, bg::cs::cartesian> point;
typedef bg::model::box<point> box;

struct inst_value{
    odb::dbInst* inst_;
    box box_;
};

struct wire_value{
    unsigned int layerNum_;
    char type_;
    odb::dbNet* net_;
    int from_x_, from_y_; // centric point
    int to_x_, to_y_; // centric point
    
    int lx_, ly_; // minimum boundary point
    int ux_, uy_; // maximum boundary point

    int width_;
    int spacing_;
    
    std::unordered_map<unsigned int, int> pitches_;

    box box_;

    void setBox(int fx, int fy, int tx, int ty);
};

struct via_value{
    unsigned int bottomLayer_;
    odb::dbNet* net_;
    int x_, y_; // centric point

    int lx_, ly_; // minimum boundary point
    int ux_, uy_; // maximum boundary point

    int width_;
    box box_;
    
    void setBox(int x, int y);
};

struct pin_value{
    std::string name_;
    int x_, y_; // centric point

    box box_;
    void setBox(int x, int y);
};

struct drc_value{
    std::string type_;

	std::string from_;
	std::string to_;

    int lx_, ly_; // minimum boundary point
    int ux_, uy_; // maximum boundary point

    box box_;

    void setBox(int lx, int ly, int ux, int uy);
};


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
    Vertex(int id, int lx, int ly, int ux, int uy, int maxLayer,  
            std::vector<odb::dbInst*> insts, 
            std::vector<wire_value> wireValues, 
            std::vector<via_value> viaValues,
            std::vector<pin_value> pinValues);

    std::vector<odb::dbInst*> getInsts();
    std::vector<Edge*>  getInEdges();
    std::vector<Edge*>  getOutEdges();

    void addInst(odb::dbInst* inst);
    std::vector<wire_value> getWireValues();
    std::vector<via_value> getViaValues();
    std::vector<pin_value> getPinValues();
    std::vector<drc_value> getDrcValues();

    void addInst(odb::dbInst* inst);
    void addWireValue(wire_value wireValue);
    void addViaValue(via_value viaValue);
    void addPinValue(pin_value pinValue);
    void addDrcValue(drc_value drcValue);
    void addInEdge(Edge* edge);
    void addOutEdge(Edge* edge);

    void setLabel(int label);
    int getLabel();
    
    // for node feature (.x)
    double getUtilization() const;
    double getRoutingCongestion(char) const;
    std::unordered_map<unsigned int, std::pair<char, unsigned int> > getEachWireLength() const;
    unsigned int getTotalWireLength(std::unordered_map<unsigned int, std::pair<char, unsigned int> >) const;
    unsigned int getWireLength(std::unordered_map<unsigned int, std::pair<char, unsigned int> >, char type) const;
    std::unordered_map<unsigned int, unsigned int> getTotalWireCapacities() const;
    double getWireCongestion(char type) const;
    double getViaUtilization() const;
    unsigned int getNumOfPin() const;
    unsigned int getNumOfDrc() const;
    double getStdOfPins(char type) const;
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
    std::vector<wire_value> wireValues_;
    std::vector<via_value> viaValues_;
    std::vector<pin_value> pinValues_;
    std::vector<drc_value> drcValues_;
    std::vector<Edge*> inEdges_;
    std::vector<Edge*> outEdges_;
    int lx_, ly_, ux_, uy_;
    int maxLayer_;
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
    void addVertex(int lx, int ly, int ux, int uy, int maxLayer, 
            std::vector<odb::dbInst*> insts, 
            std::vector<wire_value> wireValues, 
            std::vector<via_value> viaValues,
            std::vector<pin_value> pinValues);

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
