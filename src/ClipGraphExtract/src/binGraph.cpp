
#include "opendb/db.h"
#include "binGraph.h"
#include <fstream>
#include <iostream>

namespace bingraph {

using namespace odb;
using namespace std;

Vertex::Vertex(int id, int lx, int ly, int ux, int uy, vector<dbInst*> insts) :
    id_(id), lx_(lx), ly_(ly), ux_(ux), uy_(uy), insts_(insts) {
    label_ = 0;
}


vector<dbInst*> Vertex::getInsts() {
    return insts_;
}

vector<Edge*> Vertex::getInEdges() {
    return inEdges_;
}

vector<Edge*> Vertex::getOutEdges() {
    return outEdges_;
}


void
Vertex::addInst(dbInst* inst) {
    insts_.push_back(inst);
}


void
Vertex::addInEdge(Edge* edge) {
    inEdges_.push_back(edge);
}

void
Vertex::addOutEdge(Edge* edge) {
    outEdges_.push_back(edge);
}

void 
Vertex::setLabel(int label) {
    label_ = label;
}

int 
Vertex::getLabel() {
    return label_;
}



int Vertex::getId() const { return id_; }
int Vertex::getLx() const { return lx_; }
int Vertex::getLy() const { return ly_; }
int Vertex::getUx() const { return ux_; }
int Vertex::getUy() const { return uy_; }

double Vertex::getUtilization() const {
    long int overlaps = 0;
    for(dbInst* inst : insts_) {
        int overlapLx = max( lx_, inst->getBBox()->xMin() );
        int overlapLy = max( ly_, inst->getBBox()->yMin() );
        int overlapUx = min( ux_, inst->getBBox()->xMax() );
        int overlapUy = min( uy_, inst->getBBox()->yMax() );
        overlaps += (overlapUx - overlapLx) * (overlapUy - overlapLy);
    }

    long int totalArea = (ux_ - lx_) * (uy_ - ly_);

    return 1.0* overlaps / totalArea;
}

double Vertex::getAvgInEdges() const {
    if (insts_.size() == 0)
        return 0;
    else
        return 1.0 * inEdges_.size() / insts_.size();
}

double Vertex::getAvgOutEdges() const {
    if (insts_.size() == 0)
        return 0;
    else
        return 1.0 * outEdges_.size() / insts_.size();
}

double Vertex::getSequentialRatio() const {
    if (insts_.size() == 0)
        return 0;
    else {
        int numSequential = 0;
        for(dbInst* inst : insts_) {
            if(inst->getMaster()->isSequential())
                numSequential++;
        }
        return 1.0 * numSequential / insts_.size();
    }
}

Edge::Edge() {

}

Edge::Edge(Vertex* from, Vertex* to, float weight) {
    from_ = from;
    to_ = to;
    weight_ = weight;
}

Edge::~Edge() {

}

void
Edge::incrWeight(float weight) {
    weight_ += weight;
}


set<dbInst*> getSinks(dbInst* inst) {
    dbSet<dbITerm> inst_iterms = inst->getITerms();
    set<dbInst*> sinks;
    for(auto inst_iterm_itr=inst_iterms.begin(); inst_iterm_itr != inst_iterms.end(); 
            ++inst_iterm_itr) {
        dbITerm* inst_iterm = *inst_iterm_itr;
        if(inst_iterm->getNet() == NULL)
            continue;
        if(inst_iterm->getIoType() != dbIoType::OUTPUT) 
            continue;
        dbNet* net = inst_iterm->getNet();
        dbSet<dbITerm> net_iterms = net->getITerms();
        for(auto net_iterm_itr = net_iterms.begin(); net_iterm_itr != net_iterms.end(); 
                ++net_iterm_itr) {
            dbITerm* net_iterm = *net_iterm_itr;
            dbInst* sink_inst = net_iterm->getInst();
            sinks.insert(sink_inst);
        }
    }
    return sinks; 
}

Graph::Graph() {

}

Graph::~Graph() {

}

vector<Vertex*> Graph::getVertices() {
    vector<Vertex*> vertices;
    for(int i=0; i < vertices_.size(); i++)
        vertices.push_back(&vertices_[i]);
    return vertices;
}


void 
Graph::addVertex( int lx, int ly, int ux, int uy, vector<dbInst*> insts ) {
    //cout << "start addvertex" << endl;

    Vertex vert(vertices_.size(), lx, ly, ux, uy, insts);
    //for( dbInst* inst : insts ) {
    //    vert.addInst(inst);
    //}

    vertices_.push_back(vert);
   // vertices_.push_back( Vertex(vertices_.size(), lx, ly, ux, uy) );
    
    //Vertex* vertex = &vertices_.back();
    //for( dbInst* inst : insts ) {
    //    vertex->addInst(inst);
    //    inst2vertex_[inst] = vertex;
    //}

    //cout << "end addvertex" << endl;
}


void
Graph::initEdges() {
    cout << "start init edges" << endl;
    unordered_map<edge_key_type, Edge, pair_hash> pair2edge;

    for(int i=0; i < vertices_.size(); i++ ) {
        Vertex* vertex = &vertices_[i];

        for(dbInst* inst : vertex->getInsts()) {
            inst2vertex_[inst] = vertex;
        }
    }

    for(int i=0; i < vertices_.size(); i++) {
        Vertex* from = &vertices_[i];
        for(dbInst* i1 : from->getInsts()) {
            for(dbInst* i2 : getSinks(i1)) {
                Vertex* to = inst2vertex_[i2];
                int fromId = from->getId();
                int toId = to->getId();
                edge_key_type key(fromId, toId);

                if( pair2edge.find(key) == pair2edge.end() ) {
                    pair2edge[key] = Edge(from,to,1.0);
                    //edges_.push_back( Edge(from, to, 1.0) );
                    //Edge* edge = &edges_.back();
                    //pair2edge_[key] = edge;
                    //from->addOutEdge( edge );
                    //to->addInEdge( edge );
                } else {
                    pair2edge[key].incrWeight(1.0);

                    //Edge* edge = pair2edge_[key];
                    //edge->incrWeight(1.0);
                }
            }
        }
    }

    // copy
    for(auto it : pair2edge) {
        edges_.push_back(it.second);
    }

    // key mapping
    for(int i=0; i < edges_.size(); i++) {
        Edge* edge = &edges_[i];
        Vertex* from = edge->getFrom();
        Vertex* to = edge->getTo();
        edge_key_type key(from->getId(), to->getId());
        pair2edge_[key] = edge;
        from->addOutEdge( edge );
        to->addInEdge( edge );
    }


    cout << "end init edges" << endl;
}



Vertex* Edge::getFrom() const {
    return from_;
}

Vertex* Edge::getTo() const {
    return to_;
}

float Edge::getWeight() const {
    return weight_;
}

void
Graph::saveFile(const char* prefix) {
    ofstream nodeAttr;
    ofstream nodeLabel;
    ofstream edgeIndex;
    string attrFileName = string(prefix) + ".x";
    string labelFileName = string(prefix) + ".y";
    string edgeIndexFileName = string(prefix) + ".edge_index"; 

    nodeAttr.open(attrFileName, std::ios_base::out);
    nodeLabel.open(labelFileName, std::ios_base::out);
    for(auto& vert : vertices_) {
        nodeAttr << vert.getId() << " " 
                 << vert.getUtilization() << " " 
                 << vert.getAvgInEdges() << " " 
                 << vert.getAvgOutEdges() << " " 
                 << vert.getSequentialRatio() << endl;
    
        nodeLabel   << vert.getId() << " " 
                    << vert.getLabel() << endl; 
    }
    nodeAttr.close();
    nodeLabel.close();
    edgeIndex.open(edgeIndexFileName, std::ios_base::out);

    for(auto& edge : edges_) {
        Vertex* from = edge.getFrom();
        Vertex* to = edge.getTo();
        float weight = edge.getWeight();
        for(int i=0; i < weight; i++) {
            edgeIndex << from->getId() << " " << to->getId() << endl;
        }
    }
    edgeIndex.close();
}

};

