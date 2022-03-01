#include "opendb/db.h"
#include "instGraph.h"
#include <set>
#include <cmath>
#include <fstream>
#include <queue>
#include <string>

using std::set;
using std::vector;
using std::cout;
using std::endl;
using std::make_pair;
using std::queue;

namespace ClipGraphExtract {

static float 
getEdgeWeight(int fanout, EdgeWeightModel eModel);

Vertex::Vertex() : 
  inst_(nullptr), 
  weight_(0),
  id_(0) {};

Vertex::Vertex(odb::dbInst* inst, 
    int id, float weight) 
  : Vertex() {
  inst_ = inst;
  id_ = id;
  weight_ = weight;
} 

void Vertex::setInst(odb::dbInst* inst) {
  inst_ = inst;
}

void Vertex::setWeight(float weight) {
  weight_ = weight;
}

void Vertex::addInEdge(Edge* edge) {
  inEdges_.push_back(edge);
}

void Vertex::addOutEdge(Edge* edge) { 
  outEdges_.push_back(edge);
}


Edge::Edge() : 
  from_(nullptr), to_(nullptr), 
  weight_(0) {}

Edge::Edge(Vertex* from, 
    Vertex* to,
    float weight)
  : from_(from), to_(to), weight_(weight) {}

void Edge::setFrom(Vertex* vertex) {
  from_ = vertex;
}

void Edge::setTo(Vertex* vertex) {
  to_ = vertex;
}

void Edge::setWeight(float weight) {
  weight_ = weight;
}

Graph::Graph() : db_(nullptr) {}
void Graph::setDb(odb::dbDatabase* db) {
  db_ = db;
}

Graph::~Graph() {
  db_ = nullptr;
  vector<Vertex>().swap(vertices_);
  vector<Edge>().swap(edges_);
  vertexMap_.clear();
}


// with given insts.
void Graph::init(std::set<odb::dbInst*> & insts, 
    GraphModel gModel, EdgeWeightModel eModel) {
  // extract iTerm
  set<odb::dbInst*> instSet;
  set<odb::dbITerm*> iTermSet;

  // init instSet, vertexMap_, and vertices_
  int vertexId = 0;
  vertices_.reserve(insts.size()); 
  for(auto& inst : insts) {
    instSet.insert(inst);
    vertices_.push_back(Vertex(inst, vertexId, 0.0f));

    // this is safe because of "reserve"
    vertexMap_.insert(std::make_pair(inst, 
          &vertices_[vertices_.size()-1]));

    for(odb::dbITerm* iTerm : inst->getITerms()) {
      if( iTerm->getSigType() == odb::dbSigType::POWER ||
          iTerm->getSigType() == odb::dbSigType::GROUND ) {
        continue;
      }
      iTermSet.insert(iTerm); 
    }
    vertexId++;
  } 
  // extract Net
  set<odb::dbNet*> netSet;
  for(odb::dbITerm* iTerm : iTermSet) {
    odb::dbNet* net = iTerm->getNet();
    if (!net) { 
      continue;
    }
    if( net->getSigType() == odb::dbSigType::POWER ||
        net->getSigType() == odb::dbSigType::GROUND ||
        net->getSigType() == odb::dbSigType::CLOCK ) {
      continue;
    }
    netSet.insert(net);
  }

  // STAR 
  if( gModel == GraphModel::Star ) {
    for(odb::dbNet* net : netSet) {
      // nets' source port is IO-port
      if( !net->getFirstOutput() ) {
        continue;
      }

      odb::dbInst* sourceInst = net->getFirstOutput()->getInst();
      Vertex* sourceVert = dbToGraph(sourceInst);
        
      // not exists in given insts pool, then escape
      auto instPtr = insts.find(sourceInst);
      if( instPtr == insts.end() ) {
        continue;
      }

      set<odb::dbInst*> sinkInstSet;
      for(odb::dbITerm* iTerm : net->getITerms()) {
        if( iTerm->getSigType() == odb::dbSigType::POWER ||
            iTerm->getSigType() == odb::dbSigType::GROUND ) {
          continue;
        }

        // not exists in given insts pool, then escape
        auto instPtr = insts.find(iTerm->getInst());
        if( instPtr == insts.end() ) {
          continue;
        }
        sinkInstSet.insert(iTerm->getInst());
      }
     
      // get fanout for star model 
      int fanout = 1;
      for(auto& sinkInst : sinkInstSet) {
        if( sinkInst == sourceInst) {
          continue;
        }
        fanout ++;
      }
      
      const float edgeWeight = getEdgeWeight(fanout, eModel);

      // for all sink instances
      for(auto& sinkInst : sinkInstSet) {
        if( sinkInst == sourceInst) {
          continue;
        }

        Vertex* sinkVert = dbToGraph(sinkInst);
        edges_.push_back(Edge(sourceVert, sinkVert, edgeWeight));
      }
    }
  }
  // CLIQUE
  else if( gModel == GraphModel::Clique ) {
    for(odb::dbNet* net : netSet) {
      set<odb::dbInst*> netInstSet;
      for(odb::dbITerm* iTerm : net->getITerms()) {
        if( iTerm->getSigType() == odb::dbSigType::POWER ||
            iTerm->getSigType() == odb::dbSigType::GROUND ) {
          continue;
        }
        
        // not exists in given insts pool, then escape
        auto instPtr = insts.find(iTerm->getInst());
        if( instPtr == insts.end() ) {
          continue;
        }

        netInstSet.insert(iTerm->getInst());
      }

      vector<odb::dbInst*> netInstStor(netInstSet.begin(), 
          netInstSet.end());
      int fanout = netInstSet.size();
      const float edgeWeight = getEdgeWeight(fanout, eModel);

      for(int i=0; i<netInstStor.size(); i++) {
        for(int j=i+1; j<netInstStor.size(); j++) {
          odb::dbInst* inst1 = netInstStor[i];
          odb::dbInst* inst2 = netInstStor[j];

          Vertex* vert1 = dbToGraph(inst1);
          Vertex* vert2 = dbToGraph(inst2);

          edges_.push_back(Edge(vert1, vert2, edgeWeight));
          edges_.push_back(Edge(vert2, vert1, edgeWeight));
        }
      }
    }
  }

  std::string netModelName = "";
  if( gModel == GraphModel::Clique ) {
    netModelName = "Clique";
  } 
  else if( gModel == GraphModel::Star ) {
    netModelName = "Star";
  }

  cout << "TotalVertices: " << vertices_.size() << endl;
  cout << "NetModel: " << netModelName << endl; 
  cout << "TotalEdges: " << edges_.size() << endl;
  // vertex' inEdge/outEdge update
  updateVertsFromEdges();
}

// Need for BFS  search
void Graph::updateVertsFromEdges() {
  for(auto& edge : edges_) {
    Vertex* fromVert = edge.from();
    Vertex* toVert = edge.to();
    
    if( !fromVert || !toVert ) {
      cout << "ERROR: Vertex not existed!!" << endl;
      exit(1);
    }

    fromVert->addOutEdge(&edge);
    toVert->addInEdge(&edge);
  }
}

void Graph::printEdgeList() {
  cout << "edges: " << edges_.size() << endl;
  for(auto& edge: edges_) {
    cout << edge.from()->inst()->getConstName() << " ";
    cout << edge.to()->inst()->getConstName() << " " ;
    cout << edge.weight() << endl;
  }
}

void Graph::saveFile(std::string fileName) {
  std::ofstream outFile;
  outFile.open(fileName, std::ios_base::out); // overwrite
	for(auto& edge: edges_) {
		outFile << edge.from()->inst()->getConstName() << " ";
		outFile << edge.to()->inst()->getConstName() << " ";
		outFile << edge.weight() << endl;
	}
	outFile.close();
}

Vertex* Graph::dbToGraph(odb::dbInst* inst) {
  auto vertPtr = vertexMap_.find(inst);
  if( vertPtr == vertexMap_.end() ) {
    return nullptr; 
  }
  else {
    return vertPtr->second;
  }
}

static float 
getEdgeWeight(int fanout, EdgeWeightModel eModel) {
  switch( eModel ) {
    case A:
      return 4.0/(fanout*(fanout-1));
      break;
    case B:
      return 2.0/fanout;
      break;
    case C:
      return 8.0/(fanout*fanout*fanout);
      break;
    case D:
      return 2.0/(std::pow(fanout,1.5));
      break;
    case E:
      return 1.0/(fanout-1);
      break;
    // default is the same as case A
    default:
      return 4.0/(fanout*(fanout-1));
      break;
  }
}

}
