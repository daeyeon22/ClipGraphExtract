#include "opendb/db.h"
#include "sta/Sta.hh"
#include "sta/Network.hh"
#include "db_sta/dbSta.hh"
#include "db_sta/dbNetwork.hh"
#
#include "instGraph.h"
#include <set>
#include <cmath>
#include <fstream>
#include <queue>
#include <string>
#include <unordered_map>


//using std::set;
//using std::vector;
//using std::cout;
//using std::endl;
//using std::make_pair;
//using std::queue;

namespace ClipGraphExtract {

using namespace std;
using namespace odb;

static float 
getEdgeWeight(int fanout, EdgeWeightModel eModel);




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


void Graph::setSta(sta::dbSta* sta) {
    sta_ = sta;
}


void Graph::setGraphModel(GraphModel graphModel) {
    graphModel_ = graphModel;
}

void Graph::setEdgeWeightModel(EdgeWeightModel edgeWeightModel) {
    edgeWeightModel_ = edgeWeightModel;
}

void Graph::init(std::set<odb::dbInst*> &insts) {
    init(insts, graphModel_, edgeWeightModel_);
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

  //cout << "TotalVertices: " << vertices_.size() << endl;
  //cout << "NetModel: " << netModelName << endl; 
  //cout << "TotalEdges: " << edges_.size() << endl;
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



void Graph::saveNodeFeaFile(std::string fileName) {

    
    int dbu = db_->getChip()->getBlock()->getDbUnitsPerMicron();

    std::ofstream outFile;
    outFile.open(fileName, std::ios_base::out);
    for(auto &tarVertex : vertices_) {
        outFile << tarVertex.id() << " " 
                << tarVertex.getRelPosX() << " " 
                << tarVertex.getRelPosY() << " " 
                << tarVertex.getAbsSlack() << " "
                << tarVertex.getRelSlack() << " "
                << tarVertex.getNumAccPoints() << " "
                << tarVertex.getNumBlkPoints() << " " 
                << tarVertex.getNumBndPoints() << " " 
                << tarVertex.getWhiteSpaceL() << " "
                << tarVertex.getWhiteSpaceR() << " "
                << tarVertex.getWhiteSpaceT() << " "
                << tarVertex.getWhiteSpaceD() << " "
                << 1.0* tarVertex.getSize() / (dbu*dbu) << " "
                << tarVertex.getDegree() << " "
                << tarVertex.getNumInEdges() << " "
                << tarVertex.getNumOutEdges() << " "
                << tarVertex.isClocked() << " "
                << tarVertex.getSWireOverlap() << " " 
                << tarVertex.getBBoxSize() << " " 
                << tarVertex.getCellType() << " "
                << tarVertex.getCutEdges() << " " 
                << tarVertex.getIsCrit() << endl;
    }
    outFile.close();
}

void Graph::saveEdgeIdxFile(std::string fileName) {

    std::ofstream outFile;
    outFile.open(fileName, std::ios_base::out);
	for(auto& edge: edges_) {
		outFile << edge.from()->id() << " " 
                << edge.to()->id() << endl;
	}
	outFile.close();

}

void Graph::saveEdgeAttFile(std::string fileName) {
    std::ofstream outFile;
    outFile.open(fileName, std::ios_base::out);
	for(auto& edge: edges_) {
		outFile << edge.weight() << endl;
	}
	outFile.close();
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

void Graph::setNumDrvs( unordered_map<dbInst*, int> &numDrvs) {
    for(auto& val : vertexMap_) {
        Vertex* tarVertex = val.second;
        dbInst* tarInst = val.first;
        tarVertex->setNumDrvs(numDrvs[tarInst]);
    }
}



void Graph::setRelPos(  unordered_map<dbInst*, double> &relPosX, 
                        unordered_map<dbInst*, double> &relPosY ) {
    for(auto& val : vertexMap_) {
        Vertex* tarVertex = val.second;
        dbInst* tarInst = val.first;
        tarVertex->setRelPos(relPosX[tarInst], relPosY[tarInst]);
    }
}



void Graph::setSlack(double clockPeriod, std::unordered_map<dbInst*, double> &absoluteSlack) {
    for(auto& val : vertexMap_) {
        Vertex* tarVertex = val.second;
        dbInst* tarInst = val.first;
        double slack = absoluteSlack[tarInst];

        if(tarVertex == NULL) {
            cout << tarInst->getName() << " ??? nullptr" << endl;
            exit(0);
        }
        tarVertex->setSlack(clockPeriod, slack);
    }
}

void Graph::setIsCrit(std::unordered_map<dbInst*, bool> &isCrit) {
    for(auto& val : vertexMap_) {
        Vertex* tarVert = val.second;
        dbInst* tarInst = val.first;
        tarVert->setIsCrit(isCrit[tarInst]);
    }
}


void Graph::setSWireOverlap(std::unordered_map<dbInst*, double> &sWireOverlap) {
    for(auto& val : vertexMap_) {
        Vertex* tarVert = val.second;
        dbInst* tarInst = val.first;
        tarVert->setSWireOverlap(sWireOverlap[tarInst]);
    }
}


void Graph::setCutEdges(std::unordered_map<dbInst*, int> &cutEdges) {
    for(auto& val : vertexMap_) {
        Vertex* tarVert = val.second;
        dbInst* tarInst = val.first;
        tarVert->setCutEdges(cutEdges[tarInst]);
    }
}


void Graph::setCellType(std::unordered_map<dbInst*, int> &cellType) {
    for(auto& val : vertexMap_) {
        Vertex* tarVert = val.second;
        dbInst* tarInst = val.first;
        tarVert->setCellType(cellType[tarInst]);
    }
}

void Graph::setBBoxSize(std::unordered_map<dbInst*, double> &bboxSize) {
     for(auto& val : vertexMap_) {
        Vertex* tarVert = val.second;
        dbInst* tarInst = val.first;
        tarVert->setBBoxSize(bboxSize[tarInst]);
    }   
}


void Graph::setNumPoints(   std::unordered_map<dbInst*, int> &numAccPoints,
                            std::unordered_map<dbInst*, int> &numBlkPoints,
                            std::unordered_map<dbInst*, int> &numBndPoints) {

    // set properties
    for(auto& val : vertexMap_) {
        Vertex* tarVertex = val.second;
        dbInst* tarInst = val.first;
        tarVertex->setNumAccPoints(numAccPoints[tarInst]);
        tarVertex->setNumBlkPoints(numBlkPoints[tarInst]);
        tarVertex->setNumBndPoints(numBndPoints[tarInst]);

    
        

    }
}


void Graph::setViaFeature( std::unordered_map<dbInst*, double> &powerViaDistance) {
    for(auto& val : vertexMap_) {
        Vertex* tarVertex = val.second;
        dbInst* tarInst = val.first;
        tarVertex->setViaFeature(powerViaDistance[tarInst]);
    }
}
void Graph::setWhiteSpace(  std::unordered_map<dbInst*,double> &whiteSpaceL,
                            std::unordered_map<dbInst*,double> &whiteSpaceR,
                            std::unordered_map<dbInst*,double> &whiteSpaceT,
                            std::unordered_map<dbInst*,double> &whiteSpaceD ) {
    for(auto& val : vertexMap_) {
        Vertex* tarVertex = val.second;
        dbInst* tarInst = val.first;
        tarVertex->setWhiteSpaceL(whiteSpaceL[tarInst]);
        tarVertex->setWhiteSpaceR(whiteSpaceR[tarInst]);
        tarVertex->setWhiteSpaceT(whiteSpaceT[tarInst]);
        tarVertex->setWhiteSpaceD(whiteSpaceD[tarInst]);
    }
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


void Graph::print() {

    double dbu2 = db_->getChip()->getBlock()->getDbUnitsPerMicron();
    dbu2 = dbu2*dbu2;

    for(auto elem : vertexMap_) {
        dbInst* tarInst = elem.first;
        Vertex* tarVert = elem.second;

        cout << tarVert->id() << " " << tarInst->getName() << endl;
        cout << "   - relpos (x,y) : " << tarVert->getRelPosX() << " " << tarVert->getRelPosY() << endl;
        cout << "   - # acc points : " << tarVert->getNumAccPoints() << endl;
        cout << "   - # blk points : " << tarVert->getNumBlkPoints() << endl;
        cout << "   - # bnd points : " << tarVert->getNumBndPoints() << endl;
        cout << "   - # cut edges  : " << tarVert->getCutEdges() << endl;
        cout << "   - wire overlap : " << tarVert->getSWireOverlap() << endl;
        cout << "   - wspace (L/R) : " << tarVert->getWhiteSpaceL() << " " 
                                       << tarVert->getWhiteSpaceR() << endl;
        cout << "   - wspace (T/D) : " << tarVert->getWhiteSpaceT() << " " 
                                       << tarVert->getWhiteSpaceD() << endl;
        cout << "   - encoded type : " << tarVert->getCellType() << endl;
        cout << "   - size of bbox : " << tarVert->getBBoxSize() << endl;
        cout << "   - # outedges   : " << tarVert->getNumOutEdges() << endl;
        cout << "   - # inedges    : " << tarVert->getNumInEdges() << endl;
        cout << "   - degree       : " << tarVert->getDegree() << endl;
        cout << "   - cell size    : " << 1.0*tarVert->getSize() / dbu2 << endl;
        cout << "   - abs slack    : " << tarVert->getAbsSlack() << endl;
        cout << "   - rel slack    : " << tarVert->getRelSlack() << endl;
        cout << "   - is critical  : " << tarVert->getIsCrit() << endl;
        cout << endl; 
    }

}


}
