#include "math.h"
#include "opendb/db.h"
#include "binGraph.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <unordered_map>

#include <iostream>

void wire_value::setBox(int fx, int fy, int tx, int ty){
	from_x_ = fx;
    from_y_ = fy;
    to_x_ = tx;
    to_y_ = ty;

    lx_ = fx - width_/2;
    ly_ = fy - width_/2;
    ux_ = tx + width_/2;
    uy_ = ty + width_/2;

    box b(point(lx_, ly_), 
           point(ux_, uy_));
    box_ = b;
}

void via_value::setBox(int x, int y){
    x_ = x;
    y_ = y;

    lx_ = x - width_/2;
    ly_ = y - width_/2;
    ux_ = x + width_/2;
    uy_ = y + width_/2;
    
    box b(point(lx_, ly_), 
            point(ux_, uy_));
    box_ = b;
}

void pin_value::setBox(int x, int y){
    x_ = x;
    y_ = y;

    int lx_ = x - 1;
    int ly_ = y - 1;
    int ux_ = x + 1;
    int uy_ = y + 1;
    
    box b(point(lx_, ly_), 
            point(ux_, uy_));
    box_ = b;
}

void drc_value::setBox(int lx, int ly, int ux, int uy){
    lx_ = lx;
    ly_ = ly;
    ux_ = ux;
    uy_ = uy;

    box b(point(lx_, ly_), 
           point(ux_, uy_));
    box_ = b;
}

namespace bingraph {

using namespace odb;
using namespace std;

Vertex::Vertex(int id, int lx, int ly, int ux, int uy, int maxLayer, 
        vector<dbInst*> insts, vector<wire_value> wireValues,
        vector<via_value> viaValues, vector<pin_value> pinValues) :
    id_(id), lx_(lx), ly_(ly), ux_(ux), uy_(uy), maxLayer_(maxLayer), insts_(insts), 
    wireValues_(wireValues), viaValues_(viaValues), pinValues_(pinValues) {
    label_ = 0;
}


vector<dbInst*> Vertex::getInsts() {
    return insts_;
}

vector<wire_value> Vertex::getWireValues() {
    return wireValues_;
}

vector<via_value> Vertex::getViaValues() {
    return viaValues_;
}

vector<pin_value> Vertex::getPinValues() {
    return pinValues_;
}

vector<drc_value> Vertex::getDrcValues() {
    return drcValues_;
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
Vertex::addWireValue(wire_value wireValue) {
    wireValues_.push_back(wireValue);
}

void
Vertex::addViaValue(via_value viaValue) {
    viaValues_.push_back(viaValue);
}

void
Vertex::addPinValue(pin_value pinValue) {
    pinValues_.push_back(pinValue);
}

void
Vertex::addDrcValue(drc_value drcValue) {
    drcValues_.push_back(drcValue);
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

void Vertex::getNets() {
	for(wire_value wireValue : wireValues_){
		nets_.insert(wireValue.net_);

		if((wireValue.lx_ < lx_) || (wireValue.ly_ < ly_) || (wireValue.ux_ > ux_) || (wireValue.uy_ > uy_))
			globalNets_.insert(wireValue.net_);
	}

	set_difference(nets_.begin(), nets_.end(), globalNets_.begin(), globalNets_.end(), inserter(localNets_, localNets_.end()));
	
	cout << "total" << endl;
	for(auto net : nets_) cout << net->getName() << " ";
	cout << endl;

	cout << "global" << endl;
	for(auto net : globalNets_) cout << net->getName() << " ";
	cout << endl;

	cout << "local" << endl;
	for(auto net : localNets_) cout << net->getName() << " ";
	cout << endl;

	
}


double Vertex::getRoutingCongestion(char type) const {
	
	unsigned int boxLength = ux_-lx_; // Box is a square.

	unordered_map<char, double> layerDemands;	
	unordered_map<char, double> layerCaps;
	unordered_map<unsigned int, bool> tf;
	
	char type_;

	// Wire distinction
	for(wire_value wireValue : wireValues_){
		double value;
		
		if(wireValue.type_ == 'H'){
			if(!tf[wireValue.layerNum_] && wireValue.layerNum_ <= maxLayer_){
				layerCaps['L'] += boxLength;
				layerCaps['R'] += boxLength;
				tf[wireValue.layerNum_] = true;
			}

			// Select type
			if(wireValue.lx_ < lx_) type_ = 'L'; // Left
			else if(ux_ < wireValue.ux_) type_ = 'R'; // Right

			if(uy_ < wireValue.uy_)
				value = uy_-wireValue.uy_+wireValue.width_;
			
			else if(wireValue.ly_ < ly_)
				value = wireValue.ly_-ly_+wireValue.width_;

			else value = wireValue.width_;
		}
		// Up, down
		else if(wireValue.type_ == 'V'){
			if(!tf[wireValue.layerNum_] && wireValue.layerNum_ <= maxLayer_){
				layerCaps['D'] += boxLength;
				layerCaps['U'] += boxLength;
				tf[wireValue.layerNum_] = true;
			}

			// Select type
			if(wireValue.ly_ < ly_) type_ = 'D'; // Down(Lower)
			else if(uy_ < wireValue.uy_) type_ = 'U'; // Up(Upper)
		
			if(ux_ < wireValue.ux_)
				value = ux_-wireValue.ux_+wireValue.width_;
			
			else if(wireValue.lx_ < lx_)
				value = wireValue.lx_-lx_+wireValue.width_;

			else value = wireValue.width_;		
		}

		else continue; // Inside
		
		value /= wireValue.width_;
		value *= wireValue.width_+wireValue.spacing_;

		layerDemands[type_] += value;
	}

		
	unordered_map<char, double> congestion;
	if(layerCaps['L']) congestion['L'] = layerDemands['L']/layerCaps['L'];
	else congestion['L'] = 0;
	
	if(layerCaps['R']) congestion['R'] = layerDemands['R']/layerCaps['R'];
	else congestion['R'] = 0;

	if(layerCaps['D']) congestion['D'] = layerDemands['D']/layerCaps['D'];
	else congestion['D'] = 0;
	
	if(layerCaps['U']) congestion['U'] = layerDemands['U']/layerCaps['U'];
	else congestion['U'] = 0;

	congestion['T'] = (congestion['L']+congestion['R']+congestion['D']+congestion['U']) / 4;

	return congestion[type];		
}

unordered_map<unsigned int, pair<char, unsigned int> > Vertex::getEachWireLength() const { 
    unordered_map<unsigned int, pair<char, unsigned int> > wireDemands; // Track * wire length = Actual wire length

   
    for(wire_value wireValue : wireValues_){
        unsigned int overlapLx = max( lx_, wireValue.lx_ );
        unsigned int overlapLy = max( ly_, wireValue.ly_ );
        unsigned int overlapUx = min( ux_, wireValue.ux_ );
        unsigned int overlapUy = min( uy_, wireValue.uy_ );

		wireDemands[wireValue.layerNum_].first = wireValue.type_;
        
		if(wireValue.type_ == 'H')
            wireDemands[wireValue.layerNum_].second += (overlapUx - overlapLx);
		else if(wireValue.type_ == 'V')
            wireDemands[wireValue.layerNum_].second += (overlapUy - overlapLy);

    }

    return wireDemands;
}

unsigned int Vertex::getTotalWireLength(unordered_map<unsigned int, pair<char, unsigned int> > each) const {
    unsigned int sum = 0;
    for(auto e : each)
        sum += e.second.second;

    return sum;
}

unsigned int Vertex::getWireLength(unordered_map<unsigned int, pair<char, unsigned int> > each, char type) const {
    unsigned int sum = 0;
    
	for(auto e : each){
        if(e.second.first == type || type == 'T')
            sum += e.second.second;
    }
    return sum;
}

double Vertex::getWireCongestion(char type) const {
    if(wireValues_.size() == 0) return 0;

    unordered_map<unsigned int, double> wireCapacities; // Track * box length = Total wire length
    unordered_map<unsigned int, pair<char, unsigned int> > wireDemands; // Track * wire length = Actual wire length
  
    for(int n = 2; n <= 8; n++)
		wireDemands[n].second = 0;


    for(wire_value wireValue : wireValues_){
        unsigned int overlapLx = max( lx_, wireValue.lx_ );
        unsigned int overlapLy = max( ly_, wireValue.ly_ );
        unsigned int overlapUx = min( ux_, wireValue.ux_ );
        unsigned int overlapUy = min( uy_, wireValue.uy_ );
		
		wireDemands[wireValue.layerNum_].first = wireValue.type_;

        if(wireValue.type_ == 'H')
            wireDemands[wireValue.layerNum_].second += (overlapUx - overlapLx);
        else if(wireValue.type_ == 'V')
            wireDemands[wireValue.layerNum_].second += (overlapUy - overlapLy);

    }
   
    bool cap = false;
    for(wire_value wireValue : wireValues_){
        if(!cap){
            for(int n = 2; n <= 8; n++){
                if(wireValue.pitches_[n] == 0)
                    continue;
                wireCapacities[n] = (ux_ - lx_) * ((ux_ - lx_) / (double)(wireValue.pitches_[n]));
				// length * track
                cap = true;
            }
        }
    }
    
    unordered_map<unsigned int, pair<char, double> > wireUtils;
    for(int n = 2; n <= 8; n++){
		wireUtils[n].first = wireDemands[n].first; // write type => 'H', 'V'
        if(wireDemands[n].second == 0 || wireCapacities[n] == 0) wireUtils[n].second = 0;
        else wireUtils[n].second = (double) wireDemands[n].second / (double) wireCapacities[n];
    }
    double avgWireUtil = 0;
    double sum = 0;
    int cnt = 0;

    for(pair<unsigned int, pair<char, double> > wireUtil : wireUtils){
        // if(maxLayer_ < wireUtil.first) continue;
        if(type == 'T' || wireUtil.second.first == type) {
			sum += wireUtil.second.second;
        	cnt++;
		}
    } 
    
    // avgWireUtil = sum / cnt;
    avgWireUtil = sum / maxLayer_;
    return avgWireUtil;
}

double Vertex::getViaUtilization() const {
    long int overlaps = 0;
    for(via_value viaValue : viaValues_) {
        int overlapLx = max( lx_, viaValue.lx_ );
        int overlapLy = max( ly_, viaValue.ly_ );
        int overlapUx = min( ux_, viaValue.ux_ );
        int overlapUy = min( uy_, viaValue.uy_ );
        overlaps += (overlapUx - overlapLx) * (overlapUy - overlapLy);
    }

    long int totalArea = (ux_ - lx_) * (uy_ - ly_);

    return 1.0* overlaps / totalArea;
}

unsigned int Vertex::getNumOfPin() const {
    return pinValues_.size();
}

double Vertex::getStdOfPins(char type) const {
    if(getNumOfPin() == 0) return 0;

    double variance = 0;
    double variance_x = 0;
    double variance_y = 0;
    double ave_x = 0;
    double ave_y = 0;


    if(type == 'H' || type == 'T'){
        int sum = 0;
        int cnt = 0;
    
        for(pin_value pinValue : pinValues_){
            sum += pinValue.x_;
            cnt++;
        }
        ave_x = sum/(double)cnt;
        for(pin_value pinValue : pinValues_)
            variance_x += pow((pinValue.x_ - ave_x), 2);
        variance_x = variance_x/cnt;
        variance = variance_x;
    }
    
    if(type == 'V' || type == 'T'){
        int sum = 0;
        int cnt = 0;

        for(pin_value pinValue : pinValues_){
            sum += pinValue.y_;
            cnt++;
        }
        ave_y = sum/(double)cnt;
        for(pin_value pinValue : pinValues_)
            variance_y += pow((pinValue.y_ - ave_y), 2);
        variance_y = variance_y/cnt;
        variance = variance_y;
    }

    if(type == 'T')
        variance = sqrt(pow(variance_x, 2) + pow(variance_y, 2));

    return sqrt(variance);
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

unsigned int Vertex::getNumOfDrc() const {
    return drcValues_.size();
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
Graph::addVertex( int lx, int ly, int ux, int uy, int maxLayer,
                vector<dbInst*> insts,
                vector<wire_value> wireValues,
                vector<via_value> viaValues,
                vector<pin_value> pinValues) {

    Vertex vert(vertices_.size(), lx, ly, ux, uy, maxLayer, 
                insts, wireValues, viaValues, pinValues);

    vertices_.push_back(vert);
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
        unordered_map<unsigned int, pair<char, unsigned int> > each = vert.getEachWireLength();
		
		vert.getNets();
		
        nodeAttr << vert.getId() << ","
				 <<	vert.getNumOfDrc() << ","
                 << vert.getWireCongestion('T') << "," 
                 << vert.getWireCongestion('H') << "," 
                 << vert.getWireCongestion('V') << "," 
                 << vert.getUtilization() << "," 
                 << vert.getRoutingCongestion('T') << "," 
                 << vert.getRoutingCongestion('U') << "," 
                 << vert.getRoutingCongestion('D') << "," 
                 << vert.getRoutingCongestion('L') << "," 
                 << vert.getRoutingCongestion('R') << "," 
                 << vert.getNumOfDrc() << ","
                 << vert.getViaUtilization() << "," 
                 << vert.getNumOfPin() << ","
                 << vert.getStdOfPins('T') << ","
                 << vert.getStdOfPins('H') << ","
                 << vert.getStdOfPins('V') << ","
                 << vert.getTotalWireLength(each) << ","
                 << vert.getWireLength(each, 'H') << ","
                 << vert.getWireLength(each, 'V') << endl;
      
        //int cnt = 0; 
        //for(pair<string, unsigned int> e : each){
        //    cnt++;
        //    nodeAttr << e.second << ","; 
        //}
        //for(; cnt < 7; cnt++){
        //    nodeAttr << "0";
        //    if (cnt == 6) break;
        //    nodeAttr << ",";
        //}
        //nodeAttr << endl;

                 // << vert.getAvgInEdges() << " " 
                 // << vert.getAvgOutEdges() << " " 
                 // << vert.getSequentialRatio() << endl;
    
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

