
#include "opendb/geom.h"
#include "grid.h"
#include "clip_graph_ext/clipGraphExtractor.h"
#include "opendb/db.h"
#include "opendb/dbWireCodec.h"
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <fstream>

#include <regex>
#include <sstream>


using namespace std;
using namespace odb;
using namespace feature_extractor;
using namespace ClipGraphExtract;

vector<string> splitAsTokens(string str, string delim){
    vector<string> _tokens;
    size_t start, end=0;
    while(end < str.size()){
        start = end;
        while(start < str.size() && (delim.find(str[start]) != string::npos)){
            start++;
        }
        end = start;
        while(end < str.size() && (delim.find(str[end]) == string::npos)){
            end++;
        }
        if(end-start != 0){
            _tokens.push_back(string(str, start, end-start));
        }
    }
    return _tokens;
}

void ClipGraphExtractor::readRoutingReport(const char* fileName) {


    cout << "Start to read routing report (" << fileName << ")" << endl;

    ifstream inFile(fileName);
    const std::regex colon(":");
    string line;

	dbBlock* block = db_->getChip()->getBlock();
    int dbUnitMicron = block->getDbUnitsPerMicron();

	int lineNum = 0;
	string type = "0";
	string detailed = "0";
	string toNet = "0";
	string fromNet = "0";
	string toInst = "0";
    int layer=0;


    BoxRtree<Marker*> rtree;

    Grid* grid = (Grid*)grid_;

    // NEED TO BE GENERALIZED!
    // TO JKLEE
    while(getline(inFile, line)) {
        lineNum++;
        if(lineNum < 10) continue;

        std::smatch match;
        if(regex_search(line, match, colon)) {
            string head = match.prefix();
            string tail = match.suffix();

            if(head != "Bounds "){
                if(head == "  Total Violations ") continue;

                type = head;
                string delim = "()";
                vector<string> tokens = splitAsTokens(tail, delim);
                ZASSERT(tokens.size() == 4);
                for(int i = 0; i < tokens.size(); i++)
                    tokens[i] = tokens[i].substr(1, tokens[i].size()-2);

                detailed = tokens[1];
                layer = atoi(tokens[3].substr(1).c_str());

                delim = "&";
                tokens = splitAsTokens(tokens[2], delim);
                ZASSERT(tokens.size() < 3);

                if(tokens[0].substr(0, 19) == "Regular Wire of Net")
                    toNet = tokens[0].substr(20, tokens[0].size()-21);

                else if(tokens[0].substr(0, 11) == "Pin of Cell")
                    toInst = tokens[0].substr(12, tokens[0].size()-12);

                else cout << "outlier: " << line << endl;

                if(tokens.size() > 1){
                    tokens[1] = tokens[1].substr(1, tokens[1].size()-2);
                    if(tokens[1].substr(0, 19) == "Regular Wire of Net")
                        fromNet = tokens[1].substr(20, tokens[1].size()-20);
                }
            } else{
                string delim = " (),";
                vector<string> tokens = splitAsTokens(tail, delim);
                if(tokens.size() != 4) {			
                    continue;
                }

                int lx = dbUnitMicron * atof(tokens[0].c_str());
                int ly = dbUnitMicron * atof(tokens[1].c_str());
                int ux = dbUnitMicron * atof(tokens[2].c_str());
                int uy = dbUnitMicron * atof(tokens[3].c_str());

                //cout << type << " " << detailed << " " << toNet << " " << fromNet << " " << toInst << " " << layer << " ";
                //cout << "(" << lx << " " << ly << ") (" << ux << " " << uy << ")" << endl;

                // TODO
                Marker* mark = grid->createMarker(lx,ly,ux,uy);
                mark->setType(type);
                mark->setRule(detailed);
                //mark->setBoundary(Rect(lx, ly, ux, uy));
                mark->setFromNet(block->findNet(fromNet.c_str()));
                mark->setToNet(block->findNet(toNet.c_str()));
                mark->setToInst(block->findInst(toInst.c_str()));

                if(mark->getFromNet() != NULL && mark->getToNet() != NULL) {
                    mark->setTag(Marker::Tag::N2N);
                }

                if(mark->getFromNet() != NULL && mark->getToInst() != NULL) {
                    mark->setTag(Marker::Tag::N2I);
                }

                rtree.insert(make_pair(mark->getQueryBox(), mark));
                type = "0";
                detailed = "0";
                toNet = "0";
                fromNet = "0";
                toInst = "0";
                layer = 0;
            }
        }
    }
    // labeling
    for(Gcell* gcell : grid->getGcells()) {
        gcell->annotateLabel(rtree);

        if(gcell->getNumMarkers() > 0)
            gcell->print();
    
    }
}



namespace feature_extractor {

Marker* Grid::createMarker(int x1, int y1, int x2, int y2) {
    Marker* mark = new Marker();
    mark->setBoundary(Rect(x1,y1,x2,y2));
    markers_.push_back(mark);
    return mark;
}


void Gcell::annotateLabel(BoxRtree<Marker*> &rtree) {
    vector<pair<bgBox, Marker*>> queryResults;
    rtree.query(bgi::intersects(getQueryBox()), back_inserter(queryResults));
    for(auto& val : queryResults) {
        Marker* mark = val.second;
        markers_.push_back(mark);
    }
}












};

