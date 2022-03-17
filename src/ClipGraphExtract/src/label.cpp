
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
using namespace ClipGraphExtract;
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



string parseType(string substr) {
    return regex_replace(substr, regex("[:\\s]"), "");
}


string parseRule(string substr) {
    return regex_replace(substr, regex("[\\s\\(\\)]"), "");
}

string parseInstName(string substr) {
    string str =substr;
    str = regex_replace(str, regex("Blockage of Cell "), "");
    str = regex_replace(str, regex("Pin of Cell "), "");
    str = regex_replace(str, regex("\\s"), "");
    return str;
}

string parseDrv(string substr) {
    string str =substr;
    str = regex_replace(str, regex("\\s+Total Violations\\s:\\s"), "");
    str = regex_replace(str, regex("\\sViols\\."), "");
    return str;
}

string parseNetName(string substr) {
    string str =substr;
    str = regex_replace(str, regex("Regular Wire of Net "), "");
    str = regex_replace(str, regex("Special Wire of Net "), "");
    str = regex_replace(str, regex("\\s"), "");
	return str;
}

string parseLayerName(string substr) {
    return regex_replace(substr, regex("[\\(\\)\\s]"), "");
}




void ClipGraphExtractor::parseDrcReport(const char* fileName) {

    cout << "Start to read routing report (" << fileName << ")" << endl;

    ifstream inFile(fileName);
    const std::regex colon(":");
    string line;

	dbBlock* block = db_->getChip()->getBlock();
    int dbUnitMicron = block->getDbUnitsPerMicron();

    BoxRtree<Marker*> rtree;

    regex lyrRex("\\( Metal[0-9]+ \\)");
    regex startRex("[\\w]+: \\( [\\w\\s\\d-\\.]+ \\)");
    regex DrvRex("\\s+Total Violations : \\d+ Viols\\.");
    regex typeRex("[\\w]+:");
    regex ruleRex("\\( [\\w\\s\\d-]+ \\)");
    regex objRex1("Blockage of Cell [\\w\\d\\[\\]]+");
    regex objRex2("Pin of Cell [\\w\\d\\[\\]]+");
    regex objRex3("Regular Wire of Net [\\w\\d\\[\\]]+");
    regex objRex4("Special Wire of Net [\\w\\d\\[\\]]+");
    regex boxRex("\\( \\-?[0-9]+\\.[0-9]+, \\-?[0-9]+\\.[0-9]+ \\) \\( \\-?[0-9]+\\.[0-9]+, \\-?[0-9]+\\.[0-9]+ \\)");

    string typeName ="";
    string ruleName ="";
    string lyrName ="";
    
	string fromPrefix = "";
    string toPrefix = "";
 
 	string fromInst ="";
    string toInst ="";
    string fromNet ="";
    string toNet ="";
    
	Grid* grid = (Grid*)grid_;
	
	uint drvNum = 0;

	while(getline(inFile, line)) {
        smatch matStr; 
        string str = line;
        smatch m;
		// Detect parsing start pattern
        if(regex_search(str, m, startRex)) {
			
			// Detect type and delete the corresponding part
            if(regex_search(str, m, typeRex)) {
                //cout << "1" << str << endl;
                typeName = parseType(m[0].str());
                str = regex_replace(str, typeRex, "");
                //cout << str << endl;
            }
		
			// Detect layer and delete the corresponding part
			if(regex_search(str, m, lyrRex)) {
                lyrName = parseLayerName(m[0].str());
                str = regex_replace(str, lyrRex, "");
                //cout << str << endl;
            }
		
			// Detect rule and delete the corresponding part
            if(regex_search(str, m, ruleRex)) {
                ruleName = parseRule(m[0].str());
                //for(int i=0; i < m.size(); i++) {
                //    cout << m[i].str() << endl;
                //}
                str = regex_replace(str, ruleRex, "");
                //cout << str << endl;
            }
			
			// Split object1 and object2
			string delim = "&";
            vector<string> tokens = splitAsTokens(str, delim);
		
            // Parse object and delete the corresponding part
            if(regex_search(tokens[0], m, objRex1)) {
                fromInst = parseInstName(m[0].str());
                fromPrefix = "Blockage of Cell";
                tokens[0] = regex_replace(tokens[0], objRex1, "");
                //cout << tokens[0] << endl;
            } else if(regex_search(tokens[0], m, objRex2)) {
                fromInst = parseInstName(m[0].str());
                fromPrefix = "Pin of Cell";
                tokens[0] = regex_replace(tokens[0], objRex2, "");
                //cout << tokens[0] << endl;
            } else if(regex_search(tokens[0], m, objRex3)) {
                fromNet = parseNetName(m[0].str());
			    fromPrefix = "Regular Wire of Net";
				if(fromNet == "null") {
					fromNet = ""; // In case there is no any inst/net/pin...
			    	fromPrefix = "";
				}

                tokens[0] = regex_replace(tokens[0], objRex3, "");
                //cout << tokens[0] << endl;
            } else if(regex_search(tokens[0], m, objRex4)) {
                //fromNet = parseNetName(m[0].str()); // Special net cannot convert RSMT.
                fromPrefix = "Special Wire of Net";
                tokens[0] = regex_replace(tokens[0], objRex4, "");
                //cout << tokens[0] << endl;
			} else {
                cout << "exception case! here!" << endl;
                cout << tokens[0] << endl;
                exit(0);
            }
			
			if(tokens.size() > 1) {
				// parse object2
				if(regex_search(tokens[1], m, objRex1)) {
					toInst = parseInstName(m[0].str());
					toPrefix = "Blockage of Cell";
					tokens[1] = regex_replace(tokens[1], objRex1, "");
					//cout << tokens[1] << endl;
				} else if(regex_search(tokens[1], m, objRex2)) {
					toInst = parseInstName(m[0].str());
					toPrefix = "Pin of Cell";
					tokens[1] = regex_replace(tokens[1], objRex2, "");
					//cout << tokens[1] << endl;
				} else if(regex_search(tokens[1], m, objRex3)) {
					toNet = parseNetName(m[0].str());
					toPrefix = "Regular Wire of Net";
					tokens[1] = regex_replace(tokens[1], objRex3, "");
					//cout << tokens[1] << endl;
				} else if(regex_search(tokens[0], m, objRex4)) {
					//toNet = parseNetName(m[0].str()); // Special net cannot convert RSMT.
					toPrefix = "Special Wire of Net";
					tokens[0] = regex_replace(tokens[0], objRex4, "");
					//cout << tokens[0] << endl;
				} else {
					//cout << "There is only object1" << endl;
				}
			}
        } else if (regex_search(str, m, boxRex)) {
            string delim = " (),";
            vector<string> tokens = splitAsTokens(m[0].str(), delim);
            if(tokens.size() != 4) {
            	cout << "exception case!" << endl;
				cout << str << endl;
				exit(0);
            }

            int lx = dbUnitMicron * atof(tokens[0].c_str());
            int ly = dbUnitMicron * atof(tokens[1].c_str());
            int ux = dbUnitMicron * atof(tokens[2].c_str());
            int uy = dbUnitMicron * atof(tokens[3].c_str());

            //cout << typeName << " " << ruleName << " " << lyrName 
            //    << " " << fromPrefix << " " << fromInst << fromNet 
            //    << " " << toPrefix << " " << toInst << toNet << endl;

            //cout << "BOUNDS (" << lx << " " << ly << ") (" << ux << " " << uy <<")" << endl;

            
            // TODO
            Marker* mark = grid->createMarker(lx,ly,ux,uy);
            mark->setType(typeName);
            mark->setRule(ruleName);
            //mark->setLayer(
            dbNet* net1 = block->findNet(fromNet.c_str());
            dbNet* net2 = block->findNet(toNet.c_str());
            dbInst* inst1 = block->findInst(fromInst.c_str());
            dbInst* inst2 = block->findInst(toInst.c_str());

            ///////////////////////////////////////////////////////////
            if(fromNet != "") {
                if(net1 == NULL) {
                    cout << fromNet << " is NULL ptr" << endl;
                    exit(0);
                }else{
                    if(grid->getRSMT(net1)==NULL) {
                        cout << fromNet << " has no RSMT" << endl;
                        exit(0);
                    }
                }
            }
            ///////////////////////////////////////////////////////////

            mark->setFromNet(grid->getRSMT(net1));
            mark->setToNet(grid->getRSMT(net2));
            mark->setFromInst(inst1);
            mark->setToInst(inst2);

            if(fromPrefix == "Blockage of Cell") {
                mark->setFromTag(Marker::Tag::BoC);
            } else if (fromPrefix == "Pin of Cell") {
                mark->setFromTag(Marker::Tag::PoC);
            } else if (fromPrefix == "Regular Wire of Net") {
                mark->setFromTag(Marker::Tag::RWoN);
            } else {
                mark->setFromTag(Marker::Tag::NONE);
            }

            if(toPrefix == "Blockage of Cell") {
                mark->setToTag(Marker::Tag::BoC);
            } else if (toPrefix == "Pin of Cell") {
                mark->setToTag(Marker::Tag::PoC);
            } else if (toPrefix == "Regular Wire of Net") {
                mark->setToTag(Marker::Tag::RWoN);
            } else {
                mark->setToTag(Marker::Tag::NONE);
            }

            rtree.insert(make_pair(mark->getQueryBox(), mark));
/*	
			cout << typeName << ":";
			cout << ruleName << ":";
			cout << fromPrefix << ":";
			cout << fromInst << ":";
			cout << fromNet << ":";
			cout << toPrefix << ":";
			cout << toInst << ":";
			cout << toNet << ":";
			cout << lyrName << ":";
			cout << endl;
            cout << "(" << tokens[0] << " " << tokens[1] << ") (" << tokens[2] << " " << tokens[3] <<")" << endl;
			cout << endl;
*/		
			typeName ="";
			ruleName ="";
			lyrName ="";
			fromPrefix = "";
			toPrefix = "";
			fromInst ="";
			toInst ="";
			fromNet ="";
			toNet ="";
			
			drvNum++;

        } else if (regex_search(str, m, DrvRex)) {
			if(stoi(parseDrv(m[0].str())) != drvNum){
				cout << "The number of DRVs is different." << endl;
				cout << "rptNum: " << parseDrv(m[0].str()) << " parseNum: " << drvNum << endl;
				exit(0);
			}	
		} else {
            //cout << "exception case!" << endl;
			//cout << str << endl;
			//exit(0);
        }
    }
    // labeling
    for(Gcell* gcell : grid->getGcells()) {
        gcell->annotateLabel(rtree);

        if(gcell->getNumMarkers() > 0)
            gcell->print();
    
    }
}



namespace ClipGraphExtract {

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

void Grid::reportDRC() {
    int nL2L=0;
    int nL2I=0;
    int nL2G=0;
    int nG2I=0;
    int nG2G=0;
    int nI2I=0;
    int nSELF=0;
    int nERR=0;
    unordered_map<string,int> type2count;

    for(Marker* mark : markers_) {

        Marker::Category ctgy = mark->getCategory();

        switch(ctgy) {
            case Marker::Category::L2L:
                nL2L++;	break;
            case Marker::Category::L2G:
                nL2G++;	break;
            case Marker::Category::L2I:
                nL2I++;	break;
            case Marker::Category::G2I:
                nG2I++;	break;
            case Marker::Category::G2G:
                nG2G++;break;
            case Marker::Category::I2I:
                nI2I++;	break;
            case Marker::Category::SELF:
                nSELF++;break;
            default:
				nERR++;	break;
        }
        type2count[mark->getType()]++;
    }

    uint lnet, gnet, inst;



    cout << "= = = = Report DRC = = = =" << endl;
    cout << " Total #DRVs is " << markers_.size() << endl;
    cout << "   - L2L : " << nL2L << endl;
    cout << "   - L2G : " << nL2G << endl;
    cout << "   - L2I : " << nL2I << endl;
    cout << "   - G2I : " << nG2I << endl;
    cout << "   - G2G : " << nG2G << endl;
    cout << "   - I2I : " << nI2I << endl;
    cout << "   - SELF : " << nSELF << endl;
    cout << " # of exception case =" << nERR << endl;
    cout << " #DRVs due to Global =" << nL2G + nG2I + nG2G << endl;
    cout << " #DRVs due to Local =" << nL2L + nL2G + nL2I << endl;
    cout << " #DRVs due to Inst =" << nL2I + nG2I << endl;

    for(auto elem : type2count)
    cout << " #DRVs due to " << elem.first << " = " << elem.second << endl;
    cout << "= = = = = = = = = = = = = =" << endl;

}






};
