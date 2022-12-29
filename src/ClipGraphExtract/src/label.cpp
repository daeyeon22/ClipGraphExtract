
#include "opendb/geom.h"
#include "grid.h"
#include "clip_graph_ext/clipGraphExtractor.h"
#include "opendb/db.h"
#include "opendb/dbWireCodec.h"
#include "opendb/dbTransform.h"
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

string parseRegularNetName(string substr) {
    string str =substr;
    str = regex_replace(str, regex("Regular Wire of Net "), "");
    str = regex_replace(str, regex("Special Wire of Net "), "");
    str = regex_replace(str, regex("\\s"), "");
	return str;
}

string parseSpecialNetName(string substr) {
    string str =substr;
    str = regex_replace(substr, regex("Special Wire of Net "), "");
    str = regex_replace(str, regex("\\s"), "");
	return str;
}


string parseLayerName(string substr) {
    return regex_replace(substr, regex("[\\(\\)\\s]"), "");
}





void ClipGraphExtractor::parseDrcReport_(const char* fileName) {
    cout << "Start to read routing report (" << fileName << ")" << endl;
    ifstream inFile(fileName);
    const std::regex colon(":");
    string line;
	dbBlock* block = db_->getChip()->getBlock();
    int dbu = block->getDbUnitsPerMicron();


    BoxRtree<Marker*> rtree;
    regex ex("Bounds \\([0-9]+\\.[0-9]+, [0-9]+\\.[0-9]+\\) \\([0-9]+\\.[0-9]+, [0-9]+\\.[0-9]+\\)");
    
    int numDrvs = 0;
	Grid* grid = (Grid*)grid_;
    
    while(getline(inFile, line)) {
        //cout << line << endl;    
        smatch m;
        if(regex_search(line, m, ex)) {
            cout << m[0] << endl;
            numDrvs++;
            string delim = " (),";
            vector<string> tokens = splitAsTokens(m[0].str(), delim);
            int lx = dbu * atof(tokens[0].c_str());
            int ly = dbu * atof(tokens[1].c_str());
            int ux = dbu * atof(tokens[2].c_str());
            int uy = dbu * atof(tokens[3].c_str());
            // TODO
            Marker* mark = grid->createMarker(lx,ly,ux,uy);
            rtree.insert(make_pair(mark->getQueryBox(), mark));
        }
    }

    cout << "num drvs  = " << numDrvs << endl;
    for(Gcell* gcell : grid->getGcells()) {
        gcell->annotateLabel(rtree);
    }
}







void ClipGraphExtractor::parseDrcReport(const char* fileName) {

    cout << "Start to read routing report (" << fileName << ")" << endl;

    ifstream inFile(fileName);
    const std::regex colon(":");
    string line;

	dbBlock* block = db_->getChip()->getBlock();
    int dbu = block->getDbUnitsPerMicron();

    BoxRtree<Marker*> rtree;

    regex lyrRex("\\( m[0-9]+ \\)");
    regex startRex("[\\w]+: \\( [\\w\\s\\d-\\.]+ \\)");
    regex DrvRex("\\s+Total Violations : \\d+ Viols\\.");
    regex typeRex("[\\w]+:");
    regex ruleRex("\\( [\\w\\s\\d-]+ \\)");
    regex objRex1("Blockage of Cell [\\w\\d\\[\\]]+");
    regex objRex2("Pin of Cell [\\w\\d\\[\\]]+");
    regex objRex3("Regular Wire of Net [\\w\\d\\[\\]]+");
    regex objRex4("Special Wire of Net [\\w\\d\\[\\]]+");
    regex boxRex("\\( [0-9]+\\.[0-9]+, [0-9]+\\.[0-9]+ \\) \\( [0-9]+\\.[0-9]+, [0-9]+\\.[0-9]+ \\)");

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
	
	int drvNum = 0;

    typedef pair<dbNet*, dbInst*> dbValue;
    unordered_map<int, BoxRtree<dbValue>> rtrees;

    //unordered_map<dbInst*, int> numDrvs_;
    for(dbInst* tarInst : block->getInsts()) {
        numDrvs_[tarInst] = 0;

        int xOrig, yOrig;
        tarInst->getOrigin(xOrig, yOrig);
        dbTransform transform;
        transform.setOrient(tarInst->getOrient());
        transform.setOffset(Point(xOrig, yOrig));

        for(dbITerm* tarIterm : tarInst->getITerms()) {
            dbMTerm* tarMterm = tarIterm->getMTerm();
            dbNet* tarNet = tarIterm->getNet();
            for(dbMPin* tarMpin : tarMterm->getMPins()) {
                for(dbBox* tarBox : tarMpin->getGeometry()) {
                    int routingLevel = tarBox->getTechLayer()->getRoutingLevel();
                    if(routingLevel == 1 || routingLevel == 2) {
                        Rect tarBBox;
                        tarBox->getBox(tarBBox);
                        transform.apply(tarBBox);
                        int xMin = tarBBox.xMin();
                        int yMin = tarBBox.yMin();
                        int xMax = tarBBox.xMax();
                        int yMax = tarBBox.yMax();

                        rtrees[routingLevel].insert(
                                make_pair(bgBox(bgPoint(xMin, yMin), bgPoint(xMax,yMax)), dbValue(tarNet, tarInst))
                        );
                    
               
                    }
                }
            }
        }
    }


    int x, y, ext;
    for(dbNet* net : block->getNets()) {

        unordered_map<int, vector<pair<bgBox, dbValue>>> accessMetals;
        
        dbWire* wire = net->getWire();
        if( wire && wire->length() ) {
            int wl = wire->getLength();
            int wl_ = 0;
            dbWireDecoder decoder;
            decoder.begin(wire);
            vector<odb::Point> points;
            dbWireDecoder::OpCode opcode = decoder.peek();
            dbTechLayer* tarLayer;
            dbTechVia* techVia;
            while(opcode != dbWireDecoder::END_DECODE) {
                bool hasPath=false;
                switch(opcode) {
                    case dbWireDecoder::PATH:
                        tarLayer = decoder.getLayer();
                        //cout << "PATH begin " << endl;
                        //if(tarLayer != NULL)
                        //    cout << tarLayer->getName() << endl;
                        points.clear(); break;
                    case dbWireDecoder::JUNCTION: 
                        //cout << "JUNCTION" << endl;
                        break;

                    case dbWireDecoder::SHORT: 
                        //cout << "SHORT" << endl;
                        break;
                    case dbWireDecoder::TECH_VIA: {
                        techVia = decoder.getTechVia();
                        decoder.getPoint(x,y);
                        points.push_back(odb::Point(x,y));
                        //cout << "TECHVIA " << techVia->getName() << " (" << x << " " << y << ")" << endl;
                        bool isAccessed = false;
                        dbInst* tarInst;
                        for(dbBox* tarBox : techVia->getBoxes()) {
                            int xMin = x + tarBox->xMin();
                            int yMin = y + tarBox->yMin();
                            int xMax = x + tarBox->xMax();
                            int yMax = y + tarBox->yMax();
                
                            int routingLevel = tarBox->getTechLayer()->getRoutingLevel();

                            if(routingLevel == 1 || routingLevel == 2) {
                                //cout << tarBox->getTechLayer()->getName() 
                                //    << " (" << xMin << " " << yMin << ")"
                                //    << " (" << xMax << " " << yMax << ")" << endl;
                            
                                bgBox tBox(bgPoint(xMin, yMin), bgPoint(xMax, yMax));
                                if(!isAccessed) {
                                    vector<pair<bgBox,dbValue>> queryResults;
                                    rtrees[routingLevel].query(bgi::intersects(tBox), back_inserter(queryResults));
                                    for(auto it : queryResults) {
                                        dbValue val = it.second;
                                        dbNet* curNet = val.first;
                                        dbInst* curInst = val.second;

                                        if(curNet == net) {
                                            isAccessed = true;
                                            tarInst = curInst;
                                            break;
                                        }
                                    }  
                                }
                                if(isAccessed) {
                                    accessMetals[routingLevel].push_back(make_pair(tBox,dbValue(net,tarInst)));
                                }
                            }
                        }

                        hasPath=true;
                        break;
                    }
                    case dbWireDecoder::VIA: 
                        //cout << "VIA" << endl;                          
                        break;
                    case dbWireDecoder::VWIRE: 
                        //cout << "VWIRE" << endl;
                        break;
                    case dbWireDecoder::POINT_EXT:
                        decoder.getPoint(x,y,ext);
                        points.push_back(odb::Point(x,y));
                        //cout << "POINT_EXT (" << x << " " << y << ")" << endl;
                        hasPath=true;
                        break;
                    case dbWireDecoder::BTERM: 
                        //cout << "BTERM" << endl;
                        break;
                    case dbWireDecoder::ITERM: 
                        //cout << "ITERM" << endl;
                        break;
                    case dbWireDecoder::RULE: 
                        //cout << "RULE" << endl;
                        break;
                    case dbWireDecoder::POINT: {
                        decoder.getPoint(x,y);
                        points.push_back(odb::Point(x,y));
                        //cout << "POINT (" << x << " " << y << ")" <<  endl;
                        hasPath=true;

                        
                        
                        break;
                    }
                    default: break;
                }
                if(points.size() > 1) {
                    odb::Point pt1 = points[points.size()-2];
                    odb::Point pt2 = points[points.size()-1];
                    int xMin = min(pt1.getX(), pt2.getX());// - layerMinWidth[decoder.getLayer()]/2;
                    int xMax = max(pt1.getX(), pt2.getX());// + layerMinWidth[decoder.getLayer()]/2;
                    int yMin = min(pt1.getY(), pt2.getY());// - layerMinWidth[decoder.getLayer()]/2;
                    int yMax = max(pt1.getY(), pt2.getY());// + layerMinWidth[decoder.getLayer()]/2;
                    int dist = (xMax-xMin) + (yMax-yMin);
                    int width = tarLayer->getWidth();
                    xMin -= width/2;
                    xMax += width/2;
                    yMin -= width/2;
                    yMax += width/2;
                    int routingLevel = tarLayer->getRoutingLevel();

                    if( routingLevel == 1 || routingLevel == 2 ) {

                        //cout << net->getName() << " " << tarLayer->getName() << " (" << xMin <<  " " << yMin << ") (" << xMax << " " << yMax <<")" << endl;
                        vector<pair<bgBox,dbValue>> queryResults;
                        bgBox tarBox(bgPoint(xMin, yMin), bgPoint(xMax, yMax));
                        rtrees[routingLevel].query(bgi::intersects(tarBox), back_inserter(queryResults));
                        for(auto it : queryResults) {
                            dbValue val = it.second;
                            dbNet* curNet = val.first;
                            dbInst* curInst = val.second;
                            if(curNet == net) {
                                accessMetals[routingLevel].push_back(make_pair(tarBox,dbValue(net, curInst)));
                                break;
                            }
                        }
                    }
                }
                /*

                */
                opcode = decoder.next();
            }



            for(int i=1; i<=2; i++)
                for(auto val : accessMetals[i])
                    rtrees[i].insert(val);
        }

        

        /*
        for(dbSWire* swire : net->getSWires()) {
            for(dbSBox* sbox : swire->getWires()) {
                int xMin = sbox->xMin();
                int yMin = sbox->yMin();
                int xMax = sbox->xMax();
                int yMax = sbox->yMax();
                bgBox wireBox( bgPoint(xMin, yMin), bgPoint(xMax, yMax) );
                swireRtree.insert( make_pair( wireBox, net ) );
            }
        }
        */
    }
 

    //cout << "Rtree size (1) : " << rtrees[1].size() << endl;
    //cout << "Rtree size (2) : " << rtrees[2].size() << endl;

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
                fromNet = parseRegularNetName(m[0].str());
                fromPrefix = "Regular Wire of Net";
                tokens[0] = regex_replace(tokens[0], objRex3, "");
                //cout << tokens[0] << endl;
            } else if(regex_search(tokens[0], m, objRex4)) {
                //cout << "1 " <<  m[0] << endl;
                fromNet = parseSpecialNetName(m[0].str());
                fromPrefix = "Special Wire of Net";
                tokens[0] = regex_replace(tokens[0], objRex4, "");
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
					toNet = parseRegularNetName(m[0].str());
					toPrefix = "Regular Wire of Net";
					tokens[1] = regex_replace(tokens[1], objRex3, "");
					//cout << tokens[1] << endl;
                } else if(regex_search(tokens[0], m, objRex4)) {
                    //cout << "2 " << m[0] << endl;
                    toNet = parseSpecialNetName(m[0].str());
                    toPrefix = "Special Wire of Net";
                    tokens[0] = regex_replace(tokens[0], objRex4, "");
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

            int lx = dbu * atof(tokens[0].c_str());
            int ly = dbu * atof(tokens[1].c_str());
            int ux = dbu * atof(tokens[2].c_str());
            int uy = dbu * atof(tokens[3].c_str());
            // TODO
            Marker* mark = grid->createMarker(lx,ly,ux,uy);
            mark->setType(typeName);
            mark->setRule(ruleName);
            //mark->setLayer(
            dbNet* net1 = block->findNet(fromNet.c_str());
            dbNet* net2 = block->findNet(toNet.c_str());
            dbInst* inst1 = block->findInst(fromInst.c_str());
            dbInst* inst2 = block->findInst(toInst.c_str());
            dbTechLayer* tarLayer = db_->getTech()->findLayer(lyrName.c_str());
            mark->setLayer(tarLayer);
            set<dbInst*> violateInsts;
            //if(inst1 != NULL)
            //    violateInsts.insert(inst1);
            //if(inst2 != NULL)
            //    violateInsts.insert(inst2);
            if( tarLayer != NULL ) {
                int routingLevel = tarLayer->getRoutingLevel();
                if( routingLevel ==1 || routingLevel ==2 ) {
                    //bgBox tarBox(bgPoint(lx-10,ly-10), bgPoint(ux+10, uy+10));
                    bgBox tarBox(bgPoint(lx,ly), bgPoint(ux, uy));
                    vector<pair<bgBox,dbValue>> queryResults;
                    rtrees[routingLevel].query(bgi::intersects(tarBox), back_inserter(queryResults));
                    for(auto it : queryResults) {
                        dbValue val = it.second;
                        dbNet* tarNet = val.first;
                        dbInst* tarInst = val.second;
                        violateInsts.insert(tarInst);
                    }
                }
            }
            //cout << typeName << " " << ruleName << " " << lyrName 
            //    << " " << fromPrefix << " " << fromInst << fromNet 
            //    << " " << toPrefix << " " << toInst << toNet << endl;
            //cout << "BOUNDS (" << 1.0*lx/dbu << " " << 1.0*ly/dbu << ") (" << 1.0*ux/dbu << " " << 1.0*uy/dbu <<")" << endl;


            for(dbInst* tarInst : violateInsts) {
                //cout << "   - " << tarInst->getName() << endl;
                numDrvs_[tarInst]++;
            }

            
            //if(inst1 != NULL) {
            //    numDrvs_[inst1]++;
            //}

            //if(inst2 != NULL) {
            //    numDrvs_[inst2]++;
            //}
            ///////////////////////////////////////////////////////////
            if(fromNet != "") {
                if(net1 == NULL) {
                    cout << fromNet << " is NULL ptr" << endl;
                    //exit(0);
                }else{
                    if(grid->getRSMT(net1)==NULL) {
                        cout << fromNet << " has no RSMT" << endl;
                        //exit(0);
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
            } else if (fromPrefix == "Special Wire of Net") {
                mark->setFromTag(Marker::Tag::SWoN);
            } else {
                mark->setFromTag(Marker::Tag::NONE);
            }

            if(toPrefix == "Blockage of Cell") {
                mark->setToTag(Marker::Tag::BoC);
            } else if (toPrefix == "Pin of Cell") {
                mark->setToTag(Marker::Tag::PoC);
            } else if (toPrefix == "Regular Wire of Net") {
                mark->setToTag(Marker::Tag::RWoN);
            } else if (toPrefix == "Special Wire of Net") {
                mark->setToTag(Marker::Tag::SWoN);
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
				cout << "parseNum: " << parseDrv(m[0].str()) << " OriginalNum: " << drvNum << endl;
                //exit(0);
			}	
		} else {
            //cout << "exception case!" << endl;
			//cout << str << endl;
			//exit(0);
        }
    }

    cout << "# of Markers : " <<  rtree.size() << endl;

    // labeling
    for(Gcell* gcell : grid->getGcells()) {
        gcell->annotateLabel(rtree);
        //gcell->getGraph()->setNumDrvs(numDrvs_);
        //if(gcell->getNumMarkers() > 0)
        //    gcell->print();
    
    }
    //int total = 0;
    //for(auto it : numDrvs_) {
    //    dbInst* tarInst = it.first;
    //    int count = it.second;
    //    cout << tarInst->getName() << " -> " << count << endl;
    //    total += count;
    //}

    //cout << "# of total DRVs : " << total << endl;
    
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
        Marker* tarMark = val.second;
        markers_.push_back(tarMark);
   
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

    int lnet, gnet, inst;



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
