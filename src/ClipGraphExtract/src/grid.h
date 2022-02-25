#include "opendb/geom.h"



namespace Flute {
    class Tree;
};

namespace odb {
    class dbDatabase;
};



struct MyPoint {
    int x, y;
};

struct MyBox {
    MyPoint p1, p2;
};

struct MyLine {
    MyPoint p1, p2;
};


enum ObjectType {
    NET, INSTANCE, BLOCKAGE, SNET
};

class DrcMarker {
    std::string type_;
    std::string rule_;
    ObjectType objType1_, objType2_;

    void* obj1_;
    void* obj2_;
    odb::Rect bbox_;
};


class RSMT {
  private:
    odb::dbNet* net_;
    odb::Rect bbox_;
    std::vector<odb::Point> terminals_;
    Flute::Tree rsmt_;

    std::vector<Gcell*> rsmtOverlaps_;
    std::vector<Gcell*> bboxOverlaps_;
    std::vector<DrcMarker*> markers_;
    
  public:
    RSMT(odb::dbNet* net) : net_(net) {}

    std::vector<odb::Rect> decomposeRSMT();

    bgBox getQueryBox();
    odb::Rect getBBox();
    void addTerminal(int x, int y);

    bool isLocalNet();
    bool isGlobalNet();
    bool hasDRV();
    void createRSMT();
    int getNetDegree() { return terminals_.size(); }
    int getWireLengthRSMT();
    int getWireLengthHPWL();
    double getWireUniformDensity();


    
};

enum Orient {
    LEFT=0,
    RIGHT=1,
    TOP=2,
    BOTTOM=3
};

struct ResourceModel {
    int trackSupply[4] = {0};
    int trackDemand[4] = {0};
    int wireCapacity;
    int wireLength;
    
    double getWireDensity();
    double getTrackSupply(Orient type);
    double getTrackDemand(Orient type);
    double getWireCapaicty();
    double getWireLength();
};


class Gcell {
  private:   
    // Gcell features
    odb::Rect bbox_;

    ResourceModel rmGR; // using GR results
    ResourceModel rmDR; // using DR results
    ResourceModel rmRSMT; // using PLACE results

    // using placement, RSMT results
    int numInstances_;
    int numTerminals_;
    int numLocalNets_;
    int numGlobalNets_;

    int totalCellArea_;
    int totalPinArea_;

    double cellDensity_;
    double pinDensity_;
    double RUDY_;

    std::vector<odb::dbInst*> insts_;


  public:
    bgBox getQueryBox();
    odb::Rect getBBox(){ return rect_; }
    void updateResourceModelRSMT(odb::Rect seg);
    void updateResourceModelGR(odb::Rect seg);
    void addInst(odb::dbInst* inst);


    //
    void extractFeaturePL(BoxRtree<odb::dbInst*> &rtree);
    void extractFeatureEGR(SegRtree<odb::dbNet*> &rtree);
    void extractFeatureRSMT(SegRtree<RSMT*> &rtree);

    // helper
    int getArea();
    int getCellArea();
    int getPinArea();

};

class Grid {

    odb::Rect coreBBox_;
    int numCols_, numRows_;
    int gridWidth_, gridHeight_;
    int wireCapacity_;
    int trackSupply_;

    odb::dbDatabase* db_;
    std::vector<Gcell> gcells_;
    std::vector<RSMT> rsmts_;
    // After read drc.rpt
    std::vector<DrcMarker> markers_;
    

  public:
    void setCoreArea(odb::Rect& rect);
    void setGcellWidth(int width);
    void setGcellHeight(int height);
    void setWireCapacity(int wCap);
    void setTrackSupply(int tSup);



};




