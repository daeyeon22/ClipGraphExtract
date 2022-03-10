#include "grid.h"
#include "CImg.h"
#include <iostream>
#include <tinycolormap.hpp>


namespace feature_extractor {
using namespace std;
using namespace cimg_library;
using namespace odb;
using namespace tinycolormap;
typedef cimg_library::CImg< unsigned char > CImgObj;
static const unsigned char yellow[] = {255, 255, 0}, white[] = {255, 255, 255},
                           green[] = {0, 255, 0}, blue[] = {120, 200, 255},
                           darkblue[] = {69, 66, 244},
                           purple[] = {255, 100, 255}, black[] = {0, 0, 0},
                           red[] = {255, 0, 0};


enum ValueType {
  RUDY,
  WIRE_DEN,
  CELL_DEN,
  PIN_DEN,
  CHAN_DEN,
  MARKER,
  LNET_DEN,
  GNET_DEN
};

double getValue(Gcell* gcell, ValueType valType) {
    switch(valType) {
        case ValueType::RUDY:
            return min(1.0, gcell->getRUDY());
        case ValueType::WIRE_DEN:
            return min(1.0, gcell->getWireDensity(ModelType::PL));
        case ValueType::CELL_DEN:
            return min(1.0, gcell->getCellDensity());
        case ValueType::PIN_DEN:
            return min(1.0, gcell->getPinDensity());
        case ValueType::CHAN_DEN:
            return min(1.0, gcell->getChannelDensity(ModelType::PL));
        case ValueType::MARKER:
            return min(1.0, 0.3 * gcell->getNumMarkers());
        case ValueType::LNET_DEN:
            return min(1.0, gcell->getLNetDensity());
        case ValueType::GNET_DEN:
            return min(1.0, gcell->getGNetDensity());
        default:
            return 0.0;
    }
}

double getDenom(Grid* grid, ValueType valType) {
    switch(valType) {
        case ValueType::RUDY:
            return grid->getMaxRUDY();
        case ValueType::WIRE_DEN:
            return 1.0;
        case ValueType::CELL_DEN:
            return 1.0;
        case ValueType::PIN_DEN:
            return 1.0;
        case ValueType::CHAN_DEN:
            return 1.0;
        case ValueType::MARKER:
            return 1.0;
        case ValueType::LNET_DEN:
            return 1.0;
        case ValueType::GNET_DEN:
            return 1.0;
        default:
            return 1.0;
    }
}



void drawGcell(CImgObj* img, Gcell* gcell, Point origin, double scale, double opacity, double value) {
    int x1 = gcell->getBBox().xMin();
    int y1 = gcell->getBBox().yMin();
    int x2 = gcell->getBBox().xMax();
    int y2 = gcell->getBBox().yMax();
    x1 = (int) ( scale * (x1 - origin.getX()) );
    x2 = (int) ( scale * (x2 - origin.getX()) );
    y1 = (int) ( scale * (y1 - origin.getY()) );
    y2 = (int) ( scale * (y2 - origin.getY()) );

    Color tcol = GetColor(value, ColormapType::Heat);
    const unsigned char color[] = { tcol.ri(), tcol.gi(), tcol.bi() };
    img->draw_rectangle(x1, y1, x2, y2, color, opacity);
}




void saveMapImage(Grid* grid, ValueType vtype, string fileName, string dirPath) {

    // img scaling factor
    float sf = 0.001;
    int imgWidth = (int) ( sf * grid->getBoundary().dx() );
    int imgHeight = (int) ( sf * grid->getBoundary().dy() );
    int imgDepth = 1;   
    int imgChannel = 3; // RGB
    int imgSpectrum = 255; // RGBi
    Point origin = grid->getBoundary().ll();
    float opacity = 1.0;
    CImgObj img(imgWidth, imgHeight, imgDepth, imgChannel, imgSpectrum);


    //double denom = getDenom(grid, vtype);

    for(Gcell* gcell : grid->getGcells()) {
    
        double val = getValue(gcell, vtype);
        //val = val / denom;
        
        drawGcell(&img, gcell, origin, sf, opacity, val);
    }

    string imgPath = dirPath + "/" + fileName + ".jpg";
    img.save_jpeg(imgPath.c_str(), 200);

    cout << "[DONE] Save Map Image (" << imgPath << ")" << endl;
}



void Grid::saveMapImages(string dirPath) {

   
    saveMapImage(this, ValueType::RUDY, "RUDY", dirPath);
    saveMapImage(this, ValueType::WIRE_DEN, "WireDen", dirPath);
    saveMapImage(this, ValueType::CELL_DEN, "CellDen", dirPath);
    saveMapImage(this, ValueType::PIN_DEN, "PinDen", dirPath);
    saveMapImage(this, ValueType::CHAN_DEN, "ChanDen", dirPath);
    saveMapImage(this, ValueType::MARKER, "DRV", dirPath);
    saveMapImage(this, ValueType::LNET_DEN, "LNetDen", dirPath);
    saveMapImage(this, ValueType::GNET_DEN, "GNetDen", dirPath);



    /*
    int dbUnitMicron = getDb()->getChip()->getBlock()->getDbUnitsPerMicron();

    // img scaling factor
    float sf = 0.01;
    int imgWidth = (int) ( sf * getBoundary().dx() );
    int imgHeight = (int) ( sf * getBoundary().dy() );
    int imgDepth = 1;   
    int imgChannel = 3; // RGB
    int imgSpectrum = 255; // RGB
    int xOrigin = getBoundary().xMin();
    int yOrigin = getBoundary().yMin();
    float opacity = 1.0;
    string imgName = "";
    string imgPath = "";
 
    double RUDY, PinDen, CellDen, ChanDen, WireDen;
    double HasMarker;

    cout << "imageSize (" << imgWidth << " " << imgHeight << ")" << endl;
    CImgObj img1(imgWidth, imgHeight, imgDepth, imgChannel, imgSpectrum); // RUDY Map
    CImgObj img2(imgWidth, imgHeight, imgDepth, imgChannel, imgSpectrum); // PinDen Map
    CImgObj img3(imgWidth, imgHeight, imgDepth, imgChannel, imgSpectrum); // CellDen Map
    CImgObj img4(imgWidth, imgHeight, imgDepth, imgChannel, imgSpectrum); // ChanDen(Avg) Map
    CImgObj img5(imgWidth, imgHeight, imgDepth, imgChannel, imgSpectrum); // WireDen Map
    CImgObj img6(imgWidth, imgHeight, imgDepth, imgChannel, imgSpectrum); // WireDen Map

    for(Gcell* gcell : gcells_) {

        //gcell->getBBox().print();
        //cout << "   - cell density  : " << gcell->getCellDensity() << endl;
        //cout << "   - pin density   : " << gcell->getPinDensity() << endl;
        //cout << "   - RUDY          : " << gcell->getRUDY() << endl;

        RUDY = min(1.0, gcell->getRUDY());
        WireDen = min(1.0, gcell->getWireDensity(ModelType::PL));
        CellDen = min(1.0, gcell->getCellDensity());
        PinDen = min(1.0, gcell->getPinDensity());
        ChanDen = min(1.0, gcell->getChannelDensity(ModelType::PL));
        HasMarker = (gcell->getNumMarkers() > 0) ? 1.0 : 0.0;
        
        //const unsigned char *color = RLEV[lev];
        int x1 = gcell->getBBox().xMin();
        int y1 = gcell->getBBox().yMin();
        int x2 = gcell->getBBox().xMax();
        int y2 = gcell->getBBox().yMax();
        x1 = (int) ( sf * (x1 - xOrigin) );
        x2 = (int) ( sf * (x2 - xOrigin) );
        y1 = (int) ( sf * (y1 - yOrigin) );
        y2 = (int) ( sf * (y2 - yOrigin) );
        
        // RUDY
        Color tinyColor1 = GetColor(RUDY, ColormapType::Heat);
        const unsigned char color1[] = {tinyColor1.ri(), tinyColor1.gi(), tinyColor1.bi()};
        img1.draw_rectangle(x1,y1,x2,y2,color1,opacity);

        // PinDen
        Color tinyColor2 = GetColor(PinDen, ColormapType::Heat);
        const unsigned char color2[] = {tinyColor2.ri(), tinyColor2.gi(), tinyColor2.bi()};
        img2.draw_rectangle(x1,y1,x2,y2,color2,opacity);

        // CellDen
        Color tinyColor3 = GetColor(CellDen, ColormapType::Heat);
        const unsigned char color3[] = {tinyColor3.ri(), tinyColor3.gi(), tinyColor3.bi()};
        img3.draw_rectangle(x1,y1,x2,y2,color3,opacity);

        // ChanDen
        Color tinyColor4 = GetColor(ChanDen, ColormapType::Heat);
        const unsigned char color4[] = {tinyColor4.ri(), tinyColor4.gi(), tinyColor4.bi()};
        img4.draw_rectangle(x1,y1,x2,y2,color4,opacity);
        
        // WireDen
        Color tinyColor5 = GetColor(WireDen, ColormapType::Heat);
        const unsigned char color5[] = {tinyColor5.ri(), tinyColor5.gi(), tinyColor5.bi()};
        img5.draw_rectangle(x1,y1,x2,y2,color5,opacity);
        
        // WireDen
        Color tinyColor6 = GetColor(HasMarker, ColormapType::Heat);
        const unsigned char color6[] = {tinyColor6.ri(), tinyColor6.gi(), tinyColor6.bi()};
        img6.draw_rectangle(x1,y1,x2,y2,color6,opacity);
    }
    imgPath = dirPath + "/RUDY.jpg";
    img1.save_jpeg(imgPath.c_str(), 200);

    imgPath = dirPath + "/PinDen.jpg";
    img2.save_jpeg(imgPath.c_str(), 200);

    imgPath = dirPath + "/CellDen.jpg";
    img3.save_jpeg(imgPath.c_str(), 200);

    imgPath = dirPath + "/ChanDen.jpg";
    img4.save_jpeg(imgPath.c_str(), 200);

    imgPath = dirPath + "/WireDen.jpg";
    img5.save_jpeg(imgPath.c_str(), 200);
    
    imgPath = dirPath + "/DRV.jpg";
    img6.save_jpeg(imgPath.c_str(), 200);
    */
}




};
//static const unsigned char  RLEV[9][3] = {{255, 245, 235}, {254, 230, 206},{253, 208, 162},
//                                        {253, 174, 107}, {253, 141, 60}, {241, 105, 19},
//                                        {217, 72, 1},  {166, 54, 3},{127, 39, 4}},
//                            BLEV[9][3] = {{252, 251, 253}, {239, 237, 245}, {218, 218, 235},
//                                        {188, 189, 220}, {158, 154, 200}, {128, 125, 186},
//                                        {106, 81, 163}, {84, 39, 143}, {63, 0, 125}};
//

