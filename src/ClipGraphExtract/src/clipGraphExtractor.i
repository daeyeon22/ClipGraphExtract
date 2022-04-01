%module ClipGraphExtractor

%{
#include "openroad/OpenRoad.hh"
#include "clip_graph_ext/clipGraphExtractor.h"

namespace ord {
ClipGraphExtract::ClipGraphExtractor*
getClipGraphExtractor(); 
odb::dbDatabase*
getDb();
sta::dbSta*
getSta();
}


namespace ClipGraphExtract { 
enum GraphModel;
enum EdgeWeightModel;
}

using ord::getClipGraphExtractor;
using ord::getDb;
using ord::getSta;
using ClipGraphExtract::ClipGraphExtractor;
using ClipGraphExtract::GraphModel;
using ClipGraphExtract::EdgeWeightModel;

%}

%inline %{

void
save_grid_images_cmd(const char* imgDir) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->saveGridImages(imgDir);
}

void
parse_drc_report_cmd(const char* file_name) {

    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->parseDrcReport(file_name);

}

void
set_gcell_size_cmd(int numRows) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->setGcellSize(numRows);
}

void
set_max_route_layer_cmd(int maxRouteLayer) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->setMaxRouteLayer(maxRouteLayer);
}


void
set_graph_model_cmd(const char* model) 
{
  ClipGraphExtractor* graphExt = getClipGraphExtractor();
  graphExt->setGraphModel(model);
}

void
set_edge_weight_model_cmd(const char* model)
{
  ClipGraphExtractor* graphExt = getClipGraphExtractor();
  graphExt->setEdgeWeightModel(model);
}

void
set_graph_extract_save_file_name_cmd(const char* file)
{
  ClipGraphExtractor* graphExt = getClipGraphExtractor();
  graphExt->setSaveFileName(file);
}

void
set_graph_extract_save_file_prefix_cmd(const char* prefix)
{
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->setSaveFilePrefix(prefix);
}



void 
graph_extract_init_cmd()
{
  ClipGraphExtractor* graphExt = getClipGraphExtractor();
  graphExt->setDb(getDb());
  graphExt->setSta(getSta());
  graphExt->init();
}


void 
save_features_cmd(const char* dirPath, int numHops) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->saveFeatures(dirPath, numHops);
}

void
save_features_cmd(const char* dirPath) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->saveFeatures(dirPath);
}

void
save_graphs_cmd(const char* dirPath) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->saveGraphs(dirPath);
}

void
save_labels_cmd(const char* dirPath) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->saveLabels(dirPath);
}

void
graph_extract_cmd(int lx, int ly, int ux, int uy) 
{
  ClipGraphExtractor* graphExt = getClipGraphExtractor();
  graphExt->extract(lx, ly, ux, uy);
}

void
graph_extract_clear_cmd() 
{
  ClipGraphExtractor* graphExt = getClipGraphExtractor();
  graphExt->clear();
}

%}
