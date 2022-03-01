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
construct_gcell_grid_cmd(int num_rows, int max_layer) 
{
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->setDb(getDb());
    graphExt->initGcellGrid(num_rows, max_layer);
}

void
save_map_images_cmd(const char* imgDir) {
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->saveMapImages(imgDir);
}



void
analyze_congestion_cmd() 
{
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->showCongestionMap();

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
graph_extract_cmd(int lx, int ly, int ux, int uy) 
{
  ClipGraphExtractor* graphExt = getClipGraphExtractor();
  graphExt->extract(lx, ly, ux, uy);
}

void

bin_graph_extract_cmd(int num_rows, int max_layer) 
{
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->extractBinGraph(num_rows, max_layer);
}

void
bin_graph_labeling_cmd(const char* inv_rpt_file) 
{
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->labelingBinGraph(inv_rpt_file);
}

void
save_bin_graph_file_cmd()
{
    ClipGraphExtractor* graphExt = getClipGraphExtractor();
    graphExt->saveBinGraph();
}


void
graph_extract_clear_cmd() 
{
  ClipGraphExtractor* graphExt = getClipGraphExtractor();
  graphExt->clear();
}

%}
