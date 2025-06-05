#ifndef COVERAGEPLANNER_PLANNER_UTILS_H_
#define COVERAGEPLANNER_PLANNER_UTILS_H_

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

#include "coverage_planner.h"

cv::Mat preprocessMap(const cv::Mat& img);
std::vector<std::vector<cv::Point>> polygonize(const cv::Mat& bin_img);
PolygonWithHoles createPolygonWithHoles(const std::vector<std::vector<cv::Point>>& polys);
std::vector<std::vector<Point_2>> computeCellSweeps(const std::vector<Polygon_2>& cells, int sweep_step);
std::vector<Point_2> computeDenseCoveragePath(
    const Point_2& start,
    const std::vector<int>& cell_idx_path,
    std::vector<CellNode>& cell_graph,
    const std::vector<Polygon_2>& bcd_cells,
    const std::vector<std::vector<Point_2>>& cells_sweeps,
    const std::vector<std::map<int, std::list<Point_2>>>& cell_intersections);

#endif  // COVERAGEPLANNER_PLANNER_UTILS_H_
