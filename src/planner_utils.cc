#include "planner_utils.h"
#include <numeric>
#include <opencv2/highgui.hpp>

cv::Mat preprocessMap(const cv::Mat& img) {
    cv::Mat gray, bin;
    cv::cvtColor(img, gray, cv::COLOR_BGR2GRAY);
    cv::threshold(gray, bin, 250, 255, 0);
    cv::Mat erode_kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(10,10));
    cv::morphologyEx(bin, bin, cv::MORPH_ERODE, erode_kernel);
    cv::Mat open_kernel = cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(5,5));
    cv::morphologyEx(bin, bin, cv::MORPH_OPEN, open_kernel);
    return bin;
}

std::vector<std::vector<cv::Point>> polygonize(const cv::Mat& bin_img) {
    std::vector<std::vector<cv::Point>> cnts;
    std::vector<cv::Vec4i> hierarchy;
    cv::findContours(bin_img, cnts, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE);

    std::vector<int> cnt_indices(cnts.size());
    std::iota(cnt_indices.begin(), cnt_indices.end(), 0);
    std::sort(cnt_indices.begin(), cnt_indices.end(), [&cnts](int lhs, int rhs){
        return cv::contourArea(cnts[lhs]) > cv::contourArea(cnts[rhs]);
    });
    int ext_cnt_idx = cnt_indices.empty() ? -1 : cnt_indices.front();
    std::vector<std::vector<cv::Point>> contours;
    if(ext_cnt_idx >= 0){
        contours.push_back(cnts[ext_cnt_idx]);
        for(size_t i = 0; i < hierarchy.size(); ++i){
            if(hierarchy[i][3] == ext_cnt_idx){
                contours.push_back(cnts[i]);
            }
        }
    }

    std::vector<std::vector<cv::Point>> polys;
    std::vector<cv::Point> poly;
    for(auto& contour : contours){
        cv::approxPolyDP(contour, poly, 3, true);
        polys.push_back(poly);
        poly.clear();
    }
    return polys;
}

PolygonWithHoles createPolygonWithHoles(const std::vector<std::vector<cv::Point>>& polys) {
    if(polys.empty()) return PolygonWithHoles();
    std::vector<cv::Point> outer_poly = polys.front();
    std::vector<std::vector<cv::Point>> inner_polys;
    if(polys.size() > 1) inner_polys.assign(polys.begin()+1, polys.end());

    Polygon_2 outer_polygon;
    for(const auto& p : outer_poly){
        outer_polygon.push_back(Point_2(p.x, p.y));
    }

    std::vector<Polygon_2> holes(inner_polys.size());
    for(size_t i = 0; i < inner_polys.size(); ++i){
        for(const auto& pt : inner_polys[i]){
            holes[i].push_back(Point_2(pt.x, pt.y));
        }
    }
    return PolygonWithHoles(outer_polygon, holes.begin(), holes.end());
}

std::vector<std::vector<Point_2>> computeCellSweeps(const std::vector<Polygon_2>& cells, int sweep_step){
    std::vector<std::vector<Point_2>> cells_sweeps;
    for(const auto& cell : cells){
        std::vector<Point_2> cell_sweep;
        Direction_2 best_dir;
        polygon_coverage_planning::findBestSweepDir(cell, &best_dir);
        polygon_coverage_planning::visibility_graph::VisibilityGraph vis_graph(cell);
        bool counter_clockwise = true;
        polygon_coverage_planning::computeSweep(cell, vis_graph, sweep_step, best_dir, counter_clockwise, &cell_sweep);
        cells_sweeps.push_back(cell_sweep);
    }
    return cells_sweeps;
}

std::vector<Point_2> computeDenseCoveragePath(
    const Point_2& start,
    const std::vector<int>& cell_idx_path,
    std::vector<CellNode>& cell_graph,
    const std::vector<Polygon_2>& bcd_cells,
    const std::vector<std::vector<Point_2>>& cells_sweeps,
    const std::vector<std::map<int, std::list<Point_2>>>& cell_intersections)
{
    std::vector<Point_2> way_points;
    Point_2 point = start;
    std::list<Point_2> next_candidates;
    Point_2 next_point;
    std::vector<Point_2> shortest_path;

    if(doReverseNextSweep(start, cells_sweeps[cell_idx_path.front()])){
        shortest_path = getShortestPath(bcd_cells[cell_idx_path.front()], start, cells_sweeps[cell_idx_path.front()].back());
    } else {
        shortest_path = getShortestPath(bcd_cells[cell_idx_path.front()], start, cells_sweeps[cell_idx_path.front()].front());
    }
    way_points.insert(way_points.end(), shortest_path.begin(), std::prev(shortest_path.end()));
    point = way_points.back();

    for(size_t i = 0; i < cell_idx_path.size(); ++i){
        if(!cell_graph[cell_idx_path[i]].isCleaned){
            if(doReverseNextSweep(point, cells_sweeps[cell_idx_path[i]])){
                way_points.insert(way_points.end(), cells_sweeps[cell_idx_path[i]].rbegin(), cells_sweeps[cell_idx_path[i]].rend());
            }else{
                way_points.insert(way_points.end(), cells_sweeps[cell_idx_path[i]].begin(), cells_sweeps[cell_idx_path[i]].end());
            }
            cell_graph[cell_idx_path[i]].isCleaned = true;
            point = way_points.back();
            if((i+1)<cell_idx_path.size()){
                next_candidates = cell_intersections[cell_idx_path[i]].at(cell_idx_path[i+1]);
                if(doReverseNextSweep(point, cells_sweeps[cell_idx_path[i+1]])){
                    next_point = findNextGoal(point, cells_sweeps[cell_idx_path[i+1]].back(), next_candidates);
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i]], point, next_point);
                    way_points.insert(way_points.end(), std::next(shortest_path.begin()), std::prev(shortest_path.end()));
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i+1]], next_point, cells_sweeps[cell_idx_path[i+1]].back());
                }else{
                    next_point = findNextGoal(point, cells_sweeps[cell_idx_path[i+1]].front(), next_candidates);
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i]], point, next_point);
                    way_points.insert(way_points.end(), std::next(shortest_path.begin()), std::prev(shortest_path.end()));
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i+1]], next_point, cells_sweeps[cell_idx_path[i+1]].front());
                }
                way_points.insert(way_points.end(), shortest_path.begin(), std::prev(shortest_path.end()));
                point = way_points.back();
            }
        }else{
            shortest_path = getShortestPath(bcd_cells[cell_idx_path[i]],
                                            cells_sweeps[cell_idx_path[i]].front(),
                                            cells_sweeps[cell_idx_path[i]].back());
            if(doReverseNextSweep(point, cells_sweeps[cell_idx_path[i]])){
                way_points.insert(way_points.end(), shortest_path.rbegin(), shortest_path.rend());
            }else{
                way_points.insert(way_points.end(), shortest_path.begin(), shortest_path.end());
            }
            point = way_points.back();
            if((i+1)<cell_idx_path.size()){
                next_candidates = cell_intersections[cell_idx_path[i]].at(cell_idx_path[i+1]);
                if(doReverseNextSweep(point, cells_sweeps[cell_idx_path[i+1]])){
                    next_point = findNextGoal(point, cells_sweeps[cell_idx_path[i+1]].back(), next_candidates);
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i]], point, next_point);
                    way_points.insert(way_points.end(), std::next(shortest_path.begin()), std::prev(shortest_path.end()));
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i+1]], next_point, cells_sweeps[cell_idx_path[i+1]].back());
                }else{
                    next_point = findNextGoal(point, cells_sweeps[cell_idx_path[i+1]].front(), next_candidates);
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i]], point, next_point);
                    way_points.insert(way_points.end(), std::next(shortest_path.begin()), std::prev(shortest_path.end()));
                    shortest_path = getShortestPath(bcd_cells[cell_idx_path[i+1]], next_point, cells_sweeps[cell_idx_path[i+1]].front());
                }
                way_points.insert(way_points.end(), shortest_path.begin(), std::prev(shortest_path.end()));
                point = way_points.back();
            }
        }
    }
    return way_points;
}
