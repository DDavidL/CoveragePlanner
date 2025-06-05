#include "planner_utils.h"

#define DENSE_PATH

int main(){
    cv::Mat img = cv::imread("../data/basement.png");
    cv::Mat bin = preprocessMap(img);

    auto polys = polygonize(bin);
    PolygonWithHoles pwh = createPolygonWithHoles(polys);

    std::vector<Polygon_2> bcd_cells;
    polygon_coverage_planning::computeBestBCDFromPolygonWithHoles(pwh, &bcd_cells);

    auto cell_graph = calculateDecompositionAdjacency(bcd_cells);
    Point_2 start = getStartingPoint(img);
    int starting_cell_idx = getCellIndexOfPoint(bcd_cells, start);
    auto cell_idx_path = getTravellingPath(cell_graph, starting_cell_idx);

    int sweep_step = 5;
    auto cells_sweeps = computeCellSweeps(bcd_cells, sweep_step);
    auto cell_intersections = calculateCellIntersections(bcd_cells, cell_graph);

#ifdef DENSE_PATH
    auto way_points = computeDenseCoveragePath(start, cell_idx_path, cell_graph,
                                              bcd_cells, cells_sweeps,
                                              cell_intersections);
    cv::namedWindow("cover",cv::WINDOW_NORMAL);
    cv::imshow("cover", img);
    cv::waitKey();
    for(size_t i = 1; i < way_points.size(); ++i){
        cv::Point p1(CGAL::to_double(way_points[i-1].x()),CGAL::to_double(way_points[i-1].y()));
        cv::Point p2(CGAL::to_double(way_points[i].x()),CGAL::to_double(way_points[i].y()));
        cv::line(img, p1, p2, cv::Scalar(0, 64, 255));
        cv::namedWindow("cover",cv::WINDOW_NORMAL);
        cv::imshow("cover", img);
    }
    cv::waitKey(1000);
#endif
    return 0;
}
