/*
 * polygon_coverage_planning implements algorithms for coverage planning in
 * general polygons with holes. Copyright (C) 2019, Rik Bähnemann, Autonomous
 * Systems Lab, ETH Zürich
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <variant> // Required for std::get_if for point location result

#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>

#include "cgal_comm.h"
#include "visibility_polygon.h"

namespace polygon_coverage_planning {

bool computeVisibilityPolygon(const PolygonWithHoles& pwh,
                              const Point_2& query_point,
                              Polygon_2* visibility_polygon) {

  // Preconditions.
  if(!pointInPolygon(pwh, query_point)){
      std::cout<<"Query point outside of polygon."<<std::endl;
  }

  if(!isStrictlySimple(pwh)){
      std::cout<<"Polygon is not strictly simple."<<std::endl;
  }

  // Create 2D arrangement.
  typedef CGAL::Arr_segment_traits_2<K> VisibilityTraits;
  typedef CGAL::Arrangement_2<VisibilityTraits> VisibilityArrangement;
  VisibilityArrangement poly;
  CGAL::insert(poly, pwh.outer_boundary().edges_begin(),
               pwh.outer_boundary().edges_end());
  // Store main face.

  VisibilityArrangement::Face_const_handle main_face = poly.faces_begin();
  while (main_face->is_unbounded()) {
    main_face++;
  }

  for (PolygonWithHoles::Hole_const_iterator hit = pwh.holes_begin();
       hit != pwh.holes_end(); ++hit)
    CGAL::insert(poly, hit->edges_begin(), hit->edges_end());

  // Create Triangular Expansion Visibility object.
  typedef CGAL::Triangular_expansion_visibility_2<VisibilityArrangement,
                                                  CGAL::Tag_true>
      TEV;
  TEV tev(poly);

  // We need to determine the halfedge or face to which the query point
  // corresponds.
  typedef CGAL::Arr_naive_point_location<VisibilityArrangement> NaivePL;
  typedef CGAL::Arr_point_location_result<VisibilityArrangement>::Type PLResult;
  NaivePL pl(poly);
  PLResult pl_result = pl.locate(query_point);

  VisibilityArrangement::Vertex_const_handle* v = nullptr;
  VisibilityArrangement::Halfedge_const_handle* e = nullptr;
  VisibilityArrangement::Face_const_handle* f_ptr = nullptr; // Renamed to avoid confusion if f is used later
  VisibilityArrangement::Vertex_const_handle* v_ptr = nullptr; // Renamed
  VisibilityArrangement::Halfedge_const_handle* e_ptr = nullptr; // Renamed

  typedef VisibilityArrangement::Face_handle VisibilityFaceHandle;
  VisibilityFaceHandle fh;
  VisibilityArrangement visibility_arr;

  // In CGAL 6.0, Arr_point_location_result::Type is likely std::variant
  if ((f_ptr = std::get_if<VisibilityArrangement::Face_const_handle>(&pl_result))) {
    // Located in face.
    fh = tev.compute_visibility(query_point, *f_ptr, visibility_arr);
  } else if ((v_ptr = std::get_if<VisibilityArrangement::Vertex_const_handle>(&pl_result))) {
    // Located on vertex.
    // Search the incident halfedge that contains the polygon face.
    VisibilityArrangement::Halfedge_const_handle he = poly.halfedges_begin();
    // Ensure proper iteration and termination for safety
    bool found_he = false;
    for (size_t i = 0; i < poly.number_of_halfedges(); ++i, ++he) { // Iterate safely
        if (he->target() == (*v_ptr) && he->face() == main_face) {
            found_he = true;
            break;
        }
    }
    if (!found_he) {
        std::cout << "Cannot find halfedge corresponding to vertex." << std::endl;
        return false;
    }
    fh = tev.compute_visibility(query_point, he, visibility_arr);
  } else if ((e_ptr = std::get_if<VisibilityArrangement::Halfedge_const_handle>(&pl_result))) {
    // Located on halfedge.
    // Find halfedge that has polygon interior as face.
    VisibilityArrangement::Halfedge_const_handle he =
        (*e_ptr)->face() == main_face ? (*e_ptr) : (*e_ptr)->twin();
    fh = tev.compute_visibility(query_point, he, visibility_arr);
  } else {
      // This case should ideally not be reached if pl_result is a variant of the above types
      std::cout<<"Cannot locate query point on arrangement (unknown variant state)."<<std::endl;
    return false;
  }

  // Result assertion.
  if (fh->is_fictitious()) {
      std::cout<<"Visibility polygon is fictitious."<<std::endl;
    return false;
  }
  if (fh->is_unbounded()) {
      std::cout<<"Visibility polygon is unbounded."<<std::endl;
    return false;
  }

  // Convert to polygon.
  VisibilityArrangement::Ccb_halfedge_circulator curr = fh->outer_ccb();
  *visibility_polygon = Polygon_2();
  do {
    visibility_polygon->push_back(curr->source()->point());
  } while (++curr != fh->outer_ccb());

  simplifyPolygon(visibility_polygon);
  if (visibility_polygon->is_clockwise_oriented())
    visibility_polygon->reverse_orientation();

  return true;
}

}  // namespace polygon_coverage_planning
