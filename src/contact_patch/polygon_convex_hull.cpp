#include "coal/contact_patch/polygon_convex_hull.h"
#include "coal/alloca.h"

#include <algorithm>

namespace coal {
void computePolygonConvexHull(const std::vector<Vec2s>& cloud,
                              std::vector<Vec2s>& cvx_hull) {
  cvx_hull.clear();

  for (const Vec2s& point : cloud) {
    if (!point.allFinite()) return;  // no meaningful hull
  }

  if (cloud.size() <= 2) {
    // Point or segment, nothing to do.
    for (const Vec2s& point : cloud) {
      cvx_hull.emplace_back(point);
    }
    return;
  }

  // temporary copy of cloud to work on
  COAL_MAKE_ALLOCA_BOOST_SPAN(Vec2s, cloud_, cloud.size());
  for (size_t i = 0; i < cloud.size(); ++i) {
    cloud_[i] = cloud[i];
  }

  if (cloud_.size() == 3) {
    // We have a triangle, we only need to arrange it in a counter clockwise
    // fashion.
    //
    // Put the vector which has the lowest y coordinate first.
    if (cloud_[0](1) > cloud_[1](1)) {
      std::swap(cloud_[0], cloud_[1]);
    }
    if (cloud_[0](1) > cloud_[2](1)) {
      std::swap(cloud_[0], cloud_[2]);
    }
    const Vec2s& a = cloud_[0];
    const Vec2s& b = cloud_[1];
    const Vec2s& c = cloud_[2];
    const Scalar det =
        (b(0) - a(0)) * (c(1) - a(1)) - (b(1) - a(1)) * (c(0) - a(0));
    if (det < 0) {
      std::swap(cloud_[1], cloud_[2]);
    }

    for (const Vec2s& point : cloud_) {
      cvx_hull.emplace_back(point);
    }
    return;
  }

  // Andrew's monotone chain: O(n log(n)), only needs a strict lexicographic
  // order, so duplicate and colinear points are handled without special cases.
  // Colinear points are not added to the cvx-hull.
  std::sort(cloud_.begin(), cloud_.end(), [](const Vec2s& p1, const Vec2s& p2) {
    return p1(0) < p2(0) || (p1(0) == p2(0) && p1(1) < p2(1));
  });
  auto isLeftTurn = [](const Vec2s& p1, const Vec2s& p2, const Vec2s& p3) {
    const Scalar det =
        (p2(0) - p1(0)) * (p3(1) - p1(1)) - (p2(1) - p1(1)) * (p3(0) - p1(0));
    return det > Eigen::NumTraits<Scalar>::dummy_precision();
  };
  const size_t n = static_cast<size_t>(
      std::unique(cloud_.begin(), cloud_.end()) - cloud_.begin());
  if (n <= 2) {
    cvx_hull.assign(cloud_.begin(), cloud_.begin() + static_cast<long>(n));
    return;
  }
  cvx_hull.reserve(2 * n);
  // Lower hull, then upper hull; both wrap counter-clockwise.
  for (size_t i = 0; i < n; ++i) {
    while (cvx_hull.size() >= 2 && !isLeftTurn(cvx_hull[cvx_hull.size() - 2],
                                               cvx_hull.back(), cloud_[i])) {
      cvx_hull.pop_back();
    }
    cvx_hull.emplace_back(cloud_[i]);
  }
  const size_t lower_size = cvx_hull.size() + 1;
  for (size_t i = n - 1; i-- > 0;) {
    while (cvx_hull.size() >= lower_size &&
           !isLeftTurn(cvx_hull[cvx_hull.size() - 2], cvx_hull.back(),
                       cloud_[i])) {
      cvx_hull.pop_back();
    }
    cvx_hull.emplace_back(cloud_[i]);
  }
  // The last point is the first one again.
  cvx_hull.pop_back();
}
}  // namespace coal
