using Rhino.Geometry;
using System.Collections.Generic;

namespace VoussoirPlugin02.Components
{
    public static class Utils
    {
        /// <summary>
        /// Ensures that two curves (arcs) in the input list are oriented in the same parameter direction.
        /// If not, reverses the second arc in-place.
        /// </summary>
        public static void OrientArcs(List<Curve> arcs)
        {
            if (arcs == null || arcs.Count != 2)
                return;

            Point3d arc0Start = arcs[0].PointAtStart;
            Point3d arc0End = arcs[0].PointAtEnd;
            Point3d arc1Start = arcs[1].PointAtStart;
            Point3d arc1End = arcs[1].PointAtEnd;

            double dStart = arc0Start.DistanceTo(arc1Start);
            double dEnd = arc0End.DistanceTo(arc1End);
            double dCrossStartEnd = arc0Start.DistanceTo(arc1End);
            double dCrossEndStart = arc0End.DistanceTo(arc1Start);

            if (dStart + dEnd > dCrossStartEnd + dCrossEndStart)
            {
                arcs[1].Reverse();
            }
        }

        /// <summary>
        /// Try to get the four corners of a closed quad curve.
        /// </summary>
        public static bool TryGetQuadCorners(Curve crv, out Point3d p0, out Point3d p1, out Point3d p2, out Point3d p3)
        {
            p0 = p1 = p2 = p3 = Point3d.Unset;

            // 1) Try direct polyline
            if (crv.TryGetPolyline(out Polyline pl) && pl.Count >= 5)
            {
                var corners = new List<Point3d>();
                for (int i = 0; i < pl.Count - 1; i++) corners.Add(pl[i]);
                if (corners.Count != 4)
                {
                    // fallback: approximate by 4 equally spaced points along the curve
                    var t0 = crv.Domain.T0; var t1 = crv.Domain.T1;
                    var pts = new List<Point3d>();
                    for (int i = 0; i < 4; i++) pts.Add(crv.PointAt(t0 + i * (t1 - t0) / 4.0));
                    p0 = pts[0]; p1 = pts[1]; p2 = pts[2]; p3 = pts[3];
                    return true;
                }
                p0 = corners[0]; p1 = corners[1]; p2 = corners[2]; p3 = corners[3];
                return true;
            }

            // 2) Light rebuild and try again
            var dup = crv.DuplicateCurve();
            dup = dup.Rebuild(20, 1, true);
            Polyline pl2;
            if (dup.TryGetPolyline(out pl2) && pl2.Count >= 5)
            {
                p0 = pl2[0]; p1 = pl2[1]; p2 = pl2[2]; p3 = pl2[3];
                return true;
            }
            return false;
        }
        public static List<Point3d> ConvexHull2D(List<Point3d> points, Rhino.Geometry.Plane plane)
        {
            // Project points to plane
            var projected = new List<Point2d>();
            foreach (var pt in points)
            {
                double u, v;
                plane.ClosestParameter(pt, out u, out v);
                projected.Add(new Point2d(u, v));
            }

            // Compute convex hull in 2D (Graham scan)
            projected = RemoveDuplicatePoints(projected, Rhino.RhinoMath.ZeroTolerance);
            var hull2d = ComputeConvexHull2D(projected);

            // Map hull points back to 3D
            var hull3d = new List<Point3d>();
            foreach (var pt2d in hull2d)
            {
                hull3d.Add(plane.PointAt(pt2d.X, pt2d.Y));
            }
            return hull3d;
        }

        // Helper to remove duplicate points
        private static List<Point2d> RemoveDuplicatePoints(List<Point2d> points, double tolerance)
        {
            var unique = new List<Point2d>();
            foreach (var pt in points)
            {
                bool isDuplicate = false;
                foreach (var u in unique)
                {
                    if (pt.DistanceTo(u) < tolerance)
                    {
                        isDuplicate = true;
                        break;
                    }
                }
                if (!isDuplicate)
                    unique.Add(pt);
            }
            return unique;
        }

        // Graham scan convex hull algorithm
        private static List<Point2d> ComputeConvexHull2D(List<Point2d> points)
        {
            if (points.Count < 3)
                return new List<Point2d>(points);

            // Sort points by X, then Y
            points.Sort((a, b) => a.X != b.X ? a.X.CompareTo(b.X) : a.Y.CompareTo(b.Y));

            var hull = new List<Point2d>();

            // Lower hull
            foreach (var pt in points)
            {
                while (hull.Count >= 2 && Cross(hull[hull.Count - 2], hull[hull.Count - 1], pt) <= 0)
                    hull.RemoveAt(hull.Count - 1);
                hull.Add(pt);
            }

            // Upper hull
            int t = hull.Count + 1;
            for (int i = points.Count - 2; i >= 0; i--)
            {
                var pt = points[i];
                while (hull.Count >= t && Cross(hull[hull.Count - 2], hull[hull.Count - 1], pt) <= 0)
                    hull.RemoveAt(hull.Count - 1);
                hull.Add(pt);
            }

            hull.RemoveAt(hull.Count - 1); // Remove duplicate
            return hull;
        }

        private static double Cross(Point2d o, Point2d a, Point2d b)
        {
            return (a.X - o.X) * (b.Y - o.Y) - (a.Y - o.Y) * (b.X - o.X);
        }
    }
}